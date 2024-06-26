import sys
import os
import uuid
import multiprocessing
from collections import OrderedDict
from typing import List
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, CpGTabix, MotifTabix
from mpi4py import MPI


def split_cpg(filename, process=10, region=None):
    """
    This function is used to divide the cpg index into [process] parts.
    """
    total = 0
    with CpGTabix(filename, region) as tabix:
        for _ in tabix:
            total += 1
    batch, remainder = divmod(total, process)
    batches = [batch + 1 if p < remainder else batch for p in range(process)]
    cut = [sum(batches[0:(p + 1)]) for p in range(process)]
    process_regions = []
    od: OrderedDict[str, List[int]] = OrderedDict()
    with CpGTabix(filename, region) as tabix:
        cut_idx = 0
        for i, (chrom, _, start) in enumerate(tabix):
            if chrom in od:
                od[chrom][1] = start
            else:
                od[chrom] = [start, start]
            if i == cut[cut_idx] - 1:
                regions = []
                for k, v in od.items():
                    regions.append(f'{k}:{v[0]}-{v[1]}')
                process_regions.append(regions)
                od = OrderedDict()
                cut_idx += 1
        if len(od):
            regions = []
            for k, v in od.items():
                regions.append(f'{k}:{v[0]}-{v[1]}')
            process_regions.append(regions)
    return process_regions


def extract_vector(input_file, outfile=None, window: int = 4, regions=None):
    motif = Motif(window)
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        with MotifTabix(input_file, regions) as tabix:
            vector_calculator = VectorCalculator(window=window)
            while True:
                line = tabix.readline_and_parse(motif.motifs)
                if not line:
                    break
                chrom, cpg_idx, motif_count = line
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).calc()
                of.write(f"{vector_calculator}\n")


def extract_vector_multi(file_list, cpg_bed, outfile, window: int = 4, regions=None, cluster='HDBSCAN'):
    input_files = []
    groups = []
    samples = []
    with open(file_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            input_files.append(line[0])
            groups.append(line[1])
            samples.append(line[2])

    motif = Motif(window)
    tabix_arr: List[MotifTabix] = []
    lines = []
    for motif_file in input_files:
        tabix = MotifTabix(motif_file, regions)
        tabix_arr.append(tabix)
        lines.append(tabix.readline_and_parse(motif.motifs))
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        with CpGTabix(cpg_bed, regions) as cpg:
            for chrom, _, start in cpg:
                vector_calculator = VectorCalculator(window=window, cluster=cluster)
                vector_calculator.set_motif_count(chrom, start, dict())
                for i, (tabix, line) in enumerate(zip(tabix_arr, lines)):
                    if line is None:
                        continue
                    vc = VectorCalculator(window=window, cluster=cluster)
                    vc.set_motif_count(line[0], line[1], line[2], sample=samples[i], group=groups[i])
                    while vector_calculator > vc:
                        line = tabix.readline_and_parse(motif.motifs)
                        vc.set_motif_count(line[0], line[1], line[2], sample=samples[i], group=groups[i])
                    if vector_calculator == vc:
                        vector_calculator = vector_calculator + vc
                        lines[i] = tabix.readline_and_parse(motif.motifs)
                vector_calculator.calc().calc_labels_groups_samples()
                of.write(f"{vector_calculator.__str__()}\t{vector_calculator.get_labels_groups_samples_str()}\n")

    for tabix in tabix_arr:
        tabix.close()


def extract_vector_multi_process(queue, process_order, file_list, cpg_bed, outfile, window, regions, cluster='HDBSCAN'):
    try:
        extract_vector_multi(file_list, cpg_bed, outfile, window, regions, cluster)
        queue.put((process_order, 'success', outfile))
    except Exception as e:
        queue.put((process_order, 'failure', str(e)))


def get_label_from_regions(regions):
    start = ""
    end = ""
    for region in regions:
        if len(start) == 0:
            start, end = region.split(':')[1].split('-')
        else:
            _, end = region.split(':')[1].split('-')
    return f"{start}_{end}"


def get_filename_by_split_regions(split_regions, outfile):
    if outfile is None:
        outfile = f'./merge.{uuid.uuid4()}.motif.gz'
    filenames = []
    for i, regions in enumerate(split_regions):
        tag = get_label_from_regions(regions)
        filenames.append(f'{outfile}.{i}.{tag}.tmp')
    return outfile, filenames


def merge_split_filenames(outfile, filenames):
    with Output(filename=outfile, file_format='motif', bgzip=True) as out:
        for file in filenames:
            with open(file) as infile:
                for line in infile:
                    out.write(line)
    for file in filenames:
        sys.stderr.write(f"Remove {file}\n")
        os.remove(file)


def extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window: int = 4, process=1, region=None,
                                         cluster='HDBSCAN'):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    task_count = comm.Get_size()
    if task_count > 1:
        if rank == 0:
            sys.stderr.write(f"MPI Size: {task_count}\n")
            split_regions = split_cpg(cpg_bed, task_count, region)
            outfile, split_filenames = get_filename_by_split_regions(split_regions, outfile)
            split_regions_and_filenames = zip(split_regions, split_filenames)
        else:
            split_regions_and_filenames = None
            outfile = None
        regions_and_filename = comm.scatter(split_regions_and_filenames, root=0)
        extract_vector_multi(file_list, cpg_bed, regions_and_filename[1], window, regions_and_filename[0],
                             cluster=cluster)
        tmp_files = comm.gather(regions_and_filename[1], root=0)
        if rank == 0:
            merge_split_filenames(outfile, tmp_files)
    else:
        sys.stderr.write(f"Process: {process}\n")
        if process == 1:
            extract_vector_multi(file_list, cpg_bed, outfile, window, region, cluster=cluster)
        else:
            if outfile is None:
                outfile = f'./merge.{uuid.uuid4()}.motif.gz'
            split_regions = split_cpg(cpg_bed, process, region)
            outfile, split_filenames = get_filename_by_split_regions(split_regions, outfile)
            process_jobs = []
            queue = multiprocessing.Queue()
            for i, regions in enumerate(split_regions):
                sys.stderr.write(f"process {i} generating {split_filenames[i]}\n")
                p = multiprocessing.Process(target=extract_vector_multi_process,
                                            args=(
                                                queue, i, file_list, cpg_bed, split_filenames[i], window, regions,
                                                cluster))
                process_jobs.append(p)
                p.start()

            for i, job in enumerate(process_jobs):
                job.join()

            success = True
            process_files = []
            while not queue.empty():
                process_order, status, ret = queue.get()
                if status == 'success':
                    process_files.append([process_order, ret])
                else:
                    success = False
            process_files = sorted(process_files, key=lambda x: x[0])
            if success:
                merge_split_filenames(outfile, [x[1] for x in process_files])
            else:
                sys.stderr.write(f"Process Failed\n")
