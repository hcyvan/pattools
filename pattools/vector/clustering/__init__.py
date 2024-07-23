import sys
import os
import uuid
import multiprocessing
from collections import OrderedDict
from typing import List
from pattools.vector.clustering.core import do_clustering
from pattools.io import Output, CpG2Tabix
from mpi4py import MPI

__all__ = ['methylation_vector_cluster']


def split_cpg(filename, process=10, region=None):
    """
    This function is used to divide the cpg index into [process] parts.
    """
    total = 0
    with CpG2Tabix(filename, region) as tabix:
        for _ in tabix:
            total += 1
    batch, remainder = divmod(total, process)
    batches = [batch + 1 if p < remainder else batch for p in range(process)]
    cut = [sum(batches[0:(p + 1)]) for p in range(process)]
    process_regions = []
    od: OrderedDict[str, List[int]] = OrderedDict()
    with CpG2Tabix(filename, region) as tabix:
        cut_idx = 0
        for i, (chrom, _, _, start) in enumerate(tabix):
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


def do_mvc_multi_process(queue, process_order, file_list, cpg_bed, outfile, window, regions, cluster='MRESC',
                         groups=None):
    try:
        do_clustering(file_list, cpg_bed, outfile, window, regions, cluster, target_groups=groups)
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
        with open(filenames[0]) as infile:
            for line in infile:
                out.write(line)
        for file in filenames[1:]:
            with open(file) as infile:
                for line in infile:
                    if not line.startswith("#"):
                        out.write(line)
    for file in filenames:
        sys.stderr.write(f"Remove {file}\n")
        os.remove(file)


def methylation_vector_cluster(file_list, cpg_bed, outfile, window: int = 4, process=1, region=None, cluster='HDBSCAN',
                               groups=None):
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
        do_clustering(file_list, cpg_bed, regions_and_filename[1], window, regions_and_filename[0],
                      cluster=cluster, target_groups=groups)
        tmp_files = comm.gather(regions_and_filename[1], root=0)
        if rank == 0:
            merge_split_filenames(outfile, tmp_files)
    else:
        sys.stderr.write(f"Process: {process}\n")
        if process == 1:
            do_clustering(file_list, cpg_bed, outfile, window, region, cluster=cluster, target_groups=groups,
                          out_gzip=True)
        else:
            if outfile is None:
                outfile = f'./merge.{uuid.uuid4()}.motif.gz'
            split_regions = split_cpg(cpg_bed, process, region)
            outfile, split_filenames = get_filename_by_split_regions(split_regions, outfile)
            process_jobs = []
            queue = multiprocessing.Queue()
            for i, regions in enumerate(split_regions):
                sys.stderr.write(f"process {i} generating {split_filenames[i]}\n")
                p = multiprocessing.Process(target=do_mvc_multi_process,
                                            args=(
                                                queue, i, file_list, cpg_bed, split_filenames[i], window, regions,
                                                cluster, groups))
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
