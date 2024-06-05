import sys
import os
import gzip
import math
import uuid
import multiprocessing
from collections import OrderedDict
from typing import List
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, CpGTabix, MotifTabix


def split_cpg(filename, process=10):
    total = 0
    with gzip.open(filename, 'rt') as f:
        for _ in f:
            total += 1
    batch = math.floor((total + process) / process)
    process_regions = []
    od: OrderedDict[str, List[int]] = OrderedDict()
    with gzip.open(filename, 'rt') as f:
        for i, line in enumerate(f):
            line = line.strip().split('\t')
            chrom = line[0]
            start = int(line[2])
            if chrom in od:
                od[chrom][1] = start
            else:
                od[chrom] = [start, -1]
            if i > 0 and i % batch == 0:
                regions = []
                for k, v in od.items():
                    regions.append(f'{k}:{v[0]}-{v[1]}')
                process_regions.append(regions)
                od = OrderedDict()
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
            vector_calculator = VectorCalculator(window=window, min_vector_proportion_in_cluster=0.33)
            while True:
                line = tabix.readline_and_parse(motif.motifs)
                if not line:
                    break
                chrom, cpg_idx, motif_count = line
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).calc()
                of.write(f"{vector_calculator}\n")


def _extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window: int = 4, regions=None):
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
                vector_calculator = VectorCalculator()
                vector_calculator.set_motif_count(chrom, start, dict())
                for i, (tabix, line) in enumerate(zip(tabix_arr, lines)):
                    if line is None:
                        continue
                    vc = VectorCalculator()
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


def _extract_vector_from_multi_motif_file_process_wrapper(queue, process_order, file_list, cpg_bed, outfile,
                                                          window, regions):
    try:
        _extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window, regions)
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


def extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window: int = 4, process=1):
    sys.stderr.write(f"Process: {process}\n")
    if process == 1:
        _extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window)
    else:
        if outfile is None:
            outfile = f'./merge.{uuid.uuid4()}.motif.gz'
        split_regions = split_cpg(cpg_bed, process)
        process_jobs = []
        queue = multiprocessing.Queue()
        for i, regions in enumerate(split_regions):
            tag = get_label_from_regions(regions)
            process_outfile = f'{outfile}.{i}.{tag}.tmp'
            sys.stderr.write(f"process {i} generating {process_outfile}\n")
            p = multiprocessing.Process(target=_extract_vector_from_multi_motif_file_process_wrapper,
                                        args=(queue, i, file_list, cpg_bed, process_outfile, window, regions))
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
            with Output(filename=outfile, file_format='motif', bgzip=True) as out:
                for i, file in process_files:
                    with open(file) as infile:
                        for line in infile:
                            out.write(line)
        else:
            sys.stderr.write(f"Process Failed\n")
            for i, file in process_files:
                sys.stderr.write(f"Remove {file}\n")
                os.remove(file)
