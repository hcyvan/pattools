from typing import List
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, CpGTabix, MotifTabix


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


def extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window: int = 4, regions=None):
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
