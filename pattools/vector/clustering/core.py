import sys
from typing import List
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, CpGTabix, MotifTabix
from pattools.vector.utils import parse_file_list


def do_clustering(file_list, cpg_bed, outfile, window: int = 4, regions=None, cluster='HDBSCAN',
                  out_version='v1', out_gzip=False):
    input_files, groups, samples = parse_file_list(file_list)
    motif = Motif(window)
    tabix_arr: List[MotifTabix] = []
    lines = []
    for motif_file in input_files:
        tabix = MotifTabix(motif_file, regions)
        tabix_arr.append(tabix)
        lines.append(tabix.readline_and_parse(motif.motifs))
    with Output(filename=outfile, file_format='mvc', bgzip=out_gzip) as of:
        group_order = None
        if out_version == 'v2':
            group_order = sorted(list(set(groups)))
            sample_order = sorted(list(set(samples)))
            of.write(f"##FORMAT: mvc\n")
            of.write(f"##WINDOW: {window}\n")
            of.write(f"##COMMAND: {' '.join(sys.argv)}\n")
            of.write(f"##GROUP: {','.join(group_order)}\n")
            of.write(f"##SAMPLE: {','.join(sample_order)}\n")
            for group in group_order:
                group_sample = []
                for s, g in zip(samples, groups):
                    if group == g:
                        group_sample.append(s)
                of.write(f"##GROUP_SAMPLE: {group}_{','.join(group_sample)}\n")
            of.write(
                f"#chr\tcpg\tstart\tend\tmvs\tc_num\tc_center\tc_group_mvs_num\tc_group_samples_num\tc_group_samples\n")
        with CpGTabix(cpg_bed, regions) as cpg:
            for chrom, genome_idx, cpg_idx in cpg:
                vector_calculator = VectorCalculator(window=window, cluster=cluster)
                vector_calculator.set_motif_count(chrom, cpg_idx, dict())
                for i, (tabix, line) in enumerate(zip(tabix_arr, lines)):
                    if line is None:
                        continue
                    vc = VectorCalculator(window=window, cluster=cluster)
                    while vector_calculator > vc:
                        line = tabix.readline_and_parse(motif.motifs)
                    vc.set_motif_count(line[0], line[1], line[2], sample=samples[i], group=groups[i])
                    if vector_calculator == vc:
                        vector_calculator = vector_calculator + vc
                        lines[i] = tabix.readline_and_parse(motif.motifs)
                vector_calculator.calc()
                vector_calculator.calc_labels_groups_samples()
                if out_version == 'v1':
                    of.write(f"{vector_calculator.__str__()}\t{vector_calculator.get_labels_groups_samples_str()}\n")
                else:
                    of.write(f"{vector_calculator.info_v2(group_order, genome_idx)}\n")

    for tabix in tabix_arr:
        tabix.close()
