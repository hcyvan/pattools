import sys
from typing import List
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, CpG2Tabix, PatTabix
from pattools.vector.utils import parse_mv_group_sample_file
from pattools.vector.format import MvFormat
from pattools.log import logger


def do_clustering(file_list, cpg_bed, outfile, window: int = None, regions=None, cluster='MRESC', out_gzip=False,
                  target_groups=None):
    """
    This is the core function to do methylation vectors clustering
    @param file_list:
    @param cpg_bed:
    @param outfile:
    @param window:[DEPRECATED] The window size of methylation vectors
    @param regions:
    @param cluster:
    @param out_gzip:
    @param target_groups: The groups selected from the MV list file should be separated by commas ','. If not set, all
                        groups will be used.
    """
    input_files, groups, samples = parse_mv_group_sample_file(file_list, target_groups)
    _window = window
    for motif_file in input_files:
        mvf = MvFormat(motif_file)
        if _window is None:
            _window = mvf.header.window
        else:
            if _window != mvf.header.window:
                logger.error(f'Window size is not the same: {_window} [others] and {mvf.header.window} [{motif_file}]')
                raise Exception('Window size error')
    window = _window
    logger.info(f"Window size: {window}")
    motif = Motif(window)
    tabix_arr: List[PatTabix] = []
    lines = []
    for motif_file in input_files:
        tabix = PatTabix(motif_file, regions)
        tabix_arr.append(tabix)
        lines.append(tabix.readline_and_parse(motif.motifs))
    with Output(filename=outfile, file_format='mvc', bgzip=out_gzip) as of:
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
        with CpG2Tabix(cpg_bed, regions, window=window) as cpg:
            for chrom, genome_start, genome_end, cpg_idx in cpg:
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
                vector_calculator.cluster()
                of.write(f"{vector_calculator.get_mvc(group_order, genome_start, genome_end)}\n")

    for tabix in tabix_arr:
        tabix.close()
