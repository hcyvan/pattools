import sys
from pattools.vector.utils import *
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, PatTabix, CpG2Tabix
from pattools.vector.format import MvcFormat, MvFormat
from pattools.utils import is_gzip_file
from pattools.log import logger


def extract_mvs(file_list, region, outfile=None):
    input_files, groups, samples = parse_mv_group_sample_file(file_list)
    with Output(filename=outfile, bgzip=False) as of:
        of.write(f'#chrom\tCpG_index\tgroup\tsample\tmotif\n')
        for i, motif_file in enumerate(input_files):
            with PatTabix(motif_file, region) as tabix:
                for line in tabix:
                    chrom, cpg_idx, motif_count = line
                    of.write(
                        f'{chrom}\t{cpg_idx}\t{groups[i]}\t{samples[i]}\t{"|".join([str(x) for x in motif_count])}\n')


def extract_mvc(file_list, regions, outfile=None):
    if is_gzip_file(file_list):
        mvc_files = [file_list]
        groups = ['group']
    else:
        mvc_files, groups = parse_mvc_group_file(file_list)
    mvc_file_list = []
    _window = None
    for mvc_file in mvc_files:
        mvc = MvcFormat(mvc_file)
        if _window is None:
            _window = mvc.header.window
        else:
            if _window != mvc.header.window:
                logger.error(f'Window size is not the same: {_window} [others] and {mvc.header.window} [{mvc_file}]')
                raise Exception('Window size error')
        mvc.readline()
        mvc_file_list.append(mvc)
    with Output(filename=outfile) as of:
        of.write(f"##FORMAT: mvm (methylation vector matrix)\n")
        of.write(f"##WINDOW: {_window}\n")
        of.write(f"##COMMAND: {' '.join(sys.argv)}\n")
        group_str = "\t".join(groups)
        of.write(f'#chrom\tcpg\tstart\tend\t{group_str}\n')
        for region in regions:
            chrom, cpg_idx, _ = parse_region_string(region)
            start = ""
            end = ""
            mvs_list = []
            for i, mvc in enumerate(mvc_file_list):
                if mvc.mvw.cpg_idx is None or mvc.mvw.cpg_idx > cpg_idx:
                    mvs_list.append("")
                else:
                    while mvc.mvw.cpg_idx < cpg_idx:
                        mvc.readline()
                    if mvc.mvw.cpg_idx == cpg_idx:
                        mvs_list.append(mvc.mvw.mvs)
                        start = mvc.mvw.genome_start
                        end = mvc.mvw.genome_end
                        mvc.readline()
            mvs_group_str = "\t".join(mvs_list)
            of.write(f'{chrom}\t{cpg_idx}\t{start}\t{end}\t{mvs_group_str}\n')


def extract_vector(input_file, outfile=None, window: int = 4, regions=None):
    motif = Motif(window)
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        with PatTabix(input_file, regions) as tabix:
            vector_calculator = VectorCalculator(window=window)
            while True:
                line = tabix.readline_and_parse(motif.motifs)
                if not line:
                    break
                chrom, cpg_idx, motif_count = line
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).cluster()
                of.write(f"{vector_calculator}\n")


def fix_mvc(mvc_file, cpg_bed=None, out=None):
    with Output(filename=out, file_format='cgs', bgzip=True) as of:
        mvc = MvcFormat(mvc_file)
        of.write(mvc.header.encode() + "\n")
        with CpG2Tabix(cpg_bed, window=mvc.header.window) as cpg:
            _line = mvc.readline()
            for chrom, genome_start, genome_end, cpg_idx in cpg:
                if _line is None:
                    break
                if mvc.mvw.cpg_idx == cpg_idx:
                    of.write(mvc.update_and_encode(genome_start - 1, genome_end) + "\n")
                    _line = mvc.readline()
                elif mvc.mvw.cpg_idx > cpg_idx:
                    continue
                else:
                    _line = mvc.readline()


def fix_mv(mv_file, out=None):
    with Output(filename=out, file_format='cgs', bgzip=True) as of:
        mv = MvFormat(mv_file)
        mv.header.fix_header()
        of.write(mv.header.encode() + "\n")
        for _ in mv:
            of.write(mv.mvw.encode() + '\n')
