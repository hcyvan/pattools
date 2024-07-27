import sys
from pattools.vector.utils import *
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, MvTabix, CpG2Tabix
from pattools.vector.format import MvcFormat, MvFormat, BaseHeader
from pattools.utils import is_gzip_file
from pattools.log import logger


def extract_mvs(file_list, regions, outfile=None):
    if is_gzip_file(file_list):
        mvc_files = [file_list]
        col_names = ['mvs']
    else:
        mvc_files, col_names = parse_mvc_group_file(file_list)
    logger.info(f"Check the format of {len(mvc_files)} files")
    header_common = BaseHeader.check_file_list_headers(mvc_files)
    logger.info(f"File info: format={header_common.format}, window={header_common.window}")
    motif = Motif(header_common.window)
    mvc_file_list = []
    for mvc_file in mvc_files:
        if header_common.format == 'mvc':
            mvc = MvcFormat(mvc_file)
        else:
            mvc = MvFormat(mvc_file)
        mvc.readline()
        mvc_file_list.append(mvc)
    with Output(filename=outfile) as of:
        of.write(f"##FORMAT: mvm (methylation vector matrix)\n")
        of.write(f"##WINDOW: {header_common.window}\n")
        of.write(f"##COMMAND: {' '.join(sys.argv)}\n")
        _col_names = "\t".join(col_names)
        of.write(f'#chrom\tcpg\t{_col_names}\n')
        for region in regions:
            chrom, cpg_idx, _ = parse_region_string(region)
            mvs_list = []
            for i, mvc in enumerate(mvc_file_list):
                if mvc.mvw.cpg_idx is None or mvc.mvw.cpg_idx > cpg_idx:
                    mvs_list.append(motif.get_mvs())
                else:
                    while mvc.mvw.cpg_idx < cpg_idx:
                        mvc.readline()
                    if mvc.mvw.cpg_idx == cpg_idx:
                        mvs_list.append(mvc.mvw.mvs)
                        mvc.readline()
                    else:
                        mvs_list.append(motif.get_mvs())
            mvs_group_str = "\t".join(mvs_list)
            of.write(f'{chrom}\t{cpg_idx}\t{mvs_group_str}\n')


def single_cluster(input_file, outfile=None, window: int = 4, regions=None):
    motif = Motif(window)
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        with MvTabix(input_file, regions) as tabix:
            vector_calculator = VectorCalculator(window=window)
            while True:
                line = tabix.readline_and_parse(motif.motifs)
                if not line:
                    break
                chrom, cpg_idx, motif_count = line
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).cluster()
                of.write(f"{vector_calculator.get_mvc_base()}\n")


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
