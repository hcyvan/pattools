import sys
import numpy as np
from pattools.vector.utils import *
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, MvTabix, CpG2Tabix
from pattools.vector.format import MvcFormat, MvFormat, BaseHeader
from pattools.log import logger
from pattools.vector.utils import get_cpg_index_regions


def extract_mvs(file_list, region_file, outfile=None):
    regions = get_cpg_index_regions(region_file)
    mv_files, col_names = parse_mv_sample_file(file_list)
    logger.info(f"Check the format of {len(mv_files)} files")
    header_common = BaseHeader.check_file_list_headers(mv_files)
    logger.info(f"File info: format={header_common.format}, window={header_common.window}")
    motif = Motif(header_common.window)
    mvf_arr = []
    for mv_file in mv_files:
        if header_common.format == 'mvc':
            mvf = MvcFormat(mv_file)
        else:
            mvf = MvFormat(mv_file)
        mvf.readline()
        mvf_arr.append(mvf)
    with Output(filename=outfile) as of:
        of.write(f"##FORMAT: mvm (methylation vector matrix)\n")
        of.write(f"##WINDOW: {header_common.window}\n")
        of.write(f"##COMMAND: {' '.join(sys.argv)}\n")
        _col_names = "\t".join(col_names)
        of.write(f'#chrom\tcpg\t{_col_names}\n')
        for region in regions:
            chrom, cpg_idx, _ = parse_region_string(region)
            mvs_list = []
            for i, mvc in enumerate(mvf_arr):
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


def find_motifs(input_file, mvc_file, outfile=None):
    mv_files, col_names = parse_mv_sample_file(input_file, default_col_name='hint')
    logger.info(f"Check the format of {len(mv_files)} files")
    header = BaseHeader.check_file_list_headers(mv_files)
    logger.info(f"File info: format={header.format}, window={header.window}")
    mvf_arr = []
    for mv_file in mv_files:
        mvf = MvFormat(mv_file)
        mvf.readline()
        mvf_arr.append(mvf)
    mvcf = MvcFormat(mvc_file)
    with Output(filename=outfile) as of:
        of.write(f"##FORMAT: mvh (methylation vector hint)\n")
        of.write(f"##WINDOW: {header.window}\n")
        of.write(f"##COMMAND: {' '.join(sys.argv)}\n")
        _col_names = "\t".join(col_names)
        of.write(f'#chrom\tcpg\t{_col_names}\n')
        for _ in mvcf:
            # TODO: The specific MV cluster should be marked in `mv-separating`. However, this feature was not implemented
            #  in previous development versions, so the current default is to check the second cluster. This issue will
            #  be addressed in future updates.
            # TODO: The search range for cluster centers should be calculated; currently, it defaults to 1
            center = mvcf.mvw.get_cluster_centers()[1]
            chrom = mvcf.mvw.chrom
            cpg_idx = mvcf.mvw.cpg_idx
            if len(center) != header.window:
                err = f'The dimension of the center, {len(center)}, does not match the size of the window, {header.window}'
                logger.error(err)
                raise Exception(err)

            mvs_list = []
            for i, _mvf in enumerate(mvf_arr):
                if _mvf.mvw.cpg_idx is None or _mvf.mvw.cpg_idx > cpg_idx:
                    mvs_list.append(0)
                else:
                    while _mvf.mvw.cpg_idx < cpg_idx:
                        _mvf.readline()
                    if _mvf.mvw.cpg_idx == cpg_idx:
                        motif = Motif(header.window)
                        motif_count = motif.mvs2motif_count(_mvf.mvw.mvs, remove_0=True)
                        hint = 0
                        for k, v in motif_count.items():
                            obj = motif.motif2vector(k)
                            dist = np.linalg.norm(np.array(center) - np.array(obj))
                            if dist < 1:
                                # print(mvcf.mvw.chrom, mvcf.mvw.cpg_idx, dist, v, obj, center)
                                hint += v
                        _mvf.readline()
                        mvs_list.append(hint)
                    else:
                        mvs_list.append(0)
            hint_str = "\t".join([str(x) for x in mvs_list])
            of.write(f'{chrom}\t{cpg_idx}\t{hint_str}\n')


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
