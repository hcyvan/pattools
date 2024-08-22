import sys
import gzip
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
            # TODO: The search range for cluster centers should be calculated; currently, it defaults to 1
            center = mvcf.mvw.get_cluster_centers()[mvcf.mvw.get_specific_mvc()]
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
                                hint += v
                        _mvf.readline()
                        mvs_list.append(hint)
                    else:
                        mvs_list.append(0)
            hint_str = "\t".join([str(x) for x in mvs_list])
            of.write(f'{chrom}\t{cpg_idx}\t{hint_str}\n')


def is_no_mvs(s):
    pattern = r'^[0|]*$'
    ret = bool(re.match(pattern, s))
    return ret


def mvc_qc(input_file, outfile):
    with Output(outfile) as of:
        mvc = MvcFormat(input_file)
        window = mvc.header.window
        mvs_total = np.zeros(2 ** window, dtype=np.int64)
        window_not_empty_count = 0
        window_count = 0
        window_size = 0
        chr_count = 0
        chr_current = ""
        with gzip.open(input_file, 'rt') as gzf:
            i = 0
            for line in gzf:
                i += 1
                if i % 10000 == 0:
                    logger.info(f"parse {i} line")
                if line.startswith("#"):
                    continue
                window_count += 1
                line = line.strip()
                items = line.split('\t')
                mvs = items[4]
                chr = items[0]
                genome_start = int(items[2])
                genome_end = int(items[3])
                if chr != chr_current:
                    chr_count += 1

                if window_size == 0:
                    window_size = genome_end - genome_start
                else:
                    window_size = (window_size + genome_end - genome_start) / 2
                if is_no_mvs(mvs):
                    continue
                mvs_np = np.array([int(x) for x in mvs.split('|')], dtype=np.int64)
                mvs_total = mvs_total + mvs_np
                window_not_empty_count += 1
        mvs_agv = mvs_total / window_not_empty_count
        total_str = "|".join([str(x) for x in mvs_total])
        agv_str = "|".join([str(round(x, 2)) for x in mvs_agv])
        of.write("#window\ttotalMVs\tcount\tcountNotEmpty\tavgMVs\twindowBpAvg\n")
        of.write(f"{window}\t{total_str}\t{window_count}\t{window_not_empty_count}\t{agv_str}\t{window_size:.2f}\n")


def _mv_qc(mv_file, sample):
    mvf = MvFormat(mv_file)
    window = mvf.header.window
    mvs_total = np.zeros(2 ** window, dtype=np.int64)
    mvs_total_all = 0
    window_count = 0
    window_not_empty_count = 0
    for line in mvf:
        items = line.split('\t')
        mvs_np = np.array([int(x) for x in items[2].split('|')], dtype=np.int64)
        mvs_np_all = np.sum(mvs_np)
        if mvs_np_all != 0:
            window_not_empty_count += 1
        window_count += 1
        mvs_total = mvs_total + mvs_np
        mvs_total_all += mvs_np_all
    total_str = "|".join([str(x) for x in mvs_total])
    return sample, window, total_str, window_count, window_not_empty_count, mvs_total_all


def mv_qc(input_file, outfile):
    mv_files, col_names = parse_mv_sample_file(input_file, default_col_name='hint')
    logger.info(f"MV QC for  {len(mv_files)} files")
    i = 0
    with Output(outfile) as of:
        of.write("#sample\twindow\ttotalMVs\tcount\tcountNotEmpty\ttotalMVsSum\n")
        for mv_file, col_name in zip(mv_files, col_names):
            i += 1
            logger.info(f"... handling {i}/{len(mv_files)} {mv_file}")
            items = _mv_qc(mv_file, col_name)
            of.write('\t'.join([str(x) for x in items]) + '\n')


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
