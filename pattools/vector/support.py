from pattools.vector.utils import parse_mv_group_sample_file, parse_region_string
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, PatTabix, CpG2Tabix
from pattools.vector.format import MvcFormat, MvFormat


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


def extract_mvc(file_list, region, outfile=None):
    # mvc_files, _ = parse_mvc_group_file(file_list)
    mvc_files = [file_list]
    for _, mvc_file in enumerate(mvc_files):
        mvc = MvcFormat(mvc_file)
        with Output(filename=outfile) as of:
            # of.write(mvc.header.encode() + "\n")
            region_string = region.pop(0)
            chrom, start, end = parse_region_string(region_string)
            for _ in mvc:
                if start == mvc.mvw.cpg_idx:
                    of.write(mvc.mvw.encode() + "\n")
                    if len(region) == 0:
                        break
                    region_string = region.pop(0)
                    chrom, start, end = parse_region_string(region_string)
                elif start > mvc.mvw.cpg_idx:
                    continue


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
