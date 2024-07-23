from pattools.vector.utils import parse_file_list
from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output, MotifTabix


def extract_mvs(file_list, region, outfile=None):
    input_files, groups, samples = parse_file_list(file_list)
    with Output(filename=outfile, bgzip=False) as of:
        of.write(f'#chrom\tCpG_index\tgroup\tsample\tmotif\n')
        for i, motif_file in enumerate(input_files):
            with MotifTabix(motif_file, region) as tabix:
                for line in tabix:
                    chrom, cpg_idx, motif_count = line
                    of.write(
                        f'{chrom}\t{cpg_idx}\t{groups[i]}\t{samples[i]}\t{"|".join([str(x) for x in motif_count])}\n')


def extract_vector(input_file, outfile=None, window: int = 4, regions=None):
    motif = Motif(window)
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        with MotifTabix(input_file, regions) as tabix:
            vector_calculator = VectorCalculator(window=window)
            while True:
                line = tabix.readline_and_parse(motif.motifs)
                if not line:
                    break
                chrom, cpg_idx, motif_count = line
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).cluster()
                of.write(f"{vector_calculator}\n")
