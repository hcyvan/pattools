from pattools.io import Output, MotifTabix
from pattools.vector.utils import parse_file_list


def extract_mvs(file_list, region, outfile=None):
    input_files, groups, samples = parse_file_list(file_list)
    with Output(filename=outfile, bgzip=False) as of:
        of.write(f'#chrom\tCpG_index\tgroup\tsample\tmotif\n')
        for i, motif_file in enumerate(input_files):
            with MotifTabix(motif_file, region) as tabix:
                for line in tabix:
                    chrom, cpg_idx, motif_count = line
                    of.write(f'{chrom}\t{cpg_idx}\t{groups[i]}\t{samples[i]}\t{"|".join([str(x) for x in motif_count])}\n')
