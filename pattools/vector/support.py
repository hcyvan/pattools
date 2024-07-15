from pattools.io import Output, MotifTabix
from pattools.vector.utils import parse_file_list


def extract_motif_from_region(file_list, region, outfile=None):
    input_files, groups, samples = parse_file_list(file_list)
    with Output(filename=outfile, bgzip=False) as of:
        of.write(f'#group\tsample\tmotif\n')
        for i, motif_file in enumerate(input_files):
            with MotifTabix(motif_file, region) as tabix:
                line = tabix.readline()
                chrom, cpg_idx, motif_count = line
                of.write(f'{groups[i]}\t{samples[i]}\t{"|".join([str(x) for x in motif_count])}\n')
