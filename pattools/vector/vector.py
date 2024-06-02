import gzip
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output


def extract_vector(input_file, outfile=None):
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        motif_array = None
        with gzip.open(input_file, 'rb') as f:
            vector_calculator = VectorCalculator()
            vector_calculator.set_param_min_vector_proportion_in_cluster(0.33)
            for line in f:
                if not line:
                    continue
                line = line.decode().strip()
                if line.startswith('##'):
                    if line.startswith('##WINDOW'):
                        vector_calculator.set_param_window(int(line.split(' ')[1]))
                    continue
                if line.startswith('#chr'):
                    motif_array = line.split('\t')[2:]
                    continue
                chrom = line.split('\t')[0]
                cpg_idx = line.split('\t')[1]
                motif_count_arr = [int(x) for x in line.split('\t')[2:]]
                motif_count = dict(zip(motif_array, motif_count_arr))
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).calc()
                of.write(f"{vector_calculator}\n")
                # if vector_calculator.get_clusters_number() >= 2:
                #     if vector_calculator.distance_between_top_two_cluster() >= 0.90:
                #         of.write(f"{vector_calculator}\n")
