import gzip

from pattools.motif import Motif
from pattools.vector.calculator import VectorCalculator
from pattools.io import Output


def parse_line(line: str):
    chrom = line.split('\t')[0]
    cpg_idx = int(line.split('\t')[1])
    motif_count_arr = [int(x) for x in line.split('\t')[2:]]
    return chrom, cpg_idx, motif_count_arr


def read_line(fin, motif: Motif):
    while True:
        line = fin.readline()
        if line:
            line = line.decode().strip()
            if line.startswith('#'):
                continue
            else:
                chrom, cpg_idx, motif_count_arr = parse_line(line)
                motif_count = dict(zip(motif.motifs, motif_count_arr))
                return (chrom, cpg_idx, motif_count)
        else:
            return None
    return None


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
                chrom, cpg_idx, motif_count_arr = parse_line(line)
                motif_count = dict(zip(motif_array, motif_count_arr))
                vector_calculator.set_motif_count(chrom, cpg_idx, motif_count).calc()
                of.write(f"{vector_calculator}\n")
                # if vector_calculator.get_clusters_number() >= 2:
                #     if vector_calculator.distance_between_top_two_cluster() >= 0.90:
                #         of.write(f"{vector_calculator}\n")


def extract_vector_from_multi_motif_file(file_list, cpg_bed, outfile, window: int = 4):
    input_files = []
    groups = []
    samples = []
    with open(file_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            input_files.append(line[0])
            groups.append(line[1])
            samples.append(line[2])

    motif = Motif(window)
    fins = []
    lines = []

    for motif_file in input_files:
        fin = gzip.open(motif_file)
        fins.append(fin)
        line = read_line(fin, motif)
        lines.append(line)
    with Output(filename=outfile, file_format='motif', bgzip=False) as of:
        with gzip.open(cpg_bed, 'rt') as fc:
            for bed in fc:
                bed = bed.strip()
                chrom = bed.split('\t')[0]
                start = int(bed.split('\t')[2])
                vector_calculator = VectorCalculator()
                vector_calculator.set_motif_count(chrom, start, dict())
                for i, (fin, line) in enumerate(zip(fins, lines)):
                    if line is None:
                        continue
                    vc = VectorCalculator()
                    vc.set_motif_count(line[0], line[1], line[2], sample=samples[i], group=groups[i])
                    while vector_calculator > vc:
                        line = read_line(fin, motif)
                        vc.set_motif_count(line[0], line[1], line[2], group=groups[i])
                    if vector_calculator == vc:
                        vector_calculator = vector_calculator + vc
                        lines[i] = read_line(fin, motif)
                vector_calculator.calc().calc_labels_groups_samples()
                of.write(f"{vector_calculator.__str__()}\t{vector_calculator.get_labels_groups_samples_str()}\n")

    for fin in fins:
        fin.close()
