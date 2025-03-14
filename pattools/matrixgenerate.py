import os
import gzip
from pattools.io import Output, CpGBedSequential


def matrix_generate(input_dir, coordinate, depth, exclude_mode, output_file,is_gzip=False):
    files = []
    for x in open(input_dir, 'r'):
        if len(x.strip()) > 0:
            files.append(x.strip())
    header = ["#chrom", "start", "end", "index"]
    header.extend([os.path.basename(x).split('.')[0] for x in files])

    fs = []
    for file_name in files:
        fs.append(gzip.open(file_name, 'rb'))
    top = []
    for i, f in enumerate(fs):
        line = f.readline()
        decoded_line = line.decode('utf-8')
        top.append(decoded_line.split('\t'))
    with CpGBedSequential(str(coordinate)) as fin, Output(filename=output_file, bgzip=is_gzip) as fo:
        fo.write("{}\n".format('\t'.join(header)))
        for chrom, cpg_idx, genome_idx in fin:
            ratio = []
            ratio_plus = []
            ratio_minus = []
            minus1 = 0
            for i, t in enumerate(top):
                if t is None:
                    ratio.append('-1')
                    minus1 += 1
                    ratio_plus.append('-1')
                    ratio_minus.append('-1')
                else:
                    totalC = int(t[3].strip())
                    if chrom == t[0] and cpg_idx == int(t[1]):
                        if totalC >= depth:
                            ratio.append(t[2].strip())
                        else:
                            ratio.append('-1')
                            minus1 += 1
                        line = fs[i].readline().decode('utf-8')
                        if line:
                            top[i] = line.split('\t')
                        else:
                            top[i] = None
                    else:
                        ratio.append('-1')
                        minus1 += 1
            if (exclude_mode == 'all' and minus1 == len(fs)) or (exclude_mode == 'one' and minus1 >= 1):
                pass
            else:
                o_str = "{}\t{}\n".format(f"{chrom}\t{genome_idx}\t{genome_idx + 2}\t{cpg_idx}", '\t'.join(ratio))
                fo.write(o_str)
