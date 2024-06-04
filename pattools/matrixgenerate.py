import os
import gzip

def matrix_generate(input_dir, coordinate, depth, exclude_mode, output_file):
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

    merge_ratio_bed = open(output_file, 'w')
    merge_ratio_bed.write("{}\n".format('\t'.join(header)))

    idx = 0
    for line in open(coordinate, 'r'):
        idx += 1
        prefix = line.strip()
        row = prefix.split("\t")
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
                if row[0] == t[0] and int(row[3]) == int(t[1]):
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

        if merge_ratio_bed:
            if (exclude_mode == 'all' and minus1 == len(fs)) or (exclude_mode == 'one' and minus1 >= 1):
                pass
            else:
                merge_ratio_bed.write("{}\t{}\n".format(prefix, '\t'.join(ratio)))
