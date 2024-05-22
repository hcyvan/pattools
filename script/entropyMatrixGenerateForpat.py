import time
import os
import argparse

parser = argparse.ArgumentParser(description='Generate entropies matrix from *.pat.entropy.bed file of entropyExtractForpat.')
parser.add_argument('-d', '--depth', default='close', type=int, help='The totalC less than depth will bed filter out')
parser.add_argument('-c', '--coordinate', required=True, help='The coordinate bed file')
parser.add_argument('-i', '--input', required=True, help='The input file list to merge')
parser.add_argument('-o', '--output', default='./merge.ratio.bed', help='The output file')
parser.add_argument('-e', '--exclude', default='close',
                    help='exclude -1 mode: all - exclude if all sample is -1; one - exclude if contain one -1;close - '
                         'close exclude mode')
parser.add_argument('-r', '--ratio', default='B',
                    help='''
                    B: double strand ratio; +: plus strand ratio; -: minus strand ratio.
                    eg: -r + for both strand;
                        -r b+ for both and plus strand''')
parser.add_argument('-v', '--version', action='version', version='0.1')

args = parser.parse_args()

depth = args.depth
input_dir = args.input
coordinate = args.coordinate
output_file = args.output
exclude_mode = args.exclude
ratio_type = list(args.ratio)

files = []
for x in open(input_dir, 'r'):
    if len(x.strip()) > 0:
        files.append(x.strip())
header = ["#chrom", "mapinfo", "end"]
header.extend([os.path.basename(x).split('.')[0] for x in files])

fs = []
for file_name in files:
    fs.append(open(file_name, 'r'))
top = []
for i, f in enumerate(fs):
    line = f.readline()
    if line.startswith("#"):
        line = f.readline()
    top.append(line.split('\t'))
merge_ratio_bed = None
merge_ratio_bed_plus = None
merge_ratio_bed_minus = None
if 'B' in ratio_type:
    merge_ratio_bed = open(output_file, 'w')
    merge_ratio_bed.write("{}\n".format('\t'.join(header)))
if '+' in ratio_type:
    merge_ratio_bed_plus = open(output_file, 'w')
if '-' in ratio_type:
    merge_ratio_bed_minus = open(output_file, 'w')
if not any([merge_ratio_bed, merge_ratio_bed_plus, merge_ratio_bed_minus]):
    print("you should set right -r/--ratio")
    exit()

start = time.time()
idx = 0
for line in open(coordinate, 'r'):
    if idx % 1000000 == 0:
        print("index {}, pass {} min".format(idx, round((time.time() - start) / 60), 4))
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
            totalC = int(t[1].strip()) 
            if row[0] == t[0] and int(row[2]) == int(t[1]):
                if merge_ratio_bed:
                    if totalC >= depth:
                        ratio.append(t[2].strip())
                    else:
                        ratio.append('-1')
                        minus1 += 1
                if merge_ratio_bed_plus:
                    if totalC >= depth and int(t[9]) != 0:
                        ratio_plus.append(str(round(int(t[10]) / int(t[9]), 4)))
                    else:
                        ratio_plus.append('-1')
                if totalC >= depth and merge_ratio_bed_minus:
                    if int(t[12]) != 0:
                        ratio_minus.append(str(round(int(t[13]) / int(t[12]), 4)))
                    else:
                        ratio_minus.append('-1')
                line = fs[i].readline()
                if line:
                    top[i] = line.split('\t')
                else:
                    top[i] = None
            else:
                ratio.append('-1')
                minus1 += 1
                ratio_plus.append('-1')
                ratio_minus.append('-1')
    if merge_ratio_bed:
        if (exclude_mode == 'all' and minus1 == len(fs)) or (exclude_mode == 'one' and minus1 >= 1):
            pass
        else:
            merge_ratio_bed.write("{}\t{}\n".format(prefix, '\t'.join(ratio)))
    if merge_ratio_bed_plus:
        merge_ratio_bed_plus.write("{}\t{}\n".format(prefix, '\t'.join(ratio_plus)))
    if merge_ratio_bed_minus:
        merge_ratio_bed_minus.write("{}\t{}\n".format(prefix, '\t'.join(ratio_minus)))
end = time.time()
print()
print("use {} min !!".format((end - start) / 60))
