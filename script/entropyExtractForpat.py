import math
import argparse
import os
from pattools.pat import PatWindow

def calculate_entropy(patWin, deep):
    """
    @parm
    patWin = ({'CCCC': 1, 'CCCT': 1, 'CCTC': 2, 'CT.C': 1, 'TCCC': 1, 'TCCT': 1, 'TCTC': 1, 'TTCT': 2, 'TTTC': 1, 'CTCT': 1})

    @return 
    entropy value
    """
    substring_counts = {'CCCC': 0, 'CCCT': 0, 'CCTC': 0, 'CCTT': 0, 'CTCC': 0, 'CTCT': 0, 'CTTC': 0, 'CTTT': 0,
                        'TCCC': 0, 'TCCT': 0, 'TCTC': 0, 'TCTT': 0, 'TTCC': 0, 'TTCT': 0, 'TTTC': 0, 'TTTT': 0}

    total_count = 0
    for substring, count in patWin.items():
        if substring in substring_counts:
            substring_counts[substring] += count
            total_count += count   
    if total_count >= deep:
        entropies = {key: (value / total_count) * math.log2((value / total_count)) for key, value in substring_counts.items() if value != 0}
        entropy = -sum(entropies.values()) * 0.25
        entropy = round(entropy, 4)
    else:
        entropy = -1
    return entropy


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', help='Input file path')
parser.add_argument('-o', dest='output_dir', help='Output dir')
parser.add_argument('-d', dest='deep', help='total counts for pattern')

args = parser.parse_args()

patWindow = PatWindow(args.input)

base_name = os.path.splitext(os.path.basename(args.input))[0]
output_file = f"{args.output_dir}/{base_name}.entropy.bed"

with open(output_file, 'w') as f:
    for pat in patWindow:
        entropy = calculate_entropy(pat[2],int(args.deep))
        f.write(f"{pat[0]}\t{pat[1]}\t{entropy}\n")