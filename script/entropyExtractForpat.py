import math
import argparse
import os
from pattools.pat import PatWindow

def calculate_entropy(patWin, depth):
    """
    Calculate the entropy of given patterns.
    
    Parameters:
    - patWin: dict, a dictionary with substring patterns as keys and their counts as values
    - deep: int, the minimum total count required to calculate entropy
    
    Returns:
    - float, entropy value, or -1 if total count is less than `deep`
    """
    given_substrings = {'CCCC': 0, 'CCCT': 0, 'CCTC': 0, 'CCTT': 0, 'CTCC': 0, 'CTCT': 0, 'CTTC': 0, 'CTTT': 0,
                        'TCCC': 0, 'TCCT': 0, 'TCTC': 0, 'TCTT': 0, 'TTCC': 0, 'TTCT': 0, 'TTTC': 0, 'TTTT': 0}

    total_count = 0
    for substring, count in patWin.items():
        if substring in given_substrings:
            given_substrings[substring] += count
            total_count += count   
    if total_count >= depth:
        entropies = {key: (value / total_count) * math.log2((value / total_count)) for key, value in given_substrings.items() if value != 0}
        entropy = -sum(entropies.values()) * 0.25
        if entropy == -0.0:
            entropy = 0.0
    else:
        entropy = -1
    return round(entropy, 4)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input', required=True, help='Input file path')
    parser.add_argument('-o', dest='output_dir', required=True, help='Output directory')
    parser.add_argument('-d', dest='deep', type=int, required=True, help='Total counts for pattern')

    args = parser.parse_args()

    patWindow = PatWindow(args.input)
    base_name = os.path.basename(args.input)
    base_name = os.path.splitext(base_name)[0]
    if base_name.endswith('.pat'):
        base_name = base_name[:-4]
    output_file = os.path.join(args.output_dir, f"{base_name}_entropy.bed")

    with open(output_file, 'w') as f:
        for pat in patWindow:
            entropy = calculate_entropy(pat[2], args.deep)
            f.write(f"{pat[0]}\t{pat[1]}\t{entropy}\n")

if __name__ == '__main__':
    main()