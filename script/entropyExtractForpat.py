import math
import argparse
import os
from pattools.pat import PatWindow

def calculate_entropy(patWin, deep):
    """
    Calculate the entropy of given patterns.
    
    Parameters:
    - patWin: dict, a dictionary with substring patterns as keys and their counts as values
    - deep: int, the minimum total count required to calculate entropy
    
    Returns:
    - float, entropy value, or -1 if total count is less than `deep`
    """
    total_count = sum(patWin.values())

    if total_count < deep:
        return -1

    entropy = -sum((count / total_count) * math.log2(count / total_count) for count in patWin.values() if count > 0) * 0.25
    if entropy == -0.0:
        entropy = 0.0
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
