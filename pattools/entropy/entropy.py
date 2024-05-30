import math
import gzip
from pattools.pat import PatWindow
from pattools.motif import Motif
def calculate_entropy(patWin, depth, window):
    """
    Calculate the entropy of given patterns.
    
    Parameters:
    - patWin: dict, a dictionary with substring patterns as keys and their counts as values
    - depth: int, the minimum total count required to calculate entropy
    
    Returns:
    - float, entropy value, or -1 if total count is less than `depth`
    """

    motif = Motif(window)
    given_substrings = motif.motif2order_counter()

    total_count = 0
    for substring, count in patWin.items():
        if substring in given_substrings:
            given_substrings[substring] += count
            total_count += count   
    if total_count >= int(depth):
        entropies = {key: (value / total_count) * math.log2((value / total_count)) for key, value in given_substrings.items() if value != 0}
        entropy = -sum(entropies.values()) * 0.25
        if entropy == -0.0:
            entropy = 0.0
    else:
        entropy = -1
    return round(entropy, 4), total_count

def extract_entropy(input, depth, window, out):
    if not out.endswith('.gz'):
        out += '.gz'
        
    patWindow = PatWindow(input, window)
    with gzip.open(out, 'wt') as f:
        for pat in patWindow:
            entropy, count = calculate_entropy(pat[2], depth, window)
            f.write(f"{pat[0]}\t{pat[1]}\t{entropy}\t{count}\n")

if __name__ == '__main__':
    pass