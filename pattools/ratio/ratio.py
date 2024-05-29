import gzip
from pattools.pat import PatWindow

def calculate_methlevel(patWin, depth):
    """
    Calculate the methylation level of given pat.
    
    Parameters:
    - patWin: dict, a dictionary with substring patterns as keys and their counts as values
    - deep: int, the minimum total count required to calculate methylation level
    
    Returns:
    - float, methylation level value, or -1 if total count is less than `deep`
    """
    cpgsite = {'C': 0, 'T': 0}
    total_count = 0
    for substring, count in patWin.items():
        if substring in cpgsite:
            cpgsite[substring] += count
            total_count += count 
    if total_count >= int(depth):
        meth_ratio = cpgsite['C'] / total_count
        return round(meth_ratio, 4)
    else:
        return -1

def extract_ratio(input, depth, out):
    if not out.endswith('.gz'):
        out += '.gz'

    patWindow = PatWindow(input, window=1)
    with gzip.open(out, 'wt') as f:
        for pat in patWindow:
            ratio = calculate_methlevel(pat[2], depth)
            f.write(f"{pat[0]}\t{pat[1]}\t{ratio}\n")

if __name__ == '__main__':
    pass
