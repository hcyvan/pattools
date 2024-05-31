import gzip
from pattools.pat import PatWindow
from pattools.io import Output

def calculate_methlevel(patWin, depth):
    """
    Calculate the methylation level of given pat.
    
    Parameters:
    - patWin: dict, a dictionary with substring patterns as keys and their counts as values
    - depth: int, the minimum total count required to calculate methylation level
    
    Returns:
    - float, methylation level value, or -1 if total count is less than `depth`
    """
    cpgsite = {'C': 0, 'T': 0}
    total_count = 0
    for substring, count in patWin.items():
        if substring in cpgsite:
            cpgsite[substring] += count
            total_count += count 
    if total_count >= depth:
        meth_ratio = cpgsite['C'] / total_count   
    else:
        meth_ratio = -1
    return round(meth_ratio, 4), total_count

def extract_ratio(input, depth, outfile, bgzip: bool = True):

    patWindow = PatWindow(input, window=1)
    with Output(filename=outfile, bgzip=bgzip) as f:
        for pat in patWindow:
            ratio, count = calculate_methlevel(pat[2], depth)
            f.write(f"{pat[0]}\t{pat[1]}\t{ratio}\t{count}\n")

if __name__ == '__main__':
    pass
