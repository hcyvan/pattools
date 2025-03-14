from pathlib import Path
from pattools.pat import PatWindow
from pattools.io import Output, CpGBedSequential


def calculate_methlevel(patWin, depth):
    """
    Calculate the methylation level of given pat.
    
    Parameters:
    - patWin: dict, a dictionary with substring patterns as keys and their counts as values
    - depth: int, the minimum total count required to calculate methylation level
    
    Returns:
    - float, methylation level value, or -1 if total count is less than `depth`
    """
    cpg_site = {'C': 0, 'T': 0}
    total_count = 0
    for substring, count in patWin.items():
        if substring in cpg_site:
            cpg_site[substring] += count
            total_count += count
    if total_count >= depth:
        meth_ratio = round(cpg_site['C'] / total_count, 4)
    else:
        meth_ratio = -1
    return meth_ratio, cpg_site['C'], cpg_site['T']


def extract_cov(pat_file, genome_path, outfile=None, bgzip: bool = True):
    genome_path = Path(genome_path)
    cpg_bed = genome_path / 'CpG.bed.gz'
    pat_window = PatWindow(pat_file, window=1)
    pat_chr, pat_cpg_idx, pat_count_dict = pat_window.readline()
    with CpGBedSequential(cpg_bed.__str__()) as fin, Output(filename=outfile, bgzip=bgzip) as fo:
        for _, cpg_idx, genome_idx in fin:
            if pat_cpg_idx is None:
                break
            if cpg_idx < pat_cpg_idx:
                continue
            elif cpg_idx == pat_cpg_idx:
                ratio, methylated_cpg, unmethylated_cpg = calculate_methlevel(pat_count_dict, 1)
                if ratio!=-1:
                    fo.write(f"{pat_chr}\t{genome_idx}\t{genome_idx}\t{ratio}\t{methylated_cpg}\t{unmethylated_cpg}\n")
                pat_chr, pat_cpg_idx, pat_count_dict = pat_window.readline()


if __name__ == '__main__':
    pass
