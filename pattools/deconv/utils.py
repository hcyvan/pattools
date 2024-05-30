import pandas as pd
from collections import OrderedDict,Counter
from pattools.pat import PatRegion


def get_methylation_density_from_pat_by_cpg_idx(pat_file, genome_cpg_regions: OrderedDict):
    """
    Get methylation density from pat file.
    Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing for noninvasive prenatal, cancer,
     and transplantation assessments.

    Note: If the region represents 1 CpG, it is equivalent to calculating the methylation level of the CpG,
     which is called the beta value.

    :param pat_file:
    :param genome_cpg_regions: ordered dict of Genomic index => CpG index
    :return: ordered dict of Genomic index => methylation density
    """
    genome_cpg_regions_df = pd.DataFrame([[k, v] for k, v in genome_cpg_regions.items()])
    genome_cpg_regions_df.columns = ['genome', 'cpg']
    cpg_methylation_density = []
    pat = PatRegion(pat_file, genome_cpg_regions_df['cpg'].to_list())
    for _, line in enumerate(pat):
        counter = Counter()
        for item in list(line[1].elements()):
            counter += Counter(item)
        total = counter["T"] + counter["C"]
        if total:
            methy_density = (counter["C"] / total)
        else:
            methy_density = -1
        cpg_methylation_density.append([line[0], methy_density])
    cpg_methylation_density = pd.DataFrame(cpg_methylation_density)
    cpg_methylation_density.columns = ['cpg', 'methylation_density']
    ret = pd.merge(genome_cpg_regions_df, cpg_methylation_density, on=['cpg'], how='left')
    return ret[ret['methylation_density'] != -1]

def get_uxm_ratio_from_pat_by_cpg_idx(pat_file, genome_cpg_regions: OrderedDict,rlen=4):
    """
    Get uxm_ratio from pat file.
    Loyfer et al., 2023, A DNA methylation atlas of normal human cell types, Nature

    :param pat_file:
    :param genome_cpg_regions: ordered dict of Genomic index => CpG index
    :return: ordered dict of Genomic index => methylation density
    """
    genome_cpg_regions_df = pd.DataFrame([[k, v] for k, v in genome_cpg_regions.items()])
    genome_cpg_regions_df.columns = ['genome', 'cpg']
    cpg_uxm_ratio = []
    pat = PatRegion(pat_file, genome_cpg_regions_df['cpg'].to_list())
    for _, line in enumerate(pat):
        uxm = {'U':0,'M':0}
        counter = Counter({k: v for k, v in line[1].items() if k.count('C') + k.count('T') >= rlen})
        if counter :
            for k,v in counter.items():
                cont = Counter(k)
                if cont['C']/(cont['C']+cont['T']) <= 0.25:
                    uxm['U'] = uxm['U'] + v
                else:
                    uxm['M'] = uxm['M'] + v
            depth = uxm['U']+uxm['M'] 
            uxm_ratio = uxm['U']/depth
        else:
            uxm_ratio = -1
            depth = 0
        cpg_uxm_ratio.append([line[0], uxm_ratio,depth])

    cpg_uxm_ratio = pd.DataFrame(cpg_uxm_ratio)
    cpg_uxm_ratio.columns = ['cpg', 'uxm_ratio','depth']
    ret = pd.merge(genome_cpg_regions_df, cpg_uxm_ratio, on=['cpg'], how='left')
    return ret[ret['uxm_ratio'] != -1]

