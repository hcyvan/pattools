import pandas as pd
from collections import OrderedDict
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
        counter = line[1]
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
