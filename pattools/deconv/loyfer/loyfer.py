import pandas as pd
import numpy as np
from collections import Counter 
from pattools.region import GenomicRegion
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.utils import get_uxm_ratio_from_pat_by_cpg_idx
from importlib import resources

def get_marker(markerfile="Atlas.U25.l4.hg38.tsv"):
    marker_path = resources.path('pattools.deconv.loyfer', markerfile)
    marker = pd.read_csv(str(marker_path), sep='\t')
    marker = marker[~marker['name'].isna()]
    return marker

def ge_tissue_matrix_and_methylation_density(pat_file, markerfile, genome_version, cpg_bed):
    marker = get_marker(markerfile)
    TISSUE = list(Counter(marker['target']).keys())
    TISSUE = sorted(TISSUE,reverse=False)
    marker = marker[~np.any(marker[TISSUE].isna(),axis=1)]

    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker['name'].to_list())

    methy = get_uxm_ratio_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = ['name', 'cpg', 'uxm_ratio','depth']
    data = pd.merge(marker, methy, on='name', how='left')
    data = data[~data['uxm_ratio'].isna()]
    return data.loc[:, TISSUE]*data['depth'].values[:,np.newaxis], data.loc[:, 'uxm_ratio']*data['depth']


def deconvolution_loyfer(pat_file, out_file, markerfile, genome_version, cpg_bed, optimization='nnls'):
    tissue_matrix, methy_density = ge_tissue_matrix_and_methylation_density(pat_file, markerfile, genome_version, cpg_bed)
    if optimization == 'QP':
        opt_func = opt_qp
    else:
        opt_func = opt_nnls
    proportion = opt_func(tissue_matrix.values, methy_density.values)
    out = pd.DataFrame({'tissue': tissue_matrix.columns, 'proportion': proportion})
    out.to_csv(out_file, sep='\t', index=False)

if __name__ == '__main__':
    pass
