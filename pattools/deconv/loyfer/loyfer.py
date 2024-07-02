import pandas as pd
import numpy as np
from collections import Counter
from pattools.region import GenomicRegion
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.utils import get_uxm_ratio_from_pat_by_cpg_idx
from importlib import resources


def validate_ref_tissues(marker, tissue_list):
    for col in tissue_list:
        if col not in marker.columns:
            print('Invalid cell type (not in markerfile):', col)
            exit()


def get_marker(genome_version='hg38', panel='U25', marker_file=None, ignore=None, include=None):
    """
    :param genome_version: hg19 or hg38 is support
    :param panel: U25 or U250 is support
    """
    if marker_file is None:
        if genome_version in ['hg38', 'hg19'] and panel in ['U25', 'U250']:
            marker_file = f'Atlas.{panel}.l4.{genome_version}.tsv'
        else:
            raise Exception(f"Not Supported Genome Version: {genome_version} and Panel: {panel}")

    marker_path = resources.path('pattools.deconv.loyfer', marker_file)
    marker = pd.read_csv(str(marker_path), sep='\t')
    marker = marker[~marker[genome_version].isna()]
    if ignore is not None:
        validate_ref_tissues(marker, ignore)
        for col in ignore:
            del marker[col]
            marker = marker[marker.target != col]
    elif include is not None:
        validate_ref_tissues(marker, include)
        marker = marker[marker.target.isin(include)]
        keep = list(marker.columns[:8]) + include
        marker = marker[keep]
    return marker


def ge_tissue_matrix_and_methylation_density(pat_file, marker, genome_version, cpg_bed):
    TISSUE = list(Counter(marker['target']).keys())
    TISSUE = sorted(TISSUE, reverse=False)
    marker = marker[~np.any(marker[TISSUE].isna(), axis=1)]

    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker['name'].to_list())

    methy = get_uxm_ratio_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = ['name', 'cpg', 'uxm_ratio', 'depth']
    data = pd.merge(marker, methy, on='name', how='left')
    data = data[~data['uxm_ratio'].isna()]
    return data.loc[:, TISSUE] * data['depth'].values[:, np.newaxis], data.loc[:, 'uxm_ratio'] * data['depth']


def deconvolution_loyfer(pat_file, out_file, markerfile, genome_version,panel, cpg_bed, optimization='nnls', ignore=None,
                         include=None):
    marker = get_marker(genome_version, panel,markerfile, ignore, include)
    tissue_matrix, methy_density = ge_tissue_matrix_and_methylation_density(pat_file, marker, genome_version, cpg_bed)
    if optimization == 'QP':
        opt_func = opt_qp
    else:
        opt_func = opt_nnls
    proportion = opt_func(tissue_matrix.values, methy_density.values)
    out = pd.DataFrame({'tissue': tissue_matrix.columns, 'proportion': proportion})
    out.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    pass
