import pandas as pd
from pattools.region import GenomicRegion
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.utils import get_methylation_density_from_pat_by_cpg_idx
from importlib import resources

TISSUE = ['Liver', 'Lungs', 'Colon', 'SmallIntestines', 'Pancreas', 'AdrenalGlands', 'Esophagus', 'AdiposeTissues',
          'Heart', 'Brain', 'T-cells', 'B-cells', 'Neutrophils', 'Placenta']


def get_marker(genome_version='hg38'):
    marker1_path = resources.path('pattools.deconv.sun', 'markerType1.csv')
    marker2_path = resources.path('pattools.deconv.sun', 'markerType2.csv')

    marker1 = pd.read_csv(str(marker1_path), sep='\t')
    marker2 = pd.read_csv(str(marker2_path), sep='\t')
    marker_headers = ['hg19', 'hg38'] + TISSUE
    marker1 = marker1.loc[:, marker_headers]
    marker2 = marker2.loc[:, marker_headers]
    marker = pd.concat([marker1, marker2])
    marker = marker[~marker[genome_version].isna()]
    return marker


def ge_tissue_matrix_and_methylation_density(pat_file, genome_version, cpg_bed):
    marker = get_marker(genome_version)
    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker[genome_version].to_list())
    methy = get_methylation_density_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = [genome_version, 'cpg', 'methylation_density']
    methy['methylation_density'] = methy['methylation_density']*100
    data = pd.merge(marker, methy, on=genome_version, how='left')
    data = data[~data['methylation_density'].isna()]
    return data.loc[:, TISSUE], data.loc[:, 'methylation_density']


def deconvolution_sun(pat_file, out_file, genome_version, cpg_bed, optimization='nnls'):
    tissue_matrix, methy_density = ge_tissue_matrix_and_methylation_density(pat_file, genome_version, cpg_bed)
    if optimization == 'QP':
        opt_func = opt_qp
    else:
        opt_func = opt_nnls
    proportion = opt_func(tissue_matrix.values, methy_density.values)
    out = pd.DataFrame({'tissue': tissue_matrix.columns, 'proportion': proportion})
    out.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    pass
