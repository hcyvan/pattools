import pandas as pd
import gzip
import heapq
from pattools.region import GenomicRegion
from pattools.pat import PatRegion
from pattools.deconv.utils import get_methylation_density_from_pat_by_cpg_idx
from pattools.deconv.optimization import opt_qp, opt_nnls
import pysam
from collections import Counter
from collections import OrderedDict
from importlib import resources

SOURCE = ['Monocytes', 'BCells', 'Cd4tCells', 'NkCells',
          'Cd8tCells', 'Neutrophils', 'ErythrocyteProgenitors', 'Adipocytes',
          'CorticalNeurons', 'Hepatocytes', 'LungCells', 'PancreaticBetaCells',
          'PancreaticAcinarCells', 'PancreaticDuctCells',
          'VascularEndothelialCells', 'ColonEpithelialCells', 'LeftAtrium',
          'Bladder', 'Breast', 'HeadAndNeckLarynx', 'Kidney', 'Prostate',
          'Thyroid', 'UpperGI', 'UterusCervix']


def get_marker(genome_version='hg38'):
    marker_path = resources.path('pattools.deconv.moss', 'marker.csv')
    marker = pd.read_csv(str(marker_path), sep='\t')
    marker = marker[~marker[genome_version].isna()]
    return marker


def ge_source_matrix_and_methylation_density(pat_file, genome_version, cpg_bed):
    marker = get_marker(genome_version)
    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker[genome_version].to_list())
    methy = get_methylation_density_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = [genome_version, 'cpg', 'methylation_density']
    data = pd.merge(marker, methy, on=genome_version, how='left')
    data = data[~data['methylation_density'].isna()]
    return data.loc[:, SOURCE], data.loc[:, 'methylation_density']


def deconvolution_moss(pat_file, out_file, genome_version, cpg_bed, optimization='nnls'):
    source_matrix, methy_density = ge_source_matrix_and_methylation_density(pat_file, genome_version, cpg_bed)
    if optimization == 'QP':
        opt_func = opt_qp
    else:
        opt_func = opt_nnls
    proportion = opt_func(source_matrix.values, methy_density.values)
    out = pd.DataFrame({'source': source_matrix.columns, 'proportion': proportion})
    out.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    pass
