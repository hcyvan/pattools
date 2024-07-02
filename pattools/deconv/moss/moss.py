import pandas as pd
from pattools.region import GenomicRegion
from pattools.deconv.utils import get_methylation_density_from_pat_by_cpg_idx
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.base import Markers
from importlib import resources
from typing import List, Optional


class MossMarkers(Markers):
    def get_cell_types(self):
        return ['Monocytes', 'BCells', 'Cd4tCells', 'NkCells',
                'Cd8tCells', 'Neutrophils', 'ErythrocyteProgenitors', 'Adipocytes',
                'CorticalNeurons', 'Hepatocytes', 'LungCells', 'PancreaticBetaCells',
                'PancreaticAcinarCells', 'PancreaticDuctCells',
                'VascularEndothelialCells', 'ColonEpithelialCells', 'LeftAtrium',
                'Bladder', 'Breast', 'HeadAndNeckLarynx', 'Kidney', 'Prostate',
                'Thyroid', 'UpperGI', 'UterusCervix']

    def __init__(self, genome_version='hg38', include: Optional[List[str]] = None, exclude: Optional[List[str]] = None):
        super(MossMarkers, self).__init__(resources.path('pattools.deconv.moss', 'marker.csv'), genome_version, include,
                                         exclude)


def ge_source_matrix_and_methylation_density(pat_file, genome_version, cpg_bed):
    marker_obj = MossMarkers(genome_version)
    marker = marker_obj.get_markers()
    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker[genome_version].to_list())
    methy = get_methylation_density_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = [genome_version, 'cpg', 'methylation_density']
    data = pd.merge(marker, methy, on=genome_version, how='left')
    data = data[~data['methylation_density'].isna()]
    return data.loc[:, marker_obj.get_cell_types()], data.loc[:, 'methylation_density']


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
