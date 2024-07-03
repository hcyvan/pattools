import pandas as pd
from pattools.region import GenomicRegion
from pattools.deconv.utils import get_methylation_density_from_pat_by_cpg_idx
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.base import Markers
from importlib import resources
from typing import List, Optional, Dict


class MossMarkers(Markers):
    def get_cell_type_info(self) -> Dict[str, str]:
        return {
            "Monocytes": "Monocytes",
            "BCells": "B cells",
            "Cd4tCells": "T helper(CD4+) cells",
            "NkCells": "NK cells",
            "Cd8tCells": "T cytotoxic (CD8+) cells",
            "Neutrophils": "Neutrophils",
            "ErythrocyteProgenitors": "Erythrocyte progenitors",
            "Adipocytes": "Adipocytes",
            "CorticalNeurons": "Cortical neurons",
            "Hepatocytes": "Hepatocytes",
            "LungCells": "Lung cells",
            "PancreaticBetaCells": "Pancreatic beta cells",
            "PancreaticAcinarCells": "Pancreatic acinar cells",
            "PancreaticDuctCells": "Pancreatic duct cells",
            "VascularEndothelialCells": "Vascular endothelial cells",
            "ColonEpithelialCells": "Colon epithelial cells",
            "LeftAtrium": "Left atrium",
            "Bladder": "Bladder",
            "Breast": "Breast",
            "HeadAndNeckLarynx": "Head and neck larynx",
            "Kidney": "Kidney",
            "Prostate": "Prostate",
            "Thyroid": "Thyroid",
            "UpperGI": "Upper Gastrointestinal",
            "UterusCervix": "Uterus cervix",
        }

    def __init__(self, genome_version='hg38', include: Optional[List[str]] = None, exclude: Optional[List[str]] = None):
        super(MossMarkers, self).__init__(resources.path('pattools.deconv.moss', 'marker.csv'), genome_version, include,
                                          exclude)


def ge_source_matrix_and_methylation_density(pat_file, genome_version, cpg_bed, include=None, exclude=None):
    marker_obj = MossMarkers(genome_version, include=include, exclude=exclude)
    marker = marker_obj.get_markers()
    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker[genome_version].to_list())
    methy = get_methylation_density_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = [genome_version, 'cpg', 'methylation_density']
    data = pd.merge(marker, methy, on=genome_version, how='left')
    data = data[~data['methylation_density'].isna()]
    return data.loc[:, marker_obj.get_final_cell_types()], data.loc[:, 'methylation_density']


def deconvolution_moss(pat_file, out_file, genome_version, cpg_bed, optimization='nnls', exclude=None, include=None):
    source_matrix, methy_density = ge_source_matrix_and_methylation_density(pat_file, genome_version, cpg_bed,
                                                                            include=include, exclude=exclude)
    if optimization == 'QP':
        opt_func = opt_qp
    else:
        opt_func = opt_nnls
    proportion = opt_func(source_matrix.values, methy_density.values)
    out = pd.DataFrame({'tissue': source_matrix.columns, 'proportion': proportion})
    out.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    pass
