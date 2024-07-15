import pandas as pd
from pattools.region.region import GenomicRegion
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.utils import get_methylation_density_from_pat_by_cpg_idx
from pattools.deconv.base import Markers
from importlib import resources
from typing import List, Optional, Dict


class SunMarkers(Markers):
    def get_cell_type_info(self) -> Dict[str, str]:
        return {
            "Liver": "Liver",
            "Lungs": "Lungs",
            "Colon": "Colon",
            "SmallIntestines": "Small intestines",
            "Pancreas": "Pancreas",
            "AdrenalGlands": "Adrenal glands",
            "Esophagus": "Esophagus",
            "AdiposeTissues": "Adipose tissues",
            "Heart": "Heart",
            "Brain": "Brain",
            "T-cells": "T cells",
            "B-cells": "B cells",
            "Neutrophils": "Neutrophils",
            "Placenta": "Placenta",
        }

    def __init__(self, genome_version='hg38', include: Optional[List[str]] = None, exclude: Optional[List[str]] = None):
        super(SunMarkers, self).__init__(resources.path('pattools.deconv.sun', 'markerType1.csv'), genome_version,
                                         include, exclude)

    def do_fix_markers(self, marker, final_types):
        marker2_path = resources.path('pattools.deconv.sun', 'markerType2.csv')
        marker2 = pd.read_csv(str(marker2_path), sep='\t')
        select_headers = [self.genome_version] + final_types
        marker = marker.loc[:, select_headers]
        marker2 = marker2.loc[:, select_headers]
        marker = pd.concat([marker, marker2])
        marker = marker[~marker[self.genome_version].isna()]
        return marker


def ge_tissue_matrix_and_methylation_density(pat_file, genome_version, cpg_bed, include=None, exclude=None):
    marker_obj = SunMarkers(genome_version, include=include, exclude=exclude)
    marker = marker_obj.get_markers()
    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker[genome_version].to_list())
    methy = get_methylation_density_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = [genome_version, 'cpg', 'methylation_density']
    methy['methylation_density'] = methy['methylation_density'] * 100
    data = pd.merge(marker, methy, on=genome_version, how='left')
    data = data[~data['methylation_density'].isna()]
    return data.loc[:, marker_obj.get_final_cell_types()], data.loc[:, 'methylation_density']


def deconvolution_sun(pat_file, out_file, genome_version, cpg_bed, optimization='nnls', exclude=None, include=None):
    tissue_matrix, methy_density = ge_tissue_matrix_and_methylation_density(pat_file, genome_version, cpg_bed,
                                                                            include=include, exclude=exclude)
    if optimization == 'QP':
        opt_func = opt_qp
    else:
        opt_func = opt_nnls
    proportion = opt_func(tissue_matrix.values, methy_density.values)
    out = pd.DataFrame({'tissue': tissue_matrix.columns, 'proportion': proportion})
    out.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    pass
