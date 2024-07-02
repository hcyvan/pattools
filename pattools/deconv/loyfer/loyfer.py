import pandas as pd
import numpy as np
from collections import Counter
from pattools.region import GenomicRegion
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.utils import get_uxm_ratio_from_pat_by_cpg_idx
from importlib import resources
from pattools.deconv.base import Markers
from typing import List, Optional


class LoyferMarkers(Markers):
    def get_cell_types(self):
        return ['Adipocytes', 'Bladder-Ep', 'Blood-B',
                'Blood-Granul', 'Blood-Mono+Macro', 'Blood-NK', 'Blood-T',
                'Bone-Osteob', 'Breast-Basal-Ep', 'Breast-Luminal-Ep', 'Colon-Ep',
                'Colon-Fibro', 'Dermal-Fibro', 'Endothel', 'Epid-Kerat', 'Eryth-prog',
                'Fallopian-Ep', 'Gallbladder', 'Gastric-Ep', 'Head-Neck-Ep',
                'Heart-Cardio', 'Heart-Fibro', 'Kidney-Ep', 'Liver-Hep',
                'Lung-Ep-Alveo', 'Lung-Ep-Bron', 'Neuron', 'Oligodend', 'Ovary-Ep',
                'Pancreas-Acinar', 'Pancreas-Alpha', 'Pancreas-Beta', 'Pancreas-Delta',
                'Pancreas-Duct', 'Megakaryocytes', 'Prostate-Ep', 'Skeletal-Musc',
                'Small-Int-Ep', 'Smooth-Musc', 'Thyroid-Ep']

    def __init__(self, genome_version='hg38', panel='U25', marker_file=None, include: Optional[List[str]] = None,
                 exclude: Optional[List[str]] = None):
        if marker_file is None:
            marker_file = f'Atlas.{panel}.l4.{genome_version}.tsv'
        super(LoyferMarkers, self).__init__(resources.path('pattools.deconv.loyfer', marker_file), genome_version,
                                            include, exclude)


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


def deconvolution_loyfer(pat_file, out_file, marker_file, genome_version, panel, cpg_bed, optimization='nnls',
                         exclude=None,
                         include=None):
    marker_obj = LoyferMarkers(genome_version, panel, marker_file, include=include, exclude=exclude)
    marker = marker_obj.get_markers()
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
