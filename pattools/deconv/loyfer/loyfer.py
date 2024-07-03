import pandas as pd
import numpy as np
from collections import Counter
from pattools.region import GenomicRegion
from pattools.deconv.optimization import opt_qp, opt_nnls
from pattools.deconv.utils import get_uxm_ratio_from_pat_by_cpg_idx
from importlib import resources
from pattools.deconv.base import Markers
from typing import List, Optional, Dict


class LoyferMarkers(Markers):
    def get_cell_type_info(self) -> Dict[str, str]:
        return {
            "Adipocytes": "Adipocytes",
            "Bladder-Ep": "Bladder-Ep",
            "Blood-B": "Blood-B",
            "Blood-Granul": "Blood-Granul",
            "Blood-Mono+Macro": "Blood-Mono+Macro",
            "Blood-NK": "Blood-NK",
            "Blood-T": "Blood-T",
            "Bone-Osteob": "Bone-Osteob",
            "Breast-Basal-Ep": "Breast-Basal-Ep",
            "Breast-Luminal-Ep": "Breast-Luminal-Ep",
            "Colon-Ep": "Colon-Ep",
            "Colon-Fibro": "Colon-Fibro",
            "Dermal-Fibro": "Dermal-Fibro",
            "Endothel": "Endothel",
            "Epid-Kerat": "Epid-Kerat",
            "Eryth-prog": "Eryth-prog",
            "Fallopian-Ep": "Fallopian-Ep",
            "Gallbladder": "Gallbladder",
            "Gastric-Ep": "Gastric-Ep",
            "Head-Neck-Ep": "Head-Neck-Ep",
            "Heart-Cardio": "Heart-Cardio",
            "Heart-Fibro": "Heart-Fibro",
            "Kidney-Ep": "Kidney-Ep",
            "Liver-Hep": "Liver-Hep",
            "Lung-Ep-Alveo": "Lung-Ep-Alveo",
            "Lung-Ep-Bron": "Lung-Ep-Bron",
            "Neuron": "Neuron",
            "Oligodend": "Oligodend",
            "Ovary-Ep": "Ovary-Ep",
            "Pancreas-Acinar": "Pancreas-Acinar",
            "Pancreas-Alpha": "Pancreas-Alpha",
            "Pancreas-Beta": "Pancreas-Beta",
            "Pancreas-Delta": "Pancreas-Delta",
            "Pancreas-Duct": "Pancreas-Duct",
            "Megakaryocytes": "Megakaryocytes",
            "Prostate-Ep": "Prostate-Ep",
            "Skeletal-Musc": "Skeletal-Musc",
            "Small-Int-Ep": "Small-Int-Ep",
            "Smooth-Musc": "Smooth-Musc",
            "Thyroid-Ep": "Thyroid-Ep",
        }

    def __init__(self, genome_version='hg38', panel='U25', marker_file=None, include: Optional[List[str]] = None,
                 exclude: Optional[List[str]] = None):
        if marker_file is None:
            marker_file = f'Atlas.{panel}.l4.{genome_version}.tsv'
        super(LoyferMarkers, self).__init__(resources.path('pattools.deconv.loyfer', marker_file), genome_version,
                                            include, exclude)


def ge_tissue_matrix_and_methylation_density(pat_file, marker, genome_version, cpg_bed):
    TISSUE = marker.columns[1:].tolist()
    marker = marker[~np.any(marker[TISSUE].isna(), axis=1)]
    gr = GenomicRegion(cpg_bed)
    genome_cpg_idx = gr.genomic_to_cpg_idx(marker[genome_version].to_list())
    methy = get_uxm_ratio_from_pat_by_cpg_idx(pat_file, genome_cpg_idx)
    methy.columns = [genome_version, 'cpg', 'uxm_ratio', 'depth']
    data = pd.merge(marker, methy, on=genome_version, how='left')
    data = data[~data['uxm_ratio'].isna()]
    return data.loc[:, TISSUE] * data['depth'].values[:, np.newaxis], data.loc[:, 'uxm_ratio'] * data['depth']


def deconvolution_loyfer(pat_file, out_file, marker_file, genome_version, panel, cpg_bed, optimization='nnls',
                         exclude=None, include=None):
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
