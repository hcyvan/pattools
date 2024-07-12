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
            "Bladder-Ep": "Bladder epithelium cells",
            "Blood-B": "Blood B cells",
            "Blood-Granul": "Blood granulocytes",
            "Blood-Mono+Macro": "Blood monocytes and macrophages",
            "Blood-NK": "Blood NK cells",
            "Blood-T": "Blood T cells",
            "Bone-Osteob": "Bone Osteoblasts",
            "Breast-Basal-Ep": "Breast basal epithelium cells",
            "Breast-Luminal-Ep": "Breast luminal epithelium cells",
            "Colon-Ep": "Colon epithelium cells",
            "Colon-Fibro": "Colon fibroblasts",
            "Dermal-Fibro": "Dermal fibroblasts",
            "Endothel": "endothelial cells",
            "Epid-Kerat": "Epidermal keratinocytes",
            "Eryth-prog": "Erythrocyte progenitors",
            "Fallopian-Ep": "Fallopian epithelium cells",
            "Gallbladder": "Gallbladder",
            "Gastric-Ep": "Gastric epithelium cells",
            "Head-Neck-Ep": "Head and neck epithelium cells",
            "Heart-Cardio": "Heart cardiomyocytes",
            "Heart-Fibro": "Heart fibroblasts",
            "Kidney-Ep": "Kidney epithelium cells",
            "Liver-Hep": "Liver hepatocytes",
            "Lung-Ep-Alveo": "lung alveolar epithelium cells",
            "Lung-Ep-Bron": "lung bronchial epithelium cells",
            "Neuron": "Neuron",
            "Oligodend": "Oligodend",
            "Ovary-Ep": "Ovary epithelium cells",
            "Pancreas-Acinar": "Pancreas acinar cells",
            "Pancreas-Alpha": "Pancreas alpha cells",
            "Pancreas-Beta": "Pancreas beta cells",
            "Pancreas-Delta": "Pancreas delta cells",
            "Pancreas-Duct": "Pancreas duct cells",
            "Megakaryocytes": "Megakaryocytes",
            "Prostate-Ep": "Prostate epithelium cells",
            "Skeletal-Musc": "Skeletal muscle cells",
            "Small-Int-Ep": "Small intestine epithelium cells",
            "Smooth-Musc": "Smooth muscle cells",
            "Thyroid-Ep": "Thyroid epithelium cells",
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
