from abc import ABC, abstractmethod
from typing import List, Dict, Optional
import pandas as pd


class Markers(ABC):
    def __init__(self, marker_file, genome_version='hg38', include: Optional[List[str]] = None,
                 exclude: Optional[List[str]] = None):
        self.marker_file = marker_file
        self.genome_version = genome_version
        self.include = include
        self.exclude = exclude

    @abstractmethod
    def get_cell_type_info(self) -> Dict[str, str]:
        """
        The function is used to set the cell type information of this deconvolution algorithm.
        :return: A dictionary of cell type information, such as:
                    cellType1 => cellTypeInfo1
                    cellType2 => cellTypeInfo2
        """
        return dict()

    def get_cell_types(self):
        return list(self.get_cell_type_info().keys())

    def get_final_cell_types(self):
        final_types = self.get_cell_types()
        if self.include is not None:
            for cell_type in self.include:
                if cell_type not in self.get_cell_types():
                    raise Exception(f'Unsupported cell type {cell_type}')
            final_types = self.include
        if self.exclude is not None:
            final_types = list(set(final_types).difference(set(self.exclude)))
        return final_types

    def _do_get_markers(self, final_types):
        marker = pd.read_csv(str(self.marker_file), sep='\t')
        marker_headers = ['target', self.genome_version] + final_types
        marker = marker.loc[marker['target'].isin(marker_headers), marker_headers]
        marker = marker.loc[:, [self.genome_version] + final_types]
        return marker

    def do_fix_markers(self, marker, final_types):
        return marker

    def get_markers(self):
        final_types = self.get_final_cell_types()
        markers = self._do_get_markers(final_types)
        markers = self.do_fix_markers(markers, final_types)
        return markers
