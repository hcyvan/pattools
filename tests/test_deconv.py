import unittest
from pattools.deconv.sun.sun import SunMarkers
from pattools.deconv.moss.moss import MossMarkers
from pattools.deconv.loyfer.loyfer import LoyferMarkers


class TestDeconv(unittest.TestCase):

    def test_SunMarkers(self):
        marker_obj = SunMarkers(genome_version='hg38')
        markers = marker_obj.get_markers()
        select_type = markers.iloc[:, 1:].columns.tolist()
        self.assertEqual(select_type, marker_obj.get_cell_types(), f"SunMarkers cell type is right")

    def test_MossMarkers(self):
        marker_obj = MossMarkers(genome_version='hg38')
        markers = marker_obj.get_markers()
        select_type = markers.iloc[:, 1:].columns.tolist()
        self.assertEqual(select_type, marker_obj.get_cell_types(), f"SunMarkers cell type is right")

    def test_LoyferMarkers(self):
        marker_obj = LoyferMarkers(genome_version='hg38', panel='U25', include=['Smooth-Musc', 'Thyroid-Ep'])
        markers = marker_obj.get_markers()
        select_type = markers.iloc[:, 1:].columns.tolist()
        self.assertEqual(select_type, marker_obj.get_final_cell_types(), f"SunMarkers cell type is right")


if __name__ == '__main__':
    unittest.main()
