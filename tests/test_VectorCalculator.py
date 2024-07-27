import unittest
from pattools.vector.calculator import VectorCalculator
from pattools.motif import Motif


class MyTestCase(unittest.TestCase):
    def test_cluster(self):
        # chr1    6016    0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|4 2       (0.000,0.000,0.000,0.000,0.000)|(0.000,1.000,1.000,1.000,0.000)
        vc = VectorCalculator(5)
        mvs = '0|0|0|0|0|0|0|0|0|1|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|0|4'
        motif = Motif(5)
        vc.set_motif_count('chr1', 6016, motif.mvs2motif_count(mvs))
        vc.cluster()
        #print(vc.get_mvc_base())

        self.assertEqual(vc.get_clusters_number() == 2, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
