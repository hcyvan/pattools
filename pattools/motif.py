import itertools
from collections import OrderedDict


class Motif:
    """
    This class is used to generate all possible methylation motif patterns. When generating motifs
    with a CpG count of 3, we get:
          ['CCC', 'TCC', 'CTC', 'CCT', 'TTC', 'TCT', 'CTT', 'TTT']
    """

    def __init__(self, count: int = 4):
        self.count = count
        self.motifs = self._get_motif_array()

    def motif2vector(self):
        motif_vector = OrderedDict()
        for _, m in enumerate(self.motifs):
            motif_vector[m] = list(map(lambda x: dict(C=1, T=0)[x], m))
        return motif_vector

    def motif2order(self):
        motif_order = OrderedDict()
        for i, m in enumerate(self.motifs):
            motif_order[m] = i
        return motif_order

    def _get_motif_array(self):
        motifs = ["C" * self.count]
        for i in range(1, self.count + 1):
            for e in itertools.combinations(range(self.count), i):
                motif = ""
                for j in range(self.count):
                    if j in set(e):
                        motif += 'T'
                    else:
                        motif += 'C'
                motifs.append(motif)
        return motifs
