import itertools
from typing import Dict
import numpy as np
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

    def motif2order_counter(self):
        motif_order = OrderedDict()
        for _, m in enumerate(self.motifs):
            motif_order[m] = 0
        return motif_order

    def count_motifs(self, motifs: Dict[str, int]) -> OrderedDict[str, int]:
        """
        :param motifs: a dictionary of motif => motif_count
        :return: an ordered dictionary of motif => motif_count. Only motifs including C/T are retained,
                    and these motifs are sorted from CCC to TTT,
                    such as: ['CCC', 'TCC', 'CTC', 'CCT', 'TTC', 'TCT', 'CTT', 'TTT']
        """
        order_counter = self.motif2order_counter()
        for k, v in motifs.items():
            if k in order_counter:
                order_counter[k] += v
        return order_counter

    def motif_count2vectors(self, counter: dict):
        """
        Change motif count dict to motif vectors
        eg:

        {'CC': 3, 'TC': 0, 'CT': 2, 'TT': 0}
        to
        [[1,1],[1,1],[1,1], [1,0],[1,0]]

        :param counter: motif count dict.
        :return: array of motif vectors
        """
        motif2vector_map = self.motif2vector()
        vectors = []
        for k, v in counter.items():
            vectors.extend([motif2vector_map[k]] * v)
        return np.array(vectors)

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
