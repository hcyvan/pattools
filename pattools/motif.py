import itertools
from typing import Dict, List, Union
import numpy as np
from collections import OrderedDict


class Motif:
    """
    This class is used to generate all possible methylation motif patterns. When generating motifs
    with a CpG count of 3, we get:
          ['CCC', 'TCC', 'CTC', 'CCT', 'TTC', 'TCT', 'CTT', 'TTT']
    """
    cache_motif2vector = dict()
    cache_motif_array = dict()

    def __init__(self, count: int = 4):
        self.count = count
        self.motifs = self._get_motif_array()
        self.vectors = []
        self._vector_motif = dict()
        for k, v in self.motif2vector().items():
            self.vectors.append(v)
            self._vector_motif[tuple(v)] = k

    def motif2vector(self):
        m2v = self.cache_motif2vector.get(self.count)
        if m2v:
            return m2v
        m2v = OrderedDict()
        for _, m in enumerate(self.motifs):
            m2v[m] = list(map(lambda x: dict(C=1, T=0)[x], m))
        self.cache_motif2vector[self.count] = m2v
        return m2v

    def vectors2motifs(self, vectors):
        return [self._vector_motif[tuple(x)] for x in vectors]

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

    def count_motifs(self, motifs: Union[Dict[str, int], List[str]]) -> OrderedDict[str, int]:
        """
        :param motifs: a dictionary of motif => motif_count
        :return: an ordered dictionary of motif => motif_count. Only motifs including C/T are retained,
                    and these motifs are sorted from CCC to TTT,
                    such as: ['CCC', 'TCC', 'CTC', 'CCT', 'TTC', 'TCT', 'CTT', 'TTT']
        """
        order_counter = self.motif2order_counter()

        if isinstance(motifs, Dict):
            for k, v in motifs.items():
                if k in order_counter:
                    order_counter[k] += v
        if isinstance(motifs, List):
            for k in motifs:
                if k in order_counter:
                    order_counter[k] += 1
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
        return vectors

    def vectors2motif_count(self, vectors):
        motifs = self.vectors2motifs(vectors)
        return self.count_motifs(motifs)

    def _get_motif_array(self):
        motifs = self.cache_motif_array.get(self.count, None)
        if motifs is not None:
            return motifs
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
        self.cache_motif_array[self.count] = motifs
        return motifs
