import math
from typing import Tuple, Dict
from numpy.typing import NDArray
import numpy as np
from collections import Counter, OrderedDict
from sklearn.cluster import DBSCAN
from scipy.spatial.distance import pdist, squareform

from pattools.motif import Motif

Array2D = NDArray[Tuple[int, int]]
Array1D = NDArray[Tuple[int]]


class VectorCalculator(object):
    def __init__(self, window: int = 4, min_vector_proportion_in_cluster=0.33):
        self._window = None
        self._base_vector = None
        self._motif = None
        self._scale = None
        self._min_vector_proportion_in_cluster = None

        self._chr = None
        self._start = None
        self._motif_count = None

        self._vectors = np.array([])
        self._labels = np.array([])
        # _clusters_labels_count: OrderedDict of cluster_label => label_counts. The dict is ordered
        # from the most label count to the least
        self._clusters_labels_count: OrderedDict[int, int] = OrderedDict()
        # _clusters_centroid: OrderedDict of cluster_label => cluster_centroid. The dict is ordered
        # from the most label count to the least
        self._clusters_centroid: OrderedDict[int, Array1D] = OrderedDict()
        self._distance_matrix: Array2D = np.array([[0.0]])

        self.set_param_window(window)
        self.set_param_min_vector_proportion_in_cluster(min_vector_proportion_in_cluster)

    def set_param_min_vector_proportion_in_cluster(self, min_vector_proportion_in_cluster):
        self._min_vector_proportion_in_cluster = min_vector_proportion_in_cluster

    def set_param_window(self, window: int):
        self._window = window
        self._base_vector = np.zeros(self._window)
        self._motif: Motif = Motif(self._window)
        self._scale = np.sqrt(np.sum((np.ones(window)) ** 2))

    def set_motif_count(self, chrom: str, start: int, motif_count: Dict[str, int]):
        self._chr = chrom
        self._start = start
        self._motif_count = motif_count

        self._vectors = np.array([])
        self._labels = np.array([])
        self._clusters_labels_count = OrderedDict()
        self._clusters_centroid = OrderedDict()
        return self

    def calc(self):
        self._vectors = self._motif.motif_count2vectors(self._motif_count)
        if len(self._vectors):
            self._do_cluster()
            self._count_labels()
            if self.get_clusters_number():
                self._find_centroids()
                m = np.vstack([[self._base_vector], list(self._clusters_centroid.values())])
                self._distance_matrix = self._multi_dist(m, self._scale)

    def get_clusters_number(self):
        return len(self._clusters_labels_count)

    def distance_between_top_k_and_base(self, k=0):
        """
        0-based index
        k: The label of cluster with the k-th most vectors
        """
        if self.get_clusters_number():
            return self._distance_matrix[0, k + 1]
        return None

    def distance_between_top_m_and_n(self, m=0, n=1):
        """
        0-based index
        m: The label of cluster with the m-th most vectors
        n: The label of cluster with the n-th most vectors
        """
        if self.get_clusters_number():
            return self._distance_matrix[m + 1, n + 1]
        return None

    def __str__(self):
        self._motif_count.values()
        motif_count_ordered: OrderedDict[str, int] = self._motif.count_motifs(self._motif_count)
        motif_tag = '|'.join([str(x) for x in motif_count_ordered.values()])
        if self.get_clusters_number() >= 2:
            dist_top_0_1_tag = self.distance_between_top_m_and_n(0, 1)
        else:
            dist_top_0_1_tag = -1
        dist_top_0_base_tag = self.distance_between_top_k_and_base(0)
        if dist_top_0_base_tag is None:
            dist_top_0_base_tag = -1
        return f"{self._chr}\t{self._start}\t{self.get_clusters_number()}\t{dist_top_0_1_tag:.3f}\t{dist_top_0_base_tag:.3f}\t{len(self._vectors)}\t{motif_tag}"

    def _do_cluster(self):
        _min_samples = math.ceil(len(self._vectors) * self._min_vector_proportion_in_cluster)
        _dbscan = DBSCAN(eps=1, min_samples=_min_samples).fit(self._vectors)
        self._labels = np.array(_dbscan.labels_)

    def _count_labels(self):
        _labels = self._labels[self._labels != -1]
        label_counts = Counter(_labels)
        self._clusters_labels_count = OrderedDict(label_counts.most_common())

    def _find_centroids(self):
        for label, _ in self._clusters_labels_count.items():
            cluster_vectors = self._vectors[self._labels == label]
            self._clusters_centroid[label] = np.mean(cluster_vectors, axis=0)

    @staticmethod
    def _distance(vector0: Array1D, vector1: Array1D, scale: float) -> float:
        dist = np.sqrt(np.sum((vector0 - vector1) ** 2)) / scale
        return dist

    @staticmethod
    def _multi_dist(x: Array2D, scale: float = 1.0) -> Array2D:
        distances = pdist(x, metric='euclidean')
        distance_matrix = squareform(distances)
        distance_matrix = distance_matrix / scale
        return distance_matrix
