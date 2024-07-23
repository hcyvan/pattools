import math
from typing import Tuple, Dict, Any, Union, List

from numpy.typing import NDArray
import numpy as np
from collections import Counter, OrderedDict
from sklearn.cluster import DBSCAN
import hdbscan

from scipy.spatial.distance import pdist, squareform

from pattools.motif import Motif
from pattools.vector.clustering.cluster import MRESC

Array2D = np.typing.NDArray[Tuple[int, int]]
Array1D = np.typing.NDArray[Tuple[int]]

np.random.seed(1000)


class VectorCalculator(object):
    cache_motif = dict()

    def __init__(self, window: int = 4, cluster='HDBSCAN'):
        """
        :param cluster: HDBSCAN, DBSCAN, MRESC
        """
        self._window = None
        self._chr = None
        self._start = None
        self._motif_count = dict()
        self._group: List = []
        self._sample: List = []
        self._vectors: List = []
        self._labels: List = []
        # _clusters_labels_count: OrderedDict of cluster_label => label_counts. The dict is ordered
        # from the most label count to the least.
        self._clusters_labels_count: OrderedDict[int, int] = OrderedDict()
        # _clusters_centroid: OrderedDict of cluster_label => cluster_centroid. The dict is ordered
        # from the most label count to the least
        self._clusters_centroid: OrderedDict[int, Array1D] = OrderedDict()
        self._clusters_labels_group: OrderedDict[int, dict] = OrderedDict()
        self._distance_matrix: Array2D = np.array([[0.0]])
        self._window = window
        self._cluster = cluster

    def get_motif(self):
        motif = self.cache_motif.get(self._window, None)
        if motif is not None:
            return motif
        motif = Motif(self._window)
        self.cache_motif[self._window] = motif
        return motif

    def set_motif_count(self, chrom: str, start: int,
                        motif_count: Dict[str, int],
                        vectors: Union[Any, List] = None,
                        sample: Union[Any, List] = None,
                        group: Union[Any, List] = None):
        """
        :param chrom:
        :param start:
        :param motif_count: the number of each type of methylation vectors.
        :param vectors: the methylation vectors. If the vectors are not from the same sample, this parameter must be set
        :param sample: the sample of vectors. If sample is a list, it specifies the sample assignment for each
                    corresponding vector. Conversely, if the input is a scalar, it denotes that all vectors
                    are assigned to a single sample.
        :param group: the group of vectors. If group is a list, it specifies the group assignment for each
                    corresponding vector. Conversely, if the input is a scalar, it denotes that all vectors
                    are assigned to a single group.
        """
        self._chr = chrom
        self._start = start
        self._motif_count = motif_count
        total = sum(self._motif_count.values())
        if vectors is not None:
            self._vectors = vectors
        else:
            self._vectors = self.get_motif().motif_count2vectors(self._motif_count)
        if group is not None:
            if isinstance(group, list):
                self._group = group
            else:
                self._group = [group] * total
        if sample is not None:
            if isinstance(sample, list):
                self._sample = sample
            else:
                self._sample = [sample] * total
        self._labels = []
        self._clusters_labels_count = OrderedDict()
        self._clusters_centroid = OrderedDict()
        return self

    def cluster(self):
        if len(self._vectors):
            _base_vector = np.zeros(self._window)
            _scale = np.sqrt(np.sum((np.ones(self._window)) ** 2))
            self._do_cluster()
            self._count_labels()
            if self.get_clusters_number():
                self._find_centroids()
                m = np.vstack([[_base_vector], list(self._clusters_centroid.values())])
                self._distance_matrix = self._multi_dist(m, _scale)
        return self

    def __add__(self, other):
        if isinstance(other, VectorCalculator):
            _sample = self._sample + other._sample
            _group = self._group + other._group
            _vectors = self._vectors + other._vectors
            counter2 = Counter(self._motif_count) + Counter(other._motif_count)
            vc = VectorCalculator(self._window, self._cluster)
            vc.set_motif_count(self._chr, self._start, motif_count=dict(counter2), vectors=_vectors, sample=_sample,
                               group=_group)
            return vc
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, VectorCalculator):
            return self._chr == other._chr and self._start == other._start
        return NotImplemented

    def __gt__(self, other):
        if isinstance(other, VectorCalculator):
            return self._chr == other._chr and self._start > other._start
        return NotImplemented

    def __lt__(self, other):
        if isinstance(other, VectorCalculator):
            return self._chr == other._chr and self._start < other._start
        return NotImplemented

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

    def get_vectors(self, label=None):
        if label is None:
            return self._vectors
        else:
            return list(np.array(self._vectors)[np.array(self._labels) == label])

    def get_labels(self):
        return self._labels

    def _get_cluster_labels(self):
        return [x[0] for x in self._clusters_labels_count.items()]

    def distance_between_top_m_and_n(self, m=0, n=1):
        """
        0-based index
        m: The label of cluster with the m-th most vectors
        n: The label of cluster with the n-th most vectors
        """
        if self.get_clusters_number():
            return self._distance_matrix[m + 1, n + 1]
        return None

    def get_mvc(self, group_order, genome_idx):
        """
        MVC will output information about _labels, _group and _sample, such as:
        _labels:    0   0   0   0   0   0   1   1   1   1       => the marker for cluster result
        _group:     1   0   0   0   0   0   0   0   1   1       => the marker for input group
        _sample:    5   1   1   2   3   4   6   7   8   9       => the marker for each sample
        The labels show there are two clusters: 0 and 1.

        @param group_order:
        @param genome_idx:
        @return:
        """

        self._motif_count.values()
        motif_count_ordered: OrderedDict[str, int] = self.get_motif().count_motifs(self._motif_count)
        motif_tag = '|'.join([str(x) for x in motif_count_ordered.values()])
        centers = []
        for k, v in self._clusters_centroid.items():
            centers.append(f"({','.join([f'{x:.3f}' for x in v])})")
        centers = '|'.join(centers)

        _group = np.array(self._group)
        _labels = np.array(self._labels)
        _sample = np.array(self._sample)

        cluster_group_label_count = []
        for k, v in self._clusters_labels_count.items():
            idx = _labels == k
            counter = Counter(_group[idx])
            group_count = []
            for g in group_order:
                group_count.append(str(counter[g]))
            cluster_group_label_count.append(','.join(group_count))
        cluster_group_label_count = '|'.join(cluster_group_label_count)

        cluster_group_sample_count = []
        for k, v in self._clusters_labels_count.items():
            idx = _labels == k
            sample_count = []
            for g in group_order:
                idx2 = _group == g
                samples = sorted(set(_sample[idx & idx2]))
                sample_count.append(str(len(samples)))
            cluster_group_sample_count.append(','.join(sample_count))
        cluster_group_sample_count = '|'.join(cluster_group_sample_count)

        cluster_group_samples = []
        for k, v in self._clusters_labels_count.items():
            idx = _labels == k
            samples_ = []
            for g in group_order:
                idx2 = _group == g
                samples = sorted(set(_sample[idx & idx2]))
                samples_.append(':'.join([str(x) for x in samples]))
            cluster_group_samples.append(','.join(samples_))
        cluster_group_samples = '|'.join(cluster_group_samples)
        return f"{self._chr}\t{self._start}\t{genome_idx - 1}\t{genome_idx + self._window}\t{motif_tag}\t{self.get_clusters_number()}\t{centers}\t{cluster_group_label_count}\t{cluster_group_sample_count}\t{cluster_group_samples}"

    def _do_cluster(self):
        _vectors = np.array(self._vectors)
        _group = np.array(self._group)
        sample_idx = None
        if self._cluster != 'MRESC':
            sample_idx = self._get_sample_idx()
            if sample_idx is not None:
                _vectors = _vectors[sample_idx, :]
                _group = _group[sample_idx]

        if self._cluster == 'HDBSCAN':
            _labels = self._do_cluster_hdbscan(_vectors, _group)
        elif self._cluster == 'DBSCAN':
            _labels = self._do_cluster_dbscan(_vectors)
        elif self._cluster == 'MRESC':
            _labels = self._do_cluster_mresc(_vectors)
        else:
            raise Exception(f"Unknown cluster method {self._cluster}")
        if sample_idx is not None:
            _label_map = dict()
            for v, l in zip(_vectors, _labels):
                _label_map[''.join([str(x) for x in v])] = l
            _labels = []
            for v in _vectors:
                label = _label_map.get(''.join([str(x) for x in v]), -1)
                _labels.append(label)
            labels = np.array(_labels)
        else:
            labels = _labels
        self._labels = list(labels)

    @staticmethod
    def _do_cluster_hdbscan(vectors, vector_group):
        if len(vectors) < 2:
            return np.array([0])
        min_group = Counter(vector_group).most_common()[-1]
        min_cluster_size = min_group[1] // 3
        if min_cluster_size < 2:
            min_cluster_size = 2
        _labels = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size).fit_predict(vectors)
        return _labels

    @staticmethod
    def _do_cluster_dbscan(vectors):
        _min_samples = math.ceil(len(vectors) * 0.33)
        _dbscan = DBSCAN(eps=1, min_samples=_min_samples).fit(vectors)
        return _dbscan.labels_

    def _do_cluster_mresc(self, vectors):
        """
        TODO: Unify `_do_cluster_xxxx` into @staticmethod or class method
        """
        return MRESC(self._window).fit(vectors)

    def _get_sample_idx(self):
        MAX_VECTOR_SIZE = 5000
        if len(self._vectors) > MAX_VECTOR_SIZE:
            select_idx = np.random.choice(len(self._vectors), MAX_VECTOR_SIZE, replace=False)
            return select_idx
        return None

    def _count_labels(self):
        _labels = np.array(self._labels)
        _labels = _labels[_labels != -1]
        label_counts = Counter(_labels)
        self._clusters_labels_count = OrderedDict(label_counts.most_common())

    def _find_centroids(self):
        _labels = np.array(self._labels)
        _vectors = np.array(self._vectors)
        for label, _ in self._clusters_labels_count.items():
            cluster_vectors = _vectors[_labels == label]
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
