import math
from typing import Tuple, Dict, Any, Union, List
from numpy.typing import NDArray
import numpy as np
from collections import Counter, OrderedDict
from sklearn.cluster import DBSCAN
import hdbscan
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt

from pattools.motif import Motif

Array2D = np.typing.NDArray[Tuple[int, int]]
Array1D = np.typing.NDArray[Tuple[int]]

np.random.seed(1000)


class VectorCalculator(object):
    def __init__(self, window: int = 4, cluster='HDBSCAN'):
        """
        :param cluster: HDBSCAN, DBSCAN
        """
        self._window = None
        self._base_vector = None
        self._motif = None
        self._scale = None

        self._chr = None
        self._start = None
        self._motif_count = dict()
        self._group = np.array([])
        self._sample = np.array([])
        self._vectors = np.array([])

        self._labels = np.array([])
        # _clusters_labels_count: OrderedDict of cluster_label => label_counts. The dict is ordered
        # from the most label count to the least.
        self._clusters_labels_count: OrderedDict[int, int] = OrderedDict()
        # _clusters_centroid: OrderedDict of cluster_label => cluster_centroid. The dict is ordered
        # from the most label count to the least
        self._clusters_centroid: OrderedDict[int, Array1D] = OrderedDict()
        self._clusters_labels_group: OrderedDict[int, dict] = OrderedDict()
        self._distance_matrix: Array2D = np.array([[0.0]])

        self.set_param_window(window)
        self._cluster = cluster

    def set_param_window(self, window: int):
        self._window = window
        self._base_vector = np.zeros(self._window)
        self._motif: Motif = Motif(self._window)
        self._scale = np.sqrt(np.sum((np.ones(window)) ** 2))

    def set_motif_count(self, chrom: str, start: int, motif_count: Dict[str, int], sample=Union[Any, List],
                        group=Union[Any, List]):
        """
        :param chrom:
        :param start:
        :param motif_count:
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
        self._vectors = self._motif.motif_count2vectors(self._motif_count)
        if isinstance(group, (np.ndarray, list, tuple)):
            self._group = np.append(self._group, group)
        else:
            self._group = np.array([group] * len(self._vectors))
        if isinstance(sample, (np.ndarray, list, tuple)):
            self._sample = np.append(self._sample, sample)
        else:
            self._sample = np.array([sample] * len(self._vectors))

        self._labels = np.array([])
        self._clusters_labels_count = OrderedDict()
        self._clusters_centroid = OrderedDict()
        return self

    def calc(self):
        if len(self._vectors):
            self._do_cluster()
            self._count_labels()
            if self.get_clusters_number():
                self._find_centroids()
                m = np.vstack([[self._base_vector], list(self._clusters_centroid.values())])
                self._distance_matrix = self._multi_dist(m, self._scale)
        return self

    def calc_labels_groups_samples(self):
        """
        labels:     0   0   0   0   0   0   1   1   1   1       => the marker for cluster result
        groups:     1   0   0   0   0   0   0   0   1   1       => the marker for input group
        samples:    5   1   1   2   3   4   6   7   8   9       => the marker for each samples

        The labels show there are two clusters: 0 and 1.
        In cluster 0:
        + target_group [group 0]: the group to which most labels belong in cluster 0
        + cluster_groups_label_count [6]: the count of labels contained in cluster 0
        + inter_groups_label_count [5]: the count of labels contained in the intersection of cluster 0 and group 0
        + group_groups_label_count [7]: the count of labels contained in group 0
        + cluster_groups_count [2,(group0,1)]: the count of groups contained in cluster 0
        + inter_groups_count [1,(always 1)]: the count of groups contained in the intersection of cluster 0 and group 0
        + group_groups_count [1,(always 1)]: the count of groups contained in group 0
        + cluster_samples_count [5,(sample1,2,3,4,5)]: the count of samples contained in cluster 0
        + inter_samples_count [4,(sample1,2,3,4)]: the count of samples contained in the intersection of cluster 0 and group 0
        + group_samples_count [6,(sample1,2,3,4,6,7)]: the count of samples contained in group 0
        """
        self._group = np.array(self._group)
        self._labels = np.array(self._labels)
        self._sample = np.array(self._sample)
        if len(self._labels) and len(self._group) == len(self._labels):
            for label in self._get_cluster_labels():
                cluster_groups = self._group[self._labels == label]
                cluster_samples = self._sample[self._labels == label]
                target_group, target_group_label_count = Counter(cluster_groups).most_common(1)[0]
                group_groups = self._group[self._group == target_group]
                group_samples = self._sample[self._group == target_group]

                self._clusters_labels_group[label] = dict(
                    target_group=target_group,
                    cluster_groups_label_count=len(cluster_groups),
                    inter_groups_label_count=target_group_label_count,
                    group_groups_label_count=len(group_groups),
                    cluster_groups_count=len(set(cluster_groups)),
                    inter_groups_count=len(set(cluster_groups) & set(group_groups)),
                    group_groups_count=len(set(group_groups)),
                    cluster_samples_count=len(set(cluster_samples)),
                    inter_samples_count=len(set(cluster_samples) & set(group_samples)),
                    group_samples_count=len(set(group_samples))
                )
        return self

    def get_labels_groups_samples_str(self):
        outs = []
        for k, v_dict in self._clusters_labels_group.items():
            v = [
                v_dict['target_group'],
                v_dict['cluster_groups_label_count'],
                v_dict['inter_groups_label_count'],
                v_dict['group_groups_label_count'],
                v_dict['cluster_groups_count'],
                v_dict['inter_groups_count'],
                v_dict['group_groups_count'],
                v_dict['cluster_samples_count'],
                v_dict['inter_samples_count'],
                v_dict['group_samples_count'],
            ]
            outs.append(f'{k}:{v[0]}:[{v[1]}:{v[2]}:{v[3]}][{v[4]}:{v[5]}:{v[6]}][{v[7]}:{v[8]}:{v[9]}]')
        return '|'.join(outs)

    def __add__(self, other):
        if isinstance(other, VectorCalculator):
            vc = VectorCalculator(self._window, self._cluster)
            counter2 = Counter(self._motif_count) + Counter(other._motif_count)
            _group = np.append(self._group, other._group)
            _sample = np.append(self._sample, other._sample)
            vc.set_motif_count(self._chr, self._start, dict(counter2), sample=_sample, group=_group)
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

    def get_vectors(self):
        return self._vectors

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
        sample_idx = self._get_sample_idx()
        if sample_idx is not None:
            _vectors = self._vectors[sample_idx, :]
            _group = self._group[sample_idx]
        else:
            _vectors = self._vectors
            _group = self._group

        if self._cluster == 'HDBSCAN':
            _labels = self._do_cluster_hdbscan(_vectors, _group)
        else:
            _labels = self._do_cluster_dbscan(_vectors)

        if sample_idx is not None:
            _label_map = dict()
            for v, l in zip(_vectors, _labels):
                _label_map[''.join([str(x) for x in v])] = l
            _labels = []
            for v in self._vectors:
                label = _label_map.get(''.join([str(x) for x in v]), -1)
                _labels.append(label)
            labels = np.array(_labels)
        else:
            labels = _labels
        self._labels = labels

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

    def _get_sample_idx(self):
        MAX_VECTOR_SIZE = 5000
        if len(self._vectors) > MAX_VECTOR_SIZE:
            select_idx = np.random.choice(self._vectors.shape[0], MAX_VECTOR_SIZE, replace=False)
            return select_idx
        return None

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


class VectorPlot:
    def __init__(self, vector_calculator: VectorCalculator):
        self._vector_calculator: VectorCalculator = vector_calculator

    def plot_vector_cluster(self):
        vectors = self._vector_calculator.get_vectors()
        labels = self._vector_calculator.get_labels()
        colorsMap = {-1: 'gray', 0: '#1f77b4', 1: '#ff7f0e', 2: '#2ca02c', 3: "#d62728", 4: "#9467bd", 5: "#8c564b",
                     6: "#e377c2", 7: "#7f7f7f", 8: "#bcbd22", 9: "#17becf"}
        labels_color = np.vectorize(lambda x: colorsMap[x])(labels)
        print(vectors)
        tsne = TSNE(n_components=2, perplexity=30)
        X_2d = tsne.fit_transform(vectors)

        plt.figure(figsize=(4.2, 4))
        plt.scatter(X_2d[:, 0], X_2d[:, 1], c=labels_color, s=15)
        plt.title('Cluster Plot after t-SNE')
        plt.xlabel('t-SNE Dimension 1')
        plt.ylabel('t-SNE Dimension 2')
        plt.show()
