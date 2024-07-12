from pattools.motif import Motif
from scipy.spatial.distance import pdist, squareform
import numpy as np
from scipy.stats import rankdata
from collections import OrderedDict


class MRESC:
    """
    Multiple repeats with equal spacing clustering
    """

    def __init__(self, n=4):
        self.motif = Motif(n)
        x = self.motif.vectors
        distances = pdist(x, metric='euclidean')
        distance_matrix = squareform(distances)
        flattened = distance_matrix.flatten()
        ranks = rankdata(flattened, method='dense')
        rank_order_matrix = ranks.reshape(distance_matrix.shape)
        self.m_dist = rank_order_matrix - 1

    def fit(self, vectors):
        motifs = self.motif.vectors2motifs(vectors)
        motif_count = self.motif.count_motifs(motifs)
        counts = list(motif_count.values())
        clusters = [0] * len(counts)
        for i in range(len(motif_count)):
            for j in self.find_neighbors(i):
                max_nb = self.find_neighbor_max(j, counts)
                if counts[i] >= max_nb:
                    self.set_cluster(i, j, clusters)
        motif_label = OrderedDict(zip(motif_count.keys(), clusters))
        labels = []
        for m in motifs:
            labels.append(motif_label[m])
        return np.array(labels)

    def find_neighbors(self, idx):
        neighbors = []
        for i in range(len(self.m_dist[0])):
            if self.m_dist[idx, i] == 1:
                neighbors.append(i)
        return neighbors

    def find_neighbor_max(self, idx, counts):
        neighbors = self.find_neighbors(idx)
        return max([counts[x] for x in neighbors])

    def set_cluster(self, i, j, cluster):
        if cluster[i] == 0:
            cluster[i] = max(set(cluster)) + 1
        cluster[j] = cluster[i]

