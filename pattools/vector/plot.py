import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA

from pattools.vector.calculator import VectorCalculator


class VectorPlot:
    def __init__(self, vector_calculator: VectorCalculator):
        self._vector_calculator: VectorCalculator = vector_calculator

    def plot_vector_cluster(self, dim='PCA', scale=1.0, out=None, width=5, height=5):
        vectors = self._vector_calculator.get_vectors()
        labels = self._vector_calculator.get_labels()
        colors_map = self._get_colors_map_from_labels(labels)
        # labels_color = np.vectorize(lambda x: colorsMap[x])(labels)
        if dim == 'PCA':
            pca = PCA(n_components=2)
            X = pca.fit_transform(vectors)
        else:
            dim = 't-SNE'
            tsne = TSNE(n_components=2, perplexity=30)
            X = tsne.fit_transform(vectors)
        counts = self._get_counts(X, scale)
        data = pd.DataFrame({
            'x': X[:, 0],
            'y': X[:, 1],
            'count': counts,
            'cluster': labels
        })
        plt.figure(figsize=(width, height))
        scatter = sns.scatterplot(data=data, x='x', y='y', size='count', hue='cluster',
                                  palette=colors_map,
                                  sizes=(min(counts), max(counts)),
                                  legend='full')
        scatter.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.grid(False)
        plt.xlabel('Dim-1')
        plt.ylabel('Dim-2')
        if out:
            plt.savefig(out, dpi=300)
        plt.show()

    def plot_vectors_cluster3d(self, dim='PCA', scale=1.0):
        vectors = self._vector_calculator.get_vectors()
        labels = self._vector_calculator.get_labels()
        colors_map = self._get_colors_map_from_labels(labels)
        labels_color = np.vectorize(lambda x: colors_map[x])(labels)
        if dim == 'PCA':
            pca = PCA(n_components=3)
            X = pca.fit_transform(vectors)
        else:
            dim = 't-SNE'
            tsne = TSNE(n_components=3, perplexity=30)
            X = tsne.fit_transform(vectors)
        counts = self._get_counts(X, scale)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(X[:, 0], X[:, 1], X[:, 2], c=labels_color, s=counts)
        ax.set_xlabel('Dimension 1')
        ax.set_ylabel('Dimension 2')
        ax.set_zlabel('Dimension 3')
        plt.show()

    def plot_vectors_cluster3d_inter(self, dim='PCA', scale=1.0):
        vectors = self._vector_calculator.get_vectors()
        labels = self._vector_calculator.get_labels()
        colors_map = self._get_colors_map_from_labels(labels)
        labels_color = np.vectorize(lambda x: colors_map[x])(labels)
        if dim == 'PCA':
            pca = PCA(n_components=3)
            X = pca.fit_transform(vectors)
        else:
            dim = 't-SNE'
            tsne = TSNE(n_components=3, perplexity=30)
            X = tsne.fit_transform(vectors)
        counts = self._get_counts(X, scale)
        fig = go.Figure(data=[go.Scatter3d(
            x=X[:, 0], y=X[:, 1], z=X[:, 2],
            mode='markers',
            marker=dict(
                size=counts,
                color=labels_color
            )
        )])
        fig.update_layout(scene=dict(
            xaxis_title='Dimension 1',
            yaxis_title='Dimension 2',
            zaxis_title='Dimension 3'),
            width=1200,
            height=900,
            title="3D Scatter Plot")
        fig.show()

    def _get_counts(self, X, scale=1.0):
        from collections import Counter
        sizes = []
        counter = Counter()
        for loc in X:
            lb = '_'.join([str(x) for x in loc])
            counter.update([lb])
        for loc in X:
            lb = '_'.join([str(x) for x in loc])
            sizes.append(counter[lb])
        return np.log2(np.array(sizes)) * scale

    def _get_colors_map_from_labels(self, labels):
        colors_map = {-1: 'gray', 0: '#1f77b4', 1: '#ff7f0e', 2: '#2ca02c', 3: "#d62728", 4: "#9467bd", 5: "#8c564b",
                      6: "#e377c2", 7: "#7f7f7f", 8: "#bcbd22", 9: "#17becf"}
        labels_map = dict()
        for l in np.unique(labels):
            labels_map[l] = colors_map[l]
        return labels_map
