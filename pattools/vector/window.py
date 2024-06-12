import re
from typing import Optional


class VectorCluster:
    _RE = re.compile(
        r"(?P<tag>\d+):(?P<group>\d+):"
        r"\[(?P<cluster_group_label>\d+):(?P<inter_group_label>\d+):(?P<group_group_label>\d+)\]"
        r"\[(?P<cluster_group>\d+):(?P<inter_group>\d+):(?P<group_group>\d+)\]"
        r"\[(?P<cluster_sample>\d+):(?P<inter_sample>\d+):(?P<group_sample>\d+)\]")

    def __init__(self, expression: str = None):
        self._expression = expression
        self.tag = None
        self.group = None
        self.cluster_group_label = None
        self.inter_group_label = None
        self.group_group_label = None
        self.cluster_group = None
        self.inter_group = None
        self.group_group = None
        self.cluster_sample = None
        self.inter_sample = None
        self.group_sample = None
        self.decode()

    def decode(self):
        if self._expression:
            m = self._RE.search(self._expression)
            if m:
                self.tag = m.group('tag')
                self.group = m.group('group')
                self.cluster_group_label = int(m.group('cluster_group_label'))
                self.inter_group_label = int(m.group('inter_group_label'))
                self.group_group_label = int(m.group('group_group_label'))
                self.cluster_group = int(m.group('cluster_group'))
                self.inter_group = int(m.group('inter_group'))
                self.group_group = int(m.group('group_group'))
                self.cluster_sample = int(m.group('cluster_sample'))
                self.inter_sample = int(m.group('inter_sample'))
                self.group_sample = int(m.group('group_sample'))

    def __str__(self):
        return self._expression

    def __repr__(self):
        return self._expression


class VectorClusterList:
    def __init__(self, expression: str = None):
        self._expression = expression
        self._clusters = []
        self.decode()

    def decode(self):
        if self._expression:
            for exp in self._expression.split('|'):
                self._clusters.append(VectorCluster(exp))

    def filter(self, group, max_cluster_group_label_vs_inter_group_label=100000, min_inter_sample_vs_group_sample=0):
        for cluster in self._clusters:
            if cluster.group == group:
                if cluster.cluster_group_label / cluster.inter_group_label <= max_cluster_group_label_vs_inter_group_label:
                    if cluster.inter_sample / cluster.group_sample >= min_inter_sample_vs_group_sample:
                        return True
        return False

    def __str__(self):
        return self._expression

    def __repr__(self):
        return self._expression


class VectorWindow:
    def __init__(self, expression: str = None):
        self._expression = expression.strip()
        self._chr = None
        self._start = None
        self.cluster_count = None
        self.vector_count = None
        self._vector_expression = None
        self._cluster_expression = None
        self.clusters: Optional[VectorClusterList, None] = None
        self.decode()

    def decode(self):
        if self._expression:
            item = self._expression.split('\t')
            self._chr = item[0]
            self._start = item[1]
            self.cluster_count = int(item[2])
            self.vector_count = int(item[5])
            self._vector_expression = item[6]
            if len(item) >= 8:
                self._cluster_expression = item[7]
                self.clusters = VectorClusterList(self._cluster_expression)

    def __str__(self):
        return self._expression

    def __repr__(self):
        return self._expression
