import re
from collections import Counter
from pattools.io import Open


class MvcWindow:
    def __init__(self, groups, group_samples):
        self._groups = groups
        self._group_samples = group_samples
        self._chrom = None
        self._cpg_idx = None
        self._genome_start = None
        self._genome_end = None
        self._mvs = None
        self._mvs_num = None
        self._cluster_num = None
        self._cluster_center = None
        self._cluster_group_mvs_num = None
        self._cluster_group_samples_num = None
        self._cluster_group_samples = None
        self._cluster_group_mvs_num_counter = []
        self._cluster_group_samples_num_counter = []

        self._group_samples_counter = Counter(dict([(k, len(v)) for k, v in group_samples.items()]))
        self._mvc_str = None

    def decode(self, mvc_str):
        self._mvc_str = mvc_str
        items = mvc_str.strip("\n").split('\t')
        self._chrom = items[0]
        self._cpg_idx = items[1]
        self._genome_start = items[2]
        self._genome_end = items[3]
        self._mvs = items[4]
        self._cluster_num = int(items[5])
        self._cluster_center = items[6]
        self._cluster_group_mvs_num = items[7]
        self._cluster_group_samples_num = items[8]
        self._cluster_group_samples = items[9]
        self._cluster_group_mvs_num_counter = []
        self._cluster_group_samples_num_counter = []
        self._mvs_num = sum([int(x) for x in self._mvs.split('|')])
        if self._mvs_num > 0:
            for cluster in self._cluster_group_mvs_num.split('|'):
                counter = Counter(dict(zip(self._groups, [int(x) for x in cluster.split(',')])))
                self._cluster_group_mvs_num_counter.append(counter)
            for cluster in self._cluster_group_samples_num.split('|'):
                counter = Counter(dict(zip(self._groups, [int(x) for x in cluster.split(',')])))
                self._cluster_group_samples_num_counter.append(counter)
        return self

    def satisfied(self, target, min_mvs_in_cluster_frac=1.0, min_samples_in_group_frac=0.9, min_mvs_num=80):
        if self._cluster_num < 2 or self._mvs_num < min_mvs_num:
            return None
        _target_group_samples_total = self._group_samples_counter[target]
        for i in range(self._cluster_num):
            _mvs_counter = self._cluster_group_mvs_num_counter[i]
            _mvs_counter_total = sum(_mvs_counter.values())
            _mvs_counter_target = _mvs_counter[target]
            _samples_counter = self._cluster_group_samples_num_counter[i]
            _samples_counter_total = sum(_samples_counter.values())
            _samples_counter_target = _samples_counter[target]
            mvs_in_cluster_frac = _mvs_counter_target / _mvs_counter_total
            samples_in_group_frac = _samples_counter_target / _target_group_samples_total
            if mvs_in_cluster_frac >= min_mvs_in_cluster_frac and samples_in_group_frac >= min_samples_in_group_frac:
                return mvs_in_cluster_frac, samples_in_group_frac
        return None

    def calc_meta(self):
        if self._cluster_num < 2 or self._mvs_num < 80:
            return ""
        return "meta"


class BaseHeader:
    def __init__(self):
        self.window = None
        self.header = None
        self.headers = []

    @staticmethod
    def parse_header_window(header_str):
        p = re.compile(r"##WINDOW:\s*(\d+)")
        m = p.match(header_str)
        if m:
            return int(m.group(1))
        return None

    def decode(self, line):
        self.headers.append(line)
        if line.startswith('##'):
            if self.window is None:
                self.window = self.parse_header_window(line)
        elif line.startswith('#'):
            self.header = line[1:].split('\t')
        else:
            raise Exception('Not Header Info')


class MvcHeader(BaseHeader):
    def __init__(self):
        super().__init__()
        self.groups = None
        self.samples = None
        self.group_samples = None

    @staticmethod
    def parse_header_group(header_str):
        p = re.compile(r"##GROUP:\s*([a-zA-Z0-9,]+)")
        m = p.match(header_str)
        if m:
            return m.group(1).split(',')
        return None

    @staticmethod
    def parse_header_sample(header_str):
        p = re.compile(r"##SAMPLE:\s*([a-zA-Z0-9,]+)")
        m = p.match(header_str)
        if m:
            return m.group(1).split(',')
        return None

    @staticmethod
    def parse_header_group_sample(header_str):
        p = re.compile(r"##GROUP_SAMPLE:\s*([a-zA-Z0-9,]+)_([a-zA-Z0-9,]+)")
        m = p.match(header_str)
        if m:
            return m.group(1), m.group(2).split(',')
        return None, None

    def decode(self, line):
        super().decode(line)
        if line.startswith('#'):
            if self.groups is None:
                self.groups = self.parse_header_group(line)
            if self.samples is None:
                self.samples = self.parse_header_sample(line)
            if self.group_samples is None:
                self.group_samples = dict()
            _group, _sample = self.parse_header_group_sample(line)
            if _group:
                self.group_samples[_group] = _sample
        else:
            raise Exception('Not MVC Header Info')


class MvcFormat:
    def __init__(self):
        self.header = MvcHeader()
        self.mvw = None

    def filter(self, mvc_file, f_out, group, min_mvs_in_cluster_frac=1, min_samples_in_group_frac=0.9, with_meta=False):
        with Open(mvc_file) as f:
            for line in f:
                line = line.strip("\n")
                if line.startswith('#'):
                    self.header.decode(line)
                else:
                    if self.mvw is None:
                        self.mvw = MvcWindow(self.header.groups, self.header.group_samples)
                    self.mvw.decode(line)
                    _satisfied = self.mvw.satisfied(group, min_mvs_in_cluster_frac, min_samples_in_group_frac)
                    if _satisfied:
                        if with_meta:
                            mvs_in_cluster_frac, samples_in_group_frac = _satisfied
                            f_out.write(f"{line}\t{mvs_in_cluster_frac}\t{samples_in_group_frac}\n")
                        else:
                            f_out.write(line + '\n')


class MvHeader(BaseHeader):
    pass


class MvFormat:
    def __init__(self):
        self.header = MvHeader()

    def parse_header(self, mv_file):
        with Open(mv_file) as f:
            for line in f:
                line = line.strip("\n")
                if line.startswith('#'):
                    self.header.decode(line)
                else:
                    break
        return self
