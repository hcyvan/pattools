import gzip
from collections import Counter


class PatIterator:
    def __init__(self, pat_file):
        self.pat = gzip.open(pat_file, 'rb')
        self.line = self.read_next_line()

    def read_next_line(self):
        _line = self.pat.readline()
        if len(_line):
            return _line.decode().strip().split('\t')
        return None

    def __iter__(self):
        return self

    def __next__(self):
        raise StopIteration


class PatWindow(PatIterator):
    """
    This class is used to split Pat file into windows
    eg:
        patFile = '/path/to/test.chr21_22.pat.gz'
        patWindow = PatWindow(patFile)
        for win in patWindow:
            print(win)
    """

    def __init__(self, pat_file, window=4, remove_empty=True):
        super(PatWindow, self).__init__(pat_file)
        self.window = window
        self.remove_empty = remove_empty

        self.patternDictChr = None
        self.patternDictStart = None
        self.patternDict = dict()

        self.patternTmpDictChr = None
        self.patternTmpDictStart = None
        self.patternTmpDict = dict()

    def read_next_pattern_group(self):
        """
        :return (bool): whether there are more pattern group. If yes, you should use *extract_pattern* to extract
         the pattern group
        """
        self.patternDict = dict()
        self.patternDict.update(self.patternTmpDict)
        if self.line is None and not self.patternDict:
            return False
        if self.line is None and self.patternDict:
            self.patternDictChr = self.patternTmpDictChr
            self.patternDictStart = self.patternTmpDictStart
            return True
        elif self.line and self.patternDict:
            if self.patternTmpDictChr != self.line[0] or self.patternTmpDictStart != int(self.line[1]):
                self.patternDictChr = self.patternTmpDictChr
                self.patternDictStart = self.patternTmpDictStart
                return True
        self.update_pattern_dict()
        self.patternDictChr = self.line[0]
        self.patternDictStart = int(self.line[1])
        self.line = self.read_next_line()
        if self.line is None:
            return True
        while self.patternDictChr == self.line[0] and self.patternDictStart == int(self.line[1]):
            self.update_pattern_dict()
            self.line = self.read_next_line()
            if self.line is None:
                break
        return True

    def update_pattern_dict(self):
        self.patternDict[self.line[2]] = self.patternDict.get(self.line[2], 0) + int(self.line[3])

    def extract_pattern(self):
        self.patternTmpDict = dict()
        self.patternTmpDictChr = self.patternDictChr
        self.patternTmpDictStart = self.patternDictStart + 1
        ret = dict()

        for k, v in self.patternDict.items():
            if len(k) >= self.window:
                ret[k[0:self.window]] = ret.get(k[0:self.window], 0) + v
                if len(k) >= self.window + 1:
                    self.patternTmpDict[k[1:]] = self.patternTmpDict.get(k[1:], 0) + v
        return self.patternDictChr, self.patternDictStart, ret

    def __next__(self):
        while True:
            more = self.read_next_pattern_group()
            if more:
                chr, start, pattern = self.extract_pattern()
                if self.remove_empty:
                    if pattern:
                        return chr, start, pattern
                else:
                    return chr, start, pattern
            else:
                raise StopIteration


class PatRegion(PatIterator):
    """
    This class is used to extract the methylation motifs of the target CpG index regions from the pat file.
    eg:
        patFile = '/path/to/test.chr21_22.pat.gz'
        patRegion = PatRegion(patFile)
        for win in patRegion:
            print(win)
    """

    def __init__(self, pat_file, cpg_idx_regions):
        """
        :param pat_file: The pat file to extract the methylation motifs
        :param cpg_idx_regions: A list of CpG index regions, such as ['chr1:10-100', 'chr1:10-200', ...]
        """
        super(PatRegion, self).__init__(pat_file)
        self.cpg_idx_map = None
        self.cpg_idx_map_tp = None
        self.cpg_idx_map_len = None
        self.cache = []

        self._init_cpg_idx_map(cpg_idx_regions)

    def _init_cpg_idx_map(self, cpg_idx_regions):
        self.cpg_idx_map = []
        for i, cpg_idx in enumerate(cpg_idx_regions):
            if cpg_idx:
                chr = cpg_idx.split(':')[0]
                pos = cpg_idx.split(':')[1]
                self.cpg_idx_map.append([chr, int(pos.split('-')[0]), int(pos.split('-')[1]), Counter()])
        self.cpg_idx_map = sorted(self.cpg_idx_map, key=lambda x: (x[1], x[2]))
        self.cpg_idx_map_tp = 0
        self.cpg_idx_map_len = len(self.cpg_idx_map)

    def __next__(self):
        if len(self.cache) > 0:
            ret = self.cache.pop(0)
            return ret
        while True:
            if self.line is None:
                if len(self.cache) > 0:
                    break
                else:
                    raise StopIteration
            start = int(self.line[1])
            end = int(self.line[1]) + len(self.line[2]) - 1
            motif = self.line[2]
            for j in range(self.cpg_idx_map_tp, self.cpg_idx_map_len):
                chr0 = self.cpg_idx_map[j][0]
                start0 = self.cpg_idx_map[j][1]
                end0 = self.cpg_idx_map[j][2]
                if end < start0:
                    break
                elif start > end0:
                    for k in range(self.cpg_idx_map_tp, j + 1):
                        self.cache.append([f'{chr0}:{start0}-{end0}', self.cpg_idx_map[k][3]])
                    self.cpg_idx_map_tp = j + 1
                    ret = self.cache.pop(0)
                    self.line = self.read_next_line()
                    return ret
                else:
                    s_idx_delta = start0 - start if start0 - start > 0 else 0
                    e_idx_delta = end - end0 if end - end0 > 0 else 0
                    motif_delta = motif[s_idx_delta: (len(motif) - e_idx_delta)]
                    self.cpg_idx_map[j][3] += Counter([motif_delta])
            self.line = self.read_next_line()


if __name__ == '__main__':
    pass
