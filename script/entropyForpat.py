import gzip
import math
import argparse
import os

class PatWindow:
    """
    This class is used to split Pat file into windows
    eg:
        patFile = 'D:/pat/test.chr21_22.pat.gz'
        patWindow = PatWindow(patFile)
        for win in patWindow:
            print(win)
    """
    def __init__(self, pat_file, window=4, remove_empty=True):
        self.window = window
        self.remove_empty = remove_empty
        self.pat = gzip.open(pat_file, 'rb')
        self.line = self.read_next_line()

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

    def read_next_line(self):
        _line = self.pat.readline()
        if len(_line):
            return _line.decode().strip().split('\t')
        return None

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
                    self.patternTmpDict[k[1:]] = v
        return self.patternDictChr, self.patternDictStart, ret

    def __iter__(self):
        return self

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

if __name__ == '__main__':
    pass

def calculate_entropy(patWin, deep):
    """
    @parm
    patWin = ({'CCCC': 1, 'CCCT': 1, 'CCTC': 2, 'CT.C': 1, 'TCCC': 1, 'TCCT': 1, 'TCTC': 1, 'TTCT': 2, 'TTTC': 1, 'CTCT': 1})

    @return 
    entropy value
    """
    substring_counts = {'CCCC': 0, 'CCCT': 0, 'CCTC': 0, 'CCTT': 0, 'CTCC': 0, 'CTCT': 0, 'CTTC': 0, 'CTTT': 0,
                        'TCCC': 0, 'TCCT': 0, 'TCTC': 0, 'TCTT': 0, 'TTCC': 0, 'TTCT': 0, 'TTTC': 0, 'TTTT': 0}

    total_count = 0
    for substring, count in patWin.items():
        if substring in substring_counts:
            substring_counts[substring] += count
            total_count += count   
    if total_count >= deep:
        entropies = {key: (value / total_count) * math.log2((value / total_count)) for key, value in substring_counts.items() if value != 0}
        entropy = -sum(entropies.values()) * 0.25
        entropy = round(entropy, 4)
    else:
        entropy = -1
    return entropy


parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', help='Input file path')
parser.add_argument('-o', dest='output_dir', help='Output dir')
parser.add_argument('-d', dest='deep', help='total counts for pattern')

args = parser.parse_args()

patWindow = PatWindow(args.input)

base_name = os.path.splitext(os.path.basename(args.input))[0]
output_file = f"{args.output_dir}/{base_name}.entropy.bed"

with open(args.output, 'w') as f:
    for pat in patWindow:
        entropy = calculate_entropy(pat[2],int(args.deep))
        f.write(f"{pat[0]}\t{pat[1]}\t{entropy}\n")