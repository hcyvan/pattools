import sys
import pysam
import gzip
from typing import List
from collections import deque
from pattools.utils import is_gzip_file


class Open:
    def __init__(self, file_name, is_zip=None):
        if is_zip is None:
            if is_gzip_file(file_name):
                self.__is_zip = True
            else:
                self.__is_zip = False
        else:
            self.__is_zip = is_zip
        if self.__is_zip:
            self._f = gzip.open(file_name, 'rt')
        else:
            self._f = open(file_name, 'rt')

    def readline(self):
        line = self._f.readline()
        if line:
            return line.strip('\n')
        return None

    def __iter__(self):
        return self

    def __next__(self):
        line = self._f.readline()
        if line:
            if self.__is_zip:
                return line
            else:
                return line
        else:
            raise StopIteration

    def close(self):
        self._f.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class Output:
    """
    This class is used to control file output
    eg:
    with Output(filename=out.__str__(), file_format='text', bgzip=True) as of:
        of.write("hello Output\n")
    """

    def __init__(self, filename=None, file_format='cgs', bgzip=False):
        """
        :param filename: output filename. if filename is None, it will write to standard output.
        :param file_format: pat, pat format. motif, motif format
        :param bgzip: Whether to use bgzip compression for output files. False by default.
        """
        self.filename = filename
        self.file_format = file_format
        self.bgzip = bgzip
        if self.filename is None:
            self.bgzip = False
        else:
            if bgzip and not self.filename.endswith('.gz'):
                self.filename += '.gz'

        if self.filename is None:
            self.of = sys.stdout
        else:
            if self.bgzip:
                self.of = pysam.BGZFile(self.filename, 'wb')
            else:
                self.of = open(self.filename, 'w')

    def write(self, line: str):
        if self.bgzip:
            self.of.write(line.encode())
        else:
            self.of.write(line)

    def close(self):
        if self.filename is not None:
            self.of.close()
            if self.bgzip:
                if self.file_format in ['mvc', 'mv', 'pat', 'cgs']:
                    pysam.tabix_index(self.filename, csi=True, seq_col=0, start_col=1, end_col=1, force=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class _CpG2TabixRandom:
    def __init__(self, filename: str, region: str | List[str] = None, window=1):
        """
        TODO: The conversion between CpG and genome indices necessitates the specification of distinct index files
         or coordinate files in conjunction with index files.
        """
        self._filename = filename
        self._filename_csi = filename + '.csi'
        self._regions = region
        self._regions_pointer = -1
        self._iterator = None
        self._tabixfile = pysam.TabixFile(self._filename, index=self._filename_csi)
        self._is_multi_region = not (self._regions is None or isinstance(self._regions, str))
        if self._is_multi_region:
            self._regions_pointer = 0
        else:
            self._iterator = self._tabixfile.fetch(region=self._regions)

        self._queue = deque()
        for i in range(window):
            self._queue.append(self._next())

    def __iter__(self):
        return self

    def _next(self):
        if self._is_multi_region:
            if self._regions_pointer < len(self._regions):
                if self._iterator is None:
                    _region = self._regions[self._regions_pointer]
                    self._iterator = self._tabixfile.fetch(region=_region)
                try:
                    return next(self._iterator)
                except StopIteration:
                    self._regions_pointer += 1
                    self._iterator = None
                    return self._next()
            else:
                return None
        else:
            try:
                return next(self._iterator)
            except StopIteration:
                return None

    def __next__(self):
        latest_line = self._queue[-1]
        if latest_line is None:
            raise StopIteration()
        item_latest = latest_line.split('\t')
        line = self._queue.popleft()
        items = line.split('\t')
        line2 = self._next()
        self._queue.append(line2)
        return items[0], int(items[1]), int(item_latest[1]), int(items[2])

    def close(self):
        if self._tabixfile:
            self._tabixfile.close()


class _TabixRandom:
    def __init__(self, filename: str, region: str | List[str] = None):
        """
        TODO: The conversion between CpG and genome indices necessitates the specification of distinct index files
         or coordinate files in conjunction with index files.
        """
        self._filename = filename
        self._filename_csi = filename + '.csi'
        self._regions = region
        self._regions_pointer = -1
        self._iterator = None
        self._tabixfile = pysam.TabixFile(self._filename, index=self._filename_csi)
        self._is_multi_region = not (self._regions is None or isinstance(self._regions, str))
        if self._is_multi_region:
            self._regions_pointer = 0
        else:
            self._iterator = self._tabixfile.fetch(region=self._regions)

    def __iter__(self):
        return self

    def __next__(self):
        if self._is_multi_region:
            if self._regions_pointer < len(self._regions):
                if self._iterator is None:
                    _region = self._regions[self._regions_pointer]
                    self._iterator = self._tabixfile.fetch(region=_region)
                try:
                    return next(self._iterator)
                except StopIteration:
                    self._regions_pointer += 1
                    self._iterator = None
                    return self.__next__()
            else:
                raise StopIteration
        else:
            return next(self._iterator)

    def close(self):
        if self._tabixfile:
            self._tabixfile.close()


class _TabixSequential:
    def __init__(self, filename: str, region: str | List[str] = None):
        self._filename = filename
        self._regions_pointer = -1
        self._gz = gzip.open(self._filename, 'rt')

        if region is not None:
            if isinstance(region, str):
                region = [region]
            self._regions = []
            for reg in region:
                chr, idx = reg.split(':')
                start, end = idx.split('-')
                start = int(start)
                end = int(end)
                self._regions.append([chr, start, end])
                self._regions_pointer = 0
        else:
            self._regions = None

    def __iter__(self):
        return self

    def __next__(self):
        if self._regions is None:
            return next(self._gz)
        else:
            line = self._gz.readline()
            while True:
                _, s, e = self._regions[self._regions_pointer]
                if line.startswith('#'):
                    line = self._gz.readline()
                    continue
                line = line.strip()
                items = line.split('\t')
                idx_motif = int(items[1])
                if idx_motif < s:
                    line = self._gz.readline()
                    continue
                elif idx_motif > e:
                    self._regions_pointer += 1
                    if self._regions_pointer == len(self._regions):
                        raise StopIteration()
                else:
                    return line

    def close(self):
        if self._gz:
            self._gz.close()


class CpG2Tabix:
    """
    [chr]\t[genome_idx]\t[cpg_idx]
    CpG1 index: tabix -C -b 2 -e 2 -s 1 xxxx.pat.gz/mv.gz/mvc.gz
    CpG2 index: tabix -C -b 3 -e 3 -s 1 xxxx.pat.gz/mv.gz/mvc.gz
    """

    def __init__(self, filename: str, region: str | List[str] = None, window=1):
        self._tabix = _CpG2TabixRandom(filename, region, window=window)

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._tabix)

    def readline(self):
        try:
            return self.__next__()
        except StopIteration:
            return None

    def close(self):
        self._tabix.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class PatTabix:
    """
    [chr]\t[cpg_idx]\t[others]
    tabix -C -b 2 -e 2 -s 1 xxxx.pat.gz/mv.gz/mvc.gz
    """

    def __init__(self, filename: str, region: str | List[str] = None):
        RANDOM_READ_LIMIT = 8000
        # TODO: Here, we conducted a preliminary comparison. In the WSL environment on Windows, sequential
        #  read performance surpasses random read when requests exceed 8000 intervals, though more precise
        #  thresholds need to be tested across different platforms and environments.
        if region is not None and isinstance(region, list) and len(region) > RANDOM_READ_LIMIT:
            self._tabix = _TabixSequential(filename, region)
        else:
            self._tabix = _TabixRandom(filename, region)

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self._tabix)
        row = line.strip().split('\t')
        chrom = row[0]
        cpg_idx = int(row[1])
        # TODO: use a more general format
        motif_count_arr = [int(x) for x in row[2:]]
        return chrom, cpg_idx, motif_count_arr

    def readline(self):
        try:
            return self.__next__()
        except StopIteration:
            return None

    def close(self):
        self._tabix.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def readline_and_parse(self, motifs):
        items = self.readline()
        if items:
            chrom, cpg_idx, motif_count_arr = items
            motif_count = dict(zip(motifs, motif_count_arr))
            return chrom, cpg_idx, motif_count
        return None
