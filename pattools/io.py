import sys
import pysam
from typing import List


class Output:
    """
    This class is used to control file output
    eg:
    with Output(filename=out.__str__(), file_format='text', bgzip=True) as of:
        of.write("hello Output\n")
    """

    def __init__(self, filename=None, file_format='pat', bgzip=False):
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
            if not self.filename.endswith('.gz'):
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
                if self.file_format == 'motif':
                    pysam.tabix_index(self.filename, csi=True, seq_col=0, start_col=1, end_col=1, force=True)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class Tabix:
    def __init__(self, filename: str, region: str | List[str] = None):
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
                    return self._parse_line(next(self._iterator))
                except StopIteration:
                    self._regions_pointer += 1
                    self._iterator = None
                    return self.__next__()
            else:
                raise StopIteration
        else:
            return self._parse_line(next(self._iterator))

    def _parse_line(self, line: str):
        return line

    def readline(self):
        try:
            return self.__next__()
        except StopIteration:
            return None

    def close(self):
        if self._tabixfile:
            self._tabixfile.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


class CpGTabix(Tabix):
    def __init__(self, filename, region=None):
        super(CpGTabix, self).__init__(filename, region)

    def _parse_line(self, line):
        row = line.split('\t')
        return row[0], int(row[1]), int(row[2])


class MotifTabix(Tabix):
    def __init__(self, filename, region=None):
        super(MotifTabix, self).__init__(filename, region)

    def _parse_line(self, line: str):
        row = line.strip().split('\t')
        chrom = row[0]
        cpg_idx = int(row[1])
        motif_count_arr = [int(x) for x in row[2:]]
        return chrom, cpg_idx, motif_count_arr

    def readline_and_parse(self, motifs):
        items = self.readline()
        if items:
            chrom, cpg_idx, motif_count_arr = items
            motif_count = dict(zip(motifs, motif_count_arr))
            return chrom, cpg_idx, motif_count
        return None
