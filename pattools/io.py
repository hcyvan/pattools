import sys
import pysam


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
