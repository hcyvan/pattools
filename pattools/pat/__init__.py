from pattools.pat.pat import PatWindow, PatRegion
from pattools.pat.generate import compress_and_tabix_pat
from pattools.cmd import command, Cmd

__all__ = ['PatWindow', 'PatRegion']


@command('pat-generate', 'Generate pat format file')
class VectorExtractCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-o', '--out', default=None, help='The output file')

    def do(self, args):
        pass


@command('pat-compress', 'Compress and index pat file')
class VectorExtractCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-o', '--out', default=None, help='The output file')

    def do(self, args):
        compress_and_tabix_pat(args.input, args.out)
