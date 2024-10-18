from pattools.pat.pat import PatWindow, PatRegion
from pattools.pat.generate import compress_and_tabix_pat, merge_pat, bam2pat
from pattools.cmd import command, Cmd

__all__ = ['PatWindow', 'PatRegion']


@command('bam2pat', 'Generate pat format file')
class Bam2PatCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-o', '--out', default=None, help='The output file')

    def do(self, args):
        bam2pat(args.input, args.out)


@command('pat-merge', 'merge pat files')
class PatMergeCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-p', '--pat', default=[], action="append")
        parser.add_argument('-o', '--out', default=None, help='The output file')

    def do(self, args):
        merge_pat(args.pat, args.out)


@command('pat-compress', 'Compress and index pat file')
class PatCompressCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-o', '--out', default=None, help='The output file')

    def do(self, args):
        compress_and_tabix_pat(args.input, args.out)
