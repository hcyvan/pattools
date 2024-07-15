from .vector import extract_vector, extract_vector_from_multi_motif_file
from .diff import vector_diff
from .support import extract_motif_from_region

from pattools.cmd import command, Cmd


@command('vector-region', 'extract vectors')
class VectorRegionCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-r', '--region', required=True, help='The motif to extract')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')

    def do(self, args):
        extract_motif_from_region(args.input, args.region, args.out)


__all__ = ['extract_vector', 'extract_vector_from_multi_motif_file', 'vector_diff']
