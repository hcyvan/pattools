from pattools.cmd import command, Cmd
from .beta import extract_beta


@command('beta', 'This command performs methylation ratio analysis on the sample')
class BetaCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file, *.pat.gz format')
        parser.add_argument('-d', '--depth', type=int, default=1,
                            help='the minimum total count required to calculate methylation ratio')
        parser.add_argument('-o', '--out', help='The output file, *.gz format. There are four columns in total,'
                                 'representing chromosome, index, methylation ratio, and total sequencing depth of loci.'
                                 'If not set, output is sent to stdout.')

    def do(self, args):
        extract_beta(args.input, args.depth, args.out)
