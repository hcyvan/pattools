from pattools.cmd import command, Cmd
from .methscan import extract_cov

@command('pat2x', 'This command is used to convert the PAT format into other formats for further analysis with different software tools.')
class Pat2XCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-g', '--genome', required=True, help='the pre build genome directory')
        parser.add_argument('-i', '--input', required=True, help='Input file, *.pat.gz format')
        parser.add_argument('-x', '--x-format',
                            required=True,
                            choices=['cov-methscan'],
                            help='The target format')
        parser.add_argument('-o', '--out', help='The output file, *.gz format. There are four columns in total,'
                                                'representing chromosome, index, methylation ratio, and total sequencing depth of loci.'
                                                'If not set, output is sent to stdout.')

    def do(self, args):
        if args.x_format == 'cov-methscan':
            extract_cov(args.input, args.genome, args.out)
