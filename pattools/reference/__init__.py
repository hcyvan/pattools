from pattools.cmd import command, Cmd
from .genome import genome_build


@command('genome-build', 'This command is used to build reference genome')
class RegionCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-f', '--fasta', required=True, help='The reference genome to build')
        parser.add_argument('-c', '--chromosomes', help='The reference genome to build')
        parser.add_argument('-o', '--output', required=True, help='The output directory of the reference genome')

    def do(self, args):
        fasta = args.fasta
        output = args.output
        chromosomes = args.chromosomes
        genome_build(fasta, output, chromosomes)
