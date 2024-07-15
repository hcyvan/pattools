from pathlib import Path
from .region import region_cpg2genome, region_genome2cpg, trans_region_file
from pattools.cmd import command, Cmd


@command('region', 'This command is used to convert a region between different coordinate systems, such '
                   'as genome index coordinates or CpG index coordinates.')
class RegionCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-t', '--transform', choices=['cpg2genome', 'genome2cpg'], default='cpg2genome',
                            help='Conversion direction of genomic coordinates')
        parser.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
        parser.add_argument('-i', '--input', required=True,
                            help='The input can be: 1) region string, eg: chr1:250-300, 2) a file where each line is a region string.')

    def do(self, args):
        regions = []
        file_path = Path(args.input)
        if file_path.exists():
            with file_path.open(mode='r') as f:
                for line in f:
                    regions.append(line.strip())
        else:
            regions.append(args.input)
        if args.transform == 'cpg2genome':
            region_cpg2genome(regions, args.cpg_bed)
        elif args.transform == 'genome2cpg':
            region_genome2cpg(regions, args.cpg_bed)


@command('region-file', 'This command is used to convert a region between different coordinate systems, '
                        'such as genome index coordinates or CpG index coordinates.')
class RegionFileCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-t', '--transform', choices=['cpg2genome', 'genome2cpg'], default='cpg2genome',
                            help='Conversion direction of genomic coordinates')
        parser.add_argument('--column', choices=['col2', 'col3'], default='col3',
                            help='col2: chrom\ts; col3: chrom\tstart\tend')
        parser.add_argument('--out-format', choices=['bed', 'none'], default='none',
                            help='bed: 0-base,end-exclude; none: not change')
        parser.add_argument('--offset-col2-start-and-end', type=int, default=0,
                            help='If --column is set to col2, end is calculated based on this parameter. end=start+offset')
        parser.add_argument('-c', '--cpg-bed', required=True,
                            help='The cpg_bed file of the selected genome.')
        parser.add_argument('-i', '--input', required=True, help='The input file')
        parser.add_argument('-o', '--out', help='The output file')

    def do(self, args):
        trans_region_file(args.input, out_put=args.out, cpg_bed=args.cpg_bed, transform=args.transform,
                          col=args.column, out_format=args.out_format, end_offset=args.offset_col2_start_and_end)
