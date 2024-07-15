from .sun.sun import deconvolution_sun, SunMarkers
from .moss.moss import deconvolution_moss, MossMarkers
from .loyfer.loyfer import deconvolution_loyfer, LoyferMarkers
from pattools.cmd import command, Cmd


@command('deconv', 'This command is used to calculate the cellular composition of the sample according '
                   'to the pat format')
class DeconvCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-p', '--pat', required=True, help='The pat file')
        parser.add_argument('-m', '--method', choices=['sun', 'moss', 'loyfer'], default='sun',
                            help='The deconvolution method.'
                                 'sun: Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing'
                                 ' for noninvasive prenatal, cancer, and transplantation assessments. '
                                 'moss: Moss et al. Comprehensive human cell-type methylation atlas reveals origins'
                                 ' of circulating cell-free DNA in health and disease. '
                                 'loyfer: Loyfer et al. A DNA methylation atlas of normal human cell types.'
                            )
        parser.add_argument('-a', '--optimization-algorithm', choices=['nnls', 'qp'], default='nnls',
                            help='The optimization algorithm for deconvolution.')
        parser.add_argument('-g', '--genome-version', choices=['hg38', 'hg19'], default='hg38',
                            help='The genome version.')
        parser.add_argument('--panel', choices=['U25', 'U250'], default='U25',
                            help='The panel of the markers. Only works when --method loyfer is used')
        parser.add_argument('-f', '--markerfile', help='markerfile for loyfer method')
        parser.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
        parser.add_argument('-o', '--out', required=True, help='The output file')
        parser.add_argument('-in', '--include', nargs='+',
                            help='Cell/Tissue types to include. Complementary to --ignore. Use `pattools '
                                 'deconv-helper` to see supported cell/tissue types.')
        parser.add_argument('-ig', '--ignore', nargs='+',
                            help='Cell/Tissue types to remove. Use `pattools deconv-helper` to see supported '
                                 'cell/tissue types.')

    def do(self, args):
        if args.method == 'sun':
            deconvolution_sun(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                              optimization=args.optimization_algorithm, include=args.include, exclude=args.ignore)
        if args.method == 'moss':
            deconvolution_moss(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                               optimization=args.optimization_algorithm, include=args.include, exclude=args.ignore)
        if args.method == 'loyfer':
            deconvolution_loyfer(args.pat, args.out, args.markerfile, args.genome_version, args.panel,
                                 cpg_bed=args.cpg_bed,
                                 optimization=args.optimization_algorithm, include=args.include, exclude=args.ignore)


@command('deconv-helper', 'print the cell types support by pattools deconv')
class DeconvHelperCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-m', '--method', choices=['sun', 'moss', 'loyfer'], default='sun',
                            help='The deconvolution method.'
                                 'sun: Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing'
                                 ' for noninvasive prenatal, cancer, and transplantation assessments. '
                                 'moss: Moss et al. Comprehensive human cell-type methylation atlas reveals origins'
                                 ' of circulating cell-free DNA in health and disease. '
                                 'loyfer: Loyfer et al. A DNA methylation atlas of normal human cell types.'
                            )

    def do(self, args):
        alg = args.method
        if alg == 'sun':
            marker = SunMarkers()
            print("The deconvolution algorithm from Sun et al")
        elif alg == 'moss':
            marker = MossMarkers()
            print("The deconvolution algorithm from Moss et al")
        elif alg == 'loyfer':
            marker = LoyferMarkers()
            print("The deconvolution algorithm from Loyfer et al")
        else:
            raise Exception('unknown deconvolution algorithm')
        print(f"Cell/Tissue\t\tDetail")
        info = marker.get_cell_type_info()
        for k, v in info.items():
            print(f"{k}\t\t{v}")
