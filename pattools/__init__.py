import argparse
from pattools.deconv import deconvolution_sun, deconvolution_moss
from pattools.entropy import extract_entropy
from pattools.ratio import extract_ratio

def main():
    parser = argparse.ArgumentParser(prog='pattools',
                                     description='pattools is a BS-seq analysis tool suite based on pat format')
    parser.add_argument('-v', '--version', action='version', version='0.0.1')
    parser.add_argument('-q', '--quiet', action='store_true', help='print run details to stderr')
    subparsers = parser.add_subparsers(dest='sub', required=True, title='command', description='The available commands',
                                       help='select a sub command to use')
    # =====================================================================
    parser_region = subparsers.add_parser('deconv',
                                          help='This command is used to calculate the cellular composition of the '
                                               'sample according to the pat format.')
    parser_region.add_argument('-p', '--pat', required=True, help='The pat file')
    parser_region.add_argument('-m', '--method', choices=['sun', 'moss', 'loyfer'], default='sun',
                               help='The deconvolution method.'
                                    'sun: Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing'
                                    ' for noninvasive prenatal, cancer, and transplantation assessments. '
                                    'moss: Moss et al. Comprehensive human cell-type methylation atlas reveals origins'
                                    ' of circulating cell-free DNA in health and disease. '
                                    'loyfer: Loyfer et al. A DNA methylation atlas of normal human cell types.'
                               )
    parser_region.add_argument('-a', '--optimization-algorithm', choices=['nnls', 'qp'], default='nnls',
                               help='The optimization algorithm for deconvolution.')
    parser_region.add_argument('-g', '--genome-version', choices=['hg38', 'hg19'], default='hg38',
                               help='The genome version.')
    parser_region.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
    parser_region.add_argument('-o', '--out', required=True, help='The output file')
    # =====================================================================
    parser_entropy = subparsers.add_parser('entropy',
                                           help='This command performs entropy analysis on the sample')
    parser_entropy.add_argument('-i', '--input', required=True, help='Input file, *.pat.gz format')
    parser_entropy.add_argument('-d', '--depth', required=True, help='the minimum total count required to calculate entropy')
    parser_entropy.add_argument('-o', '--out', required=True, help='The output file')
    # =====================================================================
    parser_ratio = subparsers.add_parser('ratio',
                                           help='This command performs methylation ratio analysis on the sample')
    parser_ratio.add_argument('-i', '--input', required=True, help='Input file, *.pat.gz format')
    parser_ratio.add_argument('-d', '--depth', required=True, help='the minimum total count required to calculate entropy')
    parser_ratio.add_argument('-o', '--out', required=True, help='The output file, *.gzip format')    

    args = parser.parse_args()
    if args.sub == 'deconv':
        if args.method == 'sun':
            deconvolution_sun(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                              optimization=args.optimization_algorithm)
        if args.method == 'moss':
            deconvolution_moss(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                               optimization=args.optimization_algorithm)
        else:
            print("This method is not complete")
    elif args.sub == 'entropy':
        extract_entropy(args.input, args.depth, args.out)
    elif args.sub == 'ratio':
        extract_ratio(args.input, args.depth, args.out)
        # print("This method is not complete")
