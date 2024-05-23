import argparse


def main():
    parser = argparse.ArgumentParser(prog='pattools',
                                     description='pattools is a BS-seq analysis tool suite based on pat format')
    parser.add_argument('-v', '--version', action='version', version='0.0.1')
    parser.add_argument('-q', '--quiet', action='store_true', help='print run details to stderr')
    subparsers = parser.add_subparsers(dest='sub', required=True, title='command', description='The available commands',
                                       help='select a sub command to use')

    parser_region = subparsers.add_parser('deconv',
                                          help='This command is used to calculate the cellular composition of the '
                                               'sample according to the pat format.')
    parser_region.add_argument('-p', '--pat', required=True, help='The pat file')
    parser_region.add_argument('-m', '--method', required=True,
                               help='The deconvolution method.'
                                    'sun: Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing'
                                    ' for noninvasive prenatal, cancer, and transplantation assessments. '
                                    'moss: Moss et al. Comprehensive human cell-type methylation atlas reveals origins'
                                    ' of circulating cell-free DNA in health and disease. '
                                    'loyfer: Loyfer et al. A DNA methylation atlas of normal human cell types.'
                               )
    parser_region.add_argument('-o', '--out', required=True, help='The output file')

    parser_entropy = subparsers.add_parser('entropy',
                                           help='This command performs entropy analysis on the sample')
    parser_region.add_argument('-o', '--out', required=True, help='The output file')

    args = parser.parse_args()
    if args.sub == 'deconv':
        print("This method is not complete")
    elif args.sub == 'entropy':
        print("This method is not complete")
