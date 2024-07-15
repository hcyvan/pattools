import argparse
from pattools.entropy import extract_entropy
from pattools.beta import extract_beta
from pattools.format import pat2motif
from pattools.matrixgenerate import matrix_generate
from pattools.cmd import CmdFactory
from pattools.vector import *
from pattools.region import *
from pattools.deconv import *


def main():
    parser = argparse.ArgumentParser(prog='pattools',
                                     description='pattools is a BS-seq analysis tool suite based on pat format')
    parser.add_argument('-v', '--version', action='version', version='0.0.17')
    parser.add_argument('-q', '--quiet', action='store_true', help='print run details to stderr')
    subparsers = parser.add_subparsers(dest='sub', required=True, title='command', description='The available commands',
                                       help='select a sub command to use')
    CmdFactory.set_subparser(subparsers)

    # =====================================================================
    parser_entropy = subparsers.add_parser('entropy',
                                           help='This command performs entropy analysis on the sample')
    parser_entropy.add_argument('-i', '--input', required=True, help='Input file, *.pat.gz format')
    parser_entropy.add_argument('-d', '--depth', type=int, default='3',
                                help='the minimum total count required to calculate entropy')
    parser_entropy.add_argument('-w', '--window', type=int, default='4',
                                help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')
    parser_entropy.add_argument('-o', '--out', required=True,
                                help='The output file, *.gz format. There are four columns in total, '
                                     'representing chromosome, index, entropy, and total sequencing depth of loci')
    # =====================================================================
    parser_beta = subparsers.add_parser('beta',
                                        help='This command performs methylation ratio analysis on the sample')
    parser_beta.add_argument('-i', '--input', required=True, help='Input file, *.pat.gz format')
    parser_beta.add_argument('-d', '--depth', type=int, default='1',
                             help='the minimum total count required to calculate methylation ratio')
    parser_beta.add_argument('-o', '--out', required=True,
                             help='The output file, *.gz format. There are four columns in total,'
                                  'representing chromosome, index, methylation ratio, and total sequencing depth of loci')
    # =====================================================================
    parser_pat2motif = subparsers.add_parser('pat2motif',
                                             help='This command is used to convert pat file to motif file')
    parser_pat2motif.add_argument('-i', '--input', required=True, help='The input file')
    parser_pat2motif.add_argument('-o', '--out', default=None,
                                  help='The output file, If not set, output is sent to standard output.')
    parser_pat2motif.add_argument('--text', action='store_true', help='If set, files are not '
                                                                      'compressed with bgzip')
    parser_pat2motif.add_argument('-w', '--window', type=int, default='4',
                                  help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')
    # ======================================================================
    parser_matrix_generate = subparsers.add_parser('matgen',
                                                   help='This command is used to generate matrix for entropy and beta')
    parser_matrix_generate.add_argument('-i', '--input', required=True, help='This is a text, with each line being the '
                                                                             'path to each entropy or beta file')
    parser_matrix_generate.add_argument('-o', '--out', required=True, help='The output is a standard bed format file')
    parser_matrix_generate.add_argument('-c', '--coordinate', required=True,
                                        help='This is a standard CpG coordinate file '
                                             'The current path is /PUBLIC/rd/lung_cac/rawdata/cpgMapinfo')
    parser_matrix_generate.add_argument('-d', '--depth', type=int, default=3, help='the lowest depth of a matrix')
    parser_matrix_generate.add_argument('-e', '--exclude_mode', default='all',
                                        help='exclude -1 mode: all - exclude if all sample is -1,'
                                             'one - exclude if contain one -1, close exclude mode')
    # ======================================================================    
    args = parser.parse_args()
    if args.sub == 'entropy':
        extract_entropy(args.input, args.depth, args.window, args.out)
    if args.sub == 'beta':
        extract_beta(args.input, args.depth, args.out)
    if args.sub == 'pat2motif':
        pat2motif(args.input, args.out, args.window, not args.text)
    if args.sub == 'matgen':
        matrix_generate(args.input, args.coordinate, args.depth, args.exclude_mode, args.out)
    CmdFactory.run(args)
