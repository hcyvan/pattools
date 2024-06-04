import argparse
from pattools.deconv import deconvolution_sun, deconvolution_moss, deconvolution_loyfer
from pattools.entropy import extract_entropy
from pattools.beta import extract_beta
from pattools.format import pat2motif
from pattools.matrixgenerate import matrix_generate

def main():
    parser = argparse.ArgumentParser(prog='pattools',
                                     description='pattools is a BS-seq analysis tool suite based on pat format')
    parser.add_argument('-v', '--version', action='version', version='0.0.11')
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
    parser_region.add_argument('-f', '--markerfile', default='Atlas.U25.l4.hg38.tsv',
                               help='markerfile for loyfer method')
    parser_region.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
    parser_region.add_argument('-o', '--out', required=True, help='The output file')
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
    parser_vector = subparsers.add_parser('vector',
                                          help='This command performs vector analysis on the sample')
    parser_vector.add_argument('-i', '--input', required=True, help='Input file, *.motif.gz')
    parser_vector.add_argument('-o', '--out', default=None,
                               help='The output file, If not set, output is sent to standard output.')
    # =====================================================================
    parser_vector_multi = subparsers.add_parser('vector-multi',
                                                help='This command performs vector analysis on the sample')
    parser_vector_multi.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
    parser_vector_multi.add_argument('-i', '--input', required=True,
                                     help='a list file in tsv format, which contains multiple sample files,'
                                          ' sample grouping, etc. eg: <MOTIF_FILE>  <GROUP_LABEL>')
    parser_vector_multi.add_argument('-w', '--window', type=int, default='4',
                                     help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')
    parser_vector_multi.add_argument('-o', '--out', default=None,
                                     help='The output file, If not set, output is sent to standard output.')
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
    parser_matrix_generate = subparsers.add_parser('matrix_generate',
                                             help='This command is used to generate matrix for entropy and beta')
    parser_matrix_generate.add_argument('-i', '--input', required=True, help='The input file')
    parser_matrix_generate.add_argument('-o', '--out', default=None, help='The output file')
    parser_matrix_generate.add_argument('-c', '--coordinate', help='If set, files are not '
                                                                      'compressed with bgzip')
    parser_matrix_generate.add_argument('-d', '--depth', type=int, help='If set, files are not '
                                                                      'compressed with bgzip')
    parser_matrix_generate.add_argument('-e', '--exclude_mode', help='If set, files are not '
                                                                      'compressed with bgzip')
    # ======================================================================    
    args = parser.parse_args()
    if args.sub == 'deconv':
        if args.method == 'sun':
            deconvolution_sun(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                              optimization=args.optimization_algorithm)
        if args.method == 'moss':
            deconvolution_moss(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                               optimization=args.optimization_algorithm)
        if args.method == 'loyfer':
            deconvolution_loyfer(args.pat, args.out, args.markerfile, args.genome_version, cpg_bed=args.cpg_bed,
                                 optimization=args.optimization_algorithm)
    if args.sub == 'entropy':
        extract_entropy(args.input, args.depth, args.window, args.out)
    if args.sub == 'beta':
        extract_beta(args.input, args.depth, args.out)
    if args.sub == 'vector':
        extract_vector(args.input, args.out)
    if args.sub == 'vector-multi':
        extract_vector_from_multi_motif_file(args.input, args.cpg_bed, args.out, window=args.window)
    if args.sub == 'pat2motif':
        pat2motif(args.input, args.out, args.window, not args.text)
    if args.sub == 'matrix_generate':
        matrix_generate(args.input, args.coordinate, args.depth, args.exclude_mode, args.out)        
