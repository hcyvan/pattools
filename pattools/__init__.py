import argparse
from pathlib import Path
from pattools.deconv import deconvolution_sun, deconvolution_moss, deconvolution_loyfer, print_cell_type_helper
from pattools.entropy import extract_entropy
from pattools.beta import extract_beta
from pattools.format import pat2motif
from pattools.vector import extract_vector, extract_vector_from_multi_motif_file, vector_diff
from pattools.matrixgenerate import matrix_generate
from pattools.region import region_cpg2genome, region_genome2cpg


def main():
    parser = argparse.ArgumentParser(prog='pattools',
                                     description='pattools is a BS-seq analysis tool suite based on pat format')
    parser.add_argument('-v', '--version', action='version', version='0.0.17')
    parser.add_argument('-q', '--quiet', action='store_true', help='print run details to stderr')
    subparsers = parser.add_subparsers(dest='sub', required=True, title='command', description='The available commands',
                                       help='select a sub command to use')
    # =====================================================================
    parser_region = subparsers.add_parser('region',
                                          help='This command is used to convert a region between different coordinate systems, such as genome index coordinates or CpG index coordinates.')
    parser_region.add_argument('-t', '--transform', choices=['cpg2genome', 'genome2cpg'], default='cpg2genome',
                               help='Conversion direction of genomic coordinates')
    parser_region.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
    parser_region.add_argument('-i', '--input', required=True,
                               help='The input can be: 1) region string, eg: chr1:250-300, 2) a file where each line is a region string.')
    # =====================================================================
    parser_deconv = subparsers.add_parser('deconv',
                                          help='This command is used to calculate the cellular composition of the '
                                               'sample according to the pat format.')
    parser_deconv.add_argument('-p', '--pat', required=True, help='The pat file')
    parser_deconv.add_argument('-m', '--method', choices=['sun', 'moss', 'loyfer'], default='sun',
                               help='The deconvolution method.'
                                    'sun: Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing'
                                    ' for noninvasive prenatal, cancer, and transplantation assessments. '
                                    'moss: Moss et al. Comprehensive human cell-type methylation atlas reveals origins'
                                    ' of circulating cell-free DNA in health and disease. '
                                    'loyfer: Loyfer et al. A DNA methylation atlas of normal human cell types.'
                               )
    parser_deconv.add_argument('-a', '--optimization-algorithm', choices=['nnls', 'qp'], default='nnls',
                               help='The optimization algorithm for deconvolution.')
    parser_deconv.add_argument('-g', '--genome-version', choices=['hg38', 'hg19'], default='hg38',
                               help='The genome version.')
    parser_deconv.add_argument('--panel', choices=['U25', 'U250'], default='U25',
                               help='The panel of the markers. Only works when --method loyfer is used')
    parser_deconv.add_argument('-f', '--markerfile', default='Atlas.U25.l4.hg38.full.tsv',
                               help='markerfile for loyfer method')
    parser_deconv.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
    parser_deconv.add_argument('-o', '--out', required=True, help='The output file')
    parser_deconv.add_argument('-in', '--include', nargs='+',
                               help='Cell/Tissue types to include. Complementary to --ignore. Use `pattools '
                                    'deconv-helper` to see supported cell/tissue types.')
    parser_deconv.add_argument('-ig', '--ignore', nargs='+',
                               help='Cell/Tissue types to remove. Use `pattools deconv-helper` to see supported '
                                    'cell/tissue types.')
    # --------------------------------------------------------------------
    parser_deconv_helper = subparsers.add_parser('deconv-helper',
                                                 help='print the cell types support by pattools deconv')
    parser_deconv_helper.add_argument('-m', '--method', choices=['sun', 'moss', 'loyfer'], default='sun',
                                      help='The deconvolution method.'
                                           'sun: Sun et al. Plasma DNA tissue mapping by genome-wide methylation sequencing'
                                           ' for noninvasive prenatal, cancer, and transplantation assessments. '
                                           'moss: Moss et al. Comprehensive human cell-type methylation atlas reveals origins'
                                           ' of circulating cell-free DNA in health and disease. '
                                           'loyfer: Loyfer et al. A DNA methylation atlas of normal human cell types.'
                                      )
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
    parser_vector.add_argument('-w', '--window', type=int, default='4',
                               help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')
    parser_vector.add_argument('-o', '--out', default=None,
                               help='The output file, If not set, output is sent to standard output.')
    # =====================================================================
    parser_vector_multi = subparsers.add_parser('vector-multi',
                                                help='Extract vectors from multiple samples from different groups and'
                                                     ' analyze them. This command supports MPI, which can accelerate'
                                                     ' calculations in HPC')
    parser_vector_multi.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
    parser_vector_multi.add_argument('-i', '--input', required=True,
                                     help='a list file in tsv format, which contains multiple sample files,'
                                          ' sample grouping, etc. eg: <MOTIF_FILE>  <GROUP_LABEL>')
    parser_vector_multi.add_argument('-w', '--window', type=int, default='4',
                                     help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')
    parser_vector_multi.add_argument('-p', '--process', type=int, default=1,
                                     help='The number of processes used for processing')
    parser_vector_multi.add_argument('-r', '--region', default=None,
                                     help='TThe region to be processed. If not set, the entire genome is processed.'
                                          ' eg: -r chr1:10000-15000')
    parser_vector_multi.add_argument('-m', '--cluster-method', choices=['HDBSCAN', 'DBSCAN'], default='HDBSCAN',
                                     help='Algorithm for classifying all motifs in a window')
    parser_vector_multi.add_argument('-o', '--out', default=None,
                                     help='The output file, If not set, output is sent to standard output.')
    # =====================================================================
    parser_vector_diff = subparsers.add_parser('vector-diff',
                                               help='Identify the window of the differential vector within'
                                                    ' the merged file. (generate by vector-multi)')
    parser_vector_diff.add_argument('-i', '--input', required=True, help='The input merged vector files.'
                                                                         ' (generate by vector-multi)')
    parser_vector_diff.add_argument('-o', '--out', default=None,
                                    help='The output file, If not set, output is sent to standard output.')
    parser_vector_diff.add_argument('-g', '--group', required=True, help='Output group-specific '
                                                                         'differential vector window')
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
    if args.sub == 'deconv':
        if args.method == 'sun':
            deconvolution_sun(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                              optimization=args.optimization_algorithm)
        if args.method == 'moss':
            deconvolution_moss(args.pat, args.out, args.genome_version, cpg_bed=args.cpg_bed,
                               optimization=args.optimization_algorithm)
        if args.method == 'loyfer':
            deconvolution_loyfer(args.pat, args.out, args.markerfile, args.genome_version, args.panel,
                                 cpg_bed=args.cpg_bed,
                                 optimization=args.optimization_algorithm, include=args.include, ignore=args.ignore)
    if args.sub == 'deconv-helper':
        print_cell_type_helper(args.method)
    if args.sub == 'entropy':
        extract_entropy(args.input, args.depth, args.window, args.out)
    if args.sub == 'beta':
        extract_beta(args.input, args.depth, args.out)
    if args.sub == 'vector':
        extract_vector(args.input, args.out, window=args.window)
    if args.sub == 'vector-multi':
        extract_vector_from_multi_motif_file(args.input, args.cpg_bed, args.out, window=args.window,
                                             process=args.process, region=args.region, cluster=args.cluster_method)
    if args.sub == 'vector-diff':
        vector_diff(args.input, args.group, args.out)
    if args.sub == 'pat2motif':
        pat2motif(args.input, args.out, args.window, not args.text)
    if args.sub == 'matgen':
        matrix_generate(args.input, args.coordinate, args.depth, args.exclude_mode, args.out)
    if args.sub == 'region':
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
