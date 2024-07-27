from pattools.vector.clustering import methylation_vector_cluster
from pattools.vector.separating import mv_separating
from pattools.vector.support import extract_mvs, single_cluster, fix_mvc, fix_mv
from pattools.vector.vectorization import pat2mv
from pattools.cmd import command, Cmd
from pattools.vector.utils import get_cpg_index_regions


@command('mv-extract', 'extract methylation vectors from .mv or .mvc files')
class VectorExtractCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-r', '--region', required=True,
                            help='This parameter can be specified in one of the following formats: (1) a genomic region'
                                 'string, for example, chr4:6909897-6909899; or (2) a CGS format file')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')

    def do(self, args):
        regions = get_cpg_index_regions(args.region)
        extract_mvs(args.input, regions, args.out)


@command('mv-single-clustering', 'This command exclusively clusters the .mv file of a single sample, '
                                 'without merging multiple samples or including any grouping information.'
                                 'relate `mv-clustering`')
class VectorSingleClusteringCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file, *.mvc.gz')
        parser.add_argument('-w', '--window', type=int, default='4',
                            help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')

    def do(self, args):
        single_cluster(args.input, args.out, window=args.window)


@command('mv-vectorization', 'Methylation vectors vectorization')
class VectorizationCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='The input file')
        parser.add_argument('--out-version', default='v2', help='The output file format version')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')
        parser.add_argument('--text', action='store_true', help='If set, files are not '
                                                                'compressed with bgzip')
        parser.add_argument('-w', '--window', type=int, default='4',
                            help='Define the length of motif, such as ''3:CCT; 4: CCTT; 5:CCTTT'' ')

    def do(self, args):
        pat2mv(args.input, args.out, window=args.window, bgzip=(not args.text), out_version=args.out_version)


@command('mv-clustering',
         'Methylation vectors clustering. This command supports MPI, which can accelerate calculations in HPC')
class VectorClusteringCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
        parser.add_argument('-i', '--input', required=True,
                            help='a list file in tsv format, which contains multiple mv files,'
                                 ' group info and sample info, etc. eg: <MV_FILE>  <GROUP_LABEL> <SAMPLE_LABEL>')
        parser.add_argument('-r', '--region', default=None,
                            help='TThe region to be processed. If not set, the entire genome is processed.'
                                 ' eg: -r chr1:10000-15000')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')
        parser.add_argument('-g', '--groups', default=None, help="The groups to be entered should be separated by "
                                                                 "commas ','. If not set, all will be entered")
        parser.add_argument('-w', '--window', type=int,
                            help='[DEPRECATED] Define the length of MV, such as ''3:CCT; 4: CCTT; 5:CCTTT''')
        parser.add_argument('-p', '--process', type=int, default=1,
                            help='The number of processes used for processing')
        parser.add_argument('-m', '--cluster-method', choices=['HDBSCAN', 'DBSCAN', 'MRESC'],
                            default='MRESC',
                            help='Algorithm for classifying all MVs in a window')

    def do(self, args):
        methylation_vector_cluster(args.input, args.cpg_bed, args.out, window=args.window,
                                   process=args.process, region=args.region, cluster=args.cluster_method,
                                   groups=args.groups)


@command('mv-separating', 'Identify and separate distinct MVs clusters. (generate by mv-clustering)')
class VectorSeparatingCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='The input merged vector files.'
                                                                 ' (generate by vector-multi)')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')
        parser.add_argument('-g', '--group', required=True, help='Output group-specific '
                                                                 'differential vector window')
        parser.add_argument('--frac-mvs', default=0.8, type=float, help='')
        parser.add_argument('--frac-samples', default=0.1, type=float, help='')
        parser.add_argument('--with-meta', action='store_true', help='')

    def do(self, args):
        mv_separating(args.input, args.group, args.frac_mvs, args.frac_samples, output_file=args.out,
                      with_meta=args.with_meta)


@command('mv-mvc-fix', 'fix mvc file')
class VectorFixCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-c', '--cpg-bed', required=True, help='The cpg_bed file of the selected genome.')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')

    def do(self, args):
        fix_mvc(args.input, args.cpg_bed, args.out)


@command('mv-mv-fix', 'fix mv file')
class VectorFixCmd(Cmd):
    def add_argument(self, parser):
        parser.add_argument('-i', '--input', required=True, help='Input file list')
        parser.add_argument('-o', '--out', default=None,
                            help='The output file, If not set, output is sent to standard output.')

    def do(self, args):
        fix_mv(args.input, args.out)
