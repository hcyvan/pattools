import gzip

from pattools.vector.window import VectorWindow, MvcFormat
from pattools.io import Output


def vector_diff(input_file, group='0', output_file=None):
    with Output(filename=output_file, file_format='motif', bgzip=False) as fo:
        with gzip.open(input_file, 'rt') as fi:
            for l in fi:
                vw = VectorWindow(l)
                if vw.cluster_count >= 2 and vw.vector_count > 80:
                    if vw.clusters.filter(group, 1.1, 0.5):
                        fo.write(vw.__str__() + '\n')


def mv_separating(input_file, group, mvs_frac=1.0, samples_frac=0.9, output_file=None):
    with Output(filename=output_file) as of:
        mvc = MvcFormat()
        mvc.filter(input_file, of, group, mvs_frac, samples_frac)
