import gzip

from pattools.vector.window import VectorWindow
from pattools.io import Output


def vector_diff(input_file, group='0', output_file=None):
    with Output(filename=output_file, file_format='motif', bgzip=False) as fo:
        with gzip.open(input_file, 'rt') as fi:
            for l in fi:
                vw = VectorWindow(l)
                if vw.cluster_count >= 2 and vw.vector_count > 80:
                    if vw.clusters.filter(group, 1, 0.9):
                        fo.write(vw.__str__()+'\n')
