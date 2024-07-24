from pattools.vector.format import MvcFormat
from pattools.io import Output


def mv_separating(input_file, group, mvs_frac=1.0, samples_frac=0.9, output_file=None, with_meta=False):
    with Output(filename=output_file) as of:
        mvc = MvcFormat(input_file)
        mvc.filter(of, group, mvs_frac, samples_frac, with_meta=with_meta)
