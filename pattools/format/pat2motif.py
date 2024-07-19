import sys
from collections import OrderedDict
from pattools.pat import PatWindow
from pattools.motif import Motif
from pattools.io import Output


def pat2motif(filename: str, outfile: str = None, window: int = 4, bgzip: bool = True):
    motif = Motif(window)
    pat_window = PatWindow(filename, window=window)
    motif_pattern = '\t'.join(motif.motifs)
    with Output(filename=outfile, file_format='motif', bgzip=bgzip) as of:
        of.write(f"##FORMAT: mv\n")
        of.write(f"##WINDOW: {window}\n")
        of.write(f"##COMMAND: {' '.join(sys.argv)}\n")
        comment_header = f'#chr\tCpG_index\t{motif_pattern}\n'
        of.write(comment_header)
        for i, win in enumerate(pat_window):
            counter: OrderedDict[str, int] = motif.count_motifs(win[2])
            counts = '\t'.join([str(x) for x in counter.values()])
            of.write(f'{win[0]}\t{win[1]}\t{counts}\n')
