import sys

from config import CONFIG
from pattools.pat import PatWindow
from pattools.motif import Motif

patFile = CONFIG.DataRaw / 'pat' / 'test.chr21_22.pat.gz'

motif = Motif(4)
motif2vector_map = motif.motif2vector()

patWindow = PatWindow(patFile, window=4)

motif_pattern = '\t'.join(motif.motifs)
header = f'#chr\tCpG_ndex\t{motif_pattern}'
print(header)
for i, win in enumerate(patWindow):
    counter = motif.count_motifs(win[2])
    # print(motif.motif_count2vectors(counter))
    # print(dict(counter))
    counts = '\t'.join([str(x) for x in counter.values()])
    print(f'{win[0]}\t{win[1]}\t{counts}')

    # if int(win[1]) == 26802539:
    #     print(win[2])
    #     break
