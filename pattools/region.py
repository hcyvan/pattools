import pysam
from collections import OrderedDict


class GenomicRegion:
    """
    Genomic index (g_idx) is the coordinate system of the genome sequence, which is 1-based index.
    CpG index (cpg_idx) is the coordinate system of CpGs, which is 1-based index.
    Region string (region_str) is a string formatted as chrom:start-end, which is 1-based and end included.
    """

    def __init__(self, cpg_bed):
        self.cpg_bed = cpg_bed
        self.cpg_bed_csi = cpg_bed + ".csi"

    def genomic_to_cpg_idx(self, g_idx_regions):
        """
        :param g_idx_regions: Genomic index region string array
        :return: Ordered dictionary from Genomic index regions to CpG index regions
        """
        gr_cpg_map = OrderedDict()
        with pysam.TabixFile(self.cpg_bed, index=self.cpg_bed_csi) as tbx:
            for region in g_idx_regions:
                tmp = []
                for record in tbx.fetch(region=region):
                    tmp.append(record)
                if len(tmp):
                    chrom = tmp[0].split('\t')[0]
                    cpg_start = tmp[0].split('\t')[2]
                    cpg_end = tmp[-1].split('\t')[2]
                    gr_cpg_map[region] = f'{chrom}:{cpg_start}-{cpg_end}'
                else:
                    gr_cpg_map[region] = None
        return gr_cpg_map
