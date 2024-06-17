import sys
import pysam
from typing import List
from collections import OrderedDict


class GenomicRegion:
    """
    Genomic index (g_idx) is the coordinate system of the genome sequence, which is 1-based index.
    CpG index (cpg_idx) is the coordinate system of CpGs, which is 1-based index.
    Region string (region_str) is a string formatted as chrom:start-end, which is 1-based and end included.
    """

    def __init__(self, cpg_bed, csi_cpg=None):
        self.cpg_bed = cpg_bed
        self.cpg_bed_csi = cpg_bed + ".csi"
        """
        Use the second and third columns as the start to create indexes.
        --begin 2: index for genomic coordinate system
        --begin 3: cpg.csi, index for CpG system
        """
        if csi_cpg:
            self.cpg_bed_cpg_csi = csi_cpg
        else:
            self.cpg_bed_cpg_csi = cpg_bed + ".cpg.csi"

    def genomic_to_cpg_idx(self, g_idx_regions: List[str]) -> OrderedDict[str, str]:
        """
        :param g_idx_regions: Genomic index region string array: such as: ['chr1:20000-20050', 'chr1:30000-30050']
        :return: Ordered dictionary from Genomic index regions to CpG index regions
        """
        gr_cpg_map = OrderedDict()
        with pysam.TabixFile(self.cpg_bed, index=self.cpg_bed_csi) as tbx:
            for region in g_idx_regions:
                tmp = []
                try:
                    fetch_regions = tbx.fetch(region=region)
                    for record in fetch_regions:
                        tmp.append(record)
                except Exception as e:
                    sys.stderr.write(f'{e}\n')

                if len(tmp):
                    chrom = tmp[0].split('\t')[0]
                    cpg_start = tmp[0].split('\t')[2]
                    cpg_end = tmp[-1].split('\t')[2]
                    gr_cpg_map[region] = f'{chrom}:{cpg_start}-{cpg_end}'
                else:
                    gr_cpg_map[region] = None
        return gr_cpg_map

    def cpg_to_genomic_idx(self, cpg_idx_regions: List[str]) -> OrderedDict[str, str]:
        cpg_gr_map = OrderedDict()
        with pysam.TabixFile(self.cpg_bed, index=self.cpg_bed_cpg_csi) as tbx:
            for region in cpg_idx_regions:
                tmp = []
                try:
                    fetch_regions = tbx.fetch(region=region)
                    for record in fetch_regions:
                        tmp.append(record)
                except Exception as e:
                    sys.stderr.write(f'{e}\n')

                if len(tmp):
                    chrom = tmp[0].split('\t')[0]
                    start = tmp[0].split('\t')[1]
                    end = tmp[-1].split('\t')[1]
                    cpg_gr_map[region] = f'{chrom}:{start}-{end}'
                else:
                    cpg_gr_map[region] = None
        return cpg_gr_map


def region_cpg2genome(regions, cpg, cpg_index=None):
    gr = GenomicRegion(cpg, csi_cpg=cpg_index)
    regions_map = gr.cpg_to_genomic_idx(regions)
    for k, v in regions_map.items():
        print(v)


def region_genome2cpg(regions, cpg):
    gr = GenomicRegion(cpg)
    regions_map = gr.genomic_to_cpg_idx(regions)
    for k, v in regions_map.items():
        print(v)
