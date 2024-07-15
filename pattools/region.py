import sys
from pathlib import Path
from collections import OrderedDict
from typing import List
import pysam
from .io import Output


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
        TODO: If the number of regions is excessively large (for instance, reaching 100,000), performance will be
         significantly degraded, necessitating the consideration of alternative approaches.
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
        """
        TODO: If the number of regions is excessively large (for instance, reaching 100,000), performance will be
         significantly degraded, necessitating the consideration of alternative approaches.
        """
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


def region_cpg2genome(regions, cpg, cpg_index=None, show=True):
    gr = GenomicRegion(cpg, csi_cpg=cpg_index)
    regions_map = gr.cpg_to_genomic_idx(regions)
    if show is True:
        for k, v in regions_map.items():
            print(v)
    return regions_map


def region_genome2cpg(regions, cpg, show=True):
    gr = GenomicRegion(cpg)
    regions_map = gr.genomic_to_cpg_idx(regions)
    if show is True:
        for k, v in regions_map.items():
            print(v)
    return regions_map


def trans_region_file(input, out_put, cpg_bed, transform="cpg2genome", col='col3', out_format='bed', end_offset=0):
    regions = []
    file_path = Path(input)
    with file_path.open(mode='r') as f:
        for line in f:
            line = line.strip()
            items = line.split()
            if col == 'col2':
                regions.append(f"{items[0]}:{items[1]}-{int(items[1]) + end_offset}")
            else:
                regions.append(f"{items[0]}:{items[1]}-{items[2]}")
    region_map = OrderedDict()
    if transform == 'cpg2genome':
        region_map = region_cpg2genome(regions, cpg_bed, show=False)
    elif transform == 'genome2cpg':
        region_map = region_genome2cpg(regions, cpg_bed, show=False)
    with file_path.open(mode='r') as f, Output(out_put, bgzip=False) as fo:
        for line in f:
            line = line.strip()
            items = line.split()
            if col == 'col2':
                region_cpg = f"{items[0]}:{items[1]}-{int(items[1]) + end_offset}"
            else:
                region_cpg = f"{items[0]}:{items[1]}-{items[2]}"
            region_genome = region_map[region_cpg]
            chrom = region_genome.split(':')[0]
            start = int(region_genome.split(':')[1].split('-')[0])
            end = int(region_genome.split(':')[1].split('-')[1])
            if col == 'col2':
                coli = 2
            else:
                coli = 3
            if out_format == 'bed':
                fo.write(f"{chrom}\t{start - 1}\t{end}\t" + "\t".join(items[coli:]) + "\n")
            else:
                if col == 'col2':
                    fo.write(f"{chrom}\t{start}\t" + "\t".join(items[coli:]) + "\n")
                else:
                    fo.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(items[coli:]) + "\n")
