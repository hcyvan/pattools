import re
from pathlib import Path


def parse_mv_group_sample_file(file_list, target_groups=None):
    input_files = []
    groups = []
    samples = []
    if target_groups is not None:
        target_groups = target_groups.split(',')
    with open(file_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            if target_groups is not None:
                if line[1] in target_groups:
                    input_files.append(line[0])
                    groups.append(line[1])
                    samples.append(line[2])
            else:
                input_files.append(line[0])
                groups.append(line[1])
                samples.append(line[2])
    return input_files, groups, samples


def parse_mvc_group_file(file_list):
    input_files = []
    groups = []
    with open(file_list, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            input_files.append(line[0])
            groups.append(line[1])
    return input_files, groups


def parse_region_string(region_string):
    pattern = re.compile(r'^([a-zA-Z0-9]+):(\d+)-(\d+)$')
    m = pattern.match(region_string)
    if m:
        chrom = m.group(1)
        start = int(m.group(2))
        end = int(m.group(3))
        return chrom, start, end
    return None, None, None


def get_cpg_index_regions(region):
    regions = []
    if Path(region).exists():
        with open(Path(region), 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                line = line.strip()
                if line:
                    items = line.split('\t')
                    regions.append(f'{items[0]}:{items[1]}-{items[1]}')
    else:
        chrom, start, end = parse_region_string(region)
        if chrom is not None:
            for i in range(start, end + 1):
                regions.append(f'{chrom}:{i}-{i}')
        else:
            raise Exception(f'Wrong region format: {region}')
    return regions
