#!/usr/bin/env python
import argparse
import gzip
import time
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(prog='methMatrixFilter.py', description='extract dmr from matrix')
    parser.add_argument('-v', '--version', action='version', version='0.0.1')
    parser.add_argument('-q', '--quiet', default=False, type=bool, help='print run details to stderr')
    parser.add_argument('-r', '--ratio-file', required=True, help='ratio file generate by methMatrixGenerate.py')
    parser.add_argument('-d', '--depth-file', required=True, help='depth file generate by methMatrixGenerate.py')

    parser.add_argument('-c', '--exclude-count', type=int, default=None,
                        help='exclude if the count of missing value is more than *--exclude-count*')
    parser.add_argument('-p', '--exclude-proportion', type=float, default=0.1,
                        help='exclude if the proportion of missing value is more than *--exclude-proportion*. This option will be overridden by *--exclude-count* ')
    parser.add_argument('-m', '--min-depth', default=1, type=int,
                        help='The totalC less than depth will be treated as missing value')

    parser.add_argument('-or', '--out-ratio', required=True, help='The output ratio file')
    parser.add_argument('-od', '--out-depth', required=True, help='The output depth file')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    excludeIfMissCount = args.exclude_count
    excludeIfMissProportion = args.exclude_proportion
    depthMinThreshold = args.min_depth
    ratioFile = args.ratio_file
    depthFile = args.depth_file
    outRatioFile = args.out_ratio
    outDepthFile = args.out_depth
    header = None
    start = time.time()
    idx = 0
    with (gzip.open(ratioFile) as ratio_file, gzip.open(depthFile) as depth_file,
          open(outRatioFile, 'w') as out_ratio_file, open(outDepthFile, 'w') as out_depth_file):
        for ratio_line, depth_line in zip(ratio_file, depth_file):
            if idx % 500000 == 0:
                print("index {}, pass {} min".format(idx, round((time.time() - start) / 60), 4))
                out_ratio_file.flush()
                out_depth_file.flush()
            idx += 1
            ratio_line = ratio_line.decode()
            depth_line = depth_line.decode()
            if header is None:
                header = ratio_line
                out_ratio_file.write(header)
                out_depth_file.write(header)
                if excludeIfMissCount is None:
                    if excludeIfMissProportion is not None:
                        excludeIfMissCount = int((len(ratio_line.split('\t')) - 3) * excludeIfMissProportion)
                    else:
                        excludeIfMissCount = 0
                continue
            if len(ratio_line) and len(depth_line):
                ratio_line = ratio_line.split('\t')
                depth_line = depth_line.split('\t')
                depth_passed = np.array([int(x) for x in depth_line[3:]]) >= depthMinThreshold
                if np.sum(~depth_passed) < excludeIfMissCount:
                    out_ratio_file.write('\t'.join(ratio_line) + "\n")
                    out_depth_file.write('\t'.join(depth_line) + "\n")
