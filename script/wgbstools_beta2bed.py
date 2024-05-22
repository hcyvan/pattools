import numpy as np
import gzip
import argparse
import os

def merge_binary_text(binary_file, mapinfo, output_path):
    # 读取二进制文件
    binary_data = np.fromfile(binary_file, dtype=np.uint8).reshape((-1, 2))
    
    # 计算浮点数数据
    with np.errstate(divide='ignore', invalid='ignore'):
        divisor = binary_data[:, 1]
        meth_rate = np.where(divisor != 0, np.round((binary_data[:, 0] / divisor), 4), -1)
    
    # 读取.gz压缩的文本文件并逐行合并
    merged_data = []
    with gzip.open(mapinfo, 'rt') as f:
        for i, line in enumerate(f):
            line = line.strip().split('\t')  # 使用'\t'作为分隔符
            line[1:1] = [str(int(line[1]) - 1)]
            merged_line = [str(item) for item in line] + [str(meth_rate[i]), str(binary_data[i, 0]), str(binary_data[i, 1])]
            merged_data.append(merged_line)
    
    # 保存为文本文件
    base_name = os.path.splitext(os.path.basename(binary_file))[0]
    output_file = f"{output_path}/{base_name}.bed"
    np.savetxt(output_file, merged_data, delimiter='\t', fmt='%s')

# 创建命令行参数解析器
parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input_path', help='Input binary path')
parser.add_argument('-f', dest='path_file', help='Input binary path file ')
parser.add_argument('-o', dest='output_dir', help='Output dir')
parser.add_argument('-c', dest='cpg_coordinate', required=True, help='cpg_corrdinate file')
args = parser.parse_args()

# 调用函数进行文件合并
if args.input_path:
    merge_binary_text(args.input_path, args.cpg_coordinate, args.output_dir)

if args.path_file:
    with open(args.path_file, 'r') as file:
        binary_files = file.read().splitlines()

    for binary_file in binary_files:
        merge_binary_text(binary_file, args.cpg_coordinate, args.output_dir)