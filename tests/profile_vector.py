from pattools.vector.vector import extract_vector_multi

file_list = '/mnt/d/data/epiLungCancer/intermediate/vector/w4/sample_group_sample.input'
cpg_bed = '/mnt/d/project/wgbs_tools/references/hg38/CpG.bed.cpg.gz'
outfile = '/mnt/d/project/pattools/tmp/aaa.txt'
extract_vector_multi(file_list, cpg_bed, outfile, window=4, regions='chr1:1-10000', cluster='MRESC', out_version='v1')