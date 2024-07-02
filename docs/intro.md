---
title: "Introduction"
---
# Introduction
The PAT format was first proposed in wgbstools ([Github](https://github.com/nloyfer/wgbs_tools)) {cite}`Loyfer2024.05.08.593132`.
This format is compressed and indexed, providing comprehensive information about methylation patterns while excluding
the subjects' private DNA sequence data. Therefore, the pat format can effectively facilitate the transmission and
sharing of methylation data.

To this end, we developed pattools, a supplementary extension of wgbstools, to support the analysis of pat format data.
Pattools currently offers the following functionalities:

+ Multiple methylation deconvolution algorithms
+ Methylation beta value analysis
+ Methylation entropy analysis
+ Methylation vector analysis


## The PAT format

The PAT format encompasses 4 columns:

1. Chromosome Number: This column specifies the chromosome number.
2. Consecutive CpG Positions: This column lists the successive genomic positions of CpG sites. The positions
   are continuous across chromosomes, maintaining a 1-based index. It indicates the position of the first
   CpG in the column 3.
3. Methylation Motif: This column denotes the methylation status of the CpG site
   - 'C': methylated CpG
   - 'T': unmethylated CpG
   - '.': unknown methylation status.
4. Motif Occurrence Count: This column records the frequency of each motif described in the third column.

```
chr1    755     CCCTCCCCTCTTCCT 1
chr1    755     TTTT    1
chr1    756     CCCCCTC 1
chr1    756     CCCCCT....CCCC  1
chr1    758     CCCCCCCCCCCC    1
chr1    758     CCTTCCCTCCC     1
```


