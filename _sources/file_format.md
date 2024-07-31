# File format

## CpG start indexed format

The CpG Start Indexed (CGS) format comprises two mandatory columns along with several informational columns:

- __Sequence name__ [chrom]: chromosome number or sequence names
- __CpG index__ [cpg_index]: This column lists the sequential positions of CpG sites., 1-based, end-inclusive
  position of the CpG index
- __Informational columns__ [optional]: column 3 to column *n*,

```
chr1	1	c3	c4
chr1	2	c3	c4
chr1	3	c3	c4
chr1	4	c3	c4
```

The file header is prefixed with a # symbol. The entire file is compressed using *bgzip* and indexed with
*tabix*. Below is an example:

```
bgzip demo.cgs
tabix -C -s 1 -b 2 -e 2 demo.cgs.gz
```

### PAT format

PAT is a CGS-based format, containing 4 columns:

- __Sequence name__
- __CpG index__
- __Methylation motif__: This column denotes the methylation status of the CpG site
    - 'C': methylated CpGs
    - 'T': unmethylated CpGs
    - '.': unknown methylation status
- __Motif count__: This column records the frequency of each motif described in column 3

```
chr1    755     CCCTCCCCTCTTCCT 1
chr1    755     TTTT    2
chr1    756     CCCCCTC 1
chr1    756     CCCCCT....CCCC  10
chr1    758     CCCCCCCCCCCC    4
chr1    758     CCTTCCCTCCC     1
```

### MV format

MV (methylation vector) format is a CGS-based format

### MVC format

MVC (methylation vector cluster) format is a CGS-based format

### MVM format

MVM (methylation vector matrix) format is a CGS-based format
