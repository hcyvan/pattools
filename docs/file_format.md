# File format

## CpG start indexed format

The CpG Start Indexed (CGS) format comprises two mandatory columns along with several informational columns:

- chrom: chromosome number or sequence names
- cpg_index: 1-based, end-inclusive position of the CpG index
- Informational columns [optional]: column 3 to column *n*, 

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

PAT is a CGS-based format

### MV format

MV (methylation vector) format is a CGS-based format

### MVC format

MVC (methylation vector cluster) format is a CGS-based format
