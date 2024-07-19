# Beta value

```beta``` This command performs methylation ratio analysis on the sample

## Usage: 
pattools beta [-h] -i INPUT [-d DEPTH] -o OUT

## Example:

```
pattools beta -d 3 \
    -i /PUBLIC/rd/lung_cac/public/shichengguo/test.pat.gz \
    -o /PUBLIC/rd/lung_cac/public/shichengguo/test.beta
```
input
```
zcat test.pat.gz | head

chr1	5930	C	1
chr1	7420	CTTTTT	2
chr1	7420	TTTTTT	6
chr1	7442	CT	1
chr1	7442	TT	14
chr1	7454	TT	2
chr1	7454	TT	1
chr1	7553	TT	1
chr1	7553	TTT	2
chr1	7574	CCCCC	2
```
output
```
zcat test.beta.gz | head

chr1	5930	-1	1
chr1	7420	0.25	8
chr1	7421	0.0	8
chr1	7422	0.0	8
chr1	7423	0.0	8
chr1	7424	0.0	8
chr1	7425	0.0	8
chr1	7442	0.0667	15
chr1	7443	0.0	15
chr1	7454	0.0	3
```

## Options:
  
-h, --help | Show this help message and exit
 
-i INPUT, --input INPUT | Input file, *.pat.gz format
  
-d DEPTH, --depth DEPTH | The minimum total count required to calculate methylation ratio
  
-o OUT, --out OUT | The output file, *.gz format. There are four columns in total, representing chromosome, index, methylation ratio, and total sequencing depth of loci