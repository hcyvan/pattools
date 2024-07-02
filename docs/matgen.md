# Matrix generate

```matgen``` This command is used to generate matrix for entropy and beta

## Usage: 
pattools matgen [-h] -i INPUT -o OUT -c COORDINATE [-d DEPTH] [-e EXCLUDE_MODE]

## Example:
```
pattools matgen -d 3 -e all \
    -i /PUBLIC/rd/lung_cac/public/shichengguo/testbeta_path.txt \
    -o /PUBLIC/rd/lung_cac/public/shichengguo/matrix.bed \
    -c /PUBLIC/rd/lung_cac/rawdata/cpgMapinfo/cpg_coordinates_index.bed
```
input
```
cat testbeta_path

/PUBLIC/rd/lung_cac/public/shichengguo/test1.beta.gz
/PUBLIC/rd/lung_cac/public/shichengguo/test2.beta.gz
/PUBLIC/rd/lung_cac/public/shichengguo/test3.beta.gz
```
test1
```
zcat test1.beta.gz | head -n 20

chr1	2728	-1	1
chr1	2729	-1	1
chr1	5961	-1	1
chr1	5962	-1	1
chr1	5963	-1	1
chr1	5964	-1	1
chr1	5965	-1	1
chr1	5966	-1	1
chr1	7395	-1	1
chr1	7396	-1	1
chr1	7397	-1	1
chr1	7398	-1	1
chr1	7399	-1	1
chr1	7400	-1	1
chr1	7442	0.0	5
chr1	7443	0.0	5
chr1	7574	0.0526	19
chr1	7575	0.0526	19
chr1	7576	0.0526	19
chr1	7577	0.0526	19
```
test2
```
zcat test2.beta.gz | head -n 20

chr1	5930	-1	1
chr1	5961	-1	1
chr1	5962	-1	1
chr1	5963	-1	1
chr1	6018	-1	1
chr1	7395	-1	1
chr1	7420	-1	1
chr1	7421	-1	1
chr1	7422	-1	1
chr1	7423	-1	1
chr1	7424	-1	1
chr1	7425	-1	1
chr1	7442	0.0	3
chr1	7443	-1	2
chr1	7553	0.1818	11
chr1	7554	0.3333	6
chr1	7555	0.5	4
chr1	7574	0.0	30
chr1	7575	0.037	27
chr1	7576	0.037	27
```
test3
```
zcat test3.beta.gz | head -n 30 | tail -n 20


chr1	5961	-1	1
chr1	5962	-1	1
chr1	5963	-1	1
chr1	5964	-1	1
chr1	5965	-1	1
chr1	5966	-1	1
chr1	6004	-1	2
chr1	6005	-1	2
chr1	7395	-1	1
chr1	7396	-1	1
chr1	7397	-1	1
chr1	7398	-1	1
chr1	7399	-1	1
chr1	7400	-1	1
chr1	7401	-1	1
chr1	7420	-1	1
chr1	7421	-1	1
chr1	7422	-1	1
chr1	7423	-1	1
chr1	7424	-1	1
chr1	7425	-1	1
chr1	7442	0.0	12
chr1	7443	0.0	12
chr1	7553	0.2353	17
chr1	7554	0.5385	13
chr1	7555	0.25	12
chr1	7574	0.0	20
chr1	7575	0.0526	19
chr1	7576	0.1053	19
chr1	7577	0.0	15
```

coordinate (Add build mode)

output
```
head matrix.bed

#chrom	start	end	index	test1	test2	test3
chr1	779047	779049	7442	0.0	0.0	0.0
chr1	779062	779064	7443	0.0	-1	0.0
chr1	788919	788921	7553	-1	0.1818	0.2353
chr1	788924	788926	7554	-1	0.3333	0.5385
chr1	788929	788931	7555	-1	0.5	0.25
chr1	789889	789891	7574	0.0526	0.0	0.0
chr1	789894	789896	7575	0.0526	0.037	0.0526
chr1	789899	789901	7576	0.0526	0.037	0.1053
chr1	789904	789906	7577	0.0526	0.0	0.0
```

## Options:
  
-h, --help | Show this help message and exit
  
-i INPUT, --input INPUT | This is a text, with each line being the path to each entropy or beta file
  
-o OUT, --out OUT | The output is a standard bed format file
 
-c COORDINATE, --coordinate COORDINATE | This is a standard CpG coordinate file. The current path is /PUBLIC/rd/lung_cac/rawdata/cpgMapinfo
 
-d DEPTH, --depth DEPTH | The lowest depth of a matrix
  
-e EXCLUDE_MODE, --exclude_mode EXCLUDE_MODE | exclude -1 mode: all - exclude if all sample is -1, one - exclude if contain one -1, close exclude mode