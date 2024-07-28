# Methylation vector

In the PAT format file, methylated CpG sites are denoted by 'C', and unmethylated CpG sites are
denoted by 'T'. In the vector representation of methylation, a '1' corresponds to a methylated
site, while a '0' corresponds to an unmethylated site. Consequently, the sequence 'CCTTCT' in the
PAT format can be represented by the vector [1, 1, 0, 0, 1, 0]. Continuous CpG sites of length $n$
are used as a window, which contains several CpG motifs of the same length. These motifs can be
converted into vector format. For instance, in a window of length 3, the following motif is included:

$$
\begin{pmatrix}
1 & 1 & 1 \\
1 & 1 & 1 \\
1 & 1 & 1 \\
1 & 1 & 0 \\
0 & 0 & 0 \\
\end{pmatrix}
$$

The vectors corresponding to motifs within a window of length $n$ are situated within an $n$-dimensional
cartesian coordinate system.

## Usage

```{image} image/vector-flow.png
:alt: fishy
:class: bg-primary mb-1
:width: 600px
:align: center
```

### mv-vectorization

mv-vectorization is used to vectorize .pat files. You can adjust the window sizes based on your specific requirements.

```
pattools mv-vectorization -w 5 -i xxx.pat.gz -o xxx.mv.gz
```

### mv-clustering

mv-clustering is a tool used to merge all samples and perform clustering on the combined MVs

```
pattools mv-clustering -i mv_group_sample.input -c /path/to/CpG.bed.cpg.gz -o ./lung.mvc.gz
```

### mv-separating

```
pattools mv-separating -i ./lung.mvc.gz -g 1 --frac-mvs 0.8 --frac-samples 0.1 -o ./LUAD.mvc
```

### mv-extract

mv-extract is a tool designed to extract MVs from .mv and .mvc files, outputting the extracted data as an MVs matrix.

```
pattools mv-extract -i mvc_group.input -r ../LUAD.mvc -o LUAD.group.mvm
pattools mv-extract -i mv_sample.input -r ../LUAD.mvc -o LUAD.sample.mvm
```

### mv-find

The mv-find tool is designed to identify methylation vectors within smvc (specific methylation vector
clusters, obtained via mv-separating) across one or more samples.

```
pattools mv-find -i mv_sample.input -m ../LUAD.mvc -o LUAD.GSE186458.mvh
pattools mv-find -i *.mv.gz -m ../LUAD.mvc -o LUAD.GSE186458.mvh
```

#### mv-single-clustering

This command exclusively clusters the .mv file of a single sample, without merging multiple samples or including any
grouping information.

```
pattools mv-single-clustering -i /mnt/d/data/epiLungCancer/intermediate/vector/w5/mv/E05A2071.mv.gz -w5
```

## Input File List

- mv-group-sample list
- mv-sample list
- mvc-group list

### mv-sample list

[*.mv.gz]\t[sample-label]

```
/path/to/GSM5652176_Adipocytes-Z000000T7.mv.gz        Adipocytes-Z000000T7
/path/to/GSM5652177_Adipocytes-Z000000T9.mv.gz        Adipocytes-Z000000T9
/path/to/GSM5652179_Aorta-Endothel-Z00000422.mv.gz    Aorta-Endothel-Z00000422
/path/to/GSM5652180_Aorta-Endothel-Z0000043G.mv.gz    Aorta-Endothel-Z0000043G
```

