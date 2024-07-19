# Methylation vector

In the PAT format file, methylated CpG sites are denoted by 'C', and unmethylated CpG sites are denoted by 'T'. In the vector representation of methylation, a '1' corresponds to a methylated site, while a '0' corresponds to an unmethylated site. Consequently, the sequence 'CCTTCT' in the PAT format can be represented by the vector [1, 1, 0, 0, 1, 0]. Continuous CpG sites of length $n$ are used as a window, which contains several CpG motifs of the same length. These motifs can be converted into vector format. For instance, in a window of length 3, the following motif is included:

$$
\begin{pmatrix}
1 & 1 & 1 \\
1 & 1 & 1 \\
1 & 1 & 1 \\
1 & 1 & 0 \\
0 & 0 & 0 \\
\end{pmatrix}
$$

The vectors corresponding to motifs within a window of length $n$ are situated within an $n$-dimensional Cartesian coordinate system.


