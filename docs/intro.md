---
title: "Introduction"
---
# Introduction

Pattools is an advanced toolkit designed for the comprehensive analysis of Bisulfite sequencing (BS-Seq) data.
This tool enables researchers to seamlessly evaluate a variety of DNA methylation metircs using BS-seq data
and to implement algorithms developed by other researchers without the need for data format conversion.

The PAT format was first proposed in wgbstools ([Github](https://github.com/nloyfer/wgbs_tools)) {cite}`Loyfer2024.05.08.593132`.
This format is compressed and indexed, providing comprehensive information about methylation patterns while excluding
the subjects' private DNA sequence data. Therefore, the pat format can effectively facilitate the transmission and
sharing of methylation data.


Therefore, we selected the PAT format as the central hub for BS-seq data and developed pattools to facilitate
the analysis of PAT format data. Pattools currently offers the following functionalities:

+ Multiple methylation deconvolution algorithms
+ Methylation beta value analysis
+ Methylation entropy analysis
+ Methylation vector analysis

```{image} image/framework.png
:alt: fishy
:class: bg-primary mb-1
:width: 600px
:align: center
```