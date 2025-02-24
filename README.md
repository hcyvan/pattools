<img src="docs/image/logo.png" alt="logo" width="600"/>

[![Documentation](https://img.shields.io/badge/docs-Github_Page-blue?style=flat-square)](https://hcyvan.github.io/pattools/intro.html)
[![PyPI](https://img.shields.io/pypi/v/pattools_methy.svg?color=brightgreen)](https://pypi.org/project/pattools-methy/)
[![Python 3.10](https://img.shields.io/badge/python-3.10+-32cd32)](https://www.python.org/)
[![Development](https://img.shields.io/badge/development-rapid%20progress-orange)](https://github.com/hcyvan/pattools/graphs/commit-activity)
[![Developers](https://img.shields.io/github/contributors/hcyvan/pattools?color=green)](https://github.com/hcyvan/pattools/graphs/contributors)
[![Packagist License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

## Introduction

Pattools is an advanced toolkit designed for the comprehensive analysis of Bisulfite sequencing (BS-Seq) data.
This tool enables researchers to seamlessly evaluate a variety of DNA methylation metircs using BS-seq data
and to implement algorithms developed by other researchers without the need for data format conversion.


<img src="docs/image/framework.png" alt="pattools framework" width="500"/>

+ [Pattools Document Website](https://hcyvan.github.io/pattools/intro.html)
+ [Github](https://github.com/hcyvan/pattools)
+ [Open Issues](https://github.com/hcyvan/pattools/issues/)

## Development

1. Create a python 3.10+ environment, you can use conda, virtualenv, docker, etc. Or a local python environment (not
   recommended)
2. Install the project's dependent libraries

```
pip install -r requirements.txt
```

3. *Deprecated* Copy config.yaml.tpl to config.yaml and complete the configuration according to the local environment

## Unit test

```shell
python -m unittest discover -s tests
```

or

```shell
python setup.py test
```

## Installation

### PyPI

```
pip install pattools-methy
```

### Source code

Download the source code from the public repository [Github](https://github.com/hcyvan/pattools)

``` 
cd /path/to/source/code
python ./setup.py install
```

### Install from private repository (Deprecated)

#### Install from the master branch

```
pip install git+https://e.coding.net/gomicsgene2/lung_cancer_yixing/pattools.git
```

#### Install from a feature branch

For example: install from dev-xxx branch

```
pip install git+https://e.coding.net/gomicsgene2/lung_cancer_yixing/pattools.git@dev-xxx
```

## Usage (Deprecated, gradually moved to the docs directory)

### pattools SDK

### pattools toolkit

#### deconv

sun

```
pattools deconv -m sun -g hg38 \
    -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz \
    -p /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
    --include  T-cells B-cells \
    -o /mnt/d/project/pattools/tmp/out.csv
```

moss

```
pattools deconv -m moss -g hg38 \
    -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz \
    -p /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
    -o /mnt/d/project/pattools/tmp/out.csv
```

loyfer

```
pattools deconv -m loyfer -g hg38 -f Atlas.U25.l4.hg38.tsv\
    -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz \
    -p /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
    -o /mnt/d/project/pattools/tmp/out.csv
```

#### entropy

```
pattools entropy -d 3 -w 4\
    -i /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
    -o /mnt/d/project/pattools/tmp/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.entropy
```

#### region

This command is used to convert a region between different coordinate systems, such as
genome index coordinates or CpG index coordinates.

```
pattools region -t cpg2genome -i chr1:266762-266762 -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz
```

#### region-file

```
pattools region-file -t cpg2genome --column col2  --out-format bed -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz -i LUAD_1.0_0.txt -o LUAD_1.0_0.t.txt
```

## License

Distributed under the MIT License. See [LICENSE](LICENSE) for more information.

## Cite Us
Cheng et al. (2024). Novel diagnosis and therapeutic strategies for small cell lung cancer based on aberrant methylation in lung cancer subtypes. (journal not determined)

## Acknowledgements

+ [wgbstools](https://github.com/nloyfer/wgbs_tools)
+ [pysam](https://github.com/pysam-developers/pysam)
+ [samtools](http://www.htslib.org/)

