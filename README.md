# pattools

Pat format toolkit

![pattools framework](notebook/framework.png)


## Usage
copy config.yaml.tpl to config.yaml and complete the configuration according to the local environment
### pattools SDK
### pattools toolkit
#### deconv
sun
```
pattools deconv -m sun -g hg38 \
    -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz \
    -p /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
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
#### ratio
```
pattools ratio -d 3\
    -i /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
    -o /mnt/d/project/pattools/tmp/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.ratio
```
####
```
pattools pat2motif -i /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
                    -o /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.motif.gz
```
## Installation

### Source code

Download the source code and execute the following commands

``` 
cd /path/to/source/code
python ./setup.py install
```

### Install from repository

Note: This project will be synchronized to github in the future and can be installed through the github git url.

#### Install from the master branch

```
pip install git+https://e.coding.net/gomicsgene2/lung_cancer_yixing/pattools.git
```

#### Install from a feature branch

For example: install from dev-xxx branch

```shell
pip install git+https://e.coding.net/gomicsgene2/lung_cancer_yixing/pattools.git@dev-xxx
```
