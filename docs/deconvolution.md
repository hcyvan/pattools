# Deconvolution

Pattools facilitates multiple deconvolution algorithms tailored to methylation data analysis.
The deconvolution algorithms implemented by `pattools deconv` include:

- `sun` {cite}`sun_plasma_2015` ： based on the methylation densities (i.e., the percentage of CpGs being 
methylated within a 500-bp unit) of the methylation biomarker
- `moss` {cite}`moss_comprehensive_2018` : based on the beta values of the Illumina Infinium Human Methylation 450K or EPIC BeadChip arrays
- `loyfer` {cite}`Loyfer2024.05.08.593132` : based on the uxm ratios of the DNA regions

## usage
take `sun` for example

```bash
# View the cell components for each method
pattools deconv-helper -m sun 

# run 
pattools deconv -m sun -g hg38 \
    -c /mnt/d/project/wgbs_tools/references/hg38/CpG.bed.gz \
    -p /mnt/d/data/cacLung/raw/pat/GSM5652289_Blood-T-Eff-CD8-Z0000041F.hg38.pat.gz \
    -o /mnt/d/project/pattools/tmp/out.csv \
    --include T-cells B-cells Neutrophils \ 
    --ignore Liver Lung 
```

