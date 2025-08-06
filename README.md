# Hashimoto thyroiditis

This repository contains scripts, configurations, and reports for analyzing the Hashimoto thyroiditis data. The project includes workflows for analysing ResolveOME's PTA data and 10X data. Analyses include somatic mutation calling, genotyping, phylogenetic reconstruction, copy number calling, and cell type annotation.

## Repository Structure

```
.
├── bin                # utility scripts and tools
├── config             # configuration files
├── data               # input data
├── log                # log files from bsub jobs
├── out                # output files
├── reports            # generated reports and visualizations, summarising the analysis
└── src                # source code for analyses, organised by {technology}/{analysis}
    ├── 00_setup
    │   ├── 00_generate_samplesheet.R
    │   └── 01_get_bams.sh
    ├── resolveome
    │   ├── basejumper
    │   ├── nf-resolveome
    │   └── ptato
    └── trencadis-seq
        ├── nf-sc-geno
        ├── phase_1p_snps
        ├── seurat
        └── vartrix
```
