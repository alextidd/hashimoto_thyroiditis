#!/bin/bash

# modules
module load ISG/conda

# create env with compatible python and snakemake
conda create -n scan2-py311 python=3.11 -c conda-forge -c bioconda -y
conda activate scan2-py311
conda install -c conda-forge -c bioconda snakemake=9.13.4 -y

# install R + r-scan2 dependencies in same env (conda R or system R; conda recommended)
conda install -c conda-forge r-base r-devtools r-reticulate r-data.table r-dplyr -y

# add gatk3
conda install -c conda-forge -c bioconda -c jluquette -c dranew -c soil gatk=3.8