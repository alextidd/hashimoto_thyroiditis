#!/bin/bash

# modules
module load ISG/conda

# build conda env
conda create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2

# download develop r-scan2
R -e 'devtools::install_github("parklab/r-scan2")'