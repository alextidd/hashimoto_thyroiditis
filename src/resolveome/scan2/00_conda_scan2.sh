#!/bin/bash

# modules
module load ISG/conda

# build conda env
conda create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2