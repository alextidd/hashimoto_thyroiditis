#!/bin/bash

# modules
module load ISG/conda
module load samtools-1.19.2/python-3.11.6

# create new conda environment with all dependencies
conda create -n hiscanner \
  -c conda-forge -c bioconda \
  python=3.8 samtools=1.15.1 bcftools=1.13 r-base "r-mgcv>=1.8"
conda activate hiscanner
pip install hiscanner --no-cache-dir

# split fasta
mkdir -p data/hiscanner/GRCh37/split
awk '{print $1}' data/scan2/GRCh37/genome.fa.fai |
xargs -I {} sh -c 'samtools faidx data/scan2/GRCh37/genome.fa {} > data/hiscanner/GRCh37/split/{}.fasta'