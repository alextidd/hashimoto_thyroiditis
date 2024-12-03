#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J driver_coverage -o log/driver_coverage_%J.out -e log/driver_coverage_%J.err 'bash src/01b_run_nf-driver_coverage.sh'

# mapped to grch37!

# r library
module load ISG/rocker/rver/4.4.0
export R_LIBS_USER=~/R-tmp-4.4

nextflow run ../nf-driver_coverage \
  --samplesheet out/driver_coverage/samplesheet.csv \
  --bulk \
  --no_chr \
  --genes out/driver_coverage/genes.txt \
  --out_dir out/driver_coverage/ \
  --gencode_gff3 ../reference/gencode/GRCh37/gencode.v39lift37.annotation.no_chr.gff3.gz \
  --refcds ../reference/dndscv/RefCDS_human_hg19_GencodeV18_newcovariates.rda \
  -w work/driver_coverage/ \
  -with-tower \
  -N at31@sanger.ac.uk \
  -resume
