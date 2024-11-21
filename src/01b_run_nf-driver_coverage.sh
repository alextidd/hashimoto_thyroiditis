#!/bin/bash

# mapped to grch37
project_dir=/nfs/cancer_ref01/nst_links/live/3464/

nextflow run ../nf-driver_coverage \
  --samplesheet data/metadata/samplesheet.csv \
  --drivers data/driver_genes/driver_genes.txt \
  --out_dir out/driver_coverage/ \
  --gencode_gff3 data/gencode/gencode.v39.annotation.gff3.gz \
  --refcds data/dndscv/RefCDS_human_GRCh38_GencodeV18_recommended.rda \
  --validate_params false \
  -with-tower \
  -N at31@sanger.ac.uk \
  -w work/driver_coverage/ \
  -resume