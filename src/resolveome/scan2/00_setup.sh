#!/bin/bash

# get fasta with corrected dict filename for gatk
ref_dir=data/scan2/GRCh37/
lustre_ref_dir=/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh37d5/
mkdir -p $ref_dir
ln -s $lustre_ref_dir/genome.fa $ref_dir/genome.fa
ln -s $lustre_ref_dir/genome.fa.fai $ref_dir/genome.fa.fai
ln -s $lustre_ref_dir/genome.fa.dict $ref_dir/genome.dict

# get bulk
bulk_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/data/bams/matched_normal/
mkdir $bulk_dir/
cp \
  /nfs/cancer_ref01/nst_links/live/3464/PD63118b_lo0001/PD63118b_lo0001.sample.dupmarked.bam* \
  $bulk_dir/