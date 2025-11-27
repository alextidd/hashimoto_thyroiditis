#!/bin/bash

# modules
module load bedtools2-2.29.0/python-3.10.10

# dirs
ref_dir=data/scan2/GRCh37/
lustre_ref_dir=/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh37d5/
bulk_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/data/bams/matched_normal/
windows_dir=out/resolveome/scan2/windows/
mkdir -p $ref_dir $bulk_dir $windows_dir

# get fasta with corrected dict filename for gatk
ln -s $lustre_ref_dir/genome.fa $ref_dir/genome.fa
ln -s $lustre_ref_dir/genome.fa.fai $ref_dir/genome.fa.fai
ln -s $lustre_ref_dir/genome.fa.dict $ref_dir/genome.dict

# get bulk
cp \
  /nfs/cancer_ref01/nst_links/live/3464/PD63118b_lo0001/PD63118b_lo0001.sample.dupmarked.bam* \
  $bulk_dir/

# make genome windows of 100kb (autosomes only, no chr prefix)
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
  "select chrom, size from hg19.chromInfo" \
  > $windows_dir/hg19.genome
bedtools makewindows -g $windows_dir/hg19.genome -w 100000 |
sed 's/^chr//' |
sort -V -k1,1 -k2,2n |
awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/' \
> $windows_dir/windows.bed