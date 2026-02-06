#!/bin/bash
# runsub src/resolveome/scan2/00_merge_non_lymphocytes.sh -M 40000 -n 8

# modules
module load samtools-1.19.2/python-3.11.6

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
ss=$lustre_dir/data/bams/samplesheet_local.csv
bulk_dir=$lustre_dir/data/bams/merged_non_lymphocytes/
mkdir -p $bulk_dir

# 3 WGS bams from non-lymphocyte cells
# check that they have dna, get bam paths
cat data/resolveome/manual_inspection/H1_PD63118_pta_additional_annotation.tsv |
  awk -F'\t' '$10=="not lymphocyte" {print $2}' |
  awk -F, 'NR==FNR{a[$1]; next} FNR==1{next} $6=="dna" && $3 in a{print $NF}' - $ss |
  head -3 \
> $bulk_dir/bams.txt

# merge bams
samtools merge -@ 8 -b $bulk_dir/bams.txt -o $bulk_dir/merged_wgs_non_lymphocytes.tmp.bam 

# update SM tag to merged_wgs_non_lymphocytes
samtools view -H $bulk_dir/merged_wgs_non_lymphocytes.tmp.bam |
  sed 's/SM:[^\t]*/SM:merged_wgs_non_lymphocytes/g' |
  samtools reheader - $bulk_dir/merged_wgs_non_lymphocytes.tmp.bam > $bulk_dir/merged_wgs_non_lymphocytes.bam
rm $bulk_dir/merged_wgs_non_lymphocytes.tmp.bam

samtools index $bulk_dir/merged_wgs_non_lymphocytes.bam