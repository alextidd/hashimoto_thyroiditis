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
cat data/resolveome/manual_inspection/pta_additional_annotation_H1.tsv |
awk -F'\t' '$10=="not lymphocyte" {print $2}' |
awk -F, 'NR==FNR{a[$1]; next} FNR==1{next} $5=="dna" && $2 in a{print $NF}' - $ss |
head -3 \
> $bulk_dir/bams.txt

# merge bams
samtools merge -@ 8 -b $bulk_dir/bams.txt -o $bulk_dir/merged_wgs_non_lymphocytes.bam 
samtools index $bulk_dir/merged_wgs_non_lymphocytes.bam