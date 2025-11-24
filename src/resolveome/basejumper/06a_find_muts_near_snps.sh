#!/bin/bash

# modules
module load bedtools2-2.29.0/python-3.10.10

# grch37
snps=out/resolveome/nf-resolveome/muts_and_snps/PD63118/caveman_snps.tsv

# grch38
muts=out/resolveome/sequoia/Patient_both_NV_tree_all.txt

# dirs
out_dir=out/resolveome/basejumper/find_muts_near_snps/
mkdir -p $out_dir

# muts bed, liftover to grch37
cat $muts | sed 1d | cut -f1 -d" " | sed 's/\"//g' | sed 's/_/\t/g' |
awk '{print $1"\t"($2-1)"\t"$2"\t"$3"\t"$4}' \
> $out_dir/muts_grch38.bed
liftOver \
  $out_dir/muts_grch38.bed \
  ../../reference/liftover/hg38ToHg19.over.chain \
  $out_dir/muts_grch37.bed \
  $out_dir/muts_unmapped.bed

# snps bed
cat $snps | sed 1d |
awk '{print "chr"$1"\t"($2-1)"\t"$2"\t"$3"\t"$4}' > $out_dir/snps.bed

# add 150 bp window around snps
cat $out_dir/snps.bed |
awk '{print $1"\t"($2-150)"\t"($3+150)"\t"$4"\t"$5"\t"$3}' > $out_dir/snps_150bp.bed

# intersect muts and 150bp snp windows
echo -e "chr\tmut_pos\tmut_ref\tmut_alt\tsnp_ref\tsnp_alt\tsnp_pos" > $out_dir/muts_near_snps.tsv
bedtools intersect -a $out_dir/muts_grch37.bed -b $out_dir/snps_150bp.bed -wa -wb |
cut -f1,3,4,5,9,10,11 >> $out_dir/muts_near_snps.tsv