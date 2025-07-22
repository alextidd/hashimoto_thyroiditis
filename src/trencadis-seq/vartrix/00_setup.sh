#!/bin/bash

# modules
module load bcftools-1.9/python-3.11.6
module load samtools-1.19.2

# files
data_dir=$LUSTRE_TEAM/projects/hashimoto_thyroiditis/data/vartrix
mkdir -p $data_dir
vcf_nfs=/nfs/irods-cgp-sr12-sdc/intproj/3464/sample/PD63118b_lo0044/PD63118b_lo0044.v1.caveman_c.snps.vcf.gz
vcf=$data_dir/PD63118b_lo0044.v1.caveman_c.snps
snps=out/resolveome/nf-resolveome/muts_and_snps/PD63118/caveman_snps_positions.tsv
fasta=/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa

# stage and index the vcf (grch37)
cp $vcf_nfs $vcf.vcf.gz
bcftools index $vcf.vcf.gz

# subset to chr1p
bcftools view -r 1:1-121500000 $vcf.vcf.gz -Oz > $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz
rm $vcf.vcf.gz

# subset to SNPs with DP>50 and 0.3<VAF<0.7
bcftools view -R $snps $vcf.1p.vcf.gz -Oz > $vcf.tmp.vcf.gz
mv $vcf.tmp.vcf.gz $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz

# add chr prefix
echo -e "1\tchr1" > $data_dir/chr_rename.txt
bcftools annotate \
  --rename-chrs $data_dir/chr_rename.txt $vcf.1p.vcf.gz -Oz \
  -o $vcf.1p_renamed.vcf.gz
mv $vcf.1p_renamed.vcf.gz $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz

# liftover to grch38
java -Xmx4g -jar /software/CASM/modules/installs/picard-tools/picard-tools-3.1.0/picard.jar \
  LiftoverVcf \
  I=$vcf.1p.vcf.gz \
  O=$vcf.1p_grch38.vcf \
  CHAIN=../../reference/liftover/hg19ToHg38.over.chain \
  REJECT=$vcf.1p_rejected_grch38.vcf \
  R=$fasta
bgzip -c $vcf.1p_grch38.vcf > $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz
rm $data_dir/*grch38.vcf*

# remove chr prefix
echo -e "chr1\t1" > $data_dir/chr_rename.txt
bcftools annotate \
  --rename-chrs $data_dir/chr_rename.txt $vcf.1p.vcf.gz -Oz \
  -o $vcf.1p_nochr.vcf.gz
bcftools index $vcf.1p_nochr.vcf.gz

# convert fasta to remove chr
sed 's/^>chr/>/' $fasta > $data_dir/genome_nochr.fa
samtools faidx $data_dir/genome_nochr.fa

# 10X: get barcodes
iget -k /seq/illumina/runs/49/49200/cellranger/cellranger720_count_49200_7613STDY14897605_GRCh38-3_0_0/raw_feature_bc_matrix/barcodes.tsv.gz $data_dir

# 10X: stage and index bam (grch38)
tx_bam=/seq/illumina/runs/49/49200/cellranger/cellranger720_count_49200_7613STDY14897605_GRCh38-3_0_0/possorted_genome_bam.bam
iget -K $tx_bam $data_dir
iget -K $tx_bam.bai $data_dir

# PacBio: stage and index bam (grch38)
pb_bam=/seq/pacbio/r84093_20241025_154925/1_D01/21999/21999.m84093_241025_215302_s4.scisoseq.mapped.bam
iget -K $pb_bam $data_dir
iget -K $pb_bam.bai $data_dir
iget -K /seq/pacbio/r84093_20241025_154925/1_D01/21999/21999.m84093_241025_215302_s4.scisoseq.seurat_info.tar.gz $data_dir
tar -xvf $data_dir/21999.m84093_241025_215302_s4.scisoseq.seurat_info.tar.gz -C $data_dir
