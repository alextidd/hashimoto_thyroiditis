#!/bin/bash

# modules
module load bcftools-1.9/python-3.11.6
module load samtools-1.19.2
module load IRODS/1.0

# dirs
data_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/data/vartrix/
ref_dir=$data_dir/reference/
snps_dir=$data_dir/snps/
tx_dir=$data_dir/TX/
pb_dir=$data_dir/PB/
out_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/out/vartrix/
mkdir -p $data_dir $ref_dir $snps_dir $tx_dir $pb_dir $out_dir

# files
vcf_nfs=/nfs/irods-cgp-sr12-sdc/intproj/3464/sample/PD63118b_lo0044/PD63118b_lo0044.v1.caveman_c.snps.vcf.gz
vcf=$snps_dir/PD63118b_lo0044.v1.caveman_c.snps
snps_pos=out/resolveome/nf-resolveome/muts_and_snps/PD63118/caveman_snps_positions.tsv
fasta=/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa

# stage and index the vcf (grch37)
cp $vcf_nfs $vcf.vcf.gz
bcftools index $vcf.vcf.gz

# subset to chr1p
bcftools view -r 1:1-121500000 $vcf.vcf.gz -Oz > $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz
rm $vcf.vcf.gz*

# subset to SNPs with DP>50 and 0.3<VAF<0.7
bcftools view -R $snps_pos $vcf.1p.vcf.gz -Oz > $vcf.1p.subset.vcf.gz
mv $vcf.1p.subset.vcf.gz $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz

# add chr prefix
echo -e "1\tchr1" > $snps_dir/chr_rename.txt
bcftools annotate \
  --rename-chrs $snps_dir/chr_rename.txt $vcf.1p.vcf.gz -Oz \
  -o $vcf.1p.renamed.vcf.gz
mv $vcf.1p.renamed.vcf.gz $vcf.1p.vcf.gz
bcftools index $vcf.1p.vcf.gz
rm $snps_dir/chr_rename.txt

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
rm $vcf*grch38*

# remove chr prefix
echo -e "chr1\t1" > $snps_dir/chr_rename.txt
bcftools annotate \
  --rename-chrs $snps_dir/chr_rename.txt $vcf.1p.vcf.gz -Oz \
  -o $vcf.1p_nochr.vcf.gz
bcftools index $vcf.1p_nochr.vcf.gz
rm $snps_dir/chr_rename.txt

# convert fasta to remove chr
sed 's/^>chr/>/' $fasta > $ref_dir/genome_nochr.fa
samtools faidx $ref_dir/genome_nochr.fa

# 10X: get barcodes and bams
tx_dirs=(cellranger720_count_49200_7613STDY14897605_GRCh38-3_0_0 cellranger720_count_49200_7613STDY14897606_GRCh38-3_0_0)
sample_ids=(PD63118b_st0001_rna PD63118b_st0002_rna)
for dir in "${tx_dirs[@]}"; do

  echo $dir

  # directories
  irods_dir=/seq/illumina/runs/49/49200/cellranger/$dir
  bc_dir=$tx_dir/$dir/filtered_feature_bc_matrix/

  # get barcodes
  mkdir -p $bc_dir
  iget -K \
    /seq/illumina/runs/49/49200/cellranger/$dir/filtered_feature_bc_matrix/barcodes.tsv.gz \
    $bc_dir
  
  # stage and index bam (grch38)
  bam=/seq/illumina/runs/49/49200/cellranger/$dir/possorted_genome_bam.bam
  iget -K $bam $tx_dir/$dir/
  iget -K $bam.bai $tx_dir/$dir/

done

# PacBio: stage and index bam (grch38)
pb_bam=/seq/pacbio/r84093_20241025_154925/1_D01/21999/21999.m84093_241025_215302_s4.scisoseq.mapped.bam
iget -K $pb_bam $pb_dir
iget -K $pb_bam.bai $pb_dir
iget -K \
  /seq/pacbio/r84093_20241025_154925/1_D01/21999/21999.m84093_241025_215302_s4.scisoseq.seurat_info.tar.gz \
  $pb_dir
tar -xvf $pb_dir/21999.m84093_241025_215302_s4.scisoseq.seurat_info.tar.gz -C $pb_dir
