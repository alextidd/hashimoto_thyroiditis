#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -n 8 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 02_run_vartrix -o log/%J_02_run_vartrix.out -e log/%J_02_run_vartrix.err 'bash src/trencadis-seq/02_run_vartrix.sh'

# # modules
# module load bcftools-1.9/python-3.11.6
# module load samtools-1.19.2

# files
data_dir=$LUSTRE_TEAM/resolveome/data/vartrix
mkdir -p $data_dir
vcf_nfs=/nfs/irods-cgp-sr12-sdc/intproj/3464/sample/PD63118b_lo0044/PD63118b_lo0044.v1.caveman_c.snps.vcf.gz
vcf=$data_dir/PD63118b_lo0044.v1.caveman_c.snps
snps=out/nf-resolveome/muts_and_snps/PD63118/caveman_snps_positions.tsv
fasta=/lustre/scratch124/casm/references/ref_tmp/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa

# # stage and index the vcf (grch37)
# cp $vcf_nfs $vcf.vcf.gz
# bcftools index $vcf.vcf.gz

# # subset to chr1p
# bcftools view -r 1:1-121500000 $vcf.vcf.gz -Oz > $vcf.1p.vcf.gz
# bcftools index $vcf.1p.vcf.gz

# # subset to SNPs with DP>50 and 0.3<VAF<0.7
# bcftools view -R $snps $vcf.1p.vcf.gz -Oz > $vcf.tmp.vcf.gz
# mv $vcf.tmp.vcf.gz $vcf.1p.vcf.gz
# bcftools index $vcf.1p.vcf.gz

# # add chr prefix
# echo -e "1\tchr1" > $data_dir/chr_rename.txt
# bcftools annotate \
#   --rename-chrs $data_dir/chr_rename.txt $vcf.1p.vcf.gz -Oz \
#   -o $vcf.1p_renamed.vcf.gz
# mv $vcf.1p_renamed.vcf.gz $vcf.1p.vcf.gz
# bcftools index $vcf.1p.vcf.gz

# # liftover to grch38
# java -Xmx4g -jar /software/CASM/modules/installs/picard-tools/picard-tools-3.1.0/picard.jar \
#   LiftoverVcf \
#   I=$vcf.1p.vcf.gz \
#   O=$vcf.1p_grch38.vcf \
#   CHAIN=../../reference/liftover/hg19ToHg38.over.chain \
#   REJECT=$vcf.1p_grch38_rejected.vcf \
#   R=$fasta
# bgzip -c $vcf.1p_grch38.vcf > $vcf.1p_grch38.vcf.gz
# mv $vcf.1p_grch38.vcf.gz $vcf.1p.vcf.gz
# bcftools index $vcf.1p.vcf.gz

# # remove chr prefix
# echo -e "chr1\t1" > $data_dir/chr_rename.txt
# bcftools annotate \
#   --rename-chrs $data_dir/chr_rename.txt $vcf.1p.vcf.gz -Oz \
#   -o $vcf.1p_renamed.vcf.gz
# mv $vcf.1p_renamed.vcf.gz $vcf.1p.vcf.gz
# bcftools index $vcf.1p.vcf.gz

# # convert fasta to remove chr
# sed 's/^>chr/>/' $fasta > $data_dir/genome.fa
# samtools faidx $data_dir/genome.fa

# # get cell barcodes
# iget -K \
#   /seq/illumina/runs/49/49200/cellranger/cellranger720_count_49200_7613STDY14897605_GRCh38-3_0_0/filtered_feature_bc_matrix/barcodes.tsv.gz \
#   $data_dir/cell_barcodes.tsv.gz

# # stage and index bam (grch3738)
# tx_bam=/seq/illumina/runs/49/49200/cellranger/cellranger720_count_49200_7613STDY14897605_GRCh38-3_0_0/possorted_genome_bam.bam
# iget -K $tx_bam $data_dir
# iget -K $tx_bam.bai $data_dir

# run vartrix
mkdir $data_dir/out/
vartrix_linux \
  --vcf $vcf.1p.vcf.gz \
  --bam $data_dir/possorted_genome_bam.bam \
  --fasta $data_dir/genome.fa \
  --cell-barcodes out/trencadis-seq/seurat/min_3_cells_min_500_genes/annotated_cell_barcodes.txt \
  --scoring-method coverage \
  --out-matrix $data_dir/out/alt.mtx \
  --ref-matrix $data_dir/out/ref.mtx \
  --out-variants $data_dir/out/variants.txt \
  --threads 8

# move results to nfs
mv $data_dir/out/* out/trencadis-seq/vartrix/