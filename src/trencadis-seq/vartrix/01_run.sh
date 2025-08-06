#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J trencadis-seq_vartrix_01_run -o log/%J_trencadis-seq_vartrix_01_run.out -e log/%J_trencadis-seq_vartrix_01_run.err 'bash src/trencadis-seq/vartrix/01_run.sh'

# modules
module load bcftools-1.9/python-3.11.6
module load samtools-1.19.2

# dirs
data_dir=$LUSTRE_TEAM/projects/hashimoto_thyroiditis/data/vartrix/
ref_dir=$data_dir/reference/
snps_dir=$data_dir/snps/
tx_dir=$data_dir/TX/
pb_dir=$data_dir/PB/
out_dir=$LUSTRE_TEAM/projects/hashimoto_thyroiditis/out/vartrix/

# files
vcf=$snps_dir/PD63118b_lo0044.v1.caveman_c.snps
fasta=/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa
fasta_nochr=$ref_dir/genome_nochr.fa

echo "run vartrix on the 10X"
sample_ids=(PD63118b_st0001 PD63118b_st0002)
tx_dirs=(cellranger720_count_49200_7613STDY14897605_GRCh38-3_0_0 cellranger720_count_49200_7613STDY14897606_GRCh38-3_0_0)
for i in ${!sample_ids[@]} ; do

  echo "sample_id: ${sample_ids[$i]}"
  echo "tx_dir: ${tx_dirs[$i]}"

  # create dir
  out_dir_i=$out_dir/TX/filtered/${sample_ids[$i]}/
  mkdir -p $out_dir_i
  echo "vartrix out dir: $out_dir_i"

  # stage barcodes
  zcat $tx_dir/${tx_dirs[$i]}/filtered_feature_bc_matrix/barcodes.tsv.gz \
    > $out_dir_i/barcodes.tsv

  # run vartrix
  vartrix_linux \
    --vcf $vcf.1p_nochr.vcf.gz \
    --bam $tx_dir/${tx_dirs[$i]}/possorted_genome_bam.bam \
    --fasta $fasta_nochr \
    --cell-barcodes $out_dir_i/barcodes.tsv \
    --out-matrix $out_dir_i/alt.mtx \
    --ref-matrix $out_dir_i/ref.mtx \
    --out-variants $out_dir_i/variants.txt \
    --scoring-method coverage \
    --threads 8 \
    --umi
  
done

echo "run vartrix on the PacBio"

# create dir
out_dir_i=$out_dir/PB/filtered/PD63118b/
mkdir -p $out_dir_i

# stage barcodes
cp $pb_dir/genes_seurat/barcodes.tsv $out_dir_i/barcodes.tsv

# run vartrix
vartrix_linux \
  --vcf $vcf.1p.vcf.gz \
  --bam $pb_dir/21999.m84093_241025_215302_s4.scisoseq.mapped.bam \
  --fasta $fasta \
  --cell-barcodes $out_dir_i/barcodes.tsv \
  --scoring-method coverage \
  --out-matrix $out_dir_i/alt.mtx \
  --ref-matrix $out_dir_i/ref.mtx \
  --out-variants $out_dir_i/variants.txt \
  --threads 8 \
  --umi