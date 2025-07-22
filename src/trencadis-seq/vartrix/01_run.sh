#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M20000 -n 8 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J 02_run_vartrix -o log/%J_02_run_vartrix.out -e log/%J_02_run_vartrix.err 'bash src/trencadis-seq/02_run_vartrix.sh'

# modules
module load bcftools-1.9/python-3.11.6
module load samtools-1.19.2

# files
data_dir=$LUSTRE_TEAM/projects/hashimoto_thyroiditis/data/vartrix
mkdir -p $data_dir
vcf=$data_dir/PD63118b_lo0044.v1.caveman_c.snps
fasta=/lustre/scratch124/casm/references/pipeline_ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa

echo "10X: run vartrix"
mkdir -p $data_dir/out/TX/
vartrix_linux \
  --vcf $vcf.1p_nochr.vcf.gz \
  --bam $data_dir/possorted_genome_bam.bam \
  --fasta $data_dir/genome_nochr.fa \
  --cell-barcodes $data_dir/barcodes.tsv.gz \
  --out-matrix $data_dir/out/TX/alt.mtx \
  --ref-matrix $data_dir/out/TX/ref.mtx \
  --out-variants $data_dir/out/TX/variants.txt \
  --scoring-method coverage \
  --threads 8 \
  --umi

echo "PacBio: run vartrix"
mkdir -p $data_dir/out/PB/
vartrix_linux \
  --vcf $vcf.1p.vcf.gz \
  --bam $data_dir/21999.m84093_241025_215302_s4.scisoseq.mapped.bam \
  --fasta $fasta \
  --cell-barcodes $data_dir/genes_seurat/barcodes.tsv  \
  --scoring-method coverage \
  --out-matrix $data_dir/out/PB/alt.mtx \
  --ref-matrix $data_dir/out/PB/ref.mtx \
  --out-variants $data_dir/out/PB/variants.txt \
  --threads 8 \
  --umi

# # move results to nfs
# mv --force $data_dir/out/* out/trencadis-seq/vartrix/