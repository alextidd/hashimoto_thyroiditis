#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/trencadis-seq ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 04c_thyroid_nf-trencadis-seq_run -o log/%J_04c_thyroid_nf-trencadis-seq_run.out -e log/%J_04c_thyroid_nf-trencadis-seq_run.err 'bash src/04c_thyroid_nf-trencadis-seq_run.sh'

# dir
wd=$(pwd)

# r library
module load ISG/rocker/rver/4.4.0
export R_LIBS_USER=~/R-tmp-4.4

# run thyroid
(
  cd out/thyroid/nf-trencadis-seq/
  nextflow run ${wd}/../nextflow/nf-trencadis-seq \
    --samplesheet samplesheet.csv \
    --location irods \
    --out_dir ./ \
    --genes ${wd}/data/thyroid/driver_genes/driver_genes.txt \
    --gencode_gff3 ${wd}/../reference/gencode/GRCh38/gencode.v39.annotation.gff3.gz \
    --refcds ${wd}/../reference/dndscv/RefCDS_human_GRCh38_GencodeV18_recommended.rda \
    -with-tower \
    -w ${wd}/work/nf-trencadis-seq/thyroid/ \
    -N at31@sanger.ac.uk \
    -resume
)