#!/bin/bash
# runsub src/00_setup/01_get_bams.sh

# module
module load IRODS/1.0

# tmp: plate 384.3
(
  cd $LUSTRE_125/projects/hashimoto_thyroiditis/data/bams/
  head -1 samplesheet_irods.csv > samplesheet_irods.tmp.csv
  cat samplesheet_irods.csv |
  grep "run51611" >> samplesheet_irods.tmp.csv
)

# run (dna / dnahyb)
(
  cd $LUSTRE_125/projects/hashimoto_thyroiditis/data/bams/
  nextflow run $NFS_TEAM/nextflow/nf-get_bam \
    --samplesheet samplesheet_irods.tmp.csv \
    --fasta $REF_PATH/Homo_sapiens/GRCh37d5/genome.fa \
    --location irods \
    --out_dir ./ \
    -resume \
    -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/get_bams/ \
    -N at31@sanger.ac.uk \
    --cram_to_bam
)

# # run (rna)
# (
#   cd $LUSTRE_125/projects/hashimoto_thyroiditis/data/bams/
#   nextflow run $NFS_TEAM/nextflow/nf-get_bam \
#     --samplesheet samplesheet_irods_rna.csv \
#     --fasta $REF_PATH/Homo_sapiens/GRCh37d5/genome.fa \
#     --location irods \
#     --out_dir ./ \
#     -resume \
#     -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/get_bams/ \
#     -N at31@sanger.ac.uk \
#     --cram_to_bam \
#     --merge_bams
# )