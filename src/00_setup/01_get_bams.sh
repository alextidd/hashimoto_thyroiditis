#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 01_get_bams -o log/%J_01_get_bams.out -e log/%J_01_get_bams.err 'bash src/00_setup/01_get_bams.sh'

# module
module load IRODS/1.0

# run (dna only for now)
(
  cd $LUSTRE_TEAM/resolveome/data/bams/
  nextflow run $NFS_TEAM/nextflow/nf-get_bam \
    --samplesheet samplesheet_irods.csv \
    --location irods \
    --out_dir ./ \
    --cram_to_bam \
    -resume \
    -w $LUSTRE_TEAM/resolveome/work/get_bams/ \
    -N at31@sanger.ac.uk
    # --merge_bams 
)
