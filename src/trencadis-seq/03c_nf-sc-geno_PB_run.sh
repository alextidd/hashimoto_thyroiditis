#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03b_nf-sc-geno_run -o log/%J_03b_nf-sc-geno_run.out -e log/%J_03b_nf-sc-geno_run.err 'bash src/trencadis-seq/03b_nf-sc-geno_run.sh'

# run
(
  cd out/trencadis-seq/nf-sc-geno/PB/
  nextflow run $NFS_TEAM/nextflow/nf-sc-geno \
    --samplesheet samplesheet.csv \
    --location irods \
    --out_dir ./ \
    --no_cell_barcodes \
    -with-tower \
    -w $LUSTRE_TEAM/resolveome/work/nf-sc-geno/PB/ \
    -N at31@sanger.ac.uk \
    -resume
)