#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J trencadis-seq_nf-sc-geno_02_PB_run -o log/%J_trencadis-seq_nf-sc-geno_02_PB_run.out -e log/%J_trencadis-seq_nf-sc-geno_02_PB_run.err 'bash src/trencadis-seq/nf-sc-geno/02_PB_run.sh'

# run
(
  cd out/trencadis-seq/nf-sc-geno/PB/
  nextflow run $NFS_TEAM/nextflow/nf-sc-geno \
    --samplesheet samplesheet.csv \
    --location irods \
    --out_dir ./ \
    --no_cell_barcodes \
    --fasta $REF_PATH/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    -with-tower \
    -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/nf-sc-geno/PB/ \
    -N at31@sanger.ac.uk \
    -resume
)