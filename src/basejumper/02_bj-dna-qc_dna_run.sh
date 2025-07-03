#!/bin/bash
# donor_id=PD63118 ; cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 02_bj-dna-qc_dna_run -o log/%J_02_bj-dna-qc_dna_run.out -e log/%J_02_bj-dna-qc_dna_run.err "bash src/basejumper/02_bj-dna-qc_dna_run.sh ${donor_id}"

# parameters
donor_id=$1

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/BaseJumper/bj-dna-qc/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd out/BaseJumper/bj-dna-qc/dna/$donor_id/
  nextflow run $NFS_TEAM/nextflow/external/BaseJumper/bj-dna-qc \
    --input_csv samplesheet.csv \
    --publish_dir $donor_id \
    --timestamp run \
    --skip_ginkgo false \
    -c $wd/config/basejumper.config \
    -w $LUSTRE_TEAM/resolveome/work/BaseJumper/bj-dna-qc/dna/$donor_id/ \
    -profile singularity \
    -N at31@sanger.ac.uk \
    -resume 
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03b_BaseJumper_dna-qc_dna_symlink -o "log/%J_03b_BaseJumper_dna-qc_dna_symlink.out" ". ~/.bashrc ; replace_symlinks out/BaseJumper/bj-dna-qc/dna/"