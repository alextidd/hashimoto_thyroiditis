#!/bin/bash
# donor_id=PD63118 ; cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J basejumper_04a_bj-somatic-variantcalling_dna_run -o log/%J_basejumper_04a_bj-somatic-variantcalling_dna_run.out -e log/%J_basejumper_04a_bj-somatic-variantcalling_dna_run.err "bash src/basejumper/04a_bj-somatic-variantcalling_dna_run.sh $donor_id"

# parameters
donor_id=$1

# dirs
wd=$(pwd)

# modules
module load singularityce-4.1.0/python-3.11.6

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/basejumper/bj-somatic-variantcalling/sentieon_eval.lic
export LSB_EXCLUSIVE=Y

# run
(
  cd out/resolveome/basejumper/bj-somatic-variantcalling/dna/$donor_id/
  nextflow run $NFS_TEAM/nextflow/external/basejumper/bj-somatic-variantcalling \
    --input_csv samplesheet.csv \
    --publish_dir $donor_id \
    --timestamp run \
    --variant_workflow_type somatic_heuristic_filter \
    --dnascope_model_selection bioskryb129 \
    --skip_variant_annotation false \
    --skip_sigprofile false \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $LUSTRE_TEAM/projects/hashimoto_thyroiditis/work/basejumper/bj-somatic-variantcalling/dna/$donor_id/ \
    -profile singularity \
    --architecture "x86_64" \
    -N at31@sanger.ac.uk \
    -resume
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03d_basejumper_somatic-variantcalling_dna_symlink -o "log/%J_03d_basejumper_somatic-variantcalling_dna_symlink.out" "source ~/.bashrc && replace_symlinks out/resolveome/basejumper/bj-somatic-variantcalling/dna/"