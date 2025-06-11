#!/bin/bash
# donor_id=PD63118 ; cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03g_BaseJumper_somatic-variantcalling-develop_dnahyb_${donor_id}_run -o log/%J_03g_BaseJumper_somatic-variantcalling-develop_dnahyb_${donor_id}_run.out -e log/%J_03g_BaseJumper_somatic-variantcalling-develop_dnahyb_${donor_id}_run.err "bash src/03g_BaseJumper_somatic-variantcalling-develop_dnahyb_run.sh ${donor_id}"

# parameters
donor_id=$1

# dirs
wd=$(pwd)

# modules
module load singularity

# run
(
  cd out/BaseJumper/bj-somatic-variantcalling-develop/dnahyb/$donor_id/
  nextflow run $wd/../nextflow/external/BaseJumper/bj-somatic-variantcalling-parabricks-develop \
    --input_csv samplesheet.csv \
    --publish_dir $donor_id \
    --timestamp run \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $LUSTRE_TEAM/resolveome/work/BaseJumper/bj-somatic-variantcalling/dnahyb/$donor_id/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume \
    -N at31@sanger.ac.uk
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03g_BaseJumper_somatic-variantcalling-develop_dnahyb_symlink -o "log/%J_03g_BaseJumper_somatic-variantcalling-develop_dnahyb_symlink.out" "source ~/.bashrc && replace_symlinks out/BaseJumper/bj-somatic-variantcalling/dna/"