#!/bin/bash
# donor_id=PD63118 ; cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03f_basejumper_somatic-variantcalling-develop_dna_${donor_id}_run -o log/%J_03f_basejumper_somatic-variantcalling-develop_dna_${donor_id}_run.out -e log/%J_03f_basejumper_somatic-variantcalling-develop_dna_${donor_id}_run.err "bash src/03f_basejumper_somatic-variantcalling-develop_dna_run.sh ${donor_id}"

# parameters
donor_id=$1

# dirs
wd=$(pwd)

# modules
module load singularity

# run
(
  cd out/resolveome/basejumper/bj-somatic-variantcalling-develop/dna/$donor_id/
  nextflow run $wd/../nextflow/external/basejumper/bj-somatic-variantcalling-parabricks-develop \
    --input_csv samplesheet.csv \
    --publish_dir $donor_id \
    --timestamp run \
    -c $wd/config/bj-somatic-variantcalling.config \
    -c $wd/config/basejumper.config \
    -w $LUSTRE_TEAM/projects/hashimoto_thyroiditis/work/basejumper/bj-somatic-variantcalling/dna/$donor_id/ \
    -profile singularity \
    --architecture "x86_64" \
    -resume \
    -N at31@sanger.ac.uk
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03d_basejumper_somatic-variantcalling_dna_symlink -o "log/%J_03d_basejumper_somatic-variantcalling_dna_symlink.out" "source ~/.bashrc && replace_symlinks out/resolveome/basejumper/bj-somatic-variantcalling/dna/"