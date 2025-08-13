#!/bin/bash
# donor_id=PD63118
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 03c_basejumper_expression_rna_run -o log/%J_03c_basejumper_expression_rna_run.out -e log/%J_03c_basejumper_expression_rna_run.err "bash src/03c_basejumper_expression_rna_run.sh ${donor_id}"

# parameters
donor_id=$1

# dirs
wd=$(pwd)

# modules
module load singularity

# sentieon license
export SENTIEON_LICENSE=$wd/../nextflow/external/basejumper/bj-expression/sentieon_eval.lic 
export LSB_EXCLUSIVE=Y

# run bj-expression
(
  cd out/resolveome/basejumper/bj-expression/rna/$donor_id/
  nextflow run $wd/../nextflow/external/basejumper/bj-expression \
    --input_csv samplesheet.csv \
    --publish_dir $donor_id \
    --timestamp run \
    --genome GRCh38 \
    --skip_subsampling \
    --tmp_dir $TMPDIR \
    -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/basejumper/bj-expression/rna/$donor_id/ \
    -c $wd/config/basejumper.config \
    -c $wd/config/bj-expression.config \
    -profile singularity \
    -resume
)

# bsub -q basement -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J 03c_basejumper_expression_rna_symlink -o "log/%J_03c_basejumper_expression_rna_symlink.out" ". ~/.bashrc ; replace_symlinks out/resolveome/basejumper/bj-expression/rna/PD63118_run/"