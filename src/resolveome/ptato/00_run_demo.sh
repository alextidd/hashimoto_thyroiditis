#!/bin/bash
# runsub src/resolveome/ptato/00_run_demo.sh

# modules
module load singularity

# dirs
wd=$(pwd)
ptato_dir=$NFS_TEAM/nextflow/external/PTATO/
out_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/out/resolveome/ptato/demo/
work_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/work/resolveome/ptato/demo/
mkdir -p $out_dir

(
  cd $out_dir

  # run
  nextflow run $ptato_dir/ptato.nf \
    --out_dir ./ \
    -w $work_dir \
    -c $wd/config/ptato_resources.config \
    -c $wd/config/ptato_demo.config \
    -resume \
    -with-tower \
    -N at31@sanger.ac.uk
)