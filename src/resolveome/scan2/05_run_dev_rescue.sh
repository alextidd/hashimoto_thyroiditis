#!/bin/bash
# runsub src/resolveome/scan2/05_run_dev_rescue.sh

# conda env
module load ISG/conda
conda activate scan2_dev
unset R_HOME

# dirs
wd=$(pwd)
scan2_dir=$NFS_TEAM/bin/repos/SCAN2
snakefile=$scan2_dir/snakemake/Snakefile
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/merged_non_lymphocytes/
out_dir=$lustre_dir/out/resolveome/scan2/PD63118_develop_merged_normal_rescue/
work_dir=$lustre_dir/work/scan2/PD63118_develop_merged_normal_rescue/
mkdir -p $out_dir $work_dir

# set paths (user lib first to prioritize develop r-scan2, then conda lib for dependencies)
export R_LIBS_USER=/software/conda/users/at31/scan2_dev/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2_dev/bin/python
export LD_LIBRARY_PATH=/software/conda/users/at31/scan2_dev/lib:$LD_LIBRARY_PATH
export TMPDIR=$work_dir
export XDG_CACHE_HOME=$work_dir

# initiate and configure outdir (only run once to avoid retriggering pipeline)
if [ ! -f $out_dir/scan.yaml ]; then

  # initiate outdir
  $scan2_dir/bin/scan2 -d $out_dir --snakefile $snakefile init

  # configure outdir
  (
    cd $out_dir
    $scan2_dir/bin/scan2 \
      --snakefile $snakefile \
      config \
      --verbose \
      --analysis rescue \
      --rescue-target-fdr 0.01 \
      --scripts $scan2_dir/scripts \
      --resources $scan2_dir/resources \
      --ref $wd/data/scan2/GRCh37/genome.fa \
      $(ls $lustre_dir/out/resolveome/scan2/PD63118_develop_merged_normal/call_mutations/*/scan2_object.rda | awk -F/ '{print "--scan2-object " $(NF-1), $0}' | tr '\n' ' ')
    $scan2_dir/bin/scan2 --snakefile $snakefile validate
  )

else
  echo "Config already exists at $out_dir/scan.yaml, skipping init/config"
fi

# run
(
  cd $out_dir
  $scan2_dir/bin/scan2 \
    --snakefile $snakefile \
    rescue \
    --joblimit 1000 \
    --snakemake-args " --profile $wd/config/snakemake/cluster-generic-lsf"
)
