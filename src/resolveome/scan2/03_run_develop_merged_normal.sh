#!/bin/bash
# runsub src/resolveome/scan2/03_run_develop_merged_normal.sh -M 20000

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
out_dir=$lustre_dir/out/resolveome/scan2/PD63118_develop_merged_normal/
work_dir=$lustre_dir/work/scan2/PD63118_develop_merged_normal/
mkdir -p $bam_dir $work_dir

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
      --scripts $scan2_dir/scripts \
      --resources $scan2_dir/resources \
      --ref $wd/data/scan2/GRCh37/genome.fa \
      --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
      --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
      --regions-file $wd/out/resolveome/scan2/windows/windows.bed \
      --sex male \
      --gatk gatk4_joint \
      --bulk-bam $bulk_dir/merged_wgs_non_lymphocytes.bam \
      $(grep "_dna_" $bam_dir/samplesheet_local.csv | cut -f11 -d, | awk '{print "--sc-bam " $0}' | tr '\n' ' ')
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
    call_mutations \
    --joblimit 1000 \
    --snakemake-args " --profile $wd/config/snakemake/cluster-generic-lsf"
)
