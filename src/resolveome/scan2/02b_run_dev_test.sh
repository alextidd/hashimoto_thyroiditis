#!/bin/bash
# runsub src/resolveome/scan2/02b_run_test_develop.sh

# conda env
module load ISG/conda
conda activate scan2_dev
unset R_HOME

# dirs
wd=$(pwd)
scan2_dir=$NFS_TEAM/bin/repos/SCAN2
snakefile=$scan2_dir/snakemake/Snakefile
scripts=$scan2_dir/scripts
resources=$scan2_dir/resources
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/merged_non_lymphocytes/
out_dir=$lustre_dir/out/resolveome/scan2/test_develop/
work_dir=$lustre_dir/work/scan2/test_develop/
log_dir=$out_dir/log/
mkdir -p $bam_dir $work_dir $log_dir

# set paths (user lib first to prioritize develop r-scan2, then conda lib for dependencies)
export R_LIBS_USER=/software/conda/users/at31/scan2_dev/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2_dev/bin/python
export LD_LIBRARY_PATH=/software/conda/users/at31/scan2_dev/lib:$LD_LIBRARY_PATH
export TMPDIR=$work_dir
export XDG_CACHE_HOME=$work_dir

# pick bams to run
# 1 non-lymphocyte from the bulk
head -1 $bulk_dir/bams.txt > $out_dir/bams.txt
# 2 B cell
cat $bam_dir/samplesheet_local.csv |
grep -E "plate10_wellH5_dna_|plate10_wellA8_dna_" |
cut -d, -f11 \
>> $out_dir/bams.txt
# 1 T cell
cat $bam_dir/samplesheet_local.csv |
grep "plate3_wellC12_dna_" |
cut -d, -f11 \
>> $out_dir/bams.txt

# initiate outdir
if [ ! -f $out_dir/scan.yaml ]; then
  $scan2_dir/bin/scan2 -d $out_dir --snakefile $snakefile init
fi

# set up
(
  cd $out_dir
  $scan2_dir/bin/scan2 \
    --snakefile $snakefile \
    config \
    --verbose \
    --scripts $scripts \
    --resources $resources \
    --ref $wd/data/scan2/GRCh37/genome.fa \
    --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
    --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
    --regions 1:1-10000000,9:1-10000000 \
    --bulk-bam $bulk_dir/merged_wgs_non_lymphocytes.bam \
    $(cat $out_dir/bams.txt | awk '{print "--sc-bam " $0}' | tr '\n' ' ') \
    --sex male \
    --gatk gatk4_joint
  $scan2_dir/bin/scan2 --snakefile $snakefile validate
)

# run
(
  cd $out_dir
  $scan2_dir/bin/scan2 \
    --snakefile $snakefile \
    call_mutations \
    --joblimit 1000 \
    --snakemake-args " --profile $wd/config/snakemake/cluster-generic-lsf"
)