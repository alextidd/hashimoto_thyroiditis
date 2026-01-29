#!/bin/bash
# runsub src/resolveome/scan2/00_run_test.sh

# conda env
module load ISG/conda
conda activate scan2

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
out_dir=$lustre_dir/out/resolveome/scan2/test/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/merged_non_lymphocytes/
work_dir=$lustre_dir/work/scan2/test/
mkdir -p $bam_dir $out_dir $work_dir $log_dir

# set paths
export R_LIBS=/software/conda/users/at31/scan2/lib/R/library
export R_LIBS_USER=/software/conda/users/at31/scan2/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2/bin/python
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
if [ ! -d $out_dir ]; then
  scan2 -d $out_dir init
fi

# set up
(
  cd $out_dir
  scan2 config \
    --verbose \
    --ref $wd/data/scan2/GRCh37/genome.fa \
    --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
    --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
    --regions 1:2480000-2500000,9:5450000-5471000 \
    --bulk-bam $bulk_dir/merged_wgs_non_lymphocytes.bam \
    $(cat $out_dir/bams.txt | awk '{print "--sc-bam " $0}' | tr '\n' ' ')
  scan2 validate
)

# run
(
  cd $out_dir
  scan2 \
    run \
    --joblimit 1000 \
    --cluster "bsub -q basement -n {threads} -M {resources.mem_mb} -R 'select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}] span[hosts=1]' -W 336:00 -oo $log_dir/{rule}.{wildcards}.out -eo $log_dir/{rule}.{wildcards}.err" \
    --snakemake-args " --jobs 500 --restart-times 1 --latency-wait 60 --keep-going --printshellcmds --max-status-checks-per-second 5 --rerun-incomplete --resources mem_mb=240000"
)