#!/bin/bash
# runsub src/resolveome/scan2/00_run_demo.sh

# conda env
module load ISG/conda
conda activate scan2

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
out_dir=$lustre_dir/out/resolveome/scan2/demo/
bam_dir=$lustre_dir/data/scan2/demo/
work_dir=$lustre_dir/work/scan2/demo/
log_dir=$out_dir/log/
mkdir -p $bam_dir $out_dir $work_dir $log_dir

# set paths
export R_LIBS=/software/conda/users/at31/scan2/lib/R/library
export R_LIBS_USER=/software/conda/users/at31/scan2/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2/bin/python
export TMPDIR=$work_dir
export XDG_CACHE_HOME=$work_dir

# download bams
# (
#   cd $bam_dir
#   wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam
#   wget http://compbio.med.harvard.edu/scan-snv/hunamp.chr22.bam.bai
#   wget http://compbio.med.harvard.edu/scan-snv/il-11.chr22.bam
#   wget http://compbio.med.harvard.edu/scan-snv/il-11.chr22.bam.bai
#   wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam
#   wget http://compbio.med.harvard.edu/scan-snv/il-12.chr22.bam.bai
# )

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
    --regions 22:30000000-30999999,22:31000000-31999999 \
    --bulk-bam $bam_dir/hunamp.chr22.bam \
    --sc-bam $bam_dir/il-11.chr22.bam \
    --sc-bam $bam_dir/il-12.chr22.bam
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