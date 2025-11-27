#!/bin/bash
# runsub src/resolveome/scan2/01_run.sh -M 20000

# conda env
module load ISG/conda
conda activate scan2

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
out_dir=$lustre_dir/out/resolveome/scan2/PD63118/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$lustre_dir/out/resolveome/scan2/PD63118/merged_non_lymphocytes/
work_dir=$lustre_dir/work/scan2/PD63118/

# set paths
export R_LIBS=/software/conda/users/at31/scan2/lib/R/library
export R_LIBS_USER=/software/conda/users/at31/scan2/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2/bin/python
export TMPDIR=$work_dir
export XDG_CACHE_HOME=$work_dir

# initiate outdir
if [ ! -d $out_dir ]; then
  scan2 -d $out_dir init
fi

# set up
(
  cd $out_dir
  scan2 config \
    --verbose \
    --genome hs37d5 \
    --gatk gatk3_joint \
    --abmodel-n-cores 10 \
    --ref $wd/data/scan2/GRCh37/genome.fa \
    --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
    --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
    --regions-file $wd/out/resolveome/scan2/windows/windows.bed \
    --sex male \
    --bulk-bam $bulk_dir/merged_non_lymphocytes.bam \
    $(grep "_dna_" $bam_dir/samplesheet_local.csv | cut -f14 -d, | awk '{print "--sc-bam " $0}' | tr '\n' ' ')
  scan2 validate
)

# run
(
  cd $out_dir
  scan2 run \
    --joblimit 800 \
    --cluster "bsub -q basement -M {resources.mem_mb} -R'span[hosts=1] select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}]' -n {threads} -o %logdir/%J.out -e %logdir/%J.err"
)