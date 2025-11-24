#!/bin/bash
# runsub src/resolveome/scan2/01_run.sh -M 20000

# modules
module load ISG/conda
module load snakemake/9.13.4

# clear system contamination
unset PYTHONPATH

# load conda env
conda activate scan2-py311

# dirs
wd=$(pwd)
scan2_dir=$NFS_TEAM/bin/SCAN2/
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
out_dir=$lustre_dir/out/resolveome/scan2/PD63118/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/matched_normal/
work_dir=$lustre_dir/work/scan2/PD63118/

# set paths
export TMPDIR=$work_dir
export R_LIBS=/software/conda/users/at31/scan2/lib/R/library
export R_LIBS_USER=/software/conda/users/at31/scan2/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2/bin/python
export XDG_CACHE_HOME=$work_dir

# initiate outdir
if [ ! -d $out_dir ]; then
  $scan2_dir/bin/scan2 -d $out_dir init
fi

# set up
(
  cd $out_dir
  $scan2_dir/bin/scan2 config \
    --verbose \
    --genome hs37d5 \
    --gatk gatk3_joint \
    --ref $wd/data/scan2/GRCh37/genome.fa \
    --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
    --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
    --scripts $scan2_dir/bin/scan2/scripts/ \
    --regions-file $wd/out/resolveome/scan2/windows/windows.bed \
    --sex male \
    --bulk-bam $bulk_dir/PD63118b_lo0001.sample.dupmarked.bam \
    $(grep "_dna_" $bam_dir/samplesheet_local.csv | cut -f14 -d, | awk '{print "--sc-bam " $0}' | tr '\n' ' ')
  $scan2_dir/bin/scan2 validate
)

# run
(
  cd $out_dir
  $scan2_dir/bin/scan2 call_mutations \
    --joblimit 800 \
    --cluster "bsub -q basement -M {resources.mem_mb} -R'span[hosts=1] select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}]' -n {threads} -o %logdir/%J.out -e %logdir/%J.err" \
    --snakemake-args " --retries 2 --notemp --keep-going --latency-wait=60 --rerun-incomplete --executor lsf " \
    --abmodel-n-cores 10
)