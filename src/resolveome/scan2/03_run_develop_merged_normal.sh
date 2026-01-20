#!/bin/bash
# runsub src/resolveome/scan2/03_run_develop_merged_normal.sh -M 20000

# modules
module load ISG/conda
module load snakemake/9.13.4
module load gatk/4.5.0.0
module load bedtools2-2.31.1/python-3.10.10
module load samtools-1.19.2/python-3.11.6
module load bcftools-1.19/python-3.11.6

# clear system contamination
unset PYTHONPATH
unset R_HOME

# load conda env
conda activate scan2-py311

# dirs
wd=$(pwd)
scan2_dir=${NFS_TEAM}bin/repos/SCAN2
lustre_dir=${LUSTRE_125}projects/hashimoto_thyroiditis
out_dir=$lustre_dir/out/resolveome/scan2/PD63118_develop_merged_normal/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/merged_non_lymphocytes/
work_dir=$lustre_dir/work/scan2/PD63118_develop_merged_normal/

# set paths
export TMPDIR=$work_dir
# export R_LIBS=/software/conda/users/at31/scan2/lib/R/library
# export R_LIBS_USER=/software/conda/users/at31/scan2/lib/R/library
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
    --scripts $scan2_dir/scripts/ \
    --resources $scan2_dir/resources/ \
    --regions-file $wd/out/resolveome/scan2/windows/windows.bed \
    --sex male \
    --bulk-bam $bulk_dir/merged_wgs_non_lymphocytes.bam \
    $(grep "_dna_" $bam_dir/samplesheet_local.csv | cut -f14 -d, | awk '{print "--sc-bam " $0}' | tr '\n' ' ')
  $scan2_dir/bin/scan2 validate
)

# run
(
  cd $out_dir
  $scan2_dir/bin/scan2 \
    --snakefile $scan2_dir/snakemake/Snakefile \
    call_mutations \
    --joblimit 800 \
    --abmodel-n-cores 10 \
    --cluster "bsub -J {name} -q basement -M {resources.mem_mb} -R'select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}]' -n {threads} -o %logdir/%J_{name}.out -e %logdir/%J_{name}.err -env 'all'" \
    --snakemake-args " --retries 2 --notemp --keep-going --latency-wait=60 --rerun-incomplete"
)