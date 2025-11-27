#!/bin/bash
# runsub src/resolveome/scan2/00_test.sh

# conda env
module load ISG/conda
conda activate scan2

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
out_dir=$lustre_dir/out/resolveome/scan2/test/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/matched_normal/
work_dir=$lustre_dir/work/scan2/test/

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
    --ref $wd/data/scan2/GRCh37/genome.fa \
    --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
    --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
    --bulk-bam $bulk_dir/PD63118b_lo0001.sample.dupmarked.bam \
    --regions-file $wd/out/resolveome/scan2/windows/windows.bed \
    --abmodel-n-cores 10 \
    --bulk-bam $bulk_dir/PD63118b_lo0001.sample.dupmarked.bam \
    --sc-bam $bam_dir/PD63118/plate3_wellA2_dna_run49882/bam/plate3_wellA2_dna_run49882.bam \
    --sc-bam $bam_dir/PD63118/plate3_wellA3_dna_run49882/bam/plate3_wellA3_dna_run49882.bam
  scan2 validate
)

# run
(
  cd $out_dir
  scan2 run \
    --cluster "bsub -q basement -M {resources.mem_mb} -R'span[hosts=1] select[mem>{resources.mem_mb}] rusage[mem={resources.mem_mb}]' -n {threads} -o %logdir/%J.out -e %logdir/%J.err"
)