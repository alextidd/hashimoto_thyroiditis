#!/bin/bash
# runsub src/resolveome/scan2/02a_run_demo_develop.sh

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
out_dir=$lustre_dir/out/resolveome/scan2/demo_develop/
bam_dir=$lustre_dir/data/scan2/demo/
work_dir=$lustre_dir/work/scan2/demo_develop/
log_dir=$out_dir/log/
mkdir -p $bam_dir $work_dir $log_dir

# set paths (user lib first to prioritize develop r-scan2, then conda lib for dependencies)
export R_LIBS_USER=/software/conda/users/at31/scan2_dev/lib/R/library
export RETICULATE_PYTHON=/software/conda/users/at31/scan2_dev/bin/python
export LD_LIBRARY_PATH=/software/conda/users/at31/scan2_dev/lib:$LD_LIBRARY_PATH
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
    --regions 22:30000000-30999999,22:31000000-31999999 \
    --bulk-bam $bam_dir/hunamp.chr22.bam \
    --sc-bam $bam_dir/il-11.chr22.bam \
    --sc-bam $bam_dir/il-12.chr22.bam \
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