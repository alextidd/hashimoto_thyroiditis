#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -n4 -M50000 -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]'     -J resolveome_scan2_01_run     -o log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_scan2_01_run.out     -e log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_scan2_01_run.err     'bash src/resolveome/scan2/01_run.sh'

# conda env
module load ISG/conda
conda activate scan2

# dirs
wd=$(pwd)
lustre_dir=$LUSTRE_125/projects/hashimoto_thyroiditis/
out_dir=$lustre_dir/out/resolveome/scan2/
bam_dir=$lustre_dir/data/bams/
bulk_dir=$bam_dir/matched_normal/

# set snakemake cache
export XDG_CACHE_HOME=$lustre_dir/work/scan2/

# out dir
scan2 -d $out_dir init

# run
(
  cd $out_dir
  scan2 config \
    --verbose \
    --ref $wd/data/scan2/GRCh37/genome.fa \
    --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
    --shapeit-refpanel $NFS_TEAM/reference/shapeit/GRCh37/ \
    --bulk-bam $bulk_dir/PD63118b_lo0001.sample.dupmarked.bam \
    --sc-bam $bam_dir/PD63118/plate3_wellA2_dna_run49882/bam/plate3_wellA2_dna_run49882.bam \
    --sc-bam $bam_dir/PD63118/plate3_wellA3_dna_run49882/bam/plate3_wellA3_dna_run49882.bam
  scan2 validate
  scan2 run --joblimit 4
)