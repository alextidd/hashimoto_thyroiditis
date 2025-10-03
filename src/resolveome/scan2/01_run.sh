
module load ISG/conda
conda activate scan2

scan2 config \
  --verbose \
  --ref ../demo_data/genome.fa \
  --dbsnp $NFS_TEAM/reference/dbsnp/GRCh37/common_all_20180423.vcf \
  --shapeit-refpanel  $NFS_TEAM/reference/shapeit/GRCh37/ \
  --regions 22:30000000-30999999,22:31000000-31999999 \
  --bulk-bam ../demo_data/hunamp.chr22.bam \
  --sc-bam ../demo_data/il-11.chr22.bam \
  --sc-bam ../demo_data/il-12.chr22.bam

scan2 validate

scan2 run --joblimit 4