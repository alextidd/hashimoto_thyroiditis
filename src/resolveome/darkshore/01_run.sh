nextflow run $NFS_TEAM/nextflow/external/darkshore/workflows/scVC.nf \
  --fq_list bamlist.txt \
  --out ./ \
  --sample_id PD63118 \
  -c $wd/config/darkshore.config 