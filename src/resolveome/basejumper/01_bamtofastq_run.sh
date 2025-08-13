
#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J 01_bamtofastq_run -o log/%J_01_bamtofastq_run.out -e log/%J_01_bamtofastq_run.err 'bash src/basejumper/01_bamtofastq_run.sh'

# modules
module load singularity/3.11.4

# run
(
  cd $LUSTRE_125/projects/hashimoto_thyroiditis/data/fastqs/
  nextflow run nf-core/bamtofastq \
    -profile singularity,sanger \
    --input samplesheet.csv \
    --outdir . \
    -w $LUSTRE_125/projects/hashimoto_thyroiditis/work/bamtofastq/ \
    -resume \
    -N at31@sanger.ac.uk
)