
#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J resolveome_basejumper_01_bamtofastq_run -o log/%J_resolveome_basejumper_01_bamtofastq_run.out -e log/%J_resolveome_basejumper_01_bamtofastq_run.err 'bash src/resolveome/basejumper/01_bamtofastq_run.sh'

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