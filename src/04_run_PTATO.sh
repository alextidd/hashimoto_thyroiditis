#!/bin/bash
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J ptato -o log/ptato_%J.out -e log/ptato_%J.err 'bash src/04_run_PTATO.sh'

# dirs
wd=$(pwd)
cd out/ptato/

# run
module load singularity
module load ISG/rocker/rver/4.4.0 
export R_LIBS_USER=$HOME/R-tmp-4.4
nextflow run $wd/../PTATO/ptato.nf \
  --run.svs false \
  --run.cnvs false \
  --shapeit.reference /lustre/scratch125/casm/team268im/at31/PTATO/resources/hg38/shapeit/Phasing_reference_no_chr/ \
  --shapeit.maps /lustre/scratch125/casm/team268im/at31/PTATO/resources/hg38/shapeit/shapeit_maps_no_chr/ \
  -c ~/.nextflow/config \
  -c $wd/../PTATO/configs/run-template.config \
  -c $wd/../PTATO/configs/nextflow.config \
  -c $wd/../PTATO/configs/process.config \
  -c $wd/../PTATO/configs/resources.config \
  -c $wd/config/run.config \
  --out_dir ./ \
  -w $wd/work/ptato/ \
  -resume \
  -N at31@sanger.ac.uk