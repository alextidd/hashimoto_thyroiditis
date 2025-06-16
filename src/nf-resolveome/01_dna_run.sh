#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M2000 -R 'span[hosts=1] select[mem>2000] rusage[mem=2000]' -J 01_dna_run -o log/%J_01_dna_run.out -e log/%J_01_dna_run.err 'bash src/nf-resolveome/01_dna_run.sh'

# dirs
wd=$(pwd)

# run
(
  cd out/nf-resolveome/dna/
  nextflow run $NFS_TEAM/nextflow/nf-resolveome \
    --samplesheet samplesheet.csv \
    --seq_type dna \
    --bait_set_hyb $wd/out/twist/Probes_merged_ok_combined_Sanger_Immune-v1_TE-91661256_hg19.bed \
    --bait_set_vdj $wd/out/vdj_coverage/regions/ig_tcr_genes.bed \
    --fasta /lustre/scratch124/casm/references/ref_tmp/Homo_sapiens/GRCh37d5/genome.fa \
    --location local \
    --baf_chrs 1,4,9 \
    --out_dir ./ \
    -w $LUSTRE_TEAM/resolveome/work/nf-resolveome/dna/ \
    -resume \
    -N at31@sanger.ac.uk \
    -with-tower
)