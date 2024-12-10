#!/bin/bash

muts=(out/genotyping/caveman_snps.tsv out/genotyping/nanoseq_mutations.tsv)
q_val=(30) # 10)

# get sample bams
sample_bams=$(ls -1 data/pta/49882/plex*/*.cram | grep -v phix | tr '\n' ',' | sed 's/,*$//g')
sample_ids=$(ls -1 data/pta/49882/plex*/*.cram | grep -v phix | cut -d/ -f4 | tr '\n' ',' | sed 's/,*$//g')

for chr in 1 ; do #{2..22} X Y ; do
  for muts_i in ${muts[@]} ; do
    for qual in ${q_val[@]} ; do
      echo $muts_i $qual
      echo "Sample BAMs: $sample_bams"
      echo "Sample IDs: $sample_ids"
      echo

      muts_name=$(basename "${muts_i}") 
      muts_name=${muts_name%.*}
      name=genotype_pta_${muts_name}_Q${qual}_${chr}

      # Safely quoting variables in the command
      cmd="module load ISG/rocker/rver/4.4.0 ; \
           export R_LIBS_USER=\$HOME/R-tmp-4.4 ; \
           Rscript src/03b_genotype_mutations_by_chr.R \
            --pdid PD63118 \
            --chr ${chr} \
            --mutations \"${muts_i}\" \
            --out_dir \"out/genotyping/49882/${muts_name}_genotyped_by_chr/Q${qual}/\" \
            --sample_bams \"${sample_bams}\" \
            --sample_ids \"${sample_ids}\"" \

      # submit
      bsub \
        -q long -M10000 \
        -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' \
        -J $name \
        -o log/${name}_%J.out -e log/${name}_%J.err \
        "$cmd"
      echo
    done
  done
done