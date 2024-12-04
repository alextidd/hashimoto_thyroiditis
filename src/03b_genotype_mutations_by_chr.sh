#!/bin/bash

muts=(out/genotyping/caveman_snps.tsv out/genotyping/nanoseq_mutations.tsv)
q_val=(10) # 30)


for chr in {1..22} X Y ; do
  for muts_i in ${muts[@]} ; do
    for qual in ${q_val[@]} ; do
      echo $muts_i ${q_out_dir[$i]} ${q_val[$i]}
      muts_name=$(basename ${muts_i}) ; muts_name=${muts_name%.*}
      name=genotype_pta_${muts_name}_Q${qual}_${chr}
      cmd="module load ISG/rocker/rver/4.4.0 ; \
           export R_LIBS_USER=$HOME/R-tmp-4.4 ; \
           Rscript src/03b_genotype_mutations_by_chr.R \
            --pdid PD63118 \
            --chr ${chr} \
            --mutations ${muts_i}
            --out_dir out/genotyping/${muts_name}_genotyped_by_chr/Q${qual}/ \
            --q ${qual} --mask 3844 --mq 30"
      bsub \
        -q long -M10000 \
        -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' \
        -J $name \
        -o log/${name}_%J.out -e log/${name}_%J.err \
        $cmd
      echo
    done
  done
done