#!/usr/bin/env Rscript
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q long -M20000 -R 'span[hosts=1] select[mem>20000] rusage[mem=20000]' -J genotype_pta -o log/genotype_pta_%J.out -e log/genotype_pta_%J.err 'module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=$HOME/R-tmp-4.4; Rscript src/03_genotype_mutations.R'

# libraries
library(magrittr)

# dirs
dir.create("out/genotyping/")

# id of interest (mutations in the florid case)
curr_pdid <- "PD63118"

# get mutations from exome and targeted nanoseq
nanoseq_dir <- "/lustre/scratch126/casm/team268im/al28/targeted_nanoseq/"
nanoseq_muts <-
  c("exome", "immune") %>%
  purrr::set_names() %>%
  purrr::map(function(x) {
    system(paste0("ls ", nanoseq_dir, "plate_062*", x,
           "*/Analysis/plate_062*final_muts.tsv"), intern = TRUE) %>%
      readr::read_tsv() %>%
      dplyr::mutate(pdid = stringr::str_sub(sampleID, 1, 7))
  }) %>%
  dplyr::bind_rows(.id = "seq_type") %>%
  dplyr::filter(pdid == curr_pdid)
  
# run dndscv to annotate mutations
dndsout <-
  nanoseq_muts %>%
  dplyr::distinct(sampleID, chr, pos, ref, mut) %>%
  dndscv::dndscv()
nanoseq_muts <-
  nanoseq_muts %>%
  dplyr::left_join(dndsout$annotmuts, relationship = "many-to-many")

# we only want to look at mutations that are in the exome
gene_muts <-
  nanoseq_muts %>%
  dplyr::filter(!is.na(gene))

# genotype these sites in the PTA bams
ss <- readr::read_csv("out/caveman/samplesheet.csv")

# genotype all sites
pta_muts <-
  gene_muts %>%
  dplyr::select(chr, pos, ref, mut, gene) %>%
  purrr::pmap(function(chr, pos, ref, mut, gene) {
    paste(chr, pos, ref, mut, gene, "\n") %>% cat()
    ss %>%
      dplyr::filter()
    split(ss$tumour_bam, ss$sample_id) %>%
      purrr::map(function(bam) {
        calls <- deepSNV::bam2R(bam, chr, pos, pos,
                                q = 30, mask = 3844, mq = 30)
        total_depth <- sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t",
                                     "DEL", "INS", "del", "ins")], na.rm = TRUE)
        mut_depth <- sum(calls[, c(mut, tolower(mut))], na.rm = TRUE)
        tibble::tibble(chr = chr, pos = pos, ref = ref, mut = mut, gene = gene,
                       total_depth = total_depth, mut_depth = mut_depth)
      }) %>%
      dplyr::bind_rows(.id = "sample_id")
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(nanoseq_muts)

# write mut calls to out
pta_muts %>%
  readr::write_tsv("out/genotyping/pta_muts.tsv")
