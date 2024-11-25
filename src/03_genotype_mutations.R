#!/usr/bin/env Rscript
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; for chr in {1..22} X Y ; do cmd="module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=$HOME/R-tmp-4.4; Rscript src/03_genotype_mutations.R PD63118 $chr" ; echo $cmd ; bsub -q long -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J genotype_pta_${chr} -o log/genotype_pta_${chr}_%J.out -e log/genotype_pta_${chr}_%J.err $cmd ; done

# libraries
library(magrittr)

# dirs
dir.create("out/genotyping/", showWarnings = FALSE)

# parallelise by chr / pdid
args <- commandArgs(trailingOnly = TRUE)
curr_pdid <- as.character(args[1])
curr_chr <- as.character(args[2])
print(paste("pdid:", curr_pdid, ", chr:", curr_chr))

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
  dplyr::filter(pdid == curr_pdid, chr == curr_chr)

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

# expand dnv/mnv mutations to all positions
gene_muts_dnv_mnv <-
  gene_muts %>%
  dplyr::filter(type %in% c("dnv", "mnv")) %>%
  dplyr::group_by(dplyr::across(-c(pos, mut, ref))) %>%
  dplyr::reframe(pos = pos:(pos + nchar(mut) - 1),
                 ref = strsplit(ref, split = "") %>% unlist(),
                 mut = strsplit(mut, split = "") %>% unlist())
gene_muts <-
  gene_muts %>%
  dplyr::filter(!(type %in% c("dnv", "mnv"))) %>%
  dplyr::bind_rows(gene_muts_dnv_mnv)

# genotype all sites
pta_muts <-
  gene_muts %>%
  dplyr::select(chr, pos, ref, mut, gene, type) %>%
  purrr::pmap(function(chr, pos, ref, mut, gene, type) {
    paste(chr, pos, ref, mut, gene, type, "\n") %>% cat()
    split(ss$tumour_bam, ss$sample_id) %>%
      purrr::map(function(bam) {

        # query bam
        calls <- deepSNV::bam2R(bam, chr, pos, pos,
                                q = 30, mask = 3844, mq = 30)

        # count all reads at site
        total_depth <- sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t",
                                     "DEL", "INS", "del", "ins")], na.rm = TRUE)
        
        # count mutant reads at site
        if (type %in% c("snv", "dnv", "mnv")) {
          # count mutant reads at site
          mut_depth <- sum(calls[, c(mut, tolower(mut))], na.rm = TRUE)
        } else if (type %in% c("ins", "del")) {
          # count ins or del reads at site (don't check sequence)
          mut_depth <- sum(calls[, c(type, toupper(type))], na.rm = TRUE)
        } else {
          stop("mut type not recognised!")
        }
        
        tibble::tibble(chr = chr, pos = pos, ref = ref, mut = mut, gene = gene,
                       total_depth = total_depth, mut_depth = mut_depth)
      }) %>%
      dplyr::bind_rows(.id = "sample_id")
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(nanoseq_muts)

# write mut calls to out
pta_muts %>%
  readr::write_tsv(
    paste0("out/genotyping/pta_muts_", curr_pdid, "_", curr_chr, ".tsv"))