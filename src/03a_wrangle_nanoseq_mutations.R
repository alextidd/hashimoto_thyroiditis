#!/usr/bin/env Rscript

# pdid
curr_pdid <- "PD63118"

# libraries
library(magrittr)

# dirs
dir.create("out/genotyping/", showWarnings = FALSE)

# nanoseq samplesheet
ss_muts <- readr::read_tsv("data/nanoseq/samplesheet.tsv")

# get mutations from exome and targeted nanoseq
nanoseq_muts <-
  ss_muts %>%
  purrr::pmap(function(muts_file, seq_type) {
    readr::read_tsv(muts_file) %>%
      dplyr::mutate(seq_type = seq_type)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(pdid = stringr::str_sub(sampleID, 1, 7)) %>%
  dplyr::filter(pdid == curr_pdid) %>%
  dplyr::distinct(pdid, chr, pos, ref, mut)

# write mutations
nanoseq_muts %>%
  readr::write_tsv("out/genotyping/nanoseq_mutations.tsv")

# get snps from caveman
caveman_dir <- "/nfs/cancer_ref01/nst_links/live/3438/"
caveman_snps <-
  system(paste0("ls ", caveman_dir, "/", curr_pdid, "*/*.caveman_c.snps.vcf.gz"),
         intern = TRUE) %>%
  purrr::map(function(snps_vcf) {
    readr::read_tsv(snps_vcf, comment = "##") %>%
      dplyr::select(chr = `#CHROM`, pos = POS, ref = REF, mut = ALT)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(pdid = curr_pdid) %>%
  dplyr::distinct()

# get mutations types
muts_and_snps <-
  list("nanoseq_muts" = nanoseq_muts, "caveman_snps" = caveman_snps) %>%
  purrr::map(function(df) {
    df %>%
      dplyr::mutate(type = dplyr::case_when(
        nchar(ref) == 1 & nchar(mut) == 1 ~ "snv",
        nchar(ref) == 1 & nchar(mut) > 1 ~ "ins",
        nchar(ref) > 1 & nchar(mut) == 1 ~ "del",
        nchar(ref) == 2 & nchar(mut) == 2 ~ "dnv",
        nchar(ref) > 1 & nchar(mut) > 1 ~ "mnv"
      ))
  }) %>%
  dplyr::bind_rows(.id = "source")

# expand dnv/mnv mutations to all positions
if (any(muts_and_snps$type %in% c("dnv", "mnv"))) {
  muts_and_snps_dnv_mnv <-
    muts_and_snps %>%
    dplyr::filter(type %in% c("dnv", "mnv")) %>%
    dplyr::mutate(mut_id = paste(chr, pos, ref, mut, type)) %>%
    dplyr::group_by(dplyr::across(-c(pos, mut, ref))) %>%
    dplyr::reframe(pos = pos:(pos + nchar(mut) - 1),
                   ref = strsplit(ref, split = "") %>% unlist(),
                   mut = strsplit(mut, split = "") %>% unlist()) %>%
    dplyr::select(-mut_id)
  muts_and_snps <-
    muts_and_snps %>%
    dplyr::filter(!(type %in% c("dnv", "mnv"))) %>%
    dplyr::bind_rows(muts_and_snps_dnv_mnv)
}

# write 
muts_and_snps %>%
  readr::write_tsv("out/genotyping/muts_and_snps.tsv")




# # load `scores`
# load("data/dndscv/covariates_20pc_GRCh37-38.epi_strict_outliers.Rdat")

# # load `RefCDS`
# load("data/dndscv/RefCDS_GRCh37_vF.Rdat") 

# # run dndscv to annotate mutations and snps
# annots <-
#   list("nanoseq_muts" = nanoseq_muts, "caveman_snps" = caveman_snps) %>%
#   purrr::map(function(df) {
#     dndsout <-
#       df %>%
#       dplyr::transmute(sampleID = pdid, chr, pos, ref, mut) %>%
#       dplyr::distinct() %>%
#       dndscv::dndscv(max_muts_per_gene_per_sample = Inf,
#                      max_coding_muts_per_sample = Inf,
#                      outmats = TRUE,
#                      refdb = RefCDS,
#                      cv = scores,
#                      mingenecovs = 0,
#                      onesided = FALSE)    
#     dndsout$annotmuts
#   })

# # expand dnv/mnv mutations to all positions
# if (any(nanoseq_muts_annot$type %in% c("dnv", "mnv"))) {
#   nanoseq_muts_annot_dnv_mnv <-
#     nanoseq_muts_annot %>%
#     dplyr::filter(type %in% c("dnv", "mnv")) %>%
#     dplyr::group_by(dplyr::across(-c(pos, mut, ref))) %>%
#     dplyr::reframe(pos = pos:(pos + nchar(mut) - 1),
#                    ref = strsplit(ref, split = "") %>% unlist(),
#                    mut = strsplit(mut, split = "") %>% unlist())
#   nanoseq_muts_annot <-
#     nanoseq_muts_annot %>%
#     dplyr::filter(!(type %in% c("dnv", "mnv"))) %>%
#     dplyr::bind_rows(nanoseq_muts_annot_dnv_mnv)
# }

# # write these mutations to file
# nanoseq_muts_annot %>%
#   readr::write_tsv("out/genotyping/nanoseq_mutations.tsv")

# # write these snps to file
# snps %>%
#   readr::write_tsv("out/genotyping/caveman_snps.tsv")
