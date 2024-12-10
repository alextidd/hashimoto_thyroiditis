#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(optparse)

# options
option_list <- list(
  make_option("--pdid", type = "character"),
  make_option("--chr", type = "character"),
  make_option("--mutations", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--sample_bams", type = "character"),
  make_option("--sample_ids", type = "character"),
  make_option("--q", type = "numeric"),
  make_option("--mask", type = "numeric"),
  make_option("--mq", type = "numeric"))
opts <- parse_args(OptionParser(option_list = option_list))
dput(opts)

# opts <- list()
# opts$pdid <- "PD63118"
# opts$chr <- "Y"
# opts$mutations <- "out/genotyping/nanoseq_mutations.tsv"
# opts$out_dir <- "out/genotyping/49882/nanoseq_mutations_genotyped_by_chr/Q30/"
# opts$sample_bams <- "data/pta/49882/plex10/49882#10.cram,data/pta/49882/plex11/49882#11.cram,data/pta/49882/plex12/49882#12.cram,data/pta/49882/plex13/49882#13.cram,data/pta/49882/plex14/49882#14.cram,data/pta/49882/plex1/49882#1.cram,data/pta/49882/plex15/49882#15.cram,data/pta/49882/plex16/49882#16.cram,data/pta/49882/plex17/49882#17.cram,data/pta/49882/plex18/49882#18.cram,data/pta/49882/plex19/49882#19.cram,data/pta/49882/plex20/49882#20.cram,data/pta/49882/plex21/49882#21.cram,data/pta/49882/plex22/49882#22.cram,data/pta/49882/plex23/49882#23.cram,data/pta/49882/plex24/49882#24.cram,data/pta/49882/plex2/49882#2.cram,data/pta/49882/plex25/49882#25.cram,data/pta/49882/plex26/49882#26.cram,data/pta/49882/plex27/49882#27.cram,data/pta/49882/plex28/49882#28.cram,data/pta/49882/plex29/49882#29.cram,data/pta/49882/plex30/49882#30.cram,data/pta/49882/plex31/49882#31.cram,data/pta/49882/plex32/49882#32.cram,data/pta/49882/plex33/49882#33.cram,data/pta/49882/plex34/49882#34.cram,data/pta/49882/plex3/49882#3.cram,data/pta/49882/plex35/49882#35.cram,data/pta/49882/plex36/49882#36.cram,data/pta/49882/plex37/49882#37.cram,data/pta/49882/plex38/49882#38.cram,data/pta/49882/plex39/49882#39.cram,data/pta/49882/plex40/49882#40.cram,data/pta/49882/plex41/49882#41.cram,data/pta/49882/plex42/49882#42.cram,data/pta/49882/plex43/49882#43.cram,data/pta/49882/plex44/49882#44.cram,data/pta/49882/plex4/49882#4.cram,data/pta/49882/plex45/49882#45.cram,data/pta/49882/plex46/49882#46.cram,data/pta/49882/plex47/49882#47.cram,data/pta/49882/plex48/49882#48.cram,data/pta/49882/plex49/49882#49.cram,data/pta/49882/plex50/49882#50.cram,data/pta/49882/plex51/49882#51.cram,data/pta/49882/plex52/49882#52.cram,data/pta/49882/plex53/49882#53.cram,data/pta/49882/plex54/49882#54.cram,data/pta/49882/plex5/49882#5.cram,data/pta/49882/plex55/49882#55.cram,data/pta/49882/plex56/49882#56.cram,data/pta/49882/plex57/49882#57.cram,data/pta/49882/plex58/49882#58.cram,data/pta/49882/plex59/49882#59.cram,data/pta/49882/plex60/49882#60.cram,data/pta/49882/plex61/49882#61.cram,data/pta/49882/plex62/49882#62.cram,data/pta/49882/plex63/49882#63.cram,data/pta/49882/plex64/49882#64.cram,data/pta/49882/plex6/49882#6.cram,data/pta/49882/plex65/49882#65.cram,data/pta/49882/plex66/49882#66.cram,data/pta/49882/plex67/49882#67.cram,data/pta/49882/plex68/49882#68.cram,data/pta/49882/plex69/49882#69.cram,data/pta/49882/plex70/49882#70.cram,data/pta/49882/plex71/49882#71.cram,data/pta/49882/plex72/49882#72.cram,data/pta/49882/plex73/49882#73.cram,data/pta/49882/plex74/49882#74.cram,data/pta/49882/plex7/49882#7.cram,data/pta/49882/plex75/49882#75.cram,data/pta/49882/plex76/49882#76.cram,data/pta/49882/plex77/49882#77.cram,data/pta/49882/plex78/49882#78.cram,data/pta/49882/plex79/49882#79.cram,data/pta/49882/plex80/49882#80.cram,data/pta/49882/plex8/49882#8.cram,data/pta/49882/plex9/49882#9.cram"
# opts$sample_ids <- "plex10,plex11,plex12,plex13,plex14,plex1,plex15,plex16,plex17,plex18,plex19,plex20,plex21,plex22,plex23,plex24,plex2,plex25,plex26,plex27,plex28,plex29,plex30,plex31,plex32,plex33,plex34,plex3,plex35,plex36,plex37,plex38,plex39,plex40,plex41,plex42,plex43,plex44,plex4,plex45,plex46,plex47,plex48,plex49,plex50,plex51,plex52,plex53,plex54,plex5,plex55,plex56,plex57,plex58,plex59,plex60,plex61,plex62,plex63,plex64,plex6,plex65,plex66,plex67,plex68,plex69,plex70,plex71,plex72,plex73,plex74,plex7,plex75,plex76,plex77,plex78,plex79,plex80,plex8,plex9"

# dirs
dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)

# get mutations
muts <-
  readr::read_tsv(opts$mutations) %>%
  dplyr::filter(chr == opts$chr, pdid == opts$pdid)

# genotype these sites in the PTA bams
bams <-
  strsplit(opts$sample_bams, ",")[[1]] %>%
  setNames(strsplit(opts$sample_ids, ",")[[1]])

if (nrow(muts)) {
  # genotype all sites
  pta_muts <-
    muts %>%
    dplyr::distinct(chr, pos, ref, mut, type) %>%
    purrr::pmap(function(chr, pos, ref, mut, type) {
      paste(chr, pos, ref, mut, type, "\n") %>% cat()
      bams %>%
        purrr::map(function(bam) {

          # query bam
          calls <- deepSNV::bam2R(bam, chr, pos, pos, q = opts$q,
                                  mask = opts$mask, mq = opts$mq)

          # count all reads at site
          total_depth <- sum(calls[, c("A", "C", "G", "T", "a", "c", "g", "t",
                                      "DEL", "INS", "del", "ins")],
                            na.rm = TRUE)
          
          # count ref reads at site
          # take first character (in case it is a deletion)
          ref_1 <- substr(ref, 1, 1)
          ref_depth <- sum(calls[, c(ref_1, tolower(ref_1))], na.rm = TRUE)

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
          
          tibble::tibble(chr = chr, pos = pos, ref = ref, mut = mut,
                        total_depth = total_depth, ref_depth = ref_depth,
                        mut_depth = mut_depth)
        }) %>%
        dplyr::bind_rows(.id = "sample_id")
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::left_join(muts)
} else {
  pta_muts <- tibble::tibble(chr = character(), pos = integer(),
                             ref = character(), mut = character(),
                             total_depth = integer(), ref_depth = integer(),
                             mut_depth = integer(), pdid = character(),
                             type = character())
}


# write mut calls to out
pta_muts %>%
  readr::write_tsv(
    paste0(opts$out_dir, "/", opts$pdid, "_", opts$chr, ".tsv"))
