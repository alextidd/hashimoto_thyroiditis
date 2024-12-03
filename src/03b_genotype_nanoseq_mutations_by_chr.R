#!/usr/bin/env Rscript
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; for chr in {1..22} X Y ; do cmd="module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=$HOME/R-tmp-4.4; Rscript src/03b_genotype_nanoseq_mutations_by_chr.R PD63118 $chr" ; echo $cmd ; bsub -q long -M10000 -R 'span[hosts=1] select[mem>10000] rusage[mem=10000]' -J genotype_pta_${chr} -o log/genotype_pta_${chr}_%J.out -e log/genotype_pta_${chr}_%J.err $cmd ; done

# libraries
library(magrittr)
library(optparse)

# options
option_list <- list(
  make_option("--pdid", type = "character"),
  make_option("--chr", type = "character"),
  make_option("--mutations", type = "character"),
  make_option("--out_dir", type = "character"),
  make_option("--q", type = "numeric"),
  make_option("--mask", type = "numeric"),
  make_option("--mq", type = "numeric"))
opts <- parse_args(OptionParser(option_list = option_list))
print(paste0("pdid: ", opts$pdid, ", chr: ", opts$chr))

# dirs
dir.create(opts$out_dir, showWarnings = FALSE)

# get mutations
muts <-
  readr::read_tsv(opts$mutations) %>%
  dplyr::filter(chr == opts$chr, pdid == opts$pdid)

# genotype these sites in the PTA bams
ss <-
  readr::read_csv("out/caveman/samplesheet.csv") %>%
  dplyr::filter(donor_id == opts$pdid)

# genotype all sites
pta_muts <-
  muts %>%
  dplyr::distinct(chr, pos, ref, mut, type) %>%
  purrr::pmap(function(chr, pos, ref, mut, type) {
    paste(chr, pos, ref, mut, type, "\n") %>% cat()
    split(ss$tumour_bam, ss$sample_id) %>%
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

# write mut calls to out
pta_muts %>%
  readr::write_tsv(
    paste0(opts$out_dir, "/pta_muts_", opts$pdid, "_", opts$chr, ".tsv"))
