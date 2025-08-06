#!/usr/bin/env Rscript

# libraries
library(magrittr)

# function: define mutation type based on ref and alt, split up mnvs and dnvs
type_mutations <- function(df) {
  typed_df <-
    df %>%
    dplyr::mutate(type = dplyr::case_when(
      nchar(ref) == 1 & nchar(alt) == 1 ~ "snv",
      nchar(ref) == 1 & nchar(alt) > 1 ~ "ins",
      nchar(ref) > 1 & nchar(alt) == 1 ~ "del",
      nchar(ref) == 2 & nchar(alt) == 2 ~ "dnv",
      nchar(ref) > 2 & nchar(alt) > 2 ~ "mnv"
    ))

  # expand dnv/mnv mutations to all positions
  if (any(typed_df$type %in% c("dnv", "mnv"))) {
    typed_df_dnv_mnv <-
      typed_df %>%
      dplyr::filter(type %in% c("dnv", "mnv")) %>%
      dplyr::mutate(mut_id = paste(chr, pos, ref, alt, type)) %>%
      dplyr::group_by(dplyr::across(-c(pos, alt, ref))) %>%
      dplyr::reframe(pos = pos:(pos + nchar(alt) - 1),
                     ref = strsplit(ref, split = "") %>% unlist(),
                     alt = strsplit(alt, split = "") %>% unlist()) %>%
      dplyr::select(-mut_id)
    typed_df <-
      typed_df %>%
      dplyr::filter(!(type %in% c("dnv", "mnv"))) %>%
      dplyr::bind_rows(typed_df_dnv_mnv)
  }

  # return
  typed_df
}

# dirs
wd <- getwd()
data_dir <- file.path(Sys.getenv("LUSTRE_TEAM"), "projects/hashimoto_thyroiditis/data/bams/")
out_dir <- file.path(wd, "out/resolveome/nf-resolveome/muts_and_snps/")
dir.create(out_dir, recursive = TRUE)

# samplesheet
ss_bams <-
  readr::read_csv(file.path(data_dir, "samplesheet_local.csv")) %>%
  dplyr::filter(seq_type %in% c("dna", "dnahyb"))

# get common snp sites
common_snps <-
  "../../reference/nanoseq/genome_masks/GRCh37_WGNS/SNP_GRCh37.wgns.bed.gz" %>%
  gzfile() %>%
  readr::read_tsv(col_names = c("#CHROM", "START", "POS")) %>%
  dplyr::select(`#CHROM`, POS)

# get a caveman snp file from a sample with high coverage
# extract common snps that are heterozygous
# 0.3 < VAF < 0.7 and DP > 50
# PD63118b_lo0044 has the highest coverage according to picard (68X)
# PD66718b_lo0041 has the highest coverage according to CanApps calculated cov
caveman_snps <-
  c("PD66718" = "/lustre/scratch124/casm/references/nst_links/live/3438/PD66718b_lo0041/PD66718b_lo0041.caveman_c.snps.vcf.gz",
    "PD63118" = "/nfs/irods-cgp-sr12-sdc/intproj/3464/sample/PD63118b_lo0044/PD63118b_lo0044.v1.caveman_c.snps.vcf.gz") %>%
  purrr::map(function(x) {
    x %>%
      readr::read_tsv(comment = "##") %>%
      dplyr::mutate(
        DP = strsplit(INFO, ";") %>% purrr::map_chr(~ .x[grepl("^DP=", .x)]) %>%
          strsplit("=") %>% purrr::map_chr(~ .x[2]) %>% as.integer(),
        VAF = gsub(".*:", "", TUMOUR) %>% as.numeric()) %>%
      dplyr::filter(DP > 50, VAF > 0.3, VAF < 0.7) %>%
      # get those at common snp sites
      dplyr::inner_join(common_snps) %>%
      dplyr::transmute(chr = `#CHROM`, pos = POS, ref = REF, alt = ALT) %>%
      # type the mutations
      type_mutations() %>%
      dplyr::distinct()
  })

# get mutations
nanoseq_muts <-
  readr::read_tsv("data/nanoseq/hashimoto_exome_targeted_combined_muts.tsv") %>%
  dplyr::transmute(chr, pos, ref, alt = mut,
                   donor_id = substr(sampleID, 1, 7)) %>%
  dplyr::distinct() %>%
  split(.$donor_id)

# write mutations and snps
purrr::walk2(names(nanoseq_muts), nanoseq_muts, function(donor_id_i, muts_i) {
  dir.create(file.path(out_dir, donor_id_i))
  muts_i %>%
    readr::write_tsv(file.path(out_dir, donor_id_i, "nanoseq_mutations.tsv"))
})
purrr::walk2(names(caveman_snps), caveman_snps, function(donor_id_i, snps_i) {
  snps_i %>%
    readr::write_tsv(file.path(out_dir, donor_id_i, "caveman_snps.tsv"))
  snps_i %>%
    dplyr::ungroup() %>%
    dplyr::select(chr, pos) %>%
    readr::write_tsv(file.path(out_dir, donor_id_i, "caveman_snps_positions.tsv"),
                     col_names = FALSE)
})

# write samplesheets
ss <-
  ss_bams %>%
  dplyr::mutate(
    mutations = file.path(out_dir, donor_id, "nanoseq_mutations.tsv"),
    mutations = ifelse(file.exists(mutations), mutations, NA),
    snps = file.path(out_dir, donor_id, "caveman_snps.tsv"),
    snps = ifelse(file.exists(snps), snps, NA)) %>%
  {split(., .$seq_type)}
purrr::walk2(names(ss), ss, function(seq_type_i, ss_i) {
  dir.create(file.path("out/resolveome/nf-resolveome", seq_type_i))
  ss_i %>%
    readr::write_csv(file.path("out/resolveome/nf-resolveome", seq_type_i,
                               "samplesheet.csv"))
})
