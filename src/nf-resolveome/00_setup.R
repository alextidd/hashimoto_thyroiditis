#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
data_dir <- file.path(Sys.getenv("LUSTRE_TEAM"), "resolveome/data/bams/")
out_dir <- file.path(wd, "out/nf-resolveome/muts_and_snps/")
dir.create(out_dir, recursive = TRUE)

# get clean cells only (and cells that have not yet been assessed)
# (plate 10 cells are being reassessed for doublet status with the dna)
clean_cell_ids <-
  system("ls data/resolveome/manual_inspection/PD*.tsv", intern = TRUE) %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(!suspected_doublet | is.na(suspected_doublet) | plate == 10,
                !chr_dropout | is.na(chr_dropout)) %>%
  dplyr::pull(cell_id)

# samplesheet
ss_bams <-
  readr::read_csv(file.path(data_dir, "samplesheet_local.csv")) %>%
  dplyr::filter(seq_type %in% c("dna", "dnahyb"), cell_id %in% clean_cell_ids)

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

# get common snp sites
common_snps <-
  "../../reference/nanoseq/genome_masks/GRCh37_WGNS/SNP_GRCh37.wgns.bed.gz" %>%
  gzfile() %>%
  readr::read_tsv(col_names = c("#CHROM", "START", "POS")) %>%
  dplyr::select(`#CHROM`, POS)

# get a caveman snp file from a sample with high coverage
# extract common snps that are heterozygous
# 0.3 < VAF < 0.7 and DP > 50
# PD63118b_lo0044 has the highest coverage at 68X according to picard
caveman_snps <-
  "/nfs/irods-cgp-sr12-sdc/intproj/3464/sample/PD63118b_lo0044/PD63118b_lo0044.v1.caveman_c.snps.vcf.gz" %>%
  readr::read_tsv(comment = "##") %>%
  dplyr::mutate(
    DP = strsplit(INFO, ";") %>% purrr::map_chr(~ .x[grepl("^DP=", .x)]) %>%
      strsplit("=") %>% purrr::map_chr(~ .x[2]) %>% as.integer(),
    VAF = gsub(".*:", "", TUMOUR) %>% as.numeric()) %>%
  dplyr::filter(DP > 50, VAF > 0.3, VAF < 0.7) %>%
  # get those at common snp sites
  dplyr::inner_join(common_snps) %>%
  dplyr::transmute(donor_id = "PD63118", chr = `#CHROM`, pos = POS, ref = REF,
                   alt = ALT) %>%
  # type the mutations
  type_mutations() %>%
  dplyr::distinct() %>%
  split(.$donor_id)

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
  dir.create(file.path("out/nf-resolveome", seq_type_i))
  ss_i %>%
    readr::write_csv(file.path("out/nf-resolveome", seq_type_i,
                               "samplesheet.csv"))
})