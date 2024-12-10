# libraries
library(magrittr)

# dirs
dir.create("out/nf-resolveome/")
wd <- getwd()

# genotype muts and snps
muts_and_snps <-
  c("caveman_snps", "nanoseq_mutations") %>%
  purrr::set_names() %>%
  purrr::map(function(source) {
    paste0("out/genotyping/", source, ".tsv") %>%
      readr::read_tsv()
  }) %>%
  dplyr::bind_rows(.id = "source")
muts_and_snps %>%
  readr::write_tsv("out/nf-resolveome/mutations.tsv")

# samplesheet 49686
dir.create("out/nf-resolveome/49686/")
bams_dir <- "/lustre/scratch125/casm/team268im/fa8/117/PTA_49686/lane4-5/"
tibble::tibble(
  donor_id = "PD63118", run = "49686",
  bam = list.files(bams_dir, recursive = TRUE, pattern = ".cram$",
                   full.names = TRUE),
  mutations = paste0(wd, "/out/nf-resolveome/mutations.tsv")) %>%
  dplyr::mutate(id = stringr::str_extract(bam, "plex\\d+"),
                plex_n = gsub("plex", "", id)) %>%
  # remove spike ins, get plex 1-19
  dplyr::filter(!grepl("_phix.cram", bam), id %in% paste0("plex", 1:19)) %>%
  readr::write_tsv("out/nf-resolveome/49686/samplesheet.tsv")

# samplesheet 49882 (80 cells)
dir.create("out/nf-resolveome/49882/")
run <- "49882"
tibble::tibble(
  donor_id = "PD63118",
  plex_n = seq(1, 80)) %>%
  dplyr::transmute(
    donor_id = "PD63118", id = paste0("plex", plex_n), run, plex_n,
    bam = paste0("/seq/illumina/runs/", substr(run, 1, 2), "/", run, "/plex",
                 plex_n, "/", run, "#", plex_n, ".cram"),
    mutations = paste0(wd, "/out/nf-resolveome/mutations.tsv")) %>%
  readr::write_tsv("out/nf-resolveome/49882/samplesheet.tsv")
