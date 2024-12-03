# RNAseq - (CP) single cell and spatial analysis of clonal architecture in spatial tissue (3335)

# dirs
dir.create("out/driver_coverage/")

# get genes of interest
immune_panel_genes <-
  readr::read_tsv("data/immune_panel/Sanger_Immune-v1_TE-91661256_hg19_gene_list.tsv") %>%
  dplyr::pull(gene)

# get mutations
readr::read_tsv("out/genotyping/nanoseq_mutations.tsv", guess_max = Inf) %>%
  dplyr::rename("alt" = "mut") %>%
  readr::write_tsv("out/driver_coverage/nanoseq_mutations.tsv")

# get bams
bams_dir <- "/lustre/scratch125/casm/team268im/fa8/117/PTA_49686/lane4-5/"

# create a samplesheet
tibble::tibble(
  pdid = "PD63118",
  bam = list.files(bams_dir, recursive = TRUE, pattern = ".cram$", full.names = TRUE),
  mutations = "out/driver_coverage/nanoseq_mutations.tsv"
  ) %>%
  dplyr::mutate(id = stringr::str_extract(bam, "plex\\d+")) %>%
  # remove spike ins, get plex 1-19
  dplyr::filter(!grepl("_phix.cram", bam), id %in% paste0("plex", 1:19)) %>%
  readr::write_csv("out/driver_coverage/samplesheet.csv")
