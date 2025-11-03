# libraries
library(magrittr)

# 

# load bam paths
bams_dir <- file.path(Sys.getenv("LUSTRE_125"),
                      "projects/hashimoto_thyroiditis/data/bams/")
bulk <-
  tibble::tibble(donor = "PD63118b",
                 sample = "PD63118b_lo0001",
                 amp = "bulk")
ss <-
  file.path(bams_dir, "samplesheet_local.csv") %>%
  readr::read_csv() %>%
  dplyr::filter(seq_type == "dna") %>%
  dplyr::transmute(donor = donor_id, sample = id, amp = "PTA") %>%
  dplyr::bind_rows(bulk)