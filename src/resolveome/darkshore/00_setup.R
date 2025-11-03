# libraries
library(magrittr)

# dirs
lustre_dir <- file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis/")
fastqs_dir <- file.path(lustre_dir, "data/fastqs/reads")
out_dir <- file.path(lustre_dir, "out/resolveome/darkshore/")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# load samplesheet
ss <-
  readr::read_csv(file.path(lustre_dir, "data/bams/samplesheet_local.csv")) %>%
  dplyr::filter(seq_type == "dna") %>%
  dplyr::transmute(
    r1 = file.path(fastqs_dir, paste0(id, "_1.merged.fastq.gz")),
    r2 = file.path(fastqs_dir, paste0(id, "_2.merged.fastq.gz"))) %>%
  dplyr::filter(file.exists(r1) & file.exists(r2))

# save
ss %>%
  readr::write_tsv(file.path(out_dir, "fq_list.txt"), col_names = FALSE)