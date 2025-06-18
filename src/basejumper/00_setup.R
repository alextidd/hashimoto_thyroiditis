#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(purrr)
library(dplyr)

# dirs
data_dir <- file.path(Sys.getenv("LUSTRE_TEAM"), "resolveome/data/")
fastq_dir <- file.path(data_dir, "/fastqs/reads")

# read samplesheet
ss <- readr::read_csv(file.path(data_dir, "bams/samplesheet_local.csv"))

# write samplesheet for fastqs
ss %>%
  dplyr::transmute(sample_id = id, mapped = bam, index = paste0(bam, ".bai"),
                   file_type = "bam") %>%
  #check if bam exists
  dplyr::filter(file.exists(mapped)) %>%
  readr::write_csv(file.path(data_dir, "fastqs/samplesheet.csv"))

# read manual inspection results, get clean cell ids (include all of plate 10)
clean_cell_ids <-
  system("ls data/resolveome/manual_inspection/PD*.tsv", intern = TRUE) %>%
  map(readr::read_tsv) %>%
  bind_rows() %>%
  filter(!suspected_doublet | is.na(suspected_doublet),
         !chr_dropout | is.na(chr_dropout)) %>%
  pull(cell_id)

# create dnahyb/dna bj-somatic-variantcalling samplesheets
# columns: biosampleName,read1,read2,groups,isbulk,bam
ss_fastq <-
  ss %>%
  mutate(filter_lvl = "all") %>%
  bind_rows(
    ss %>%
      filter(cell_id %in% clean_cell_ids) %>%
      mutate(filter_lvl = "filter")) %>%
  transmute(
     donor_id, biosampleName = id,
    read1 = file.path(fastq_dir, paste0(id, "_1.merged.fastq.gz")),
    read2 = file.path(fastq_dir, paste0(id, "_2.merged.fastq.gz")),
    groups = donor_id, isbulk = FALSE, bam = "", seq_type, filter_lvl) %>%
  {split(., .$filter_lvl)} %>%
  map(function(df) {
    split(df %>% select(-seq_type, -filter_lvl), df$seq_type) %>%
      map(~ split(.x %>% select(-donor_id), .x$donor_id))
  })

# set up all runs
bj_runs <-
  list(list(pipeline = "bj-dna-qc", lvl = "all", seq_type = "dna"),
       list(pipeline = "bj-expression", lvl = "all", seq_type = "rna"),
       list(pipeline = "bj-somatic-variantcalling", lvl = "filter", seq_type = "dna"),
       list(pipeline = "bj-somatic-variantcalling", lvl = "filter", seq_type = "dnahyb"),
       list(pipeline = "bj-somatic-variantcalling-develop", lvl = "filter", seq_type = "dna"),
       list(pipeline = "bj-somatic-variantcalling-develop", lvl = "filter", seq_type = "dnahyb"))

# save samplesheets
bj_runs %>%
  walk(function(run) {
    lis <- ss_fastq[[run$lvl]][[run$seq_type]]
    walk2(names(lis), lis, function(i, df) {
      out_dir <- file.path("out/BaseJumper", run$pipeline, run$seq_type, i)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      df %>%
        select(biosampleName, read1, read2) %>%
        readr::write_csv(file.path(out_dir, "samplesheet.csv"))
    })
  })
