#!/usr/bin/env Rscript

# libraries
library(stringr)
library(dplyr)
library(magrittr)

# args
file <- commandArgs(trailingOnly = TRUE)[1]
out_dir <- commandArgs(trailingOnly = TRUE)[2]

# read and process input data
msig <- read.table(file, sep = "\t", header = TRUE)

# change structure of trinucleotide counts
msig <- msig %>%
  mutate(MutationType = paste0(
    substr(MutationType, 1, 1),
    substr(MutationType, 3, 3),
    substr(MutationType, 7, 7),
    substr(MutationType, 5, 5)
  ))

# reformat hierarchy: 1) patient and 2) healthy/ibd state
rownames(msig) <- NULL
msig <- tibble::column_to_rownames(msig, "MutationType")

# exclude specific samples
exclude <- c("PD67979c_lo0006")
msig <- msig %>% 
  dplyr::select(!contains("PD67979c_lo0006"))

# create sample key based on colnames
sample_key <- data.frame(samples = colnames(msig))
sample_key$patient <- str_extract(sample_key$samples, "PD\\d{5}")
sample_key$ibd_status <- ifelse(
  str_detect(sample_key$samples, "PD67|PD68"),
  "ibd",
  "normal"
)

# construct hierarchy vector
samples <- colnames(msig)
sample_patient_vector <- as.vector(sample_key$patient)
sample_status_vector <- as.vector(sample_key$ibd_status)
sample_key_vector_colnames <- paste0(sample_status_vector, "::", samples)
colnames(msig) <- sample_key_vector_colnames

# readjust msig dataframe
msig <- tibble::rownames_to_column(msig, "MutationType")

# write output table
write.table(msig, file.path(out_dir, "msigHDP_input.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)