#!/usr/bin/env Rscript

# libraries
library(dplyr)
library(stringr)

# args
msig_dir <- commandArgs(trailingOnly = TRUE)[1]
msig_output <- commandArgs(trailingOnly = TRUE)[2]
sigprofiler_matrix <- commandArgs(trailingOnly = TRUE)[3]

# reformat mutations x signatures matrix
df <- read.csv(msig_output, header = TRUE)

# identify columns that match the pattern 'hdp.'
hdp_cols <- grep("^hdp\\.\\d+", colnames(df), value = TRUE)
new_names <- paste0("SBS96", LETTERS[seq_along(hdp_cols)])
colnames(df)[match(hdp_cols, colnames(df))] <- new_names
colnames(df)[1:2] <- c("Mutation.type", "Trinucleotide")

df <-
  df %>%
  mutate(MutationType = paste0(str_extract(Trinucleotide, "\\w"),
                              "[", Mutation.type, "]",
                              str_extract(Trinucleotide, "(?<=\\w{2})\\w"))) %>%
  select(MutationType, all_of(colnames(df)[3:length(colnames(df))])) %>%
  arrange(MutationType)

df %>%
  write.table(file.path(msig_dir, "extracted.signatures.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

# reformat mutations x sample matrix
msig <- read.table(sigprofiler_matrix, sep = "\t", header = TRUE)
row.names(msig) <- msig$MutationType
msig <- msig[, sapply(msig, is.numeric)]
msig <- msig[, colSums(msig) != 0]
msig <- tibble::rownames_to_column(msig, "MutationType")

# save
msig %>%
  write.table(file.path(msig_dir, "sigpro_matrix.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)