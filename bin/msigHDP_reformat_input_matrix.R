#libraries
library(stringr)
library(dplyr)
library(magrittr)

file <- commandArgs(trailingOnly = TRUE)[1]
# file <- "out/resolveome/signatures/sigprofiler/mutational_signatures/msighdp_branches/input/hashimoto_thyroiditis.SBS96.all"

# Read and process input data
msig <- read.table(file, sep = "\t", header = TRUE)

# Change structure of trinucleotide counts
msig <- msig %>%
  mutate(MutationType = paste0(
    substr(MutationType, 1, 1),
    substr(MutationType, 3, 3),
    substr(MutationType, 7, 7),
    substr(MutationType, 5, 5)
  ))

# Reformat hierarchy: 1) patient and 2) healthy/ibd state
rownames(msig) <- NULL
msig <- tibble::column_to_rownames(msig, "MutationType")

# Exclude specific samples
exclude <- c("PD67979c_lo0006")
msig <- msig %>% 
  dplyr::select(!contains("PD67979c_lo0006"))

# Create sample key based on colnames
sample_key <- data.frame(samples = colnames(msig))
sample_key$patient <- str_extract(sample_key$samples, "PD\\d{5}")
sample_key$ibd_status <- ifelse(
  str_detect(sample_key$samples, "PD67|PD68"),
  "ibd",
  "normal"
)

# Construct hierarchy vector
samples <- colnames(msig)
sample_patient_vector <- as.vector(sample_key$patient)
sample_status_vector <- as.vector(sample_key$ibd_status)
sample_key_vector_colnames <- paste0(sample_status_vector, "::", samples)
colnames(msig) <- sample_key_vector_colnames

# Readjust msig dataframe
msig <- tibble::rownames_to_column(msig, "MutationType")

# Write output table
write.table(
  msig,
  file.path(
    dirname(file),
    "msigHDP_input.txt"
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)