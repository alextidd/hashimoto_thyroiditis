# libraries
library(magrittr)

# function: convert mutation type names
hdp_to_sigpro <- function(hdp_names) {
  # old format example: "C>A,A-A"
  # 1. Split "C>A" from "A-A"
  mut <- sub(",.*", "", hdp_names)     # "C>A"
  flanks <- sub(".*,", "", hdp_names)  # "A-A"
  
  # 2. Extract left and right bases
  left  <- sub("-.*", "", flanks)      # "A"
  right <- sub(".*-", "", flanks)      # "A"
  
  # 3. Construct new format: A[C>A]A
  paste0(left, "[", mut, "]", right)
}

# dirs
seq_dir <- "out/resolveome/sequoia"
out_dir <- "out/resolveome/signatures/matrices"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# read muts per branch
muts_per_branch <-
  file.path(seq_dir, "Patient_both_assigned_to_branches.txt") %>%
  read.table(header = TRUE) %>%
  janitor::clean_names() %>%
  dplyr::rename()

# get contexts
trinuc_mut_mat <-
  muts_per_branch %>%
  alexr::muts_to_96_contexts(
    fasta = "../../reference/gatk/GRCh38/Homo_sapiens_assembly38.fasta.gz")

# save hdp format
trinuc_mut_mat %>%
  write.table(file.path(out_dir, "trinuc_mut_mat_hdp.txt"))

# save sigprofiler format
trinuc_mut_mat %>%
  t() %>%
  tibble::as_tibble(rownames = "MutationType") %>%
  dplyr::transmute(MutationType = hdp_to_sigpro(MutationType)) %>%
  dplyr::arrange(MutationType) %>%
  readr::write_tsv(file.path(out_dir, "trinuc_mut_mat_sigpro.txt"))

