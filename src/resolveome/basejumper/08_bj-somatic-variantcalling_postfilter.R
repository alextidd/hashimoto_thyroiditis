# TODO: remove germline SNPs (from caveman output)
# TODO: remove variants in >50% of cells
# TODO: remove RNA editing sites
# TODO: remove variants with global VAF > 0.25
# TODO: remove variants with cell VAF < 0.25 or mutant reads < 5
# TODO: rerun Sequoia filters - exclude cells with copy number changes from contributing to the germline filter
# TODO: generate a heterozygous SNP histogram per cell to look for outliers

# libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(tibble)
library(ape)
source("bin/build_phylogeny.R")

# dirs
seq_dir <- "out/resolveome/sequoia/"
dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

# load matrices
mats <-
  c("NV", "NR") %>%
  purrr::set_names() %>%
  purrr::map(function(i) {
    file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis",
              "out/resolveome/basejumper/bj-somatic-variantcalling/dna/PD63118/",
              "PD63118_run/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA",
              paste0("Sequoia_group_null_bino-10_rhosnp0.4_rhoindel0.4_mincov10_maxcov500_both_",
                     i, "_filtered_all.txt")) %>%
      read.table()
  })

# save matrices
purrr::walk2(names(mats), mats, function(i, mat) {
  write.table(mat, file = file.path(seq_dir, paste0(i, ".tsv")),
              quote = FALSE, sep = "\t")
})

# calculate VAFs and global VAFs
vafs <- mats$NV / (mats$NR + mats$NV)
global_vafs <- rowSums(mats$NV) / (rowSums(mats$NR) + rowSums(mats$NV))

# plot global VAF distribution
tibble::enframe(global_vafs, name = "mut_id", value = "global_vaf") %>%
  ggplot(aes(x = global_vaf, fill = grepl("chr1_", mut_id))) +
  geom_histogram(bins = 100) +
  lims(x = c(NA, 1)) +
  scale_y_continuous(trans = scales::pseudo_log_trans()) +
  theme_bw() +
  geom_vline(xintercept = 0.2, linetype = "dashed", color = "red")

# count variants in >50% of cells
n_cells_w_var <- rowSums(mats$NV > 0) / ncol(mats$NV)
table(n_cells_w_var > 0.5)

# count calls with <= 5 reads
table(mats$NR <= 5)

# count calls with mut_vaf < 0.3
table(vafs > 0.3)

# apply filters
filtered_mats <- mats
filtered_mats$NV <- filtered_mats$NV[filtered_mats$NR <= ]
filtered_mats <-
  mats %>%
  purrr::map(function(i) {
    # Set values to 0 where filters fail
    i[mats$NR <= 5] <- 0
    i[vafs < 0.3] <- 0
    # keep only variants with global VAF < 0.2
    i <- i[global_vafs < 0.2, ]
    # return
    i
  })


# remove variants with no remaining calls
filtered_mats <-
  filtered_mats %>%
  purrr::map(~ .x[rowSums(filtered_mats$NR) > 0, ])

# save nr/nv
purrr::walk2(names(filtered_mats), filtered_mats, function(i, mat) {
  write.table(mat, file = file.path(seq_dir, paste0(i, "_filtered.tsv")),
              quote = FALSE, sep = "\t")
})

# run sequoia
system(paste0("Rscript bin/build_phylogeny.R",
              " --input_nv ", seq_dir, "/NV.tsv",
              " --input_nr ", seq_dir, "/NR.tsv",
              " --output_dir ", seq_dir,
              " --snv_rho 0.4",
              " --indel_rho 0.4",
              " --germline_cutoff -10",
              " --min_cov 10",
              " --max_cov 500"))

# run sequoia (with params from basejumper run)
system(paste0("Rscript bin/build_phylogeny.R",
              " --input_nv ", seq_dir, "/NV_filtered.tsv",
              " --input_nr ", seq_dir, "/NR_filtered.tsv",
              " --output_dir ", seq_dir,
              " --snv_rho 0.4",
              " --indel_rho 0.4",
              " --germline_cutoff -10",
              " --min_cov 0", # --min_cov 10
              " --max_cov 500"))
