# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M50000 -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]' -J resolveome_signatures_03_run_sigfit -o log/%J_resolveome_signatures_03_run_sigfit.out -e log/%J_resolveome_signatures_03_run_sigfit.err 'Rscript src/resolveome/signatures/03_run_sigfit.R'

# libraries
library(sigfit)
library(magrittr)
library(ggtree)
library(ape)
library(RColorBrewer)

# dirs
seq_dir <- "out/resolveome/sequoia/20250918/"
hdp_dir <- "out/resolveome/signatures/hdp/blood"
out_dir <- "out/resolveome/signatures/sigfit/blood"
dir.create(out_dir, showWarnings = FALSE)

# reload objects
trinuc_mut_mat <-
  read.table(file.path(hdp_dir, "trinuc_mut_mat.txt"))
ref <- read.table("out/resolveome/signatures/cosmic_v3.4_ScF_ScB_SBSblood.tsv")
tree <- ape::read.tree(file.path(seq_dir, "Patient_both_tree_relabelled.tree"))

# define the final sigs
final_sigs <-
  c("machado_2022_SBSblood", "SBS9", "SBS17a", "SBS17b", "lodato_2018_ScB")
final_ref <- t(as.matrix(ref[, final_sigs]))

# keep branches with >50 muts
hdp_counts <- trinuc_mut_mat[rowSums(trinuc_mut_mat) > 50, ]

# run sigfit on each branch separately
fit <- list()
sf_exposures <- list()
for (k in rownames(hdp_counts)) {
  print(k)
  sample_counts <- hdp_counts[k, , drop = FALSE]

  # fit signatures
  fit[[k]] <- fit_signatures(counts = sample_counts,
                             signatures = final_ref,
                             iter = 20000, warmup = 10000, seed = 1756,
                             model = "poisson", chains = 4)

  # extract sf_exposures
  sf_exposures[[k]] <-
    retrieve_pars(fit[[k]], par = "exposures", hpd_prob = 0.95)

  # drop signatures with <5% contribution and refit
  keep_sigs <- colnames(sf_exposures[[k]]$mean)[sf_exposures[[k]]$mean > 0.05]
  if (length(keep_sigs) < ncol(sf_exposures[[k]]$mean)) {
    fit[[k]] <- fit_signatures(counts = sample_counts,
                               signatures = final_ref[keep_sigs, ],
                               iter = 20000, warmup = 10000,
                               model = "poisson", chains = 4)
    # extract sf_exposures
    sf_exposures[[k]] <-
      retrieve_pars(fit[[k]], par = "exposures", hpd_prob = 0.95)
  }
}

# save fits
# saveRDS(fit, file.path(out_dir, "sigfit_fits_per_branch.rds"))
saveRDS(sf_exposures, file.path(out_dir, "sigfit_exposures_per_branch.rds"))

# combine exposures into matrix
sf_exp <-
  sf_exposures %>%
  purrr::map(~ .x$mean) %>%
  dplyr::bind_rows() %>%
  t()
sf_exp[is.na(sf_exp)] <- 0

# create tree plot with fitted signatures colored along branches
pdf(file.path(out_dir, "tree_with_branch_length_sigfit.pdf"),
    height = 10, width = 10)

# get exposure colours
all_cols <- brewer.pal(n = nrow(sf_exp), "Set3")
names(all_cols) <- rownames(sf_exp)

# create tree
plot(tree, cex = 0.7)

# for each sample, draw rectangles showing signature proportions
for (sample in colnames(sf_exp)) {
  n <- as.numeric(substr(sample, 9, nchar(sample)))
  x_end <- tree_df$x[n]
  x_start <- tree_df$x[tree_df$parent[n]]
  x_intv <- x_end - x_start
  y <- node.height(tree)[n]
  tipnum <- sum(tree_df$isTip)

  # stack signature exposures proportionally along branch length
  for (s in rownames(sf_exp)) {
    x_end <- x_start + sf_exp[s, sample] * x_intv
    rect(ybottom = y - min(0.015 * tipnum, 0.3),
         ytop = y + min(0.015 * tipnum, 0.3),
         xleft = x_start, xright = x_end, col = all_cols[s], lwd = 0.25)
    x_start <- x_end
  }
}

# axis and legend
axisPhylo(side = 1, backward = FALSE)
legend("topright", title = "signatures", legend = names(all_cols),
       fill = all_cols, bty = "n", cex = 0.8, ncol = 1, xjust = 0.5)

dev.off()