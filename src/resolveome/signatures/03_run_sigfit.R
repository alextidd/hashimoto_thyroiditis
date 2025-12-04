# runsub src/resolveome/signatures/03_run_sigfit.R -R -M 50000

# libraries
library(sigfit)
library(magrittr)
library(ggtree)
library(ape)
library(RColorBrewer)

# dirs
seq_dir <- "out/resolveome/sequoia/"
hdp_dir <- "out/resolveome/signatures/hdp/1,5,blood,luquette,8"
out_dir <- "out/resolveome/signatures/sigfit/1,5,blood,luquette,7b,18"
dir.create(out_dir, showWarnings = FALSE)

# reload objects
trinuc_mut_mat <-
  read.table(file.path(hdp_dir, "trinuc_mut_mat.txt"))
ref <- read.table("out/resolveome/signatures/reference_signatures.tsv")
tree <- ape::read.tree(file.path(seq_dir, "Patient_both_tree_relabelled.tree"))
tree_df <- as.data.frame(ggtree::fortify(tree))

# define the final sigs
final_sigs <-
  c("machado_2022_SBSblood",
    "SBS1", "SBS5", "SBS9", "SBS17a", "SBS17b", "SBS18", "SBS7b",
    "luquette_2022_PTA_artefact")
final_ref <- t(as.matrix(ref[, final_sigs]))
writeLines(final_sigs, file.path(out_dir, "final_signatures.txt"))

# keep branches with >50 muts
hdp_counts <- trinuc_mut_mat[rowSums(trinuc_mut_mat) > 50, ]

# run sigfit on each branch separately
sf_exposures <- list()
for (k in rownames(hdp_counts)) {
  print(k)

  k_file <- file.path(out_dir, paste0("sigfit_exposures_", k, ".rds"))
  if (file.exists(k_file)) {
    sf_exposures[[k]] <- readRDS(k_file)
  } else {

    # fit signatures
    sample_counts <- hdp_counts[k, , drop = FALSE]
    fit <- fit_signatures(counts = sample_counts,
                          signatures = final_ref,
                          iter = 20000, warmup = 10000, seed = 1756,
                          model = "poisson", chains = 4)

    # extract sf_exposures
    sf_exposures[[k]] <-
      retrieve_pars(fit, par = "exposures", hpd_prob = 0.95)

    # drop signatures with <5% contribution and refit
    keep_sigs <- colnames(sf_exposures[[k]]$mean)[sf_exposures[[k]]$mean > 0.05]
    if (length(keep_sigs) > 1 & length(keep_sigs) < ncol(sf_exposures[[k]]$mean)) {
      fit <- fit_signatures(counts = sample_counts,
                            signatures = final_ref[keep_sigs, , drop = FALSE],
                            iter = 20000, warmup = 10000,
                            model = "poisson", chains = 4)
      # extract sf_exposures
      sf_exposures[[k]] <-
        retrieve_pars(fit, par = "exposures", hpd_prob = 0.95)
    }

    # save intermediate results
    saveRDS(sf_exposures[[k]], k_file)
  }
}

# save exposures
saveRDS(sf_exposures, file.path(out_dir, "sigfit_exposures_per_branch.rds"))

# combine exposures into matrix
sf_exp <-
  sf_exposures %>%
  purrr::map(~ .x$mean) %>%
  dplyr::bind_rows() %>%
  dplyr::select(dplyr::all_of(final_sigs)) %>%
  t()
sf_exp[is.na(sf_exp)] <- 0

# save exposures matrix
write.table(sf_exp, file.path(out_dir, "sigfit_exposures_per_branch.tsv"))

# create tree plot with fitted signatures colored along branches
pdf(file.path(out_dir, "tree_with_branch_length_sigfit.pdf"),
    height = 10, width = 10)

# get exposure colours
all_cols <- brewer.pal(n = nrow(sf_exp), "Set3")
names(all_cols) <- rownames(sf_exp)

# get tip colours by celltype
cell_annots <-
  tibble::tibble(cell_ID = as.numeric(tree$tip.label)) %>%
  dplyr::left_join(
    "data/resolveome/manual_inspection/pta_additional_annotation_H1.tsv" %>%
      readr::read_tsv())
tip_cols <- ifelse(cell_annots$celltype_VDJ_recomb == "B cell", "blue",
            ifelse(cell_annots$celltype_VDJ_recomb == "alpha-beta T cell",
                    "red", "grey"))

# create tree
plot(tree, cex = 0.7, label.offset = 0.01 * max(tree_df$x), tip.color = tip_cols)

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
legend("bottomright", title = "celltype",
       legend = c("B cell", "alpha-beta T cell", "not lymphocyte"),
       col    = c("blue", "red", "grey"),
       pch    = 19,
       pt.cex = 1,
       bty    = "n")

dev.off()