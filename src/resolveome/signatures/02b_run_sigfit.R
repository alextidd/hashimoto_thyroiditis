# runsub src/resolveome/signatures/02b_run_sigfit.R -R -M 40000

# libraries
library(sigfit)
library(magrittr)
library(ggtree)
library(ape)
library(RColorBrewer)

# dirs
seq_dir <- "out/resolveome/sequoia/"
out_dir <- "out/resolveome/signatures/sigfit/"
dir.create(out_dir, showWarnings = FALSE)

# get ordering of substitutions
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[",
                   rep(sub_vec, each = 16), "]",
                   rep(c("A", "C", "G", "T"), times = 4))

# reload objects
trinuc_mut_mat <-
  read.table("out/resolveome/signatures/matrices/trinuc_mut_mat_hdp.txt",
             check.names = FALSE)
ref <-
  readr::read_tsv("out/resolveome/signatures/ref_sigs/ref_sigs.tsv") %>%
  dplyr::mutate(Type = factor(Type, levels = full_vec)) %>%
  dplyr::arrange(Type) %>%
  tibble::column_to_rownames("Type") %>%
  as.matrix()
tree <- ape::read.tree(file.path(seq_dir, "Patient_both_tree_relabelled.tree"))
tree_df <- as.data.frame(ggtree::fortify(tree))
cell_annots <-
  "data/resolveome/manual_inspection/pta_additional_annotation_H1.tsv" %>%
  readr::read_tsv() %>%
  dplyr::transmute(
    cell_id = well_ID,
    tip_id = cell_ID,
    celltype = dplyr::case_when(celltype_SHM == "Mature B cell" ~ "B cell",
                                celltype_VDJ_recomb == "B cell" ~ "B cell",
                                TRUE ~ celltype_VDJ_recomb),
    celltype = factor(celltype,
                      levels = c("B cell",
                                 "alpha-beta T cell", "not lymphocyte")),
    tip_col = dplyr::case_when(
      celltype == "B cell" ~ "blue",
      celltype == "alpha-beta T cell" ~ "red",
      celltype == "not lymphocyte" ~ "grey"
    ))

# define the final sigs
final_sigs <-
  c("SBS1", "SBS5", "SBS9", "SBS17a", "SBS17b", "SBS85",
    "machado_2022_SBSblood", "petljak_2019_ScF")
final_ref <- t(as.matrix(ref[, final_sigs]))

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
# sf_exposures <- readRDS(file.path(out_dir, "sigfit_exposures_per_branch.rds"))

# combine exposures into matrix
sf_exp <-
  sf_exposures %>%
  purrr::map(~ .x$mean) %>%
  dplyr::bind_rows() %>%
  dplyr::select(dplyr::all_of(final_sigs)) %>%
  t()
sf_exp[is.na(sf_exp)] <- 0

# create tree plot with fitted signatures colored along branches
pdf(file.path(out_dir, "tree_with_branch_length_sigfit.pdf"),
    height = 10, width = 10)

# get exposure colours
sig_cols <- brewer.pal(n = nrow(sf_exp), "Set2")
names(sig_cols) <- rownames(sf_exp)

# get celltype colours
tip_cols <-
  tibble::tibble(tip_id = as.numeric(tree$tip.label)) %>%
  dplyr::left_join(cell_annots) %>%
  dplyr::pull(tip_col)

# create tree
plot(tree, cex = 0.5, label.offset = 0.01 * max(tree_df$x), tip.color = tip_cols)

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
         xleft = x_start, xright = x_end, col = sig_cols[s], lwd = 0.25)
    x_start <- x_end
  }
}

# axis and legend
axisPhylo(side = 1, backward = FALSE)
legend("topright", title = "signatures", legend = names(sig_cols),
       fill = sig_cols, bty = "n", cex = 0.5, ncol = 1, xjust = 0.5)
legend("bottomright", title = "celltype",
       legend = c("B cell", "alpha-beta T cell", "not lymphocyte"),
       col    = c("blue", "red", "grey"),
       pch    = 19,
       pt.cex = 1,
       bty    = "n")

dev.off()