# libraries
library(deconstructSigs)
library(magrittr)
library(dplyr)

# dirs
seq_dir <- "out/resolveome/sequoia/20250918/"
out_dir <- file.path("out/resolveome/signatures/deconstruct_sigs")
dir.create(out_dir, showWarnings = FALSE)

# load grch38
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

# load ref
ref <- read.table("out/resolveome/signatures/cosmic_v3.4_ScF_ScB_SBSblood.tsv")

# transform ref
sigs_ref <- as.data.frame(t(ref))

# load tree
tree <- ape::read.tree(file.path(seq_dir, "Patient_both_tree_relabelled.tree"))
tree_df <- as.data.frame(ggtree::fortify(tree))

# load muts per branch
muts_per_branch <-
  file.path(seq_dir, "Patient_both_assigned_to_branches.txt") %>%
  read.table(header = TRUE)

# generate sigs input from branches
sigs_input <-
  muts_per_branch %>%
  janitor::clean_names() %>%
  deconstructSigs::mut.to.sigs.input(bsg = genome, sample.id = "sample_id")

# function: deconstruct signatures of interest
deconstruct_sigs <- function(sigs_of_interest) {

  # run deconstructSigs
  test <-
    rownames(sigs_input) %>%
    purrr::set_names() %>%
    purrr::map(function(i) {
      print(i)
      deconstructSigs::whichSignatures(
        tumor.ref = sigs_input,
        signatures.ref = sigs_ref[sigs_of_interest, ],
        sample.id = i,
        contexts.needed = TRUE)
    })

  # get the exposures
  sigs_exposures <-
    test %>%
    purrr::map(~ .x$weights %>% dplyr::mutate(unknown = .x$unknown)) %>%
    dplyr::bind_rows() %>%
    t()

  # colour signatures
  if (nrow(sigs_exposures) == 2) {
    all_cols <- c("#ffcc6d", "#b56de1")
  } else {
    all_cols <- RColorBrewer::brewer.pal(n = nrow(sigs_exposures), "Set3")
  }
  names(all_cols) <- rownames(sigs_exposures)
  all_cols["unknown"] <- "grey"

  # colour tips by celltype
  cell_annots <-
    tibble::tibble(cell_ID = as.numeric(tree$tip.label)) %>%
    left_join(
      "data/resolveome/manual_inspection/pta_additional_annotation_H1.tsv" %>%
        readr::read_tsv())
  tip_cols <- ifelse(cell_annots$celltype_VDJ_recomb == "B cell", "blue",
              ifelse(cell_annots$celltype_VDJ_recomb == "alpha-beta T cell",
                     "red", "grey"))

  # plot signature attribution of branches
  plot(tree, label.offset = 0.01 * max(tree_df$x), tip.color = tip_cols, cex = 0.7)
  for (k in colnames(sigs_exposures)) {
    n <- as.numeric(substr(k, 9, nchar(k)))
    x_end <- tree_df$x[n]
    x_start <- tree_df$x[tree_df$parent[n]]
    x_intv <- x_end - x_start
    y <- node.height(tree)[n]
    tipnum <- sum(tree_df$isTip)
    for (s in rownames(sigs_exposures)) {
      x_end <- x_start + sigs_exposures[s, k] * x_intv
      rect(ybottom = y - min(0.015 * tipnum, 0.3),
           ytop = y + min(0.015 * tipnum, 0.3),
           xleft = x_start, xright = x_end, col = all_cols[s], lwd = 0.25)
      x_start <- x_end
    }
  }
  axisPhylo(side = 1, backward = FALSE)
  legend("topright", title = "Signatures", legend = names(all_cols),
         fill = all_cols, bty = "n", cex = 0.8, ncol = 1, xjust = 0.5)
  legend("bottomright", title = "Celltype",
         legend = c("B cell", "alpha-beta T cell", "not lymphocyte"),
         col    = c("blue", "red", "grey"),
         pch    = 19,
         pt.cex = 1,
         bty    = "n")

}

pdf(file.path(out_dir, "tree_with_branch_length_deconstruct_sigs.pdf"), height = 10, width = 10)
deconstruct_sigs(c("machado_2022_SBSblood", "petljak_2019_ScF"))
deconstruct_sigs(c("machado_2022_SBSblood", "lodato_2018_ScB"))
deconstruct_sigs(c("machado_2022_SBSblood", "petljak_2019_ScF", "SBS17a", "SBS17b", "SBS9"))
deconstruct_sigs(c("machado_2022_SBSblood", "lodato_2018_ScB", "SBS17a", "SBS17b", "SBS9"))
deconstruct_sigs(c("machado_2022_SBSblood", "lodato_2018_ScB", "SBS17a", "SBS17b", "SBS9"))
deconstruct_sigs(c("SBS1", "SBS5", "SBS9", "lodato_2018_ScB", "SBS17b", "SBS19", "SBS17a"))
deconstruct_sigs(c("machado_2022_SBSblood", "lodato_2018_ScB", "SBS17", "SBS9"))
deconstruct_sigs(c("machado_2022_SBSblood", "lodato_2018_ScB", "SBS17", "SBS9", "machado_2022_SignatureIg"))
dev.off()
