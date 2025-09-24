
# run deconstructSigs
library(deconstructSigs)

# load grch38
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

# load ref


# transform ref
sigs_ref <- as.data.frame(t(ref))
colnames(sigs_ref) <- colnames(signatures.nature2013)

# generate sigs input from branches
sigs_input <-
  muts_per_branch %>%
  janitor::clean_names() %>%
  deconstructSigs::mut.to.sigs.input(bsg = genome, sample.id = "sample_id")

# deconstruct signatures SBSblood (SBS1+5+19), ScF
sigs_of_interest <-
  c("SBS1", "SBS5", "SBS19", "ScF", "SBS17a", "SBS17b", "SBS9")
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

sigs_exposures <-
  test %>%
  purrr::map(~ .x$weights) %>%
  dplyr::bind_rows() %>%
  t()
samples <- colnames(sigs_exposures)

# colour signatures
all_cols <- brewer.pal(n = nrow(sigs_exposures), "Set3")
names(all_cols) <- rownames(sigs_exposures)

# colour tips by celltype
cell_annots <-
  tibble::tibble(id = tree$tip.label) %>%
  left_join(
    file.path(Sys.getenv("LUSTRE_125"),
              "projects/hashimoto_thyroiditis/out/resolveome/basejumper",
              "bj-somatic-variantcalling/dna/PD63118/samplesheet.csv") %>%
      readr::read_csv() %>%
      transmute(id = biosampleName, cell_id = gsub("_dna.*", "", id)) %>%
      # add new cells ids + additional annotations
      left_join(
        "data/resolveome/manual_inspection/pta_additional_annotation_H1.tsv" %>%
          readr::read_tsv() %>%
          mutate(cell_id_new = cell_ID, cell_id = well_ID)
      )
  )
tip_cols <- ifelse(cell_annots$celltype_VDJ_recomb == "B cell", "blue",
            ifelse(cell_annots$celltype_VDJ_recomb == "alpha-beta T cell", "red", "grey"))

pdf("test.pdf", height = 20, width = 20)
plot(tree, label.offset = 0.01 * max(tree_df$x), tip.color = tip_cols)
for (k in seq_along(samples)) {
  n <- as.numeric(substr(samples[k], 9, nchar(samples[k])))
  x_end <- tree_df$x[n]
  x_start <- tree_df$x[tree_df$parent[n]]
  x_intv <- x_end - x_start
  y <- node.height(tree)[n]
  tipnum <- sum(tree_df$isTip)
  for (s in rownames(sigs_exposures)) {
    x_end <- x_start + sigs_exposures[s, samples[k]] * x_intv
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
dev.off()