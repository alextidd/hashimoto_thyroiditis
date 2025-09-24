# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M50000 -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]' -J resolveome_sequoia_02_extract_signatures -o log/%J_resolveome_sequoia_02_extract_signatures.out -e log/%J_resolveome_sequoia_02_extract_signatures.err 'Rscript src/resolveome/sequoia/02_extract_signatures.R'

# libraries
library(magrittr)
library(dplyr)
library(data.table)
library(hdp)
library(RColorBrewer)
library(lsa)
library(lattice)
library(ape)
library(ggtree)
options(stringsAsFactors = FALSE)

# function: get 96 muts context
mutlist_to_96_contexts <- function(mutlist, genomeFile) {
  library("GenomicRanges")
  library("Rsamtools")
  library("MASS")
  samples <- unique(mutlist$SampleID)
  trinuc_mut_mat <- matrix(0, ncol = 96, nrow = length(samples))
  for (n in seq_along(samples)) {
    s <- samples[n]
    mutations <- as.data.frame(
      mutlist[mutlist$SampleID == s, c("Chr", "Pos", "Ref", "Alt")])
    colnames(mutations) <- c("chr", "pos", "ref", "mut")
    mutations$pos <- as.numeric(mutations$pos)
    mutations <- mutations[
      (mutations$ref %in% c("A", "C", "G", "T")) &
      (mutations$mut %in% c("A", "C", "G", "T")) &
      mutations$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
    mutations$trinuc_ref <- as.vector(
      scanFa(genomeFile, 
             GRanges(mutations$chr, 
                     IRanges(as.numeric(mutations$pos) - 1, 
                             as.numeric(mutations$pos) + 1))))
    ntcomp <- c(T = "A", G = "C", C = "G", A = "T")
    mutations$sub <- paste(mutations$ref, mutations$mut, sep = ">")
    mutations$trinuc_ref_py <- mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A", "G")) { # Purine base
        mutations$sub[j] <- paste(ntcomp[mutations$ref[j]], 
                                 ntcomp[mutations$mut[j]], sep = ">")
        mutations$trinuc_ref_py[j] <- paste(
          ntcomp[rev(strsplit(mutations$trinuc_ref[j], split = "")[[1]])], 
          collapse = "")
      }
    }
    freqs <- table(paste(mutations$sub, 
                       paste(substr(mutations$trinuc_ref_py, 1, 1), 
                             substr(mutations$trinuc_ref_py, 3, 3), 
                             sep = "-"), 
                       sep = ","))
    sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
    ctx_vec <- paste(rep(c("A", "C", "G", "T"), each = 4), 
                    rep(c("A", "C", "G", "T"), times = 4), sep = "-")
    full_vec <- paste(rep(sub_vec, each = 16), 
                     rep(ctx_vec, times = 6), sep = ",")
    freqs_full <- freqs[full_vec]
    freqs_full[is.na(freqs_full)] <- 0
    names(freqs_full) <- full_vec
    trinuc_mut_mat[n, ] <- freqs_full
    print(s)
  }
  colnames(trinuc_mut_mat) <- full_vec
  rownames(trinuc_mut_mat) <- samples
  return(trinuc_mut_mat)
}

# dirs
seq_dir <- "out/resolveome/sequoia/20250918/"
sig_subdir <- "lodato_2018"
out_dir <- file.path("out/resolveome/signatures/", sig_subdir)
artefact_file <- ""
dir.create(out_dir, showWarnings = FALSE)

# read muts per branch
muts_per_branch <-
  file.path(seq_dir, "Patient_both_assigned_to_branches.txt") %>%
  read.table(header = TRUE)

# get contexts
trinuc_mut_mat <-
  mutlist_to_96_contexts(
    mutlist = muts_per_branch[, c(1:4, 7)],
    genomeFile = "../../reference/gatk/GRCh38/Homo_sapiens_assembly38.fasta.gz")
samples <- rownames(trinuc_mut_mat)
key_table <- data.frame(Sample = samples,
                        Patient = substr(samples, 1, 7))

# save
write.table(trinuc_mut_mat, file.path(out_dir, "trinuc_mut_mat.txt"))
write.table(key_table, file.path(out_dir, "key_table.txt"))

# set up hdp run
freq <- table(key_table$Patient)
hdp_mut <-
  hdp_init(ppindex = c(0, rep(1, length(freq)),
                       rep(2:(length(freq) + 1), times = freq)),
           cpindex = c(1, rep(2, length(freq)),
                       rep(3:(length(freq) + 2), times = freq)),
           hh = rep(1, 96),
           alphaa = rep(1, length(freq) + 2),
           alphab = rep(1, length(freq) + 2))
hdp_mut <-
  hdp_setdata(hdp_mut,
              dpindex = (length(freq) + 2):numdp(hdp_mut),
              as.data.frame(trinuc_mut_mat))

# run hdp posterior
chlist <- vector("list", 4)
for (i in 1:4) {

  print(i)

  # activate DPs, 10 initial components
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), 
                               initcc = 10, seed = i * 200)

  # run hdp posterior
  chlist[[i]] <- hdp_posterior(hdp_activated,
                               burnin = 10000,
                               n = 100,
                               space = 200,
                               cpiter = 3,
                               seed = i * 1e3)
}
saveRDS(chlist, file.path(out_dir, "hdp_chains.rds"))
# chlist <- readRDS(file.path(out_dir, "hdp_chains.rds"))

# example multi object
mut_example_multi <- hdp_multi_chain(chlist)

# plot qc
pdf(file.path(out_dir, "QC_plots.pdf"))
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
p1 <- lapply(chains(mut_example_multi), plot_lik, bty = "L", start = 1000)
p2 <- lapply(chains(mut_example_multi), plot_numcluster, bty = "L")
p3 <- lapply(chains(mut_example_multi), plot_data_assigned, bty = "L")
dev.off()

# extract components
mut_example_multi <- hdp_extract_components(mut_example_multi)

# colour contexts
trinuc_context <- sapply(strsplit(colnames(mut_count), "\\."), `[`, 4)
group_factor <- as.factor(
  rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each = 16))
mut_colours <- c("dodgerblue", "black", "red", "grey70",
                 "olivedrab3", "plum2")

# plot hdp components
for (i in 0:mut_example_multi@numcomp) {
  pdf(file.path(out_dir, paste0("hdp_component_", i, ".pdf")),
      width = 12, height = 4)
  plot_comp_distn(mut_example_multi, cat_names = trinuc_context,
                  grouping = group_factor, col = mut_colours, comp = i,
                  col_nonsig = "grey80", show_group_labels = TRUE)
  dev.off()
}

# plot signature attribution
pdf(file.path(out_dir, "signature_attribution.pdf"), width = 10, height = 8)
plot_dp_comp_exposure(
  mut_example_multi, 
  dpindices = (length(mut_example_multi@comp_dp_counts) -
                 nrow(trinuc_mut_mat) + 1):length(mut_example_multi@comp_dp_counts),
  incl_nonsig = TRUE, ylab_exp = "Signature exposure", leg.title = "Signature",
  col = c(RColorBrewer::brewer.pal(12, "Set3"), "magenta", "firebrick",
          RColorBrewer::brewer.pal(8, "Dark2")))
dev.off()

freq <- table(key_table$Group)

dp_distn <- comp_dp_distn(mut_example_multi)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
mean_assignment <- t(dp_distn$mean[
  length(freq) + 1 + 1:nrow(trinuc_mut_mat), , drop = FALSE])
colnames(mean_assignment) <- rownames(trinuc_mut_mat)
write.table(mean_assignment, file.path(out_dir, "mean_assignment_hdp.txt"))

mean_sigs <- as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
lower <- mean_sigs <- as.data.frame(t(comp_categ_distn(mut_example_multi)$mean))
write.table(mean_sigs, file.path(out_dir, "hdp_sigs.txt"))

mut_cols <- rep(c("dodgerblue", "black", "red", "grey70",
                  "olivedrab3", "plum2"), each = 16)

# load HDP signatures
hdp_sigs <- read.table(file.path(out_dir, "hdp_sigs.txt"))

# load cosmic reference signatures
ref <-
  read.table("../../reference/cosmic/COSMIC_v3.4_SBS_GRCh38.txt", header = TRUE)
rownames(ref) <- ref$Type
ref <- ref[, -1]
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[", 
                   rep(sub_vec, each = 16), "]", 
                   rep(c("A", "C", "G", "T"), times = 4))
ref <- ref[full_vec, ]
ref <- apply(ref, 2, as.numeric)
ref[is.na(ref) | ref == 0] <- 0.00001
ref <- t(t(ref) / colSums(ref))

# load artefact signature
petljak <- readr::read_tsv("data/petljak_2019/mmc1.tsv")
ref <- cbind(ref, ScF = petljak$`SBS sc_F`)

# assess cosine similarities for all reference signatures
cosine_matrix <- data.frame(matrix(nrow = ncol(hdp_sigs), ncol = ncol(ref)))
rownames(cosine_matrix) <- colnames(hdp_sigs)
colnames(cosine_matrix) <- colnames(ref)

for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n, m] <- cosine(
      x = hdp_sigs[, rownames(cosine_matrix)[n]],
      y = ref[, colnames(cosine_matrix)[m]]
    )
  }
}

# plot output
pdf(file.path(out_dir, "cosine_similarities.pdf"), height = 5, width = 15)
color.palette <- colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1, ]), 
          col.regions = color.palette, aspect = "fill", 
          scales = list(x = list(rot = 90)))
dev.off()

colnames(hdp_sigs) <- gsub("X", "N", colnames(hdp_sigs))
gdsigs <- c("SBS1", "SBS4", "SBS5", "SBS7a", "SBS7b", "SBS7c",
            "SBS2", "SBS13", "SBS16", "SBS18", "SBS92",
            "SBS9", "SBS17a", "SBS17b", "SBS34", "SBS41",
            "SBS40a", "SBS40c", "SBS88", "ScF")

add <- "v1" # add something to titles to differentiate multiple runs [OPTIONAL]

signatures <- t(ref[, gdsigs])
sample_list <- paste0("N", c(0:(ncol(hdp_sigs) - 1)))
colnames(hdp_sigs) <- paste0("N", c(0:(ncol(hdp_sigs) - 1)))
profiles <- hdp_sigs[, sample_list]

signature_fraction <- matrix(NA, nrow = nrow(signatures),
                             ncol = length(sample_list))
rownames(signature_fraction) <- rownames(signatures)
colnames(signature_fraction) <- sample_list
maxiter <- 1000

for (j in seq_along(sample_list)) {
  freqs <- profiles[, j]
  freqs[is.na(freqs)] <- 0
  # EM algowith to estimate the signature contribution
  alpha <- runif(nrow(signatures))
  alpha <- alpha / sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr <- t(array(alpha, dim = c(nrow(signatures), 96))) * t(signatures)
    probs <- contr / array(rowSums(contr), dim = dim(contr))
    probs <- probs * freqs
    old_alpha <- alpha
    alpha <- colSums(probs) / sum(probs)
    if (sum(abs(alpha - old_alpha)) < 1e-5) {
      break
    }
  }
  # Saving the signature contributions for the sample
  print(j / length(sample_list))
  signature_fraction[, j] <- alpha
}

sigs_deconv_r2 <- list()
for (n in seq_along(sample_list)) {
  sigs_deconv_r2[[n]] <- rownames(signature_fraction)[
    signature_fraction[, n] > 0.1]
}
names(sigs_deconv_r2) <- colnames(signature_fraction)

sigs_to_deconv <- names(sigs_deconv_r2)[
  unlist(lapply(sigs_deconv_r2, length)) > 1]

ref_sigs_r2 <- sort(unique(unlist(sigs_deconv_r2)))
signature_fraction_r2 <- matrix(NA, ncol = length(sigs_to_deconv),
                             nrow = length(ref_sigs_r2))
rownames(signature_fraction_r2) <- ref_sigs_r2
colnames(signature_fraction_r2) <- sigs_to_deconv

# repeat the deconvolution with the identified constitutive signatures
n <- 1
for (s in sigs_to_deconv) {
  gdsigs <- sigs_deconv_r2[[s]]
  signatures <- t(ref[, gdsigs])

  signature_fraction <- matrix(NA, nrow = nrow(signatures), 
                             ncol = length(sample_list))
  rownames(signature_fraction) <- rownames(signatures)
  colnames(signature_fraction) <- sample_list
  maxiter <- 1000

  freqs <- profiles[, s]
  freqs[is.na(freqs)] <- 0

  alpha <- runif(nrow(signatures))
  alpha <- alpha / sum(alpha) # Random start (seems to give ~identical results)
  for (iter in 1:maxiter) {
    contr <- t(array(alpha, dim = c(nrow(signatures), 96))) * t(signatures)
    probs <- contr / array(rowSums(contr), dim = dim(contr))
    probs <- probs * freqs
    old_alpha <- alpha
    alpha <- colSums(probs) / sum(probs)
    if (sum(abs(alpha - old_alpha)) < 1e-5) {
      break
    }
  }
  # saving the signature contributions for the sample
  signature_fraction_r2[gdsigs, n] <- alpha
  n <- n + 1
  reconsbs <- rep(0, 96)
  for (g in gdsigs) {
    reconsbs <- reconsbs + (ref[, g] * alpha[g])
  }
  cosine_reconst <- cosine(x = reconsbs, y = hdp_sigs[, s])
  print(paste0(s, ": ", cosine_reconst))
  pdf(file.path(out_dir, paste0("HDP_", s, "_reconstitution_", add, ".pdf")),
      height = 10)
  par(mfrow = c(length(alpha) + 2, 1))
  par(mar = c(1, 2, 4, 1))
  barplot(hdp_sigs[, s], col = mut_cols,
          main = paste0("HDP ", s), names.arg = "")
  barplot(reconsbs, col = mut_cols, 
          main = paste0("Reconstituted ", s, " cosine similarity to original: ",
                        round(cosine_reconst, digits = 2)))
  for (g in gdsigs) {
    add_plot <- ""
    if (grepl("SBS", g)) add_plot <- "PCAWG "
    barplot(ref[, g], col = mut_cols,
            main = paste0(add_plot, g, " accounts for ",
                          round(alpha[g], digits = 2)))
  }
  dev.off()
}

exposures <- read.table(file.path(out_dir, "mean_assignment_hdp.txt"))
all_cols <- brewer.pal(n = nrow(exposures), "Set3")
tree <-
  read.tree(file.path(seq_dir, "Patient_snv_tree_with_branch_length.tree"))
tree_df <- as.data.frame(fortify(tree))
samples <- colnames(exposures)

sigs <- rownames(exposures)
names(all_cols) <- sigs
pdf(file.path(out_dir, "tree_with_branch_length_coloured_final.pdf"),
    height = 11, width = 4)
plot(tree, label.offset = 0.01 * max(tree_df$x), show.tip.label = FALSE)
for (k in seq_along(samples)) {
  n <- as.numeric(substr(samples[k], 9, nchar(samples[k])))
  x_end <- tree_df$x[n]
  x_start <- tree_df$x[tree_df$parent[n]]
  x_intv <- x_end - x_start
  y <- node.height(tree)[n]
  tipnum <- sum(tree_df$isTip)
  for (s in sigs) {
    x_end <- x_start + exposures[s, samples[k]] * x_intv
    rect(ybottom = y - min(0.015 * tipnum, 0.3),
         ytop = y + min(0.015 * tipnum, 0.3),
         xleft = x_start, xright = x_end, col = all_cols[s], lwd = 0.25)
    x_start <- x_end
  }
}
axisPhylo(side = 1, backward = FALSE)
legend("topright", title = "Signatures", legend = sigs,
       fill = all_cols, bty = "n", cex = 0.8, ncol = 1, xjust = 0.5)
dev.off()

# run deconstructSigs
library(deconstructSigs)

# load grch38
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38

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
