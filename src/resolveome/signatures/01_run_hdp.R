# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -q basement -M50000 -R 'span[hosts=1] select[mem>50000] rusage[mem=50000]' -J resolveome_signatures_01_run_hdp -o log/%J_resolveome_signatures_01_run_hdp.out -e log/%J_resolveome_signatures_01_run_hdp.err 'Rscript src/resolveome/signatures/01_run_hdp.R'

# parameters
mut_cols <- rep(c("dodgerblue", "black", "red", "grey70",
                  "olivedrab3", "plum2"), each = 16)       
sigs_to_exclude <-
  c("petljak_2019_ScF", "SBS1", "SBS5", "SBS19")

# dirs
seq_dir <- "out/resolveome/sequoia/blood/"
out_dir <- file.path("out/resolveome/signatures/hdp/blood/")
dir.create(out_dir, showWarnings = FALSE)

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
      if (mutations$ref[j] %in% c("A", "G")) { # purine base
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

#------------------------------------------------------------------
# hierarchical dirichlet process (hdp) signature extraction
#------------------------------------------------------------------
# hdp is a bayesian nonparametric method that automatically determines
# the number of mutational signatures present in the data without
# requiring prior specification. it models signatures as mixtures of
# dirichlet distributions across the 96 trinucleotide contexts.

# initialize hdp structure with hierarchical clustering
# the hierarchy allows signatures to be shared across patients while
# allowing patient-specific signature activities
freq <- table(key_table$Patient)

# hdp hierarchy setup:
# - level 0: global base distribution (shared across all patients)
# - level 1: patient-level distributions (one per patient)
# - level 2: sample-level distributions (individual samples/branches)
hdp_mut <-
  hdp_init(
    # parent-child relationships: 0->1->2 hierarchy
    ppindex = c(0, rep(1, length(freq)),
                rep(2:(length(freq) + 1), times = freq)),
    # concentration parameter indices for each dp
    cpindex = c(1, rep(2, length(freq)),
                rep(3:(length(freq) + 2), times = freq)),
    # base distribution (uniform across 96 contexts)
    hh = rep(1, 96),
    # gamma hyperparameters for concentration parameters
    alphaa = rep(1, length(freq) + 2),
    alphab = rep(1, length(freq) + 2))

# attach mutation count data to the bottom-level dps
# each sample gets its own dp with observed trinucleotide counts
hdp_mut <-
  hdp_setdata(hdp_mut,
              dpindex = (length(freq) + 2):numdp(hdp_mut),
              as.data.frame(trinuc_mut_mat))

#------------------------------------------------------------------
# run hdp mcmc sampling
#------------------------------------------------------------------
# run multiple independent mcmc chains for convergence assessment
# each chain explores the posterior distribution of signatures and
# their activities across samples

hdp_chains <- vector("list", 4)
for (i in 1:4) {

  print(paste("running chain", i, "of 4"))

  # activate dps and initialize with 10 components
  # this provides starting clusters that can merge/split during sampling
  hdp_activated <- dp_activate(hdp_mut, 1:numdp(hdp_mut), 
                               initcc = 10, seed = i * 200)

  # run mcmc sampling:
  # - burnin: discard first 10,000 iterations (chain equilibration)
  # - n: collect 100 posterior samples
  # - space: sample every 200 iterations to reduce autocorrelation
  # - cpiter: update concentration parameters every 3 iterations
  hdp_chains[[i]] <- hdp_posterior(hdp_activated,
                               burnin = 10000,
                               n = 100,
                               space = 200,
                               cpiter = 3,
                               seed = i * 1e3)
}

# save raw chains for reproducibility
saveRDS(hdp_chains, file.path(out_dir, "hdp_chains.rds"))
# hdp_chains <- readRDS(file.path(out_dir, "hdp_chains.rds"))

#------------------------------------------------------------------
# process hdp results
#------------------------------------------------------------------
# combine multiple chains and extract consensus signatures

# create multi-chain object for convergence diagnostics
hdp_multi_chains <- hdp_multi_chain(hdp_chains)

# generate quality control plots to assess chain convergence
pdf(file.path(out_dir, "QC_plots.pdf"))
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# likelihood traces (should stabilize after burnin)
p1 <- lapply(chains(hdp_multi_chains), plot_lik, bty = "L", start = 1000)
# number of active clusters over time
p2 <- lapply(chains(hdp_multi_chains), plot_numcluster, bty = "L")
# proportion of data assigned to clusters
p3 <- lapply(chains(hdp_multi_chains), plot_data_assigned, bty = "L")
dev.off()

# extract consensus signature profiles across chains
# this averages signatures with high posterior support and
# filters out those that appear in few chains (low confidence)
hdp_multi_chains <- hdp_extract_components(hdp_multi_chains)

#------------------------------------------------------------------
# visualize extracted signatures
#------------------------------------------------------------------

# set up trinucleotide context visualization
trinuc_context <- sapply(strsplit(colnames(mut_count), "\\."), `[`, 4)
group_factor <- as.factor(
  rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each = 16))
mut_colours <- c("dodgerblue", "black", "red", "grey70",
                 "olivedrab3", "plum2")

# plot each extracted signature as a 96-context profile
# component 0 is the "background" component (noise/artefacts)
for (i in 0:hdp_multi_chains@numcomp) {
  pdf(file.path(out_dir, paste0("hdp_component_", i, ".pdf")),
      width = 12, height = 4)
  plot_comp_distn(hdp_multi_chains, cat_names = trinuc_context,
                  grouping = group_factor, col = mut_colours, comp = i,
                  col_nonsig = "grey80", show_group_labels = TRUE)
  dev.off()
}

# plot signature activities across all samples
pdf(file.path(out_dir, "signature_attribution.pdf"), width = 10, height = 8)
plot_dp_comp_exposure(
  hdp_multi_chains, 
  # show only sample-level dps (exclude higher-level hierarchy)
  dpindices = (
    length(hdp_multi_chains@comp_dp_counts) -
      nrow(trinuc_mut_mat) + 1):length(hdp_multi_chains@comp_dp_counts),
  incl_nonsig = TRUE, ylab_exp = "signature exposure", leg.title = "signature",
  col = c(RColorBrewer::brewer.pal(12, "Set3"), "magenta", "firebrick",
          RColorBrewer::brewer.pal(8, "Dark2")))
dev.off()

#------------------------------------------------------------------
# extract and save signature matrices
#------------------------------------------------------------------

freq <- table(key_table$Group)

# extract signature-sample assignment matrix
# this shows which signatures are active in which samples
dp_distn <- comp_dp_distn(hdp_multi_chains)
ndp <- nrow(dp_distn$mean)
ncomp <- ncol(dp_distn$mean)
exposures <- t(dp_distn$mean[length(freq) + 1 + 1:nrow(trinuc_mut_mat), ,
                             drop = FALSE])
colnames(exposures) <- rownames(trinuc_mut_mat)
write.table(exposures, file.path(out_dir, "mean_assignment_hdp.txt"))

# extract signature profiles (96-context distributions)
hdp_sigs <- as.data.frame(t(comp_categ_distn(hdp_multi_chains)$mean))
colnames(hdp_sigs) <- paste0("N", colnames(hdp_sigs))
write.table(hdp_sigs, file.path(out_dir, "hdp_sigs.txt"))

#------------------------------------------------------------------
# signature comparison and deconvolution
#------------------------------------------------------------------
# compare extracted signatures to known cosmic signatures using
# cosine similarity and deconvolve them into known signature components

# load cosmic reference signatures + wga artefacts
ref <- read.table("out/resolveome/signatures/cosmic_v3.4_ScF_ScB_SBSblood.tsv")
ref <- ref[, !(colnames(ref) %in% sigs_to_exclude)]

# calculate cosine similarity between hdp and reference signatures
# cosine similarity ranges 0-1, with 1 = perfect match
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

# visualize cosine similarity matrix as heatmap
pdf(file.path(out_dir, "cosine_similarities.pdf"), height = 5, width = 15)
color.palette <- colorRampPalette(c("white", "orange", "purple"))
levelplot(t(cosine_matrix[dim(cosine_matrix)[1]:1, ]), 
          col.regions = color.palette, aspect = "fill", 
          scales = list(x = list(rot = 90)))
dev.off()

#------------------------------------------------------------------
# signature deconvolution using EM algorithm
#------------------------------------------------------------------
# decompose complex hdp signatures into linear combinations of
# known reference signatures using expectation-maximization

# define reference signatures of interest (normal somatic + artefacts)
gdsigs <- c("machado_2022_SBSblood", "SBS4", "SBS7a", "SBS7b", "SBS7c",
            "SBS2", "SBS13", "SBS16", "SBS18", "SBS92",
            "SBS9", "SBS17a", "SBS17b", "SBS34", "SBS41",
            "SBS40a", "SBS40c", "SBS88", "lodato_2018_ScB")

# version identifier for output files
add <- "v1"

signatures <- t(ref[, gdsigs])
sample_list <- paste0("N", c(0:(ncol(hdp_sigs) - 1)))
colnames(hdp_sigs) <- paste0("N", c(0:(ncol(hdp_sigs) - 1)))
profiles <- hdp_sigs[, sample_list]

# initialize deconvolution results matrix
signature_fraction <- matrix(NA, nrow = nrow(signatures),
                             ncol = length(sample_list))
rownames(signature_fraction) <- rownames(signatures)
colnames(signature_fraction) <- sample_list
maxiter <- 1000

# run em algorithm for each hdp signature
for (j in seq_along(sample_list)) {
  freqs <- profiles[, j]
  freqs[is.na(freqs)] <- 0

  # em algorithm to estimate the signature contribution
  # initialize with random signature weights
  alpha <- runif(nrow(signatures))
  alpha <- alpha / sum(alpha)

  # em iterations until convergence
  for (iter in 1:maxiter) {
    # e-step: calculate expected contributions
    contr <- t(array(alpha, dim = c(nrow(signatures), 96))) * t(signatures)
    probs <- contr / array(rowSums(contr), dim = dim(contr))
    probs <- probs * freqs

    # m-step: update signature weights
    old_alpha <- alpha
    alpha <- colSums(probs) / sum(probs)

    # check convergence
    if (sum(abs(alpha - old_alpha)) < 1e-5) {
      break
    }
  }

  # saving the signature contributions for the sample
  print(j / length(sample_list))
  signature_fraction[, j] <- alpha
}

# filter to signatures contributing >10% to any sample
sigs_deconv_r2 <- list()
for (n in seq_along(sample_list)) {
  sigs_deconv_r2[[n]] <-
    rownames(signature_fraction)[signature_fraction[, n] > 0.1]
}
names(sigs_deconv_r2) <- colnames(signature_fraction)

# identify samples requiring further deconvolution (multiple signatures)
sigs_to_deconv <-
  names(sigs_deconv_r2)[unlist(lapply(sigs_deconv_r2, length)) > 1]

# second round deconvolution with refined signature sets
ref_sigs_r2 <- sort(unique(unlist(sigs_deconv_r2)))
signature_fraction_r2 <- matrix(NA, ncol = length(sigs_to_deconv),
                                nrow = length(ref_sigs_r2))
rownames(signature_fraction_r2) <- ref_sigs_r2
colnames(signature_fraction_r2) <- sigs_to_deconv

# repeat the deconvolution with the identified constitutive signatures
n <- 1
for (s in sigs_to_deconv) {
  # use only signatures identified in round 1 for this sample
  gdsigs <- sigs_deconv_r2[[s]]
  signatures <- t(ref[, gdsigs])

  signature_fraction <- matrix(NA, nrow = nrow(signatures),
                               ncol = length(sample_list))
  rownames(signature_fraction) <- rownames(signatures)
  colnames(signature_fraction) <- sample_list
  maxiter <- 1000

  freqs <- profiles[, s]
  freqs[is.na(freqs)] <- 0

  # em algorithm (same as above but with reduced signature set)
  alpha <- runif(nrow(signatures))
  alpha <- alpha / sum(alpha)
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

  # reconstruct signature and calculate fit quality
  reconsbs <- rep(0, 96)
  for (g in gdsigs) {
    reconsbs <- reconsbs + (ref[, g] * alpha[g])
  }
  cosine_reconst <- cosine(x = reconsbs, y = hdp_sigs[, s])
  print(paste0(s, ": ", cosine_reconst))

  # plot deconvolution results
  pdf(file.path(out_dir, paste0("HDP_", s, "_reconstitution_", add, ".pdf")),
      height = 10)
  par(mfrow = c(length(alpha) + 2, 1))
  par(mar = c(1, 2, 4, 1))

  # original hdp signature
  barplot(hdp_sigs[, s], col = mut_cols,
          main = paste0("HDP ", s), names.arg = "")

  # reconstructed signature
  barplot(reconsbs, col = mut_cols,
          main = paste0("reconstituted ", s, " cosine similarity to original: ",
                        round(cosine_reconst, digits = 2)))

  # individual reference signature contributions
  for (g in gdsigs) {
    add_plot <- ""
    if (grepl("SBS", g)) add_plot <- "pcawg "
    barplot(ref[, g], col = mut_cols,
            main = paste0(add_plot, g, " accounts for ",
                          round(alpha[g], digits = 2)))
  }
  dev.off()
}

#------------------------------------------------------------------
# phylogenetic visualization
#------------------------------------------------------------------
# map signature activities onto the phylogenetic tree to show
# how mutational processes changed during tissue evolution

# load tree
tree <-
  ape::read.tree(file.path(seq_dir, "Patient_both_tree_relabelled.tree"))
tree_df <- as.data.frame(fortify(tree))

# get exposure colours
all_cols <- brewer.pal(n = nrow(exposures), "Set3")
names(all_cols) <- rownames(exposures)

# create tree plot with signature exposures colored along branches
pdf(file.path(out_dir, "tree_with_branch_length_hdp.pdf"),
    height = 10, width = 10)
plot(tree, label.offset = 0.01 * max(tree_df$x), cex = 0.7)

# for each sample, draw rectangles showing signature proportions
for (sample in colnames(exposures)) {
  n <- as.numeric(substr(sample, 9, nchar(sample)))
  x_end <- tree_df$x[n]
  x_start <- tree_df$x[tree_df$parent[n]]
  x_intv <- x_end - x_start
  y <- node.height(tree)[n]
  tipnum <- sum(tree_df$isTip)

  # stack signature exposures proportionally along branch length
  for (sig in rownames(exposures)) {
    x_end <- x_start + exposures[sig, sample] * x_intv
    rect(ybottom = y - min(0.015 * tipnum, 0.3),
         ytop = y + min(0.015 * tipnum, 0.3),
         xleft = x_start, xright = x_end, col = all_cols[sig], lwd = 0.25)
    x_start <- x_end
  }
}

# axis and legend
axisPhylo(side = 1, backward = FALSE)
legend("topright", title = "signatures", legend = rownames(exposures),
       fill = all_cols, bty = "n", cex = 0.8, ncol = 1, xjust = 0.5)

dev.off()