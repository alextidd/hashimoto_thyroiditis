# runsub src/resolveome/signatures/01_run_hdp.R -R -M 50000

# palette
mut_colours_96 <- rep(c("dodgerblue", "black", "red", "grey70",
                        "olivedrab3", "plum2"), each = 16)
mut_colours <- c("dodgerblue", "black", "red", "grey70",
                 "olivedrab3", "plum2")

# dirs
seq_dir <- "out/resolveome/sequoia/"
out_dir <- file.path("out/resolveome/signatures/hdp/1,5,blood,luquette,8/")
dir.create(out_dir, showWarnings = FALSE)

# libraries
source("bin/utils.R")
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
hdp_multi_chains_1 <- hdp_multi_chain(hdp_chains)

# generate quality control plots to assess chain convergence
pdf(file.path(out_dir, "QC_plots.pdf"))
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
# likelihood traces (should stabilize after burnin)
p1 <- lapply(chains(hdp_multi_chains_1), plot_lik, bty = "L", start = 1000)
# number of active clusters over time
p2 <- lapply(chains(hdp_multi_chains_1), plot_numcluster, bty = "L")
# proportion of data assigned to clusters
p3 <- lapply(chains(hdp_multi_chains_1), plot_data_assigned, bty = "L")
dev.off()

# extract consensus signature profiles across chains
# this averages signatures with high posterior support and
# filters out those that appear in few chains (low confidence)
hdp_multi_chains <- hdp_extract_components(hdp_multi_chains_1)

#------------------------------------------------------------------
# visualize extracted signatures
#------------------------------------------------------------------

# set up trinucleotide context visualization
trinuc_context <- sapply(strsplit(colnames(mut_count), "\\."), `[`, 4)
group_factor <-
  as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each = 16))

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
# hdp_sigs <- read.table(file.path(out_dir, "hdp_sigs.txt"))

#------------------------------------------------------------------
# signature comparison and deconvolution
#------------------------------------------------------------------
# compare extracted signatures to known cosmic signatures using
# cosine similarity and deconvolve them into known signature components

# define reference signatures of interest
# (normal somatic + artefact + from studies)
pta_artefact_sig <-
  "luquette_2022_PTA_artefact"
normal_somatic_sigs <-
  c("SBS1", "SBS5",             # endogenous clock-like
    "SBS4", "SBS92",            # tobacco
    "SBS7a", "SBS7b", "SBS7c",  # UV
    "SBS2", "SBS13",            # APOBEC
    "SBS88",                    # colibactin
    "SBS18",                    # oxidative
    "SBS9",                     # SHM
    "SBS8",                     # found in mem B cells in Machado et al. (2022)
    "SBS16", "SBS17a", "SBS17b", "SBS34", "SBS41", "SBS40a", "SBS40c") # unknown
study_sigs <-
  c("machado_2022_SignatureIg", "machado_2022_SBSblood")
sigs_of_interest <-
  c(pta_artefact_sig, normal_somatic_sigs, study_sigs)
writeLines(sigs_of_interest,
           file.path(out_dir, "sigs_of_interest.txt"))

# load reference signatures
ref <- read.table("out/resolveome/signatures/reference_signatures.tsv")

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

signatures <- t(ref[, sigs_of_interest])
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
cosine_reconst_ls <- list()
for (s in sigs_to_deconv) {
  # use only signatures identified in round 1 for this sample
  sigs_of_interest_i <- sigs_deconv_r2[[s]]
  signatures <- t(ref[, sigs_of_interest_i])

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
  signature_fraction_r2[sigs_of_interest_i, n] <- alpha
  n <- n + 1

  # reconstruct signature and calculate fit quality
  reconsbs <- rep(0, 96)
  for (g in sigs_of_interest_i) {
    reconsbs <- reconsbs + (ref[, g] * alpha[g])
  }
  cosine_reconst <- cosine(x = reconsbs, y = hdp_sigs[, s])
  print(paste0(s, ": ", cosine_reconst))
  cosine_reconst_ls[[s]] <- cosine_reconst[1, 1]

  # plot deconvolution results
  pdf(file.path(out_dir, paste0("HDP_", s, "_reconstitution.pdf")),
      height = 10)
  par(mfrow = c(length(alpha) + 2, 1))
  par(mar = c(1, 2, 4, 1))

  # original hdp signature
  barplot(hdp_sigs[, s], col = mut_colours_96,
          main = paste0("HDP ", s), names.arg = "")

  # reconstructed signature
  barplot(reconsbs, col = mut_colours_96,
          main = paste0("reconstituted ", s, " cosine similarity to original: ",
                        round(cosine_reconst, digits = 2)))

  # individual reference signature contributions
  for (g in sigs_of_interest_i) {
    add_plot <- ""
    if (grepl("SBS", g)) add_plot <- "pcawg "
    barplot(ref[, g], col = mut_colours_96,
            main = paste0(add_plot, g, " accounts for ",
                          round(alpha[g], digits = 2)))
  }
  dev.off()
}

# save signature fractions
signature_fraction_r2 %>%
  t() %>%
  tibble::as_tibble(rownames = "component") %>%
  tidyr::pivot_longer(cols = -c("component"),
                      names_to = "signature", values_to = "contribution") %>%
  dplyr::filter(!is.na(contribution)) %>%
  dplyr::mutate(cosine_reconstruction = purrr::map_dbl(component, ~ cosine_reconst_ls[[.x]])) %>%
  dplyr::arrange(component, -contribution) %>%
  readr::write_tsv(file.path(out_dir, "hdp_em_signature_deconvolution.tsv"))

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


# rerun with all signatures (from Luke)
subset_ref <- t(ref[, setdiff(colnames(ref), c("lodato_2018_ScB", "petljak_2019_ScF", "SBS17"))])
cosmic_fit <-
  fit_signatures(1e4 * comp_categ_distn(hdp_multi_chains)$mean, subset_ref,
                 iter = 3000, chains = 4, cores = 2, seed = 0xC0FFEE)
apply(retrieve_pars(cosmic_fit, "exposures")$mean, 1, function(x) rownames(subset_ref)[x > 0.1])

# rerun with expected signatures
exp_signatures <- signatures[setdiff(rownames(signatures), "SBS5"), ]
cosmic_fit <-
  fit_signatures(1e4 * comp_categ_distn(hdp_multi_chains)$mean, exp_signatures,
                 iter = 3000, chains = 4, cores = 2, seed = 0xC0FFEE)
apply(retrieve_pars(cosmic_fit, "exposures")$mean, 1, function(x) rownames(exp_signatures)[x > 0.1])

# rerun with sigs of interest
cosmic_fit <-
  fit_signatures(1e4 * comp_categ_distn(hdp_multi_chains)$mean, signatures,
                 iter = 3000, chains = 4, cores = 2, seed = 0xC0FFEE)
apply(retrieve_pars(cosmic_fit, "exposures")$mean, 1, function(x) rownames(signatures)[x > 0.1])

retrieve_pars(cosmic_fit, "exposures")$mean %>%
  