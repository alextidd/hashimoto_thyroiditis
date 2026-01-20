# runsub src/resolveome/signatures/02b_hdp.R -R -M 10000

# libraries
library(magrittr)
library(hdp)

# dirs
mat_dir <- "out/resolveome/signatures/matrices/"
out_dir <- "out/resolveome/signatures/hdp/"
dir.create(out_dir, showWarnings = FALSE)

# palette
mut_colours <- c("dodgerblue", "black", "red", "grey70",
                 "olivedrab3", "plum2")
mut_colours_96 <- rep(mut_colours, each = 16)

# read branch matrix, transpose
trinuc_mut_mat <-
  file.path(mat_dir, "trinuc_mut_mat_hdp.txt") %>%
  read.table(header = TRUE, row.names = 1)

# create sample-patient key (assumes pdids)
samples_to_patients <-
  tibble::tibble(sample = rownames(trinuc_mut_mat)) %>%
  dplyr::mutate(patient = substr(sample, 1, 7))
samples_per_patient <- table(samples_to_patients$patient)

# # troubleshooting ----
# source("bin/utils.R")
# seq_dir <- "out/resolveome/sequoia"
# # read muts per branch
# muts_per_branch <-
#   file.path(seq_dir, "Patient_both_assigned_to_branches.txt") %>%
#   read.table(header = TRUE)

# # get contexts
# trinuc_mut_mat <-
#   mutlist_to_96_contexts(
#     mutlist = muts_per_branch[, c(1:4, 7)],
#     genomeFile = "../../reference/gatk/GRCh38/Homo_sapiens_assembly38.fasta.gz")
# samples_to_patients <-
#   data.frame(sample = rownames(trinuc_mut_mat),
#              patient = substr(rownames(trinuc_mut_mat), 1, 7))
# samples_per_patient <- table(samples_to_patients$patient)
# # ----

# setup hdp hierarchy
hdp_in <-
  hdp_init(
    # parent-child relationships: 0->1->2 hierarchy
    ppindex = c(0, rep(1, length(samples_per_patient)),
                rep(2:(length(samples_per_patient) + 1),
                times = samples_per_patient)),
    # concentration parameter indices for each dp
    cpindex = c(1, rep(2, length(samples_per_patient)),
                rep(3:(length(samples_per_patient) + 2),
                times = samples_per_patient)),
    # base distribution (uniform across 96 contexts)
    hh = rep(1, 96),
    # gamma hyperparameters for concentration parameters
    alphaa = rep(1, length(samples_per_patient) + 2),
    alphab = rep(1, length(samples_per_patient) + 2))

# attach mutation count data to sample-level dps
# samples are nodes 2 through (number of patients + 1)
hdp_mut <-
  hdp_setdata(hdp_in,
              dpindex = (length(samples_per_patient) + 2):numdp(hdp_in),
              as.data.frame(trinuc_mut_mat))

# run hdp mcmc sampling
hdp_chains <-
  c(1:4) %>%
  purrr::map(function(i) {
    print(paste("running chain", i, "of 4"))

    hdp_activated <-
      dp_activate(hdp_mut, 1:numdp(hdp_mut), initcc = 10, seed = i * 200)

    hdp_posterior(hdp_activated,
                  burnin = 10000, n = 100, space = 200, cpiter = 3,
                  seed = i * 1e3)
  })

# save chains
saveRDS(hdp_chains, file = file.path(out_dir, "hdp_chains.rds"))
# hdp_chains <- readRDS(file.path(out_dir, "hdp_chains.rds"))

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

# plot extracted signatures
# component 0 is the "background" component (noise/artefacts)
for (i in 0:hdp_multi_chains@numcomp) {
  pdf(file.path(out_dir, paste0("hdp_component_", i, ".pdf")),
      width = 12, height = 4)
  plot_comp_distn(
    hdp_multi_chains,
    cat_names = sapply(strsplit(colnames(mut_count), "\\."), `[`, 4),
    grouping = as.factor(rep(c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"), each = 16)),
    col = mut_colours, comp = i, col_nonsig = "grey80",
    show_group_labels = TRUE)
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
  col = RColorBrewer::brewer.pal(12, "Set3"))
dev.off()

# extract signature-sample assignment matrix
# this shows which signatures are active in which samples
dp_distn <- comp_dp_distn(hdp_multi_chains)
exposures <-
  t(dp_distn$mean[length(samples_per_patient) + 1 + 1:nrow(trinuc_mut_mat), ,
                  drop = FALSE])
colnames(exposures) <- rownames(trinuc_mut_mat)
write.table(exposures, file.path(out_dir, "mean_assignment_hdp.txt"))

# extract signature profiles (96-context distributions)
hdp_sigs <- as.data.frame(t(comp_categ_distn(hdp_multi_chains)$mean))
colnames(hdp_sigs) <- paste0("N", colnames(hdp_sigs))
write.table(hdp_sigs, file.path(out_dir, "hdp_sigs.txt"))

# load reference signatures
ref <- read.table("out/resolveome/signatures/reference_signatures.tsv")

# calculate cosine similarity between hdp and reference signatures
# cosine similarity ranges 0-1, with 1 = perfect match
cosine_matrix <- data.frame(matrix(nrow = ncol(hdp_sigs), ncol = ncol(ref)))
rownames(cosine_matrix) <- colnames(hdp_sigs)
colnames(cosine_matrix) <- colnames(ref)
for (n in 1:nrow(cosine_matrix)) {
  for (m in 1:ncol(cosine_matrix)) {
    cosine_matrix[n, m] <-
      lsa::cosine(x = hdp_sigs[, rownames(cosine_matrix)[n]],
                  y = ref[, colnames(cosine_matrix)[m]])
  }
}

# visualize cosine similarity matrix as heatmap
pdf(file.path(out_dir, "cosine_similarities.pdf"), height = 5, width = 15)
color.palette <- colorRampPalette(c("white", "orange", "purple"))
t(cosine_matrix[dim(cosine_matrix)[1]:1, ]) %>%
  lattice::levelplot(col.regions = color.palette, aspect = "fill", 
                     scales = list(x = list(rot = 90)))
dev.off()
