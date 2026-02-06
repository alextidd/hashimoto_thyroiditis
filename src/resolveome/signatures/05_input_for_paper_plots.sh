#!/bin/bash

# create a folder with all inputs for the paper plots
mkdir reports/20260206_paper_plots/inputs
cp \
  out/resolveome/sequoia/Patient_both_tree_relabelled.tree \
  out/resolveome/signatures/sigfit/sigfit_exposures_per_branch_with_indels.rds \
  out/resolveome/signatures/matrices/trinuc_mut_mat_hdp.txt \
  out/resolveome/signatures/sigfit/per_branch_exposures_filled.rds \
  data/resolveome/manual_inspection/20250902_pta_additional_annotation_H1.tsv \
  data/resolveome/manual_inspection/H1_PD63118_pta_additional_annotation.tsv \
  data/resolveome/heatmap/H1_PD63118_DNAHyb_filtered_mutations.tsv \
  data/resolveome/heatmap/H1_PD63118_DNAHyb_rescued_mutations.tsv \
  reports/20260206_paper_plots/inputs/