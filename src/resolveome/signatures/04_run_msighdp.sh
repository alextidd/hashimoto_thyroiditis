#!/bin/bash
# runsub src/resolveome/signatures/02c_msighdp.sh -M 10000

# dirs
wd=$(pwd)
out_dir=out/resolveome/signatures/
mkdir -p $out_dir/msighdp

# reformat matrix for msighdp
Rscript bin/msigHDP_reformat_input_matrix.R \
  $out_dir/branch_vcfs/matrix_generator/output/SBS/PD63118.SBS96.all \
  $out_dir/msighdp/

# run msighdp
(
  cd $out_dir/msighdp/
  Rscript $wd/bin/msigHDP_master.R \
    msigHDP_input.txt \
    ./
)

# reformat the output
Rscript bin/msigHDP_reformat_output.R \
  $out_dir/msighdp/ \
  $out_dir/msighdp/extracted_signatures.csv \
  $out_dir/branch_vcfs/matrix_generator/output/SBS/PD63118.SBS96.all