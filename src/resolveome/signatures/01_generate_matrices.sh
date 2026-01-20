#!/bin/bash

# modules
module load python-3.11.6/perl-5.38.0 openjdk-11.0.20.1_1
module load sigprofiler/matrixGenerator-1.2.30

# dirs
seq_dir=out/resolveome/sequoia/
out_dir=out/resolveome/signatures/
mkdir -p $out_dir

# construct a matrix and vcf for each branch
mkdir -p $out_dir/branch_vcfs
python bin/get_branch_vcf.py \
  --vcf_path $seq_dir/Patient_snv_assigned_to_branches.txt \
  --outdir $out_dir/branch_vcfs/ \
  --prefix PD63118

# run sigprofilermatrixgenerator
mkdir -p $out_dir/sigprofiler/matrix_generator/
SigProfilerMatrixGenerator matrix_generator \
  PD63118 \
  GRCh38 \
  $out_dir/branch_vcfs/matrix_generator/

# reorder columns
Rscript -e "
  read.table('$out_dir/branch_vcfs/matrix_generator/output/SBS/PD63118.SBS96.all', header=TRUE, row.names=1)
  "

# plot mutational spectra
mkdir -p $out_dir/mutational_spectra
python3 -c "import sigProfilerPlotting as sigPlt; sigPlt.plotSBS('$out_dir/branch_vcfs/matrix_generator/output/SBS/PD63118.SBS96.all', '$out_dir/mutational_spectra/', 'PD63118', '96', percentage=False)"