#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -n 8 -q basement -M70000 -R 'span[hosts=1] select[mem>70000] rusage[mem=70000]'     -J resolveome_signatures_04_run_sigprofiler     -o log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_signatures_04_run_sigprofiler.out     -e log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_signatures_04_run_sigprofiler.err     'bash src/resolveome/signatures/04_run_sigprofiler.sh'

# modules
module load python-3.11.6/perl-5.38.0 openjdk-11.0.20.1_1
module load sigprofiler/extractor-1.1.24-GRCh38-GRCh37
module load sigprofiler/matrixGenerator-1.2.30

# dirs
seq_dir=out/resolveome/sequoia/
out_dir=out/resolveome/signatures/sigprofiler/
mkdir -p $out_dir

# get_branch_vcf.py will extract the snvs from each branch of a tree
# it will then construct a matrix and vcf for each branch
mkdir -p $out_dir/branch_vcfs
python bin/get_branch_vcf.py \
  --vcf_path $seq_dir/Patient_snv_assigned_to_branches.txt \
  --outdir $out_dir/branch_vcfs/ \
  --prefix PD63118

# move all of the output files into one directory
mkdir -p $out_dir/branch_vcfs/matrix_generator
mkdir -p $out_dir/branch_vcfs/vcf_with_header
mv $out_dir/branch_vcfs/PD*/matrix_generator/* $out_dir/branch_vcfs/matrix_generator/
mv $out_dir/branch_vcfs/PD*/vcf_with_header/* $out_dir/branch_vcfs/vcf_with_header/
rm -rf $out_dir/branch_vcfs/PD*

# run sigprofilermatrixgenerator
SigProfilerMatrixGenerator matrix_generator \
  hashimoto_thyroiditis \
  GRCh38 \
  $out_dir/branch_vcfs/matrix_generator/

# move the output to the directory annotate_trees
mkdir $out_dir/annotate_trees
mv $out_dir/branch_vcfs/matrix_generator/output/ $out_dir/annotate_trees/

# then copy the output in annotate trees to signatures directory
mkdir -p $out_dir/mutational_signatures/sigprofiler_branches/{input,output}
cp -r \
  $out_dir/annotate_trees/output/* \
  $out_dir/mutational_signatures/sigprofiler_branches/input/

# construct mutational spectra plots
mkdir $out_dir/mutational_signatures/sigprofiler_branches/sample_spectra/
python3 -c "import sigProfilerPlotting as sigPlt; sigPlt.plotSBS('$out_dir/mutational_signatures/sigprofiler_branches/input/SBS/hashimoto_thyroiditis.SBS96.all', '$out_dir/mutational_signatures/sigprofiler_branches/sample_spectra/', 'hashimoto_thyroiditis', '96', percentage=False)"

# extract mutational signatures using sigprofiler
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export PYTHONPATH="/software/CASM/modules/installs/sigprofiler/sigprofiler-1.1.24/env/lib/python3.11/site-packages:$PYTHONPATH"
export TORCH_USE_CUDA_DSA=0
export TORCH_DISABLE_GPU_DIAGNOSTICS=1
export PYTHONWARNINGS="ignore"
export PYTORCH_JIT=0
python3 -W ignore bin/run_Sigprofiler_Extractor.py \
  --input_data "$out_dir/mutational_signatures/sigprofiler_branches/input/SBS/hashimoto_thyroiditis.SBS96.all" \
  --output_dir "$out_dir/mutational_signatures/sigprofiler_branches/output/" \
  --project "hashimoto_thyroiditis" \
  --reference_genome "GRCh38"

# decompose mutational signatures (only run if you don't like the suggested solution)
mkdir -p $out_dir/mutational_signatures/sigprofiler_branches/output/SBS96/All_Solutions/SBS96_3_Signatures/decomposition/
python3 \
  bin/run_manual_decomp_fit.py \
  --samples $out_dir/mutational_signatures/sigprofiler_branches/input/output/SBS/hashimoto_thyroiditis.SBS96.all \
  --signatures $out_dir/mutational_signatures/sigprofiler_branches/output/SBS96/All_Solutions/SBS96_3_Signatures/Signatures/SBS96_S3_Signatures.txt \
  --reference GRCh38 \
  --output $out_dir/mutational_signatures/sigprofiler_branches/output/SBS96/All_Solutions/SBS96_3_Signatures/decomposition/

# decompose mutational signatures with custom reference


# plot mutational signatures on tree
Rscript bin/plot_tree.R \
  $out_dir/mutational_signatures/sigprofiler_branches/output/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt \
  out/resolveome/sequoia \
  $out_dir/mutational_signatures/sigprofiler_branches/output/sequoia/

