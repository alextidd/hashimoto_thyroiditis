#!/bin/bash
# cd /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis ; bsub -n 8 -q basement -M70000 -R 'span[hosts=1] select[mem>70000] rusage[mem=70000]'     -J resolveome_signatures_04_run_sigprofiler     -o log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_signatures_04_run_sigprofiler.out     -e log/$(date +%Y-%m-%d-%H%M)_%J_resolveome_signatures_04_run_sigprofiler.err     'bash src/resolveome/signatures/04_run_sigprofiler.sh'

# modules
module load python-3.11.6/perl-5.38.0 openjdk-11.0.20.1_1
module load sigprofiler/extractor-1.1.24-GRCh38-GRCh37
module load sigprofiler/matrixGenerator-1.2.30

# dirs
sp_dir=out/resolveome/signatures/sigprofiler/
out_dir=out/resolveome/signatures/msighdp/
mkdir -p $out_dir

# reformat matrix for msighdp
Rscript bin/msigHDP_reformat_input_matrix.R \
  $sp_dir/mutational_signatures/sigprofiler_branches/input/SBS/hashimoto_thyroiditis.SBS96.all \
  $out_dir

# run msighdp
Rscript bin/msigHDP_master.R \
  $out_dir/msigHDP_input.txt \
  $out_dir

# reformat the output
Rscript bin/msigHDP_reformat_output.R \
  $out_dir \
  $out_dir/extracted.signatures.csv \
  $sp_dir/mutational_signatures/sigprofiler_branches/input/SBS/hashimoto_thyroiditis.SBS96.all

# run manual decomposition
mkdir -p $out_dir/decomposition/
python bin/run_manual_decomp_fit.py \
  --samples $sp_dir/mutational_signatures/sigprofiler_branches/input/SBS/hashimoto_thyroiditis.SBS96.all \
  --signatures $out_dir/extracted.signatures.txt \
  --output $out_dir/decomposition/ \
  --reference_genome GRCh38

# run manual decomposition (custom ref signatures)

mkdir -p $out_dir/decomposition_custom/
python bin/run_manual_decomp_custom_ref.py \
  --samples $sp_dir/mutational_signatures/sigprofiler_branches/input/SBS/hashimoto_thyroiditis.SBS96.all \
  --signatures $out_dir/extracted.signatures.txt \
  --output $out_dir/decomposition_custom/ \
  --reference_genome GRCh38 \
  --database /lustre/scratch125/casm/teams/team294/users/hm14/lfs_panbody_project/reference_datasets/mitchell_pham_nature_2025_blood_sbs1_5_split.txt