#!/bin/bash
# runsub src/resolveome/signatures/02a_sigprofiler_extractor.sh -M 10000

# modules
module load sigprofiler/extractor-1.1.24-GRCh38-GRCh37

# dirs
out_dir=out/resolveome/signatures/
mkdir -p $out_dir/sigprofiler/extractor/

# run sigprofilerextractor
python3 -W ignore bin/run_Sigprofiler_Extractor.py \
  --input_data $out_dir/matrices/trinuc_mut_mat_sigpro.txt \
  --output_dir $out_dir/sigprofiler/extractor/ \
  --project PD63118 \
  --reference_genome GRCh38