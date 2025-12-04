#!/usr/bin/env Rscript

# libraries
library(mSigHdp)
library(ICAMS)
library(cosmicsig)

# args
input_file <- commandArgs(trailingOnly = T)[1]
out_dir <- commandArgs(trailingOnly = T)[2]

#note: the input file require the trinucleotides as format ACTA, meaning A[C>A]T 
input_data <- read.table(input_file, header = T, check.names = F, sep = "\t", quote = "", row.names = 1)

input_data <- ICAMS::as.catalog(input_data, catalog.type = "counts")

results <- mSigHdp::RunHdpxParallel(
  input.catalog        = round(input_data),
  out.dir              = out_dir,
  num.child.process    = 20, 
  CPU.cores            = 20,
  seedNumber           = 42,
  K.guess              = 10, #3 found in the original study, increased for large cohort size
  burnin               = 5000,
  burnin.multiplier    = 20,
  post.n               = 200, 
  post.space           = 100, 
  multi.types          = TRUE,
  overwrite            = TRUE,
  gamma.alpha          = 1,
  gamma.beta           = 20, 
  high.confidence.prop = 0.9,
  checkpoint           = TRUE,
  verbose              = TRUE) 

save(results, file = file.path(out_dir, "multi.chains.Rdata"))
