#!/usr/bin/env python
import sys, argparse
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerExtractor import sigpro as sig
from SigProfilerAssignment import Analyzer as Analyze

exclude=["AA_signatures"]


def main():
  parser = argparse.ArgumentParser()
  parser.add_argument("--samples", required=True, help="matrix of samples")
  parser.add_argument("--output", required=True, help="path/to/output/folder")
  parser.add_argument("--signatures", required=True, help="path/to/input/denovo/signatures/file")
  parser.add_argument("--reference_genome", required=True, help="reference genome")
  parser.add_argument("--database", required=True, help="signatures database")
  args = parser.parse_args()

  print(args) # Optional: useful for debugging

  # Call the decompose_fit function here using parsed arguments
  Analyze.decompose_fit(
    samples=args.samples,
    output=args.output,
    input_type="matrix",  
    signatures=args.signatures,
	signature_database=args.database,
	exclude_signature_subgroups=exclude,
    genome_build=args.reference_genome
  )

if __name__ == "__main__":
  main()

