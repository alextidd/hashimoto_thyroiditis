#!/usr/bin/env Rscript

# the longest 1p LOH event is observed in cells plate3_wellA2, affecting the
# entire 1p arm, and is therefore expected to be homozygous for all SNPs.
# we can therefore phase all SNPs in PD63118 chr1p.
# then, i will look for genic SNPs, as these will be most likely to be picked up
# in the snRNAseq.

# libraries
library(magrittr)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)
source("bin/utils.R")

# directories
wd <- getwd()
out_dir <- "out/trencadis-seq/nf-sc-geno/"
dir.create(file.path(out_dir, "snps"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "TX"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "PB"), recursive = TRUE, showWarnings = FALSE)

# save annotated TX barcodes
barcodes <- list()
ct_annot_dir <- "out/trencadis-seq/seurat/ncells_3_nfeature_500-2000_ncount_1000-15000_percent_mt_5/nf-trencadis-seq_thyroid_annotate_celltypes_cache/html/"
barcodes$TX <-
  list.files(ct_annot_dir, pattern = "seu_annotated", full.names = TRUE) %>%
  readRDS() %>%
  colnames()
writeLines(barcodes$TX, file.path(wd, out_dir, "TX/barcodes.txt"))

# save annotated PB barcodes
barcodes$PB <- rev_comp(barcodes$TX)
writeLines(barcodes$PB, file.path(wd, out_dir, "PB/barcodes.txt"))

# phase the snps
phased_snps_grch37 <-
  "out/resolveome/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(id == "plate3_wellA2_dna_run49882" & chr == 1, pos < 121500000 &
           (mut_vaf < 0.1 | mut_vaf > 0.9) & total_depth > 4) %>%
  # assign haplotypes
  mutate(hap_A = case_when(mut_vaf > 0.9 ~ alt,
                           mut_vaf < 0.1 ~ ref),
         hap_B = case_when(mut_vaf > 0.9 ~ ref,
                           mut_vaf < 0.1 ~ alt))

# lift over snps
phased_snps <-
  phased_snps_grch37 %>%
  {GRanges(
    seqnames = paste0("chr", .$chr),
    ranges = IRanges(start = .$pos,
                     end = .$pos),
    ref = .$ref, alt = .$alt,
    hap_A = .$hap_A, hap_B = .$hap_B
  )} %>%
  liftOver(import.chain("../../reference/liftover/hg19ToHg38.over.chain")) %>%
  unlist() %>%
  {tibble::tibble(chr = as.double(gsub("chr", "", as.character(seqnames(.)))),
                  pos = start(.), ref = .$ref, alt = .$alt,
                  hap_A = .$hap_A, hap_B = .$hap_B)} %>%
  dplyr::distinct()

# save phased snps
muts_file <- file.path(wd, out_dir, "snps/PD63118_phased_snps.tsv")
phased_snps %>%
  readr::write_tsv(muts_file)

# make samplesheet
ss <-
  readr::read_csv("data/trencadis-seq/metadata/samplesheet.csv") %>%
  select(id, bam, kit) %>%
  mutate(cell_barcodes = file.path(wd, out_dir, kit, "barcodes.txt"),
         mutations = muts_file) %>%
  {split(., .$kit)}
purrr::walk2(names(ss), ss, function(kit_i, ss_i) {
  ss_i %>%
    select(-kit) %>%
    readr::write_csv(file.path(out_dir, kit_i, "samplesheet.csv"))
})
