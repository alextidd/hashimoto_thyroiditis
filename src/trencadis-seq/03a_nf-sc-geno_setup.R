#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)

# dirs
wd <- getwd()
data_dir <- file.path(wd, "out/trencadis-seq/seurat/min_3_cells_min_500_genes/nf-trencadis-seq_thyroid_annotate_celltypes_cache/html/")
out_dir <- file.path(wd, "out/trencadis-seq/nf-sc-geno/")

# save annotated barcodes (from seu object)
barcodes_file <- file.path(out_dir, "barcodes.txt")
seu <-
  list.files(data_dir, pattern = "seu_annotated", full.names = TRUE) %>%
  readRDS()
colnames(seu) %>% writeLines(barcodes_file)

# get snps
snps_grch37 <-
  "out/nf-resolveome/muts_and_snps/PD63118/caveman_snps.tsv" %>%
  readr::read_tsv() %>%
  dplyr::filter(chr == 1, pos < 121500000)

# get coding snps
snps_grch37_coding <-
  snps_grch37 %>%
  dplyr::transmute(sampleID = "", chr, pos, ref, alt) %>%
  dplyr::distinct() %>%
  dndscv::dndscv(max_muts_per_gene_per_sample = Inf, outp = 1,
                 max_coding_muts_per_sample = Inf) %>%
  {.$annotmuts} %>%
  dplyr::select(chr, pos, ref, alt = mut)

# lift over snps
snps <-
  snps_grch37_coding %>%
  {GRanges(
    seqnames = paste0("chr", .$chr),
    ranges = IRanges(start = .$pos,
                     end = .$pos),
    ref = .$ref, alt = .$alt
  )} %>%
  liftOver(import.chain("../../reference/liftover/hg19ToHg38.over.chain")) %>%
  unlist() %>%
  {tibble::tibble(chr = gsub("chr", "", as.character(seqnames(.))),
                  pos = start(.),
                  ref = .$ref, alt = .$alt)} %>%
  dplyr::distinct()

# write snps
snps_file <- file.path(out_dir, "caveman_snps_1p.tsv")
snps %>% readr::write_tsv(snps_file)

# write the samplesheet (id, bam, mutations, cell_barcodes)
tibble::tibble(
  id = "TX_WT_PD63118",
  bam = "/lustre/scratch125/casm/teams/team268/at31//resolveome/data/vartrix/possorted_genome_bam.bam",
  mutations = snps_file,
  cell_barcodes = barcodes_file) %>%
  readr::write_csv("out/trencadis-seq/nf-sc-geno/samplesheet.csv")
