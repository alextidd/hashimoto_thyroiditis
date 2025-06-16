#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(rtracklayer)
library(GenomicRanges)
library(biomaRt)

# function: get the reverse complement of a single DNA sequence
rev_comp <- function(seq) {
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  nucleotides <- unlist(strsplit(seq, split = ""))
  complemented <- complement[nucleotides]
  reverse_complemented <- rev(complemented)
  return(paste(reverse_complemented, collapse = ""))
}

# dirs
wd <- getwd()
data_dir <- file.path(wd, "out/thyroid/seurat/min_3_cells_min_500_genes/nf-trencadis-seq_thyroid_annotate_celltypes_cache/html/")
out_dir <- file.path(wd, "out/thyroid/nf-trencadis-seq")

# load seu object
seu <-
  list.files(data_dir, pattern = "seu_annotated", full.names = TRUE) %>%
  readRDS()

# save celltypes with 10X barcodes
cts_tx <-
  seu@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  dplyr::transmute(celltype, barcode)
cts_tx %>%
  readr::write_csv(file.path(out_dir, "barcode_TX_celltypes_annotations.csv"))
colnames(seu) %>%
  writeLines(file.path(out_dir, "barcodes_TX.txt"))

# save celltypes with PB barcodes
cts_pb_file <- file.path(out_dir, "barcode_PB_celltypes_annotations.csv")
barcodes_pb_file <- file.path(out_dir, "barcodes_PB.txt")
cts_pb <-
  cts_tx %>%
  dplyr::rowwise() %>%
  dplyr::mutate(barcode = rev_comp(gsub("-1$", "", barcode)))
cts_pb %>%
  readr::write_csv(cts_pb_file)
cts_pb$barcode %>%
  writeLines(barcodes_pb_file)

# get driver gene coords
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
driver_genes <-
  getBM(attributes = c("hgnc_symbol", "chromosome_name", "start_position",
                       "end_position", "strand"),
        filters = "hgnc_symbol",
        values = readLines("data/thyroid/driver_genes/driver_genes.txt"),
        mart = ensembl) %>%
  tibble::as_tibble() %>%
  dplyr::transmute(gene = hgnc_symbol, chr = chromosome_name,
                   start = start_position, end = end_position) %>%
  # only the genes on the non-alternative chromosomes
  dplyr::filter(!grepl("^H", chr))
driver_genes_pos <-
  driver_genes %>%
  dplyr::group_by(gene, chr, start, end) %>%
  dplyr::reframe(pos = start:end) %>%
  dplyr::select(chr, pos, gene)

# get mutations
muts_file <- file.path(out_dir, "mutations_in.tsv")
muts_grch37 <-
  readr::read_tsv("../resolveome/out/nf-resolveome/PD63118/mutations.tsv") %>%
  dplyr::rename(pos_grch37 = pos) %>%
  dplyr::filter(source == "nanoseq_mutations")

# perform the liftOver
muts <-
  muts_grch37 %>%
  {GRanges(
    seqnames = paste0("chr", .$chr),
    ranges = IRanges(start = .$pos_grch37,
                     end = .$pos_grch37),
    pos_grch37 = .$pos_grch37
  )} %>%
  liftOver(import.chain("../reference/ucsc/hg19ToHg38.over.chain")) %>%
  unlist() %>%
  {tibble::tibble(chr = gsub("chr", "", as.character(seqnames(.))),
                  pos = start(.),
                  pos_grch37 = .$pos_grch37)} %>%
  dplyr::distinct() %>%
  dplyr::right_join(muts_grch37) %>%
  # get only mutations in drivers
  dplyr::inner_join(driver_genes_pos) %>%
  dplyr::mutate(chr = paste0("chr", chr))

# write the mutations
muts %>%
  readr::write_tsv(muts_file)

# write the samplesheet (id, bam, mutations, celltypes, cell_barcodes)
readr::read_csv("data/thyroid/metadata/samplesheet.csv") %>%
  dplyr::transmute(id, bam,
                   mutations = muts_file,
                   celltypes = cts_pb_file,
                   cell_barcodes = barcodes_pb_file) %>%
  readr::write_csv("out/thyroid/nf-trencadis-seq/samplesheet.csv")