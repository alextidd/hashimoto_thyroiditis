# libraries
library(Matrix)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

# dirs
data_dir <- "out/trencadis-seq/vartrix/"

# read in the cell barcodes
barcodes <- readLines("out/trencadis-seq/seurat/min_3_cells_min_500_genes/annotated_cell_barcodes.txt")

# read in the variants
snps <- readLines(file.path(data_dir, "variants.txt"))

# read in the sparse genotype matrices
mats <-
  c("alt", "ref") %>%
  purrr::set_names() %>%
  purrr::map(function(x) {
    mat <- as.data.frame(as.matrix(readMM(paste0(data_dir, x, ".mtx"))))
    colnames(mat) <- barcodes
    rownames(mat) <- snps
    mat
  })

# subset to snps with any reads
total_matrix <- mats$ref + mats$alt
total_matrix <- total_matrix[rowSums(total_matrix) > 0, ]
mats <- purrr::map(mats, ~ .x[rownames(total_matrix), ])
vaf_matrix <- mats$alt / total_matrix

# heatmap of vafs
pheatmap::pheatmap(
  vaf_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "VAFs of SNPs in VarTrix",
  filename = "out/trencadis-seq/vartrix/vafs_heatmap.png")

# histogram of vaf
tibble(vaf = unlist(as.vector(vaf_matrix))) %>%
  filter(is.finite(vaf)) %>%
  ggplot(aes(x = vaf)) +
  geom_histogram()

# phase the snps
geno_snps <-
  "out/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(id == "plate3_wellA2_dna_run49882", chr == 1, pos < 121500000,
         mut_vaf %in% c(0, 1), total_depth > 4) %>%
  # assign haplotypes
  mutate(haplotype = case_when(mut_vaf == 1 ~ "A", 
                               mut_vaf == 0 ~ "B"))

"out/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(chr == 1, pos < 121500000)

# lift over
geno_snps <-
  geno_snps %>%
  {GRanges(seqnames = paste0("chr", .$chr),
           ranges = IRanges(start = .$pos, end = .$pos),
           pos_grch37 = .$pos)} %>%
  liftOver(import.chain("../../reference/liftover/hg38ToHg19.over.chain")) %>%
  unlist() %>%
  as_tibble() %>%
  transmute(chr = gsub("chr", "", seqnames), start, end, pos_grch37) %>%
  inner_join(geno_snps %>% dplyr::rename(pos_grch37 = pos))