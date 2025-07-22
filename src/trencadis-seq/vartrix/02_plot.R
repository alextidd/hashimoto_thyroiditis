# libraries
library(Matrix)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)

# dirs
data_dir <- "out/trencadis-seq/vartrix/TX/"

# read in the cell barcodes
barcodes <- readLines(gzfile("/lustre/scratch125/casm/teams/team268/at31/projects/hashimoto_thyroiditis/data/vartrix/barcodes.tsv.gz"))

# read in the variants
snps <- readLines(file.path(data_dir, "variants.txt"))

# read in the sparse genotype matrices
mats <-
  c("alt", "ref") %>%
  purrr::set_names() %>%
  purrr::map(function(x) {
    mat <- readMM(file.path(data_dir, paste0(x, ".mtx")))
    dimnames(mat) <- list(snps, barcodes)
    mat
  })

# subset to annotated barcodes
seu <- readRDS("out/trencadis-seq/seurat/min_3_cells_min_500_genes/nf-trencadis-seq_thyroid_annotate_celltypes_cache/html/seu_annotated_cda8b5cd83213c7275398699b2692b7f.rds")
annot_barcodes <- colnames(seu)
mats <- purrr::map(mats, ~ .x[, annot_barcodes])

# subset to snps with any reads
total_matrix <- mats$ref + mats$alt
total_matrix <- total_matrix[rowSums(total_matrix) > 0, ]
mats <- purrr::map(mats, ~ .x[rownames(total_matrix), ])
vaf_matrix <- mats$alt / total_matrix

# get depths and vafs
p_dat <-
  tibble(vaf = unlist(as.vector(vaf_matrix)),
         total_depth = unlist(as.vector(total_matrix))) %>%
  filter(is.finite(vaf))

# histogram of vafs
p_dat %>%
  mutate(total_depth_bin = cut(
           total_depth, breaks = c(-Inf, 0, 5, 10, 20, 50, 100, Inf),
           labels = c("0", "<5", "5-10", "10–20", "20–50", "50–100", ">100"))) %>%
  ggplot(aes(x = vaf, fill = total_depth_bin)) +
  geom_histogram() +
  scale_fill_viridis_d()

# histogram of intermediate vafs
p_dat %>%
  filter(vaf > 0, vaf < 1) %>%
  mutate(total_depth_bin = cut(
           total_depth, breaks = c(-Inf, 0, 5, 10, 20, 50, 100, Inf),
           labels = c("0", "<5", "5-10", "10–20", "20–50", "50–100", ">100"))) %>%
  ggplot(aes(x = vaf, fill = total_depth_bin)) +
  geom_histogram() +
  scale_fill_viridis_d()

# heatmap of vafs
pheatmap::pheatmap(
  vaf_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "VAFs of SNPs in VarTrix",
  filename = "out/trencadis-seq/vartrix/vafs_heatmap.png")

# coverage
p_dat %>%
  count(cov) %>%
  filter(cov > 0) %>%
  ggplot(aes(x = cov, y = n)) +
  geom_col()
p_dat %>%
  count(cov) %>%
  ggplot(aes(x = cov, y = n)) +
  geom_col() +
  scale_y_log10()

# get celltype annotations
seu <- readRDS("out/trencadis-seq/seurat/min_3_cells_min_500_genes/nf-trencadis-seq_thyroid_annotate_celltypes_cache/html/seu_annotated_cda8b5cd83213c7275398699b2692b7f.rds")

# phase the snps
geno_snps <-
  "out/resolveome/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(id == "plate3_wellA2_dna_run49882", chr == 1, pos < 121500000,
         mut_vaf %in% c(0, 1), total_depth > 4) %>%
  # assign haplotypes
  mutate(haplotype = case_when(mut_vaf == 1 ~ "A", 
                               mut_vaf == 0 ~ "B"))

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
