#!/usr/bin/env Rscript

# the longest 1p LOH event is observed in cells plate3_wellA2, affecing the
# entire 1p arm, and is therefore expected to be homozygous all SNPs.
# we can therefore phase all SNPs in PD63118 chr1p.
# then, i will look for genic SNPs, as these will be most likely to be picked up
# in the snRNAseq in order to classify cells.

# libraries
library(magrittr)
library(dndscv)
library(dplyr)
library(ggplot2)
library(biomaRt)

# get common het snps
snps_1p <-
  readr::read_tsv("out/resolveome/nf-resolveome/muts_and_snps/PD63118/caveman_snps.tsv") %>%
  filter(chr == 1, pos < 121500000)

# get tnfsrf14 position
tnfrsf14_start <- 2487810

# prep biomart
mart <-
  useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
             host = "https://grch37.ensembl.org")

# get exonic snps
exons_1p <-
  getBM(attributes = c("external_gene_name", "ensembl_transcript_id",
                       "ensembl_exon_id", "chromosome_name", "exon_chrom_start",
                       "exon_chrom_end", "strand"),
        mart = mart) %>%
  transmute(gene = external_gene_name, chr = as.character(chromosome_name),
            exon = ensembl_exon_id,
            start = exon_chrom_start, end = exon_chrom_end) %>%
  filter(chr == "1", start < 121500000) %>%
  distinct() %>%
  group_by(gene, chr, exon) %>%
  reframe(pos = start:end) %>%
  ungroup() %>%
  distinct(gene, chr, pos)
snps_1p_in_exons <-
  snps_1p %>%
  inner_join(exons_1p)

# get genic snps
genes_1p <-
  getBM(attributes = c("external_gene_name", "chromosome_name",
                       "ensembl_gene_id",
                       "start_position", "end_position"),
        mart = mart) %>%
  transmute(gene = external_gene_name, chr = as.character(chromosome_name),
            start = start_position, end = end_position, ensembl_gene_id) %>%
  filter(chr == "1", start < 121500000) %>%
  distinct() %>%
  group_by(gene, chr, ensembl_gene_id) %>%
  reframe(pos = start:end) %>%
  ungroup() %>%
  distinct(gene, chr, pos)
snps_1p_in_genes <-
  snps_1p %>%
  inner_join(genes_1p)

# phase the snps
geno_snps <-
  "out/resolveome/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(id == "plate3_wellA2_dna_run49882", chr == 1, pos < 121500000,
         mut_vaf < 0.1 | mut_vaf > 0.9,
         total_depth > 4) %>%
  # assign haplotypes
  mutate(haplotype = ifelse(mut_vaf == 1, "A", "B"))

geno_snps %>%
  arrange(-total_depth) %>%
  mutate(row = row_number()) %>%
  ggplot(aes(x = row, y = total_depth)) +
  geom_line()

# look at baf of snps across loh region
p_dat <-
  "out/resolveome/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(id == "plate3_wellA2_dna_run49882", chr == 1, pos < 121500000)
p_dat %>%
  mutate(mut_vaf = ifelse(total_depth == 0, 0, mut_vaf),
         mut_baf = 1 - mut_vaf,
         total_depth_bin = cut(
           total_depth, breaks = c(-Inf, 0, 10, 20, 50, 100, Inf),
           labels = c("0", "<10", "10–20", "20–50", "50–100", ">100"))) %>%
  ggplot(aes(x = pos, colour = total_depth_bin)) +
  geom_point(aes(y = mut_vaf)) +
  geom_point(aes(y = mut_baf)) +
  scale_colour_viridis_d()
