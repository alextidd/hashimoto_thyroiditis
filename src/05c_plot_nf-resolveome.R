# libraries
library(magrittr)
library(ggplot2)
p <- list()

# samplesheet
ss <-
  system("ls out/nf-resolveome/*/samplesheet.tsv", intern = TRUE) %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows()
genes_pos <-
  readr::read_tsv("out/driver_coverage/genes/genes.bed",
                  col_types = list(chr = readr::col_character())) %>%
  dplyr::mutate(start2 = start, end2 = end) %>%
  dplyr::group_by(chr, start, end, gene, strand) %>%
  dplyr::reframe(pos = start2:end2)

# get summary
summ <-
  ss %>%
  purrr::pmap(function(donor_id, id, run, plex_n, ...) {
    paste0("out/nf-resolveome/", run, "/", donor_id, "/", id,
            "/mosdepth/", id, ".mosdepth.summary.txt") %>%
      readr::read_tsv() %>%
      dplyr::mutate(id = id, run = run, donor_id = donor_id, plex_n = plex_n) %>%
      # remove extra chromosomes
      dplyr::filter(!grepl("^GL", chrom), !grepl("Un", chrom),
                    !grepl("random|HLA", chrom), !grepl("alt$", chrom))
  }) %>%
  dplyr::bind_rows() 

# plot summary
p[["summ"]] <-
  summ %>%
  dplyr::filter(chrom %in% c("total", "total_region")) %>%
  dplyr::mutate(lvl = ifelse(grepl("region$", chrom), "region", "global"),
                chrom = gsub("_.*", "", chrom)) %>%
  ggplot(aes(x = tidytext::reorder_within(id, -mean, run), y = mean, fill = lvl)) +
  geom_col(position = "dodge") +
  theme_minimal() +
  labs(title = "Mean coverage",
       x = "Plex",
       y = "Mean coverage") +
  ggh4x::facet_grid2(~ run, scales ="free_x", space="free_x") +
  theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
  tidytext::scale_x_reordered()

# get cov
cov <-
  ss %>%
  purrr::pmap(function(donor_id, id, run, plex_n, ...) {
    c("region", "global") %>%
      purrr::set_names() %>%
      purrr::map(function(lvl) {
        readr::read_tsv(paste0("out/nf-resolveome/", run, "/", donor_id, "/", id,
                               "/mosdepth/", id, ".mosdepth.", lvl, ".dist.txt"),
                        col_names = c("chr", "cov", "prop")) %>%
          # remove extra chromosomes
          dplyr::filter(!grepl("^GL", chr), !grepl("Un", chr),
                        !grepl("random|HLA", chr), !grepl("alt$", chr),
                        prop >= 0.01)
      }) %>%
      dplyr::bind_rows(.id = "lvl") %>%
      dplyr::mutate(id = id, run = run, donor_id = donor_id, plex_n = plex_n)
  }) %>%
  dplyr::bind_rows()

# plot cov
p[["cov"]] <-
  cov %>%
  dplyr::filter(chr == "total") %>%
  ggplot(aes(x = cov, y = prop, color = id, group = id)) +
  geom_line() +
  labs(title = paste("Coverage distribution"),
       x = "Coverage",
       y = "Proportion of bases at coverage") +
  theme_minimal() +
  facet_grid(run ~ lvl) +
  theme(legend.position = "none")

# get mutations
muts_and_snps <-
  readr::read_tsv("out/nf-resolveome/mutations.tsv") %>%
  dplyr::distinct(source, pdid, chr, pos, ref, mut)
geno <-
  ss %>%
  purrr::pmap(function(donor_id, id, run, plex_n, ...) {
    paste0("out/nf-resolveome/", run, "/", donor_id, "/", id,
            "/genotyping/", id, "_genotyped_mutations.tsv") %>%
      readr::read_tsv(col_types = list(chr = readr::col_character())) %>%
      dplyr::mutate(id = id, run = run, donor_id = donor_id, plex_n = plex_n)
  }) %>%
  dplyr::bind_rows()

# plot VAF distribution
mut_depth_bins <- c("0", "1", ">1", ">5", ">10", ">50")
p[["vaf_dist"]] <-
  geno %>%
  dplyr::left_join(muts_and_snps) %>%
  dplyr::mutate(mut_vaf = mut_depth / total_depth,
                mut_depth_bin = dplyr::case_when(mut_depth > 50 ~ ">50",
                                                  mut_depth > 10 ~ ">10",
                                                  mut_depth > 5 ~ ">5",
                                                  mut_depth > 1 ~ ">1",
                                                  mut_depth == 1 ~ "1",
                                                  mut_depth == 0 ~ "0") %>%
                  factor(levels = mut_depth_bins)) %>%
  ggplot(aes(x = mut_vaf, fill = mut_depth_bin, colour = mut_depth_bin)) +
  geom_histogram(bins = 100) +
  ggtitle("Distribution of VAFs in PTA samples") +
  facet_grid(source ~ run, scales = "free_y") +
  theme_minimal() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)

# plot nonzero VAF distribution
p[["vaf_dist_nonzero"]] <-
  geno %>%
  dplyr::filter(mut_depth > 0) %>%
  dplyr::left_join(muts_and_snps) %>%
  dplyr::mutate(mut_vaf = mut_depth / total_depth,
                mut_depth_bin = dplyr::case_when(mut_depth > 50 ~ ">50",
                                                  mut_depth > 10 ~ ">10",
                                                  mut_depth > 5 ~ ">5",
                                                  mut_depth > 1 ~ ">1",
                                                  mut_depth == 1 ~ "1") %>%
                  factor(levels = mut_depth_bins)) %>%
  ggplot(aes(x = mut_vaf, fill = mut_depth_bin, colour = mut_depth_bin)) +
  geom_histogram(bins = 100) +
  ggtitle("Distribution of non-zero VAFs in PTA samples") +
  ggh4x::facet_grid2(run ~ source, scales = "free_y", independent = "y") +
  theme_minimal() +
  viridis::scale_fill_viridis(discrete = TRUE) +
  viridis::scale_color_viridis(discrete = TRUE)

# plot VAF distribution with increasing mut depth
seq(1,10) %>%
  purrr::walk(function(min_mut_depth) {
    p[[paste0("vaf_dist_>=", min_mut_depth)]] <<-
      geno %>%
      dplyr::filter(mut_depth >= min_mut_depth) %>%
      dplyr::left_join(muts_and_snps) %>%
      dplyr::mutate(mut_vaf = mut_depth / total_depth) %>%
      ggplot(aes(x = mut_vaf)) +
      geom_histogram(bins = 100) +
      ggtitle(paste("Distribution of VAFs in PTA samples with mut depth >= ", min_mut_depth)) +
      facet_grid(source ~ run, scales = "free_y") +
      theme_minimal()
  })

# plot mutations
pdf("test.pdf")
purrr::walk(p, print)
dev.off()

