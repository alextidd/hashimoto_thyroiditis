# libraries
library(magrittr)

# dirs
out_dir <- "out/resolveome/basejumper/find_muts_near_snps"
seq_dir <- "out/resolveome/sequoia/"

# load muts on tree
muts <-
  c("NV", "NR") %>%
  purrr::set_names() %>%
  purrr::map(function(i) {
    paste0(seq_dir, "Patient_both_", i, "_tree_all.txt") %>%
      read.table() %>%
      tibble::as_tibble(rownames = "mut_id") %>%
      tidyr::separate_wider_delim(cols = mut_id, delim = "_",
                                  names = c("chr", "pos", "ref", "alt")) %>%
      tidyr::pivot_longer(cols = -c(chr, pos, ref, alt), names_to = "id")
  }) %>%
  dplyr::bind_rows(.id = "name") %>%
  tidyr::pivot_wider() %>%
  dplyr::rename(alt_depth = NV, total_depth = NR) %>%
  dplyr::mutate(pos = as.numeric(pos)) %>%
  # filter on depth
  dplyr::filter(alt_depth > 0, total_depth >= 5) %>%
  dplyr::mutate(alt_vaf = alt_depth / total_depth)

# load muts near snps
muts_near_snps <-
  file.path(out_dir, "muts_near_snps.tsv") %>%
  readr::read_tsv()

# pinpoint muts near snps to their cells
muts_near_snps %>%
  dplyr::inner_join(muts, by = c("chr", "mut_pos" = "pos")) %>%
  dplyr::mutate(distance = abs(mut_pos - snp_pos)) %>%
  dplyr::arrange(distance)
