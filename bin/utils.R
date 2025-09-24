# function: plot VAF heatmap
plot_vaf_heatmap <- function(p_dat, p_title = "", annotations = c(), show_rownames = TRUE) {
  library(magrittr)

  # prepare data
  p_dat2 <-
    p_dat %>%
    dplyr::ungroup() %>%
    # only keep mutations with VAF > 0
    dplyr::filter(alt_vaf > 0) %>%
    dplyr::mutate(mut_id = paste(gene, chr, pos, ref, alt, sep = "_")) %>%
    # count number of cells with each mutation
    dplyr::add_count(mut_id, name = "n_cells_w_mut")

  # reshape data for heatmap
  heatmap_data <-
    p_dat2 %>%
    reshape2::dcast(mut_id + n_cells_w_mut ~ cell_id, value.var = "alt_vaf")
  rownames(heatmap_data) <- heatmap_data$mut_id
  heatmap_matrix <-
    heatmap_data %>%
    dplyr::select(-mut_id, -n_cells_w_mut) %>%
    as.matrix()
  n_total_muts <- dplyr::n_distinct(p_dat2$mut_id)
  n_total_cells <- dplyr::n_distinct(p_dat2$cell_id)

  # calculate the number of mutations per column (id) + other annots
  annotations_col <-
    p_dat2 %>%
    dplyr::count(cell_id, !!!syms(annotations)) %>%
    tibble::column_to_rownames("cell_id")

  # replace NAs with 0
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # prepare annotations for columns
  annotation_row <- data.frame(n_cells_w_mut = heatmap_data$n_cells_w_mut)
  rownames(annotation_row) <- heatmap_data$mut_id

  # plot
  pheatmap::pheatmap(
    mat = heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    annotation_row = annotation_row,
    annotation_col = annotations_col,
    color = rev(viridis::magma(100)),
    main = paste(p_title, "\nVAF heatmap\n",
                 n_total_cells, "cells,", n_total_muts, "mutations"))
}

# function: plot vaf distribution with depth information
plot_vaf_dist <- function(p_dat, p_title = "VAF distribution") {
  n_unique_muts <-
    p_dat %>%
    dplyr::distinct(chr, pos, ref, alt) %>%
    nrow() %>%
    format(big.mark = ",", scientific = FALSE)
  n_total_muts <-
    p_dat %>%
    nrow() %>%
    format(big.mark = ",", scientific = FALSE)

  p_dat_binned <-
    p_dat %>%
    dplyr::mutate(
      alt_vaf = ifelse(total_depth == 0, 0, alt_vaf),
      total_depth_bin = cut(
        total_depth, breaks = c(-Inf, 0, 10, 20, 50, 100, Inf),
        labels = c("0", "<10", "10–20", "20–50", "50–100", ">100")),
      alt_vaf_bin = cut(alt_vaf, breaks = seq(0, 1, length.out = 30 + 1),
                        include.lowest = TRUE),
      alt_vaf_min = stringr::str_extract(alt_vaf_bin, "[\\d\\.]+") %>% as.numeric(),
      alt_vaf_high = as.numeric(stringr::str_extract(alt_vaf_bin, "(?<=,)\\s*[\\d\\.]+")),
      alt_vaf_mid = (alt_vaf_min + alt_vaf_high) / 2) %>%
    dplyr::count(alt_vaf_bin, alt_vaf_mid, total_depth_bin)
  p <-
    p_dat_binned %>%
    ggplot2::ggplot(ggplot2::aes(x = alt_vaf_mid, y = n, fill = total_depth_bin)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = p_title,
                  subtitle = paste(n_total_muts, "total mutations,",
                                   n_unique_muts, "unique mutations"),
         x = "VAF bin", y = "n mutations") +
    ggplot2::lims(x = c(0, 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  print(p)
}

# function: get the reverse complement of vector of barcodes
rev_comp <- function(seq_vec) {
  seq_vec %>%
    gsub("-1$", "", .) %>%
    chartr("ATGC", "TACG", .) %>%
    stringi::stri_reverse()
}

# function: plot the trinucleotide context of the muts
plot_mut_trinucs <- function(p_dat, p_title = "") {
  sub_colours <-
    c("C>A" = "dodgerblue", "C>G" = "black", "C>T" = "red",
      "T>A" = "grey70", "T>C" = "olivedrab3", "T>G" = "plum2")
  nucs <- c("T", "C", "G", "A")
  py_nucs <- c("C", "T")
  all_subs <-
    expand.grid(do.call(paste0, expand.grid(nucs, py_nucs, nucs)), nucs) %>%
    dplyr::transmute(trinuc_ref_py = Var1, ref = substr(Var1, 2, 2), alt = Var2,
                     sub_py = paste0(ref, ">", alt)) %>%
    dplyr::filter(alt != ref)

  p_dat %>%
    dplyr::mutate(n_cells = dplyr::n_distinct(cell_id)) %>%
    dplyr::add_count(mut_id, name = "n_cells_w_mut") %>%
    dplyr::count(sub_py, trinuc_ref_py) %>%
    dplyr::right_join(all_subs) %>%
    dplyr::mutate(n = ifelse(is.na(n), 0, n)) %>%
    ggplot(aes(x = trinuc_ref_py, y = n, fill = sub_py)) +
    geom_col() +
    facet_grid(~ sub_py, scales = "free_x") +
    guides(x = guide_axis(angle = 90)) +
    theme_minimal() +
    theme(axis.text.x = element_text(family = "mono"),
          legend.position = "none") +
    scale_fill_manual(values = sub_colours) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = p_title,
         subtitle = paste(format(nrow(p_dat), big.mark = ","), "mutations"),
         x = "", y = "number of mutations")
}