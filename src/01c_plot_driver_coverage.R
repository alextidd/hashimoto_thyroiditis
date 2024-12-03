# libraries
library(magrittr)

# samplesheet
ss <- readr::read_csv("out/driver_coverage/samplesheet.csv")
genes <- readLines("out/driver_coverage/genes.txt")

# get coverage
cov <-
  ss$id %>%
  purrr::set_names() %>%
  purrr::map(function(id) {
    genes %>%
      purrr::set_names() %>%
      purrr::map(function(gene) {
        cov_file <-
          file.path("out/driver_coverage/runs/", id, gene,
                    paste0(id, "_", gene, "_coverage.tsv"))
        if (file.exists(cov_file)) {
          readr::read_tsv(cov_file, col_names = c("chr", "pos", "cov"))
        } else {
          tibble::tibble()
        }
      }) %>%
      dplyr::bind_rows(.id = "gene")
  }) %>%
  dplyr::bind_rows(.id = "id") %>%
  dplyr::mutate(id_n = as.integer(gsub("plex", "", id)),
                id = forcats::fct_reorder(id, id_n))

# plot
ggsave("TNFRRSF14_coverage.png",
       cov %>%
         dplyr::filter(gene == "TNFRSF14") %>%
         ggplot(aes(x = pos, y = cov)) +
         geom_area() +
         theme_classic() +
         facet_grid(id ~ .),
       width = 10, height = 10, units = "in", dpi = 300)
