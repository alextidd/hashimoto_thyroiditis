# libraries
library(magrittr)

# get ordering of substitutions
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[",
                   rep(sub_vec, each = 16), "]",
                   rep(c("A", "C", "G", "T"), times = 4))

# load cosmic reference signatures
refs <-
  c("v2", "v3.4") %>%
  purrr::set_names() %>%
  purrr::map(function(v) {
    readr::read_tsv(paste0("../../reference/cosmic/COSMIC_", v, "_SBS_GRCh38.txt")) %>%
      dplyr::mutate(Type = factor(Type, levels = full_vec)) %>%
      dplyr::arrange(Type) %>%
      tibble::column_to_rownames("Type") %>%
      as.matrix()
  })

# use v3.4 + Signature_17 from v2
ref <- refs$v3.4
ref <- cbind(ref, SBS17 = refs$v2[, "Signature_17"])

# replace 0 with small value and normalise
ref[is.na(ref) | ref == 0] <- 0.00001
ref <- t(t(ref) / colSums(ref))

# add SBSblood from machado 2022
machado <-
  readr::read_tsv("data/signatures/machado_2022/S8_finalsignaturetable.tsv") %>%
  tidyr::pivot_longer(cols = -c("Signature")) %>%
  dplyr::mutate(Type = factor(paste0(substr(name, 1, 1), "[",
                                     substr(name, 2, 2), ">",
                                     substr(name, 6, 6), "]",
                                     substr(name, 3, 3)),
                              levels = full_vec)) %>%
  dplyr::arrange(Type) %>%
  tidyr::pivot_wider(names_from = "Signature", values_from = "value")
ref <- cbind(ref, machado_2022_SBSblood = machado$SBSblood, machado_2022_SignatureIg = machado$Signature.Ig)

# add artefact signature ScF from petljak 2019
petljak <-
  readr::read_tsv("data/signatures/petljak_2019/mmc1.tsv") %>%
  dplyr::mutate(
    Type = factor(paste0(substr(`Mutation Subtype`, 1, 1), "[",
                         `Mutation Type`, "]",
                         substr(`Mutation Subtype`, 3, 3)),
                  levels = full_vec)) %>%
  dplyr::arrange(Type)
ref <- cbind(ref, petljak_2019_ScF = petljak$`SBS sc_F`)

# add artefact signature ScB from lodato 2018
lodato <-
  "data/signatures/lodato_2018/Lodato2018_SignatureData_Aging.csv" %>%
  readr::read_csv() %>%
  dplyr::mutate(
    Type = factor(paste0(substr(`...1`, 1, 1), "[", `...2`, "]",
                         substr(`...1`, 3, 3)), levels = full_vec)) %>%
  dplyr::arrange(Type)
ref <- cbind(ref, lodato_2018_ScB = lodato$B)

# save ref
ref %>%
  write.table(
    file = "out/resolveome/signatures/cosmic_v3.4_ScF_ScB_SBSblood.tsv",
    sep = "\t", quote = FALSE)
