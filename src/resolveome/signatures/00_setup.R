# libraries
library(magrittr)

# get ordering of substitutions
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
full_vec <- paste0(rep(c("A", "C", "G", "T"), each = 4), "[",
                   rep(sub_vec, each = 16), "]",
                   rep(c("A", "C", "G", "T"), times = 4))

# load cosmic reference signatures
ref <-
  readr::read_tsv("../../reference/cosmic/COSMIC_v3.4_SBS_GRCh38.txt") %>%
  dplyr::mutate(Type = factor(Type, levels = full_vec)) %>%
  dplyr::arrange(Type) %>%
  tibble::column_to_rownames("Type") %>%
  as.matrix()

# replace 0 with small value and normalise
ref[is.na(ref) | ref == 0] <- 0.00001
ref <- t(t(ref) / colSums(ref))

# add artefact signature ScF from petljak2019
petljak <-
  readr::read_tsv("data/signatures/petljak_2019/mmc1.tsv")
ref <- cbind(ref, petljak_2019_ScF = petljak$`SBS sc_F`)

# add artefact signature ScB from lodato2018
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
  write.table(file = "out/resolveome/signatures/cosmic_v3.4_ScF_ScB.txt",
              sep = "\t", quote = FALSE)
