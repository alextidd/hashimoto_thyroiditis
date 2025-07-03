# libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)

# phase the snps and lift over
phased_snps_grch37 <-
  "out/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  readr::read_tsv() %>%
  filter(id == "plate3_wellA2_dna_run49882" & chr == 1, pos < 121500000 &
           (mut_vaf < 0.1 | mut_vaf > 0.9) & total_depth > 4) %>%
  # assign haplotypes
  mutate(hap_A = case_when(mut_vaf > 0.9 ~ alt,
                           mut_vaf < 0.1 ~ ref),
         hap_B = case_when(mut_vaf > 0.9 ~ ref,
                           mut_vaf < 0.1 ~ alt))

# lift over snps
phased_snps <-
  phased_snps_grch37 %>%
  {GRanges(
    seqnames = paste0("chr", .$chr),
    ranges = IRanges(start = .$pos,
                     end = .$pos),
    hap_A = .$hap_A, hap_B = .$hap_B
  )} %>%
  liftOver(import.chain("../../reference/liftover/hg19ToHg38.over.chain")) %>%
  unlist() %>%
  {tibble::tibble(chr = as.double(gsub("chr", "", as.character(seqnames(.)))),
                  pos = start(.),
                  hap_A = .$hap_A, hap_B = .$hap_B)} %>%
  dplyr::distinct()

# get snps
snps <-
  readr::read_tsv("out/trencadis-seq/nf-sc-geno/TX_WT_PD63118/TX_WT_PD63118_genotyped_mutations_per_cell.tsv")

snps %>%
  left_join(geno_snps %>% select(chr, pos, ref, alt, haplotype)) %>%
  filter(!is.na(haplotype))

snps %>%
  arrange(barcode, pos) %>%
  ggplot(aes(x = pos, y = barcode, fill = alt_vaf)) +
  geom_tile()

snps %>%
  filter(alt_vaf > 0, alt_vaf < 1) %>%
  ggplot(aes(x = alt_vaf)) +
  geom_histogram()
