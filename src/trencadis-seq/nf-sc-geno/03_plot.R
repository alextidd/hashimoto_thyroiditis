# libraries
library(magrittr)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(GenomicRanges)

# phase the snps and lift over
phased_snps_grch37 <-
  readRDS("out/phase_snps/phase_snps_cache/phased_loh_genotyped_chr1_snps.rds")
  # "out/resolveome/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
  # readr::read_tsv() %>%
  # filter(id == "plate3_wellA2_dna_run49882" & chr == 1, pos < 121500000 &
  #          (mut_vaf < 0.1 | mut_vaf > 0.9) & total_depth > 4) %>%
  # # assign haplotypes
  # mutate(hap_A = case_when(mut_vaf > 0.9 ~ alt,
  #                          mut_vaf < 0.1 ~ ref),
  #        hap_B = case_when(mut_vaf > 0.9 ~ ref,
  #                          mut_vaf < 0.1 ~ alt))

# get phased snps from plate3_wellA2
phased_snps_grch37 <-
  readRDS("out/phase_snps/phase_snps_cache/plate3_wellA2_phased_snps.rds") %>%
  filter(!is.na(hap_A_type)) %>%
  mutate(hap_A = case_when(hap_A_type == "ref" ~ ref,
                           hap_A_type == "alt" ~ alt),
         hap_B = case_when(hap_B_type == "ref" ~ ref,
                           hap_B_type == "alt" ~ alt))

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
  readr::read_tsv("out/trencadis-seq/nf-sc-geno/TX_WT_PD63118/TX_WT_PD63118_genotyped_mutations_per_cell.tsv") %>%
  dplyr::filter(chr == 1, pos < 121500000)

# join snps with phased snps
snps %>%
  dplyr::left_join(phased_snps) %>%
  dplyr::mutate(haplotype)

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



# phase the snps
phased_snps_grch37 <-
  "out/resolveome/nf-resolveome/dna/PD63118/genotyping/snps/PD63118_genotyped_snps.tsv" %>%
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
    ref = .$ref, alt = .$alt,
    hap_A = .$hap_A, hap_B = .$hap_B
  )} %>%
  liftOver(import.chain("../../reference/liftover/hg19ToHg38.over.chain")) %>%
  unlist() %>%
  {tibble::tibble(chr = as.double(gsub("chr", "", as.character(seqnames(.)))),
                  pos = start(.), ref = .$ref, alt = .$alt,
                  hap_A = .$hap_A, hap_B = .$hap_B)} %>%
  dplyr::distinct()

# make samplesheet
ss <-
  readr::read_csv("data/trencadis-seq/metadata/samplesheet.csv") %>%
  {split(select(., id, bam), .$kit)}
ss$PB$barcodes <- barcodes_pb
ss$TX$barcodes <- barcodes_tx
# TX_WT_PD63118
# /lustre/scratch125/casm/teams/team268/at31/projects/hashimoto_thyroiditis/data/vartrix/possorted_genome_bam.bam
# /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/out/trencadis-seq/nf-sc-geno/caveman_snps_1p.tsv
# /nfs/casm/team268im/at31/projects/hashimoto_thyroiditis/out/trencadis-seq/nf-sc-geno/barcodes.txt

# plot phasing
p_size <- 0.8
p_alpha <- 0.2
plate3_wellA2 %>%
  inner_join(phased_snps_grch37 %>% mutate(chr = as.character(chr))) %>%
  mutate(haplotype = case_when(mut_vaf < 0.1 & ref == hap_A ~ "A",
                               mut_vaf < 0.1 & ref == hap_B ~ "B",
                               mut_vaf > 0.9 & alt == hap_A ~ "A",
                               mut_vaf > 0.9 & alt == hap_B ~ "B"),
         mut_baf = 1 - mut_vaf) %>%
  #filter(!is.na(haplotype)) %>%
  tidyr::pivot_longer(cols = c(mut_vaf, mut_baf), names_to = "vaf_type") %>%
  ggplot(aes(x = pos, y = value, colour = haplotype)) +
  # add line at vaf = 0.5
  geom_hline(yintercept = 0.5, colour = "red") +
  # add baf points
  geom_point(size = p_size, alpha = p_alpha) +
  # add baf bands
  scale_x_continuous(expand = c(0, 0)) +
  ggh4x::facet_grid2(. ~ chr, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(panel.spacing = unit(0, "lines"),
        panel.border = element_rect(color = "grey", fill = NA,
                                    linewidth = 0),
        strip.background = element_rect(color = "grey", fill = NA,
                                        linewidth = 0, linetype = "solid"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "none") +
  scale_colour_manual(values = c("A" = "#ff7be7", "B" = "blue"), na.value = "black")