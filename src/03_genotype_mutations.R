# libraries
library(magrittr)

# get mutations from exome and targeted nanoseq
nanoseq_dir <- "/lustre/scratch126/casm/team268im/al28/targeted_nanoseq/"
nanoseq_muts <-
  c("exome", "immune") %>%
  purrr::set_names() %>%
  purrr::map(function(x) {
    system(paste0("ls ", nanoseq_dir, "plate_062*", x,
           "*/Analysis/plate_062*final_muts.tsv"), intern = TRUE) %>%
      readr::read_tsv() %>%
      dplyr::mutate(pdid = stringr::str_sub(sampleID, 1, 7))
  }) %>%
  dplyr::bind_rows(.id = "seq_type") %>%
  dplyr::distinct(pdid, sampleID, chr, pos, ref, mut, seq_type, times_called)

# get mutations in the florid case
florid_muts <-
  nanoseq_muts %>%
  dplyr::filter(pdid == "PD63118")
  
# run dndscv to get annotate mutations
dndsout <- dndscv::dndscv(florid_muts)
sel_cv <- dndsout$sel_cv
print(head(sel_cv), digits = 3)
signif_genes <- sel_cv[sel_cv$qglobal_cv < 0.1, c("gene_name","qglobal_cv")]
rownames(signif_genes) <- NULL
print(signif_genes)

# genotype these sites in the PTA bams
ss <- readr::read_csv("out/caveman/samplesheet.csv")
annotmuts <-
  dplyr::left_join(florid_muts, dndsout$annotmuts) %>%
  dplyr::arrange(-times_called) %>%
  dplyr::filter(!is.na(gene)) %>%
  dplyr::select(pdid, chr, pos, ref, mut, gene, times_called)

# genotype the site
bam <- ss[1, ]$tumour_bam
site <- annotmuts[1, ]

pta_muts <-
  annotmuts %>%
  dplyr::filter(nchar(ref) == 1, nchar(mut) == 1) %>%
  purrr::pmap(function(pdid, chr, pos, ref, mut, gene, times_called) {
    split(ss$tumour_bam, ss$sample_id) %>%
      purrr::map(function(bam) {
        calls <- deepSNV::bam2R(bam, chr, pos, pos)
        tibble::tibble(chr = chr, pos = pos, ref = ref, mut = mut, gene = gene,
                      times_called = times_called, calls = calls) %>%
          dplyr::mutate(n_ref = calls[, ref] + calls[, tolower(ref)],
                        n_mut = calls[, mut] + calls[, tolower(mut)])
      }) %>%
      dplyr::bind_rows(.id = "sample_id")
  }) %>%
  dplyr::bind_rows()

pta_muts %>%
  dplyr::mutate(vaf = n_mut / (n_ref + n_mut)) %>%
  dplyr::select(vaf, n_ref, n_mut, everything()) %>%
  dplyr::arrange(-vaf) %>%
  dplyr::filter(vaf > 0) %>%
  dplyr::count(chr, pos, ref, mut, gene) %>%
  dplyr::arrange(-n)
