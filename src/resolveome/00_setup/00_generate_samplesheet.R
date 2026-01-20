#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
data_dir <- file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis/data/bams/")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# run_id-to-plate conversion
# "49686" = 1
# "49882" = 3, "49900" = 3, "49901" = 3, "50072" = 3
# "50227" = 10, "50382" = 10
# "50367" = 11
# "51611" = 384.3

# read in sequencescape pool files
seqscape <-
  list.files("data/resolveome/sequencescape/", pattern = "tsv$",
             full.names = TRUE) %>%
  {purrr::set_names(., basename(.) %>% tools::file_path_sans_ext())} %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows(.id = "id") %>%
  janitor::clean_names() %>%
  tidyr::separate_wider_delim("id", delim = "_", names = c("run_id", "lane"))

# read in manifests
manifests_og <-
  "data/resolveome/manifest/" %>%
  list.files(pattern = ".xlsx$", full.names = TRUE) %>%
  {setNames(., basename(.) %>% tools::file_path_sans_ext())} %>%
  purrr::map(function(file) {
    file %>%
      readxl::read_excel(skip = 8) %>%
      janitor::clean_names()
  })

# fix manifests
# - supplier_sample_name standardisation - {donor_id}_P{plate_number}_{seq_type}_{well}
# - filter out unwanted samples
manifests <- manifests_og

# run49686_lane3 (7761stdy_manifest_25650_031024)
# filter samples and fix supplier sample names, adding plate number
manifests[["7761stdy_manifest_25650_031024"]] <-
  manifests_og[["7761stdy_manifest_25650_031024"]] %>%
  dplyr::filter(grepl("^Hashimoto", supplier_sample_name)) %>%
  dplyr::mutate(
    donor_id = "PD63118",
    plate = 1,
    seq_type = "dnahyb",
    well = gsub(".* - ", "", supplier_sample_name))

# extract metadata for the standardised manifests
manifests <-
  manifests_og[c("7761stdy_manifest_26617_Hyb_manifest_run50227_lane7_corrected",
                 "7761stdy_manifest_26618_RNA_manifest_run50227_lane8_corrected",
                 "7761stdy_manifest_26891_150425_R_D",
                 "7761stdy_manifest_26916_230425_R_D")] %>%
  purrr::map(function(m) {
    m %>%
    tidyr::separate_wider_delim(
      supplier_sample_name, delim = "_",
      names = c("donor_id", "plate", "seq_type", "well")) %>%
    dplyr::mutate(
      plate = as.numeric(gsub("P", "", plate)),            
      seq_type = ifelse(seq_type == "HYB", "DNAHYB", seq_type) %>% tolower())
  })

# run49686_lane4-5 (7761stdy_manifest_25654_041024)
# fix supplier sample names, adding plate number
manifests[["7761stdy_manifest_25654_041024"]] <-
  manifests_og[["7761stdy_manifest_25654_041024"]] %>%
  dplyr::mutate(
    donor_id = "PD63118",
    plate = 1,
    seq_type = "dna",
    well = gsub(".* - ", "", supplier_sample_name))

# run49882_lane1-8 (7761stdy_manifest_26013_211124_DNA)
# fix supplier sample names, adding plate number
manifests[["7761stdy_manifest_26013_211124_DNA"]] <-
  manifests_og[["7761stdy_manifest_26013_211124_DNA"]] %>%
  dplyr::mutate(
    donor_id = "PD63118",
    plate = 3,
    seq_type = "dna",
    well = gsub(".*_", "", supplier_sample_name))

# run49900_lane2 (7761stdy_manifest_26014_211124_Hyb)
# fix supplier sample names, adding plate number
manifests[["7761stdy_manifest_26014_211124_Hyb"]] <-
  manifests_og[["7761stdy_manifest_26014_211124_Hyb"]] %>%
  dplyr::mutate(
    donor_id = "PD63118",
    plate = 3,
    seq_type = "dnahyb",
    well = gsub(".*_", "", supplier_sample_name))

# run49901_lane8 (7894stdy_manifest_26015_211124_RNA)
# remove control (!= "1 cell") and low-yield (quant < 2) wells
# also, there is a typo in the manifest, where 7894STDY15290419 and
# 7894STDY15290420 are both assigned as PD63118_RNA_F10, while no sample is
# assigned to PD63118_RNA_G10. i will fix this by assigning 7894STDY15290420 to
# PD63118_RNA_G10. add plate number to supplier sample name.
pre_pcr_quants <-
  "data/resolveome/plate_layout/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout_dna_pre_pcr_quants.tsv" %>%
  readr::read_tsv() %>%
  tidyr::pivot_longer(cols = -`...1`, values_to = "well_quant") %>%
  dplyr::mutate(well = paste0(`...1`, name))
plate_layout <-
  "data/resolveome/plate_layout/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout_plate_layout.tsv" %>%
  readr::read_tsv() %>%
  tidyr::pivot_longer(cols = -`...1`, values_to = "well_content") %>%
  dplyr::mutate(well = paste0(`...1`, name))
pass_samples <-
  dplyr::full_join(pre_pcr_quants, plate_layout) %>%
  dplyr::filter(well_quant >= 2, well_content == "1 cell") %>%
  dplyr::mutate(supplier_sample_name = paste0("PD63118_RNA_", well)) %>%
  dplyr::pull(supplier_sample_name)
manifests[["7894stdy_manifest_26015_211124_RNA"]] <-
  manifests_og[["7894stdy_manifest_26015_211124_RNA"]] %>%
  dplyr::mutate(
    supplier_sample_name = ifelse(sanger_sample_id == "7894STDY15290420",
                                  "PD63118_RNA_G10", supplier_sample_name)) %>%
  dplyr::filter(supplier_sample_name %in% pass_samples) %>%
  dplyr::mutate(
    donor_id = "PD63118",
    plate = 3,
    seq_type = "rna",
    well = gsub(".*_", "", supplier_sample_name))

# run51611_lane7-8 (7761stdy_manifest_28305_201125_AL28)
manifests[["7761stdy_manifest_28305_201125_AL28"]] <-
  manifests_og[["7761stdy_manifest_28305_201125_AL28"]] %>%
  # remove controls
  dplyr::filter(grepl("^PD", supplier_sample_name),
                !grepl("10_nuclei", supplier_sample_name)) %>%
  # fix supplier sample names, adding plate number
  # (now using 384 kits - must add prefix to prevent plate number clash)
  dplyr::mutate(
    cell_id = supplier_sample_name,
    donor_id = gsub("_.*", "", supplier_sample_name),
    plate = 384.3,
    seq_type = "dnahyb") %>%
  dplyr::select(-well)

# combine standardised manifests, select necessary columns only
# from the manifests, we use the following metadata...
# - supplier_sample_name - separated into pdid, plate, seq_type and well
# - donor_id_required_for_ega - used as donor_id in samplesheet
manifest <-
  manifests %>%
  dplyr::bind_rows(.id = "manifest_file") %>%
  dplyr::select(manifest_file, sanger_sample_id,
                cell_id, donor_id, plate, seq_type, well)

# get clean cells only (and cells that have not yet been assessed)
# (plate 10 cells are being reassessed for doublet status with the dna)
# PD63121 and PD63126 cells have not yet been assessed
clean_cell_ids <-
  system("ls data/resolveome/manual_inspection/PD*.tsv", intern = TRUE) %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(!suspected_doublet | is.na(suspected_doublet),
                !chr_dropout | is.na(chr_dropout)) %>%
  dplyr::pull(cell_id)

# ss irods - combine manifest and sequencescape
# generate cell_id - plate{plate}_well{well}
# generate id      - plate{plate}_well{well}_{seq_type}_run{run_id}
ss_irods <-
  dplyr::inner_join(manifest, seqscape, by = "sanger_sample_id") %>%
  dplyr::mutate(
    lane_tmp = ifelse(lane == "lane1-8", "", lane),
    lane_n_tmp = ifelse(lane == "lane1-8", "", paste0("_", gsub("lane", "",
                                                                lane)))) %>%
  dplyr::transmute(
    cell_id = dplyr::case_when(!is.na(cell_id) ~ cell_id,
                               TRUE ~ paste0("plate", plate, "_well", well)),
    # merge the two duplicate RNA runs
    id = dplyr::case_when(
      run_id %in% c(49901, 50072) & seq_type == "rna" ~
        paste0(cell_id, "_", seq_type, "_merged"),
      TRUE ~ paste0(cell_id, "_", seq_type, "_run", run_id)),
    plate, well, seq_type, run_id, lane, plex_n = npg_aliquot_index,
    donor_id, manifest_file,
    bam = paste0("/seq/illumina/runs/", substr(run_id, 1, 2), "/", run_id, "/",
                 lane_tmp, "/plex", npg_aliquot_index, "/", run_id, lane_n_tmp,
                 "#", npg_aliquot_index, ".cram")) %>%
  # filter out low quality cells
  dplyr::filter(
    cell_id %in% clean_cell_ids |
    plate == 384.3)

# write ss irods (dna/dnahyb)
ss_irods %>%
  dplyr::filter(seq_type != "rna") %>%
  readr::write_csv(file.path(data_dir, "samplesheet_irods.csv"))

# write ss irods (rna)
ss_irods %>%
  dplyr::filter(seq_type == "rna") %>%
  readr::write_csv(file.path(data_dir, "samplesheet_irods_rna.csv"))

# ss local (with merged RNA BAMs collapsed)
ss_local <-
  ss_irods %>%
  dplyr::mutate(bam = paste0(data_dir, "/", donor_id, "/", id, "/bam/", id,
                             ".bam")) %>%
  dplyr::group_by(id, donor_id) %>%
  dplyr::summarise(dplyr::across(everything(), ~ paste(unique(.),
                                                       collapse = "|")))

# write ss local
ss_local %>%
  readr::write_csv(paste0(data_dir, "/samplesheet_local.csv"))
