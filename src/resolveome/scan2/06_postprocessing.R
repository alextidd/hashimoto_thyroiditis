# runsub -M 64000 -R src/resolveome/scan2/06_postprocessing.R

# libraries
library(scan2)
library(magrittr)
library(dplyr)

# dirs
scan2_dir <- file.path(Sys.getenv("LUSTRE_125"), "projects/hashimoto_thyroiditis/out/resolveome/scan2/PD63118_develop_merged_normal_rescue/")
rescue_dir <- file.path(scan2_dir, "objects")
out_dir <- "out/resolveome/scan2/PD63118_develop_merged_normal_rescue/plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# load metadata
md <-
  readr::read_csv("out/resolveome/scan2/id_to_bam_id.csv") %>%
  dplyr::filter(grepl("^plate", id), grepl("_dna_", id))

# save raw to out_dir_i/raw and relaxed to out_dir_i/relaxed for each sample for objects and plots
purrr::pwalk(md, function(id, cell_id, bam_id) {

  # md %>% head(1) %>% as.vector() %>% list2env(envir = globalenv())

  cat(id, cell_id, bam_id, "\n")
  out_dir_i <- paste0(out_dir, "/", bam_id)
  dir.create(file.path(out_dir_i, "raw"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(out_dir_i, "relaxed"), showWarnings = FALSE, recursive = TRUE)
  # summaries <- list(); summaries$relaxed <- readRDS(file.path(out_dir_i, "relaxed/summary.rds")); summaries$raw <- readRDS(file.path(out_dir_i, "raw/summary.rds"))

  cat("loading scan2 object...\n")
  scan2_obj <- paste0(rescue_dir, "/", bam_id, "_scan2_object_rescue.rda")
  load(scan2_obj)
  summaries <- list()

  cat("making raw summary object...\n")
  summaries[["raw"]] <- make.summary.scan2(results)
  saveRDS(summaries[["raw"]], file.path(out_dir_i, "raw/summary.rds"))

  cat("writing raw vcf...\n")
  write.vcf(results, file.path(out_dir_i, "raw/variants.vcf"), overwrite = TRUE)

  cat("relaxing results...\n")
  # allow for reads in bulk so long as the probability of being a germline mutation is very low
  relaxed_params <- list(max.bulk.alt = 100, max.bulk.af = 1, max.bulk.binom.prob = 1e-5)
  relaxed_results <- update.static.filter.params(results, new.params = list(snv = relaxed_params, indel = relaxed_params))

  cat("making relaxed summary object...\n")
  summaries[["relaxed"]] <- make.summary.scan2(relaxed_results)
  saveRDS(summaries[["relaxed"]], file.path(out_dir_i, "relaxed/summary.rds"))

  cat("writing relaxed vcf...\n")
  write.vcf(relaxed_results, file.path(out_dir_i, "relaxed/variants.vcf"), overwrite = TRUE)

  purrr::walk2(summaries, names(summaries), function(summary, name) {
    cat("making plots for", name, "results...\n")
    
    pdf(file.path(out_dir_i, name, "abmodel.cov.pdf"))
    plot.abmodel.cov(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "gc.bias.pdf"))
    plot.gc.bias(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "ab.distn.pdf"))
    plot.ab.distn(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "binned.counts.pdf"))
    plot.binned.counts(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "depth.profile.pdf"))
    plot.depth.profile(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "dp.distn.pdf"))
    plot.dp.distn(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "target.fdr.effect.pdf"))
    plot.target.fdr.effect(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "mapd.pdf"))
    plot.mapd(summary)
    dev.off()
    pdf(file.path(out_dir_i, name, "mutburden.pdf"))
    plot.mutburden(summary)
    dev.off()

    if (name == "raw") {
      pdf(file.path(out_dir_i, name, "mutsig.rescue.pdf"), width = 15)
      plot.mutsig.rescue(summary)
      dev.off()    
      pdf(file.path(out_dir_i, name, "sensitivity.pdf"))
      plot.sensitivity(summary)
      dev.off()
      pdf(file.path(out_dir_i, name, "sensitivity.covs.pdf"), width = 20)
      plot.sensitivity.covs(summary)
      dev.off()
    }

  })
})