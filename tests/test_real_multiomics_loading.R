# ============================================================
# Script: test_real_multiomics_loading.R
# Purpose: Real-data checks for the matched TCGA/ICGC expression +
#          methylation cohorts. Local only: skips if data/ is absent.
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: impute, yaml; code/load_multiomics_data.R
# Run: Rscript tests/test_real_multiomics_loading.R
# ============================================================

source("tests/test_helpers.R")
source("code/preprocess_multimodal_yfb.R")
source("code/load_multiomics_data.R")

root <- yaml::read_yaml("config/globals.yml")$multiomics_data$root
if (!dir.exists(file.path(root, "tcga")) || !dir.exists(file.path(root, "icgc"))) {
  cat("SKIP: matched multi-omics data not found under", root, "\n")
  quit(status = 0)
}

d <- suppressMessages(build_multiomics_cohorts())
cohorts <- d[c("training", "validation_primary", "validation_all")]

run_test("MOReal-T1: cohort sizes match the 10/1 inventory", {
  assert_equal(d$summary$n, c(144L, 50L, 67L))
  assert_equal(d$summary$events, c(75, 30, 40))
})

run_test("MOReal-T2: every cohort passes the multimodal YFB data contract", {
  spec <- preprocess_multimodal_yfb_training(d$training$Y, d$training$time,
                                             d$training$event)$training_spec
  for (v in cohorts[-1]) {
    aligned <- align_multimodal_yfb_prediction(v$Y, spec)
    assert_equal(length(aligned$ignored_features$methylation), 0L)
  }
})

run_test("MOReal-T3: subject IDs agree across expression, methylation and survival", {
  for (x in cohorts) {
    assert_equal(rownames(x$Y$expression), rownames(x$Y$methylation))
    assert_equal(rownames(x$Y$expression), names(x$time))
    assert_equal(names(x$time), names(x$event))
  }
})

run_test("MOReal-T4: expression is on the log2 scale and methylation on [0, 1]", {
  for (x in cohorts) {
    assert_true(max(x$Y$expression) < 25)
    assert_true(min(x$Y$methylation) >= 0 && max(x$Y$methylation) <= 1)
  }
})

run_test("MOReal-T5: validation cohorts use exactly the training features", {
  for (x in cohorts[-1]) {
    assert_equal(colnames(x$Y$expression), colnames(d$training$Y$expression))
    assert_equal(colnames(x$Y$methylation), colnames(d$training$Y$methylation))
  }
})

run_test("MOReal-T6: the primary-PDAC set is the PDAC-histology primary tumours", {
  cl <- d$icgc_clinical
  primary <- cl[cl$primary_pdac, ]
  assert_true(all(primary$HistoSubtype == "Pancreatic Ductal Adenocarcinoma"))
  assert_equal(sort(names(d$validation_primary$time)), sort(primary$icgc_donor_id))
})

run_test("MOReal-T7: ICGC survival agrees with PACA_AU_seq for every RNA sample", {
  pdac_root <- yaml::read_yaml("config/globals.yml")$multiomics_data$pdac_root
  info <- readRDS(file.path(root, "icgc", "icgc_matched_data.rds"))$info.expr
  paca <- readRDS(file.path(pdac_root, "original", "PACA_AU_seq.survival_data.rds"))
  m <- merge(info, paca, by.x = "icgc_sample_id", by.y = "sampID")
  assert_equal(nrow(m), nrow(info))
  assert_true(identical(m$survival_months, m$time))
  assert_true(identical(as.numeric(m$censored == "death"), m$event))
})

report_results("Real multi-omics data loading")
