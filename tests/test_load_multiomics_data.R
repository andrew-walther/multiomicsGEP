# ============================================================
# Script: test_load_multiomics_data.R
# Purpose: Tests for the matched TCGA/ICGC multi-omics loader helpers
#          (synthetic inputs only; real-data checks are in
#          tests/test_real_multiomics_loading.R).
# Author: Andrew Walther
# Created: 2026-10-01
# Dependencies: impute; code/load_multiomics_data.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("build_multiomics_cohorts")) source("code/load_multiomics_data.R")

run_test("MOData-T1: unusable survival is dropped and reported by reason", {
  out <- suppressMessages(filter_usable_survival(
    c("a", "b", "c", "d", "e"), c(5, NA, 0, 3, 7), c(1, 0, 1, NA, 0), "toy"))
  assert_equal(names(out$time), c("a", "e"))
  assert_equal(unname(out$event), c(1, 0))
  assert_equal(out$dropped$reason, c("missing time", "time <= 0", "missing or invalid event"))
})

run_test("MOData-T2: ICGC clinical collapses duplicate RNA samples to one donor row", {
  info <- data.frame(
    icgc_donor_id = c("D1", "D1", "D2", "D3"),
    submitted_donor_id = c("S1", "S1", "S2", "S3"),
    survival_months = c(10, 10, 4, NA),
    censored = c("death", "death", "censor", NA),
    HistoSubtype = c("PDAC", "PDAC", "PDAC", "IPMN"),
    Sample.type = c("Primary tumour", "Metastatic tumour", "Primary tumour", "Primary tumour"),
    stringsAsFactors = FALSE)
  cl <- collapse_icgc_clinical(info, "PDAC", "Primary tumour")
  assert_equal(cl$icgc_donor_id, c("D1", "D2", "D3"))
  assert_equal(cl$event, c(1, 0, NA))
  # D1 has one metastatic sample, so it is not primary PDAC; D3 is IPMN
  assert_equal(cl$primary_pdac, c(FALSE, TRUE, FALSE))
})

run_test("MOData-T3: conflicting survival across a donor's samples stops", {
  info <- data.frame(icgc_donor_id = c("D1", "D1"), submitted_donor_id = "S1",
                     survival_months = c(10, 12), censored = "death",
                     HistoSubtype = "PDAC", Sample.type = "Primary tumour")
  err <- tryCatch(collapse_icgc_clinical(info, "PDAC", "Primary tumour"),
                  error = function(e) conditionMessage(e))
  assert_true(grepl("conflicting survival", err))
})

run_test("MOData-T4: KNN imputation fills only missing entries", {
  set.seed(1)
  base <- runif(10)
  # 60 CpGs that are near-copies of one profile, so a KNN estimate is accurate.
  # As in Yusha's code, only CpGs with a missing value go to impute.knn, so
  # neighbours come from those CpGs. With 30 of them, k = 10 neighbours is
  # well below the pool size, as it is on the real data (~8,000 CpGs).
  m <- t(replicate(60, base + rnorm(10, sd = 0.01)))
  na_idx <- cbind(1:30, rep(1:10, 3))
  truth <- m[na_idx]
  m[na_idx] <- NA
  filled <- suppressMessages(impute_methylation_knn(m, "toy"))
  assert_false(anyNA(filled))
  assert_equal(filled[31:60, ], m[31:60, ])
  assert_near(filled[na_idx], truth, tol = 0.05)
})

run_test("MOData-T5: screening uses training data and returns the requested sizes", {
  set.seed(2)
  expr <- matrix(rnorm(40 * 12, 5), 40, 12, dimnames = list(paste0("g", 1:40), paste0("s", 1:12)))
  expr["g1", ] <- 5                       # zero variance: must be excluded
  meth <- matrix(runif(60 * 12, 0.4, 0.6), 60, 12,
                 dimnames = list(paste0("cg", 1:60), paste0("s", 1:12)))
  meth["cg7", ] <- seq(0, 1, length.out = 12)   # largest variance
  f <- screen_multiomics_features(list(expression = expr, methylation = meth), 10, 5)
  assert_length(f$genes, 10)
  assert_length(f$cpgs, 5)
  assert_false("g1" %in% f$genes)
  assert_equal(f$cpgs[1], "cg7")
})

run_test("MOData-T6: cohorts are returned subjects x features on the screened features", {
  cohort <- list(
    expression = matrix(1:6, 3, 2, dimnames = list(c("g1", "g2", "g3"), c("s1", "s2"))),
    methylation = matrix(0.5, 2, 2, dimnames = list(c("cg1", "cg2"), c("s1", "s2"))),
    time = c(s1 = 3, s2 = 4), event = c(s1 = 1, s2 = 0))
  out <- as_multimodal_yfb_cohort(cohort, list(genes = c("g3", "g1"), cpgs = "cg2"))
  assert_equal(dimnames(out$Y$expression), list(c("s1", "s2"), c("g3", "g1")))
  assert_equal(out$Y$expression["s2", "g3"], 6L)
  assert_equal(dim(out$Y$methylation), c(2L, 1L))
  # passes the multimodal YFB data contract
  ok <- preprocess_multimodal_yfb_training(out$Y, out$time, out$event)
  assert_equal(ok$training_spec$subject_ids, c("s1", "s2"))
})
