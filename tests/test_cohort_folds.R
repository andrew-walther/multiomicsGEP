# ============================================================
# Script: test_cohort_folds.R
# Purpose: Verify cohort/event folds and their use by CV helpers.
# Author: Andrew Walther
# Created: 2026-11-11
# Dependencies: base R, survival, ebnm, tests/test_helpers.R
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
source("code/train_test_split.R")
source("code/select_K.R")
source("code/select_alpha_cv.R")
source("code/compute_cv_loglik.R")
# cv_survival_loglik checks these dependencies by name before using local mocks.
if (!exists("fit_ssbmf_yfb")) source("code/fit_ssbmf_yfb.R")
if (!exists("predict_ssbmf_yfb")) source("code/predict_ssbmf_yfb.R")

# Fixtures and fold behavior ----

.cohort_fold_data <- local({
  cohort <- rep(c("TCGA", "CPTAC"), c(39, 21))
  status <- c(rep(c(1L, 0L), length.out = 39),
              rep(c(0L, 1L), length.out = 21))
  Y <- cbind(id = seq_along(status), x = seq_along(status)^2,
             z = sin(seq_along(status)), w = cos(seq_along(status)))
  list(Y = Y, time = rev(seq_along(status)), status = status, cohort = cohort)
})

run_test("CF1: joint cells balanced, partition exact, all training folds have both cohorts", {
  d <- .cohort_fold_data
  for (seed in c(1L, 42L, 2026L)) {
    fd <- create_stratified_folds(d$status, 5L, seed, cohort_id = d$cohort)
    assert_equal(sort(unlist(fd$folds, use.names = FALSE)), seq_along(d$status))
    counts <- table(d$cohort, d$status, factor(fd$fold_id, levels = 1:5))
    assert_true(all(apply(counts, c(1, 2), function(x) max(x) - min(x)) <= 1))
    for (k in 1:5) {
      te <- fd$folds[[k]]
      tr <- setdiff(seq_along(d$status), te)
      assert_equal(sort(unique(d$cohort[tr])), sort(unique(d$cohort)))
      assert_true(all(table(factor(d$status[te], levels = 0:1)) > 0))
      assert_true(all(fd$fold_id[te] == k))
    }
  }
})

run_test("CF2: cohort-fold assignment is reproducible", {
  d <- .cohort_fold_data
  a <- create_stratified_folds(d$status, 5, 23, cohort_id = d$cohort)
  b <- create_stratified_folds(d$status, 5, 23, cohort_id = factor(d$cohort))
  assert_equal(a$fold_id, b$fold_id)
  c <- create_stratified_folds(d$status, 5, 24, cohort_id = d$cohort)
  assert_false(identical(a$fold_id, c$fold_id))
})

run_test("CF3: explicit NULL preserves the original event-only assignment", {
  status <- rep(c(1L, 0L), 15)
  set.seed(42)
  ev <- sample(which(status == 1), sum(status == 1))
  cn <- sample(which(status == 0), sum(status == 0))
  expected <- integer(length(status))
  expected[ev] <- rep(1:5, length.out = length(ev))
  expected[cn] <- rep(1:5, length.out = length(cn))
  assert_equal(create_stratified_folds(status, 5, 42)$fold_id, expected)
  assert_equal(create_stratified_folds(status, 5, 42, cohort_id = NULL)$fold_id, expected)
})

run_test("CF4: invalid cohort labels fail explicitly", {
  d <- .cohort_fold_data
  for (bad in list(d$cohort[-1], replace(d$cohort, 1, NA_character_))) {
    err <- tryCatch(create_stratified_folds(d$status, 5, cohort_id = bad), error = identity)
    assert_true(inherits(err, "error"))
    assert_true(grepl("cohort_id", conditionMessage(err)))
  }
})

run_test("CF5: a cohort confined to one held-out fold is rejected", {
  status <- rep(c(1L, 0L), 10)
  cohort <- c("rare", rep("common", 19))
  err <- tryCatch(create_stratified_folds(status, 5, cohort_id = cohort), error = identity)
  assert_true(inherits(err, "error"))
  assert_true(grepl("training fold", conditionMessage(err)))
})

run_test("CF6: sparse joint cells retain original patient indices", {
  status <- c(rep(1L, 12), rep(0L, 12))
  cohort <- rep("common", 24)
  cohort[c(12, 23, 24)] <- "small"
  fd <- create_stratified_folds(status, 3, 11, cohort_id = cohort)
  assert_equal(sort(unlist(fd$folds, use.names = FALSE)), 1:24)
  for (te in fd$folds) assert_equal(sort(unique(cohort[-te])), c("common", "small"))
})

# CV forwarding with inspected training rows ----

for (.cf_helper in c("select_K_cv_LB", "select_K_cv_YFB", "select_alpha_cv",
                     "cv_survival_loglik", "bicv_genomics_loglik")) {
  run_test(paste("CF integration:", .cf_helper), local({
    d <- .cohort_fold_data
    seen <- list()
    # Local dependency replacements inspect actual row IDs. No model fit is run.
    # Their signatures deliberately omit fold_cohort_id, so leaking it fails.
    fit_ssbmf_lb <- fit_ssbmf_yfb <- function(Y, time, status, K, verbose,
        alpha = 0.5, max_iter = 100, tol = 1e-5, prior_LF = "point_exponential",
        prior_beta = "normal", cohort_id = NULL, strata_id = NULL, beta_cohort_id = NULL) {
      seen[[length(seen) + 1L]] <<- as.integer(Y[, 1])
      list(EF = matrix(1, ncol(Y), K), EBeta = rep(1, K),
           EF_norms = rep(1, K), projection_center = NULL)
    }
    predict_ssbmf_lb <- predict_ssbmf_yfb <- function(Y_test, EF, EBeta,
        EF_norms = NULL, cohort_id_test = NULL, projection_center = NULL,
        projection_strata = NULL) list(risk_scores = Y_test[, 1])
    calc_cox_taylor_yf <- function(eta, time, status, strata = NULL)
      list(logPL = -sum(eta))
    .fit_genomics_only <- function(Y, K, max_iter, tol, prior_LF) {
      seen[[length(seen) + 1L]] <<- as.integer(Y[, 1])
      list(EL = matrix(1, nrow(Y), K), EF = matrix(1, ncol(Y), K))
    }
    if (startsWith(.cf_helper, "select_K_cv")) {
      fn <- select_K_cv; environment(fn) <- environment()
      fn(d$Y, d$time, d$status, K_grid = c(1L, 2L), n_folds = 3,
         seed = 42, model = sub("select_K_cv_", "", .cf_helper),
         fold_cohort_id = d$cohort)
    } else if (.cf_helper == "select_alpha_cv") {
      fn <- select_alpha_cv; environment(fn) <- environment()
      fn(d$Y, d$time, d$status, alpha_grid = c(0.2, 0.8), n_folds = 3,
         K_max = 1L, fold_cohort_id = d$cohort)
    } else if (.cf_helper == "cv_survival_loglik") {
      fn <- cv_survival_loglik; environment(fn) <- environment()
      fn(d$Y, d$time, d$status, K = 1L, n_folds = 3, fold_cohort_id = d$cohort)
    } else {
      fn <- bicv_genomics_loglik; environment(fn) <- environment()
      # Four columns carry independent patient IDs so either gene block reveals rows.
      x <- matrix(rep(seq_along(d$status), 4), ncol = 4)
      fn(x, d$status, K = 1L, n_row_folds = 3, n_col_folds = 2,
         fold_cohort_id = d$cohort)
    }
    expected <- create_stratified_folds(d$status, 3, 42, cohort_id = d$cohort)$folds
    assert_length(seen, if (.cf_helper == "cv_survival_loglik") 3L else 6L)
    for (ids in seen) {
      assert_equal(sort(unique(d$cohort[ids])), sort(unique(d$cohort)))
      assert_true(any(vapply(expected, function(te)
        identical(ids, setdiff(seq_along(d$status), te)), logical(1))))
    }
  }))
}

report_results("test_cohort_folds.R")
if (sys.nframe() == 0L && .test_env$failed > 0L) quit(status = 1L)
