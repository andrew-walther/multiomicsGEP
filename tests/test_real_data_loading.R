# =============================================================================
# tests/test_real_data_loading.R
# Test suite for real-data helpers in run_factor_modular_simulation.R:
#   filter_top_genes(), load_real_data(), pool_datasets()
# Plus a regression test verifying synthetic RMSE has not regressed.
#
# Usage (from project root):
#   Rscript tests/test_real_data_loading.R
#
# Data-dependent tests (T2, T3, T4) are auto-skipped when PDAC_DATA_ROOT
# is not set or does not exist — safe to run in CI.
# =============================================================================

source("tests/test_helpers.R")

# Source the helpers defined inside the runner script.
# fit_modular.R has a runner block that errors when DATA_MODE="real" and
# real_Y is NULL — we suppress it the same way the runner does.
suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
# Source the runner to get filter_top_genes / load_real_data / pool_datasets.
# DATA_MODE env var is not set, so execution falls through to the synthetic
# branch (which sets seed, runs quickly, writes nothing if table/figure dirs
# are not created — we only need the helper functions to be defined).
suppressMessages(tryCatch(
  source("results/legacy/modular_sim_factor/run_factor_modular_simulation.R"),
  error = function(e) invisible(NULL)
))

# Determine whether real-data tests can run
pdac_root <- Sys.getenv("PDAC_DATA_ROOT", unset = path.expand(
  "~/OneDrive - University of North Carolina at Chapel Hill/UNC Dissertation (Liu)/PDAC_data"
))
real_data_available <- dir.exists(pdac_root)

ALL_DATASETS_LOCAL <- c("TCGA_PAAD", "CPTAC", "Dijk", "Moffitt_GEO_array",
                         "PACA_AU_array", "PACA_AU_seq", "Puleo_array")

# =============================================================================
# T1: filter_top_genes() — pure unit tests, no data needed
# =============================================================================

cat("=== T1: filter_top_genes() ===\n")

run_test("T1.1: selects top-N most variable genes", {
  set.seed(1)
  Y <- matrix(c(rnorm(50, sd = 10),   # first 5 cols: high variance
                rnorm(50, sd = 0.01)), # last 5 cols: near-zero variance
              nrow = 10, ncol = 10)
  gnames <- paste0("gene", 1:10)
  out <- filter_top_genes(Y, gnames, top_n = 5)
  assert_equal(ncol(out$Y), 5L)
  # All returned genes should be from the high-variance block (cols 1-5)
  assert_true(all(out$gene_names %in% paste0("gene", 1:5)),
              "expected only high-variance genes to be kept")
})

run_test("T1.2: passthrough when top_n >= p", {
  set.seed(2)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10)
  gnames <- paste0("g", 1:10)
  out <- filter_top_genes(Y, gnames, top_n = 10)
  assert_equal(ncol(out$Y), 10L)
  assert_equal(out$gene_names, gnames)
})

run_test("T1.3: passthrough when top_n = NULL", {
  set.seed(3)
  Y <- matrix(rnorm(100), nrow = 10, ncol = 10)
  gnames <- paste0("g", 1:10)
  out <- filter_top_genes(Y, gnames, top_n = NULL)
  assert_equal(ncol(out$Y), 10L)
  assert_equal(out$gene_names, gnames)
})

run_test("T1.4: top_n = 1 returns single gene", {
  set.seed(4)
  Y <- matrix(rnorm(50), nrow = 10, ncol = 5)
  gnames <- paste0("g", 1:5)
  out <- filter_top_genes(Y, gnames, top_n = 1)
  assert_equal(ncol(out$Y), 1L)
  assert_length(out$gene_names, 1L)
})

run_test("T1.5: gene_names length matches Y columns after filtering", {
  set.seed(5)
  Y <- matrix(rnorm(200), nrow = 20, ncol = 10)
  gnames <- paste0("gene", 1:10)
  out <- filter_top_genes(Y, gnames, top_n = 6)
  assert_equal(length(out$gene_names), ncol(out$Y))
})

run_test("T1.6: output Y has correct nrow (unchanged)", {
  set.seed(6)
  Y <- matrix(rnorm(200), nrow = 20, ncol = 10)
  gnames <- letters[1:10]
  out <- filter_top_genes(Y, gnames, top_n = 4)
  assert_equal(nrow(out$Y), 20L)
})

# =============================================================================
# T2: load_real_data() — skipped if data unavailable
# =============================================================================

cat("\n=== T2: load_real_data() ===\n")

if (!real_data_available) {
  cat("  [SKIP] PDAC data root not found — skipping T2, T3, T4\n\n")
} else {

  for (ds in ALL_DATASETS_LOCAL) {

    run_test(sprintf("T2.1 [%s]: loads without error", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_true(!is.null(d))
    })

    run_test(sprintf("T2.2 [%s]: n > 0, p == top_n", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_true(d$n > 0, "n must be positive")
      assert_equal(d$p, 500L)
    })

    run_test(sprintf("T2.3 [%s]: no NAs in Y, time, status", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_false(anyNA(d$Y),      "Y contains NA")
      assert_false(anyNA(d$time),   "time contains NA")
      assert_false(anyNA(d$status), "status contains NA")
    })

    run_test(sprintf("T2.4 [%s]: time > 0", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_true(all(d$time > 0), "all times must be positive")
    })

    run_test(sprintf("T2.5 [%s]: status in {0, 1}", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_true(all(d$status %in% c(0L, 1L)), "status must be 0 or 1")
    })

    run_test(sprintf("T2.6 [%s]: Y is numeric matrix", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_true(is.numeric(d$Y),  "Y must be numeric")
      assert_true(is.matrix(d$Y),   "Y must be a matrix")
    })

    run_test(sprintf("T2.7 [%s]: dimensions consistent (nrow Y == length time == length status)", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_equal(nrow(d$Y), length(d$time))
      assert_equal(nrow(d$Y), length(d$status))
    })

    run_test(sprintf("T2.8 [%s]: Y is column-centred (colMeans ~ 0)", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      max_col_mean <- max(abs(colMeans(d$Y)))
      assert_true(max_col_mean < 1e-10,
                  sprintf("max |colMean| = %.2e; Y must be column-centred", max_col_mean))
    })

    run_test(sprintf("T2.9 [%s]: gene_names length matches p", ds), {
      d <- load_real_data(ds, pdac_root, top_n = 500)
      assert_equal(length(d$gene_names), d$p)
    })
  }

  # =============================================================================
  # T3: pool_datasets() — RNA-seq trio
  # =============================================================================

  cat("\n=== T3: pool_datasets() ===\n")

  rnaseq_names <- c("TCGA_PAAD", "Dijk", "PACA_AU_seq")

  run_test("T3.1: gene intersection is non-empty for RNA-seq trio", {
    ds_list <- lapply(rnaseq_names, load_real_data,
                      pdac_root = pdac_root, top_n = 500)
    names(ds_list) <- rnaseq_names
    common <- Reduce(intersect, lapply(ds_list, "[[", "gene_names"))
    assert_true(length(common) > 0, "no common genes across RNA-seq trio")
  })

  run_test("T3.2: pooled nrow == sum of individual n", {
    ds_list <- lapply(rnaseq_names, load_real_data,
                      pdac_root = pdac_root, top_n = 500)
    names(ds_list) <- rnaseq_names
    pooled   <- pool_datasets(ds_list)
    expected <- sum(sapply(ds_list, "[[", "n"))
    assert_equal(pooled$n, expected)
    assert_equal(nrow(pooled$Y), expected)
  })

  run_test("T3.3: pooled ncol == length of gene intersection", {
    ds_list <- lapply(rnaseq_names, load_real_data,
                      pdac_root = pdac_root, top_n = 500)
    names(ds_list) <- rnaseq_names
    common <- Reduce(intersect, lapply(ds_list, "[[", "gene_names"))
    pooled <- pool_datasets(ds_list)
    assert_equal(pooled$p, length(common))
    assert_equal(ncol(pooled$Y), length(common))
  })

  run_test("T3.4: pooled Y has no NAs", {
    ds_list <- lapply(rnaseq_names, load_real_data,
                      pdac_root = pdac_root, top_n = 500)
    names(ds_list) <- rnaseq_names
    pooled <- pool_datasets(ds_list)
    assert_false(anyNA(pooled$Y), "pooled Y contains NA")
  })

  run_test("T3.5: dataset_labels length == pooled n", {
    ds_list <- lapply(rnaseq_names, load_real_data,
                      pdac_root = pdac_root, top_n = 500)
    names(ds_list) <- rnaseq_names
    pooled <- pool_datasets(ds_list)
    assert_equal(length(pooled$dataset_labels), pooled$n)
  })

  run_test("T3.6: pool_datasets() errors on empty gene intersection", {
    # Build two datasets with disjoint gene names
    d1 <- list(Y = matrix(1:10, 2, 5), gene_names = letters[1:5], n = 2L)
    d2 <- list(Y = matrix(1:10, 2, 5), gene_names = letters[6:10], n = 2L)
    err <- tryCatch(pool_datasets(list(a = d1, b = d2)), error = function(e) e)
    assert_true(inherits(err, "error"),
                "expected an error when gene intersection is empty")
  })

  # =============================================================================
  # T4: Synthetic regression — RMSE and convergence unchanged
  # =============================================================================

  cat("\n=== T4: Synthetic Regression ===\n")

  run_test("T4.1: synthetic RMSE in [0.95, 1.05]", {
    # Re-run the synthetic DGP with the canonical parameters (n=250, p=1000,
    # K=5, seed=42) and verify RMSE is near 1.0 (true noise SD).
    # This guards against accidental changes to fit_modular.R or update_*.R
    # that would break the benchmark results.
    set.seed(42)
    n_s <- 250; p_s <- 1000; K_s <- 5
    L_s <- matrix(rnorm(n_s * K_s), n_s, K_s)
    F_s <- matrix(0, p_s, K_s)
    for (k in 1:K_s) {
      active <- sample(1:p_s, round(p_s * 0.05))
      F_s[active, k] <- rnorm(length(active), 0, 5)
    }
    Y_s   <- L_s %*% t(F_s) + matrix(rnorm(n_s * p_s), n_s, p_s)
    B_s   <- c(1.5, -1.2, 0.8, -0.5, 0.0)
    eta_s <- as.vector(L_s %*% B_s)
    raw_t <- (-log(runif(n_s)) / (0.01 * exp(eta_s)))^(1 / 1.5)
    cen_t <- rexp(n_s, rate = 1 / 50)
    t_s   <- pmin(raw_t, cen_t)
    ev_s  <- as.integer(raw_t <= cen_t)

    # prior_LF = "point_normal" matches the DGP (L ~ N(0,1), both signs).
    # The default "point_exponential" constrains loadings to be non-negative,
    # which prevents recovery of this DGP and inflates RMSE to ~2.4.
    res_s <- fit_supervised_mf_modular(Y_s, t_s, ev_s, K = K_s,
                                        max_iter = 300, tol = 1e-3,
                                        prior_LF = "point_normal",
                                        verbose = FALSE)
    final_rmse <- tail(res_s$history$rmse, 1)
    assert_true(final_rmse >= 0.95 && final_rmse <= 1.05,
                sprintf("RMSE = %.4f; expected in [0.95, 1.05]", final_rmse))
  })

  run_test("T4.2: synthetic run converges", {
    set.seed(42)
    n_s <- 250; p_s <- 1000; K_s <- 5
    L_s <- matrix(rnorm(n_s * K_s), n_s, K_s)
    F_s <- matrix(0, p_s, K_s)
    for (k in 1:K_s) {
      active <- sample(1:p_s, round(p_s * 0.05))
      F_s[active, k] <- rnorm(length(active), 0, 5)
    }
    Y_s   <- L_s %*% t(F_s) + matrix(rnorm(n_s * p_s), n_s, p_s)
    B_s   <- c(1.5, -1.2, 0.8, -0.5, 0.0)
    eta_s <- as.vector(L_s %*% B_s)
    raw_t <- (-log(runif(n_s)) / (0.01 * exp(eta_s)))^(1 / 1.5)
    cen_t <- rexp(n_s, rate = 1 / 50)
    t_s   <- pmin(raw_t, cen_t)
    ev_s  <- as.integer(raw_t <= cen_t)

    res_s <- fit_supervised_mf_modular(Y_s, t_s, ev_s, K = K_s,
                                        max_iter = 300, tol = 1e-3,
                                        prior_LF = "point_normal",
                                        verbose = FALSE)
    assert_true(isTRUE(res_s$history$converged),
                sprintf("expected convergence; got converged=%s after %d iters",
                        res_s$history$converged, res_s$history$n_iter))
  })

}  # end if (real_data_available)

# =============================================================================
# T5: per_platform_standardize_cohorts() — pure unit tests, no data needed
# =============================================================================

cat("\n=== T5: per_platform_standardize_cohorts() ===\n")

suppressMessages(tryCatch(
  source("code/preprocess_desurv.R"),
  error = function(e) invisible(NULL)
))

run_test("T5.1: per-platform colMeans ≈ 0 after standardization", {
  set.seed(10)
  # Two cohorts with very different location scales to verify each is centred
  Y1 <- matrix(rnorm(50 * 20, mean = 100, sd = 10), nrow = 50, ncol = 20)
  Y2 <- matrix(rnorm(30 * 20, mean = -5,  sd = 0.5), nrow = 30, ncol = 20)
  out <- per_platform_standardize_cohorts(list(cohort1 = Y1, cohort2 = Y2))
  max_mean1 <- max(abs(colMeans(out$cohort1)))
  max_mean2 <- max(abs(colMeans(out$cohort2)))
  assert_true(max_mean1 < 1e-10,
              sprintf("cohort1 max |colMean| = %.2e; expected < 1e-10", max_mean1))
  assert_true(max_mean2 < 1e-10,
              sprintf("cohort2 max |colMean| = %.2e; expected < 1e-10", max_mean2))
})

run_test("T5.2: per-platform colSDs ≈ 1 after standardization", {
  set.seed(11)
  Y1 <- matrix(rnorm(50 * 20, mean = 100, sd = 10), nrow = 50, ncol = 20)
  Y2 <- matrix(rnorm(30 * 20, mean = -5,  sd = 0.5), nrow = 30, ncol = 20)
  out <- per_platform_standardize_cohorts(list(cohort1 = Y1, cohort2 = Y2))
  sds1 <- apply(out$cohort1, 2, sd)
  sds2 <- apply(out$cohort2, 2, sd)
  assert_true(max(abs(sds1 - 1)) < 1e-10,
              sprintf("cohort1 max |colSD - 1| = %.2e", max(abs(sds1 - 1))))
  assert_true(max(abs(sds2 - 1)) < 1e-10,
              sprintf("cohort2 max |colSD - 1| = %.2e", max(abs(sds2 - 1))))
})

run_test("T5.3: constant-gene columns handled without error (SD floor)", {
  # A cohort where some columns are constant — SD = 0 before floor
  Y_const <- matrix(0, nrow = 10, ncol = 5)          # all zeros: SD = 0
  Y_norm  <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)
  out <- tryCatch(
    per_platform_standardize_cohorts(list(const = Y_const, norm = Y_norm)),
    error = function(e) NULL
  )
  assert_false(is.null(out), "function must not error on constant-column cohort")
  assert_false(anyNA(out$const), "output for constant cohort must not contain NA")
})

run_test("T5.4: output dimensions are unchanged", {
  set.seed(12)
  Y1 <- matrix(rnorm(50 * 20), nrow = 50, ncol = 20)
  Y2 <- matrix(rnorm(30 * 20), nrow = 30, ncol = 20)
  out <- per_platform_standardize_cohorts(list(a = Y1, b = Y2))
  assert_equal(nrow(out$a), 50L)
  assert_equal(ncol(out$a), 20L)
  assert_equal(nrow(out$b), 30L)
  assert_equal(ncol(out$b), 20L)
})

# =============================================================================
# T6: build_pdac_genesets() and its components — real-data (PDAC_DATA_ROOT +
# DeSurv SI appendix PDF) integration tests. Auto-skipped when either is
# unavailable, same convention as T2-T4 above.
# =============================================================================

cat("\n=== T6: build_pdac_genesets() (Moffitt/Bailey/DeSurv custom gene sets) ===\n")

source("code/pathway_enrichment.R")

.desurv_si_pdf <- file.path(dirname(pdac_root), "papers", "DeSurv", "si_appendix.pdf")
desurv_si_available <- real_data_available && file.exists(.desurv_si_pdf) &&
  nchar(Sys.which("pdftotext")) > 0

if (!real_data_available) {
  cat("  [SKIPPED] T6.1-T6.2 (load_moffitt_bailey_genesets): PDAC_DATA_ROOT not set/found\n")
} else {
  run_test("T6.1: load_moffitt_bailey_genesets returns 6 non-empty gene sets", {
    res <- load_moffitt_bailey_genesets(pdac_root)
    expected_names <- c("Moffitt_BasalLike", "Moffitt_Classical", "Bailey_Squamous",
                         "Bailey_Immunogenic", "Bailey_PancreaticProgenitor", "Bailey_ADEX")
    assert_true(all(expected_names %in% names(res)), msg = "missing expected gene sets")
    assert_true(all(vapply(res, length, integer(1)) > 0), msg = "found an empty gene set")
  })

  run_test("T6.2: Moffitt basal/classical are each exactly 25 genes (published classifier size)", {
    res <- load_moffitt_bailey_genesets(pdac_root)
    assert_length(res$Moffitt_BasalLike, 25)
    assert_length(res$Moffitt_Classical, 25)
  })
}

if (!desurv_si_available) {
  cat("  [SKIPPED] T6.3 (extract_desurv_genesets): DeSurv si_appendix.pdf or pdftotext not available\n")
} else {
  run_test("T6.3: extract_desurv_genesets returns 270 genes for each of D1/D2/D3, no cross-overlap", {
    res <- extract_desurv_genesets(.desurv_si_pdf)
    assert_length(res$D1, 270)
    assert_length(res$D2, 270)
    assert_length(res$D3, 270)
    assert_true(length(intersect(res$D1, res$D2)) == 0, msg = "D1/D2 should be disjoint top-N lists")
    assert_true(length(intersect(res$D1, res$D3)) == 0, msg = "D1/D3 should be disjoint top-N lists")
    assert_true(length(intersect(res$D2, res$D3)) == 0, msg = "D2/D3 should be disjoint top-N lists")
  })
}

if (!real_data_available || !desurv_si_available) {
  cat("  [SKIPPED] T6.4-T6.5 (build_pdac_genesets): requires both PDAC_DATA_ROOT and the DeSurv SI PDF\n")
} else {
  .pdac_genesets_test_dir <- file.path(tempdir(), "pathway_enrichment_test")

  run_test("T6.4: build_pdac_genesets assembles all 9 sets and writes rds + manifest", {
    res <- build_pdac_genesets(pdac_root, .desurv_si_pdf, .pdac_genesets_test_dir)
    assert_true(length(res) == 9, msg = sprintf("expected 9 gene sets, got %d", length(res)))
    assert_true(file.exists(file.path(.pdac_genesets_test_dir, "pdac_genesets.rds")))
    assert_true(file.exists(file.path(.pdac_genesets_test_dir, "genesets_manifest.txt")))
  })

  run_test("T6.5: manifest records a source citation line for every gene set", {
    manifest <- readLines(file.path(.pdac_genesets_test_dir, "genesets_manifest.txt"))
    res <- readRDS(file.path(.pdac_genesets_test_dir, "pdac_genesets.rds"))
    for (nm in names(res)) {
      assert_true(any(grepl(nm, manifest, fixed = TRUE)),
                  msg = paste("manifest missing citation line for", nm))
    }
  })

  unlink(.pdac_genesets_test_dir, recursive = TRUE)
}

# =============================================================================
# Summary
# =============================================================================

report_results("tests/test_real_data_loading.R")
