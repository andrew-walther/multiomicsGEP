# =============================================================================
# tests/test_fit_modular_cohort.R
# Integration tests for cohort_id extension in fit_supervised_mf_modular().
# Verifies backward compatibility, output structure, dimension invariants,
# ELBO behaviour, and shrinkage properties.
# =============================================================================

# ---------------------------------------------------------------------------
# Shared fixture: small deterministic dataset (n=40, p=30, K=3, C=2)
# ---------------------------------------------------------------------------
local({
  set.seed(99)
  n <- 40L; p <- 30L; K <- 3L
  .e <<- list(
    n = n, p = p, K = K,
    Y      = matrix(rnorm(n * p), n, p),
    time   = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    cohort2 = rep(c("A", "B"), each = n / 2),       # C=2 -> C-1=1
    cohort3 = rep(c("X", "Y", "Z"), length.out = n)  # C=3 -> C-1=2
  )
})

# ---------------------------------------------------------------------------
# T1: Backward compatibility — cohort_id=NULL matches frozen ELBO fixture
# ---------------------------------------------------------------------------
run_test("LBCohort-T1: cohort_id=NULL ELBO matches frozen baseline", {
  # Baseline regenerated 2026-07-12 (Phase 1a): elbo_full's assembly now boosts
  # the survival ELBO term by a factor of p under the default norm_convention=
  # "per_p" (see DECISIONS.md 2026-07-12). EL/EF/EBeta themselves are IDENTICAL
  # to the pre-Phase-1a values for this pure-noise fixture (verified directly --
  # all zero in both cases, correct EBNM shrinkage with no real signal); only
  # the reported elbo_full value legitimately changed.
  set.seed(99)
  n <- 40L; p <- 30L; K <- 3L
  Y      <- matrix(rnorm(n * p), n, p)
  time   <- rexp(n, 0.1)
  status <- rbinom(n, 1, 0.7)

  fit <- suppressMessages(
    fit_supervised_mf_modular(Y, time, status,
                              K = K, max_iter = 20L, alpha = 0.5,
                              verbose = FALSE, cohort_id = NULL)
  )
  baseline <- readRDS("tests/fixtures/lb_cohort_null_elbo_baseline.rds")
  # Compare only the non-zero ELBO entries (zeros pad unrun iterations)
  n_iter <- sum(fit$history$elbo_full != 0)
  assert_near(fit$history$elbo_full[seq_len(n_iter)],
              baseline[seq_len(n_iter)],
              tol = 1e-6,
              msg = "cohort_id=NULL ELBO must be bit-identical to baseline")
})

# ---------------------------------------------------------------------------
# T2: Output structure with cohort_id — required fields present
# ---------------------------------------------------------------------------
run_test("LBCohort-T2: cohort_id='A'/'B' returns cohort output fields", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L,
                              verbose = FALSE, cohort_id = .e$cohort2)
  )
  assert_true(!is.null(fit$L_cohort),   msg = "L_cohort missing from return")
  assert_true(!is.null(fit$EF_cohort),  msg = "EF_cohort missing from return")
  assert_true(!is.null(fit$EF2_cohort), msg = "EF2_cohort missing from return")
  assert_true(!is.null(fit$cohort_id),  msg = "cohort_id missing from return")
})

# ---------------------------------------------------------------------------
# T3: K-factor matrices retain K columns; cohort matrix is p x (C-1)
# ---------------------------------------------------------------------------
run_test("LBCohort-T3: EF has K cols; EF_cohort is p x (C-1)", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L,
                              verbose = FALSE, cohort_id = .e$cohort2)
  )
  # K-factor matrices unchanged in width
  assert_true(ncol(fit$EL) == .e$K, msg = "EL ncol should be K")
  assert_true(ncol(fit$EF) == .e$K, msg = "EF ncol should be K")
  # Cohort F row: p genes, C-1 = 1 cohort column
  assert_true(nrow(fit$EF_cohort) == .e$p, msg = "EF_cohort nrow should be p")
  assert_true(ncol(fit$EF_cohort) == 1L,   msg = "EF_cohort ncol should be C-1=1")
  # L_cohort: n rows, C-1 = 1 column (binary indicators)
  assert_true(nrow(fit$L_cohort) == .e$n, msg = "L_cohort nrow should be n")
  assert_true(ncol(fit$L_cohort) == 1L,   msg = "L_cohort ncol should be C-1=1")
})

# ---------------------------------------------------------------------------
# T4: Single-cohort input (nlevels == 1) degrades to no cohort fields
# ---------------------------------------------------------------------------
run_test("LBCohort-T4: nlevels(cohort_id)==1 returns no cohort fields", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L,
                              verbose = FALSE,
                              cohort_id = rep("A", .e$n))
  )
  assert_true(is.null(fit$EF_cohort),  msg = "EF_cohort should be NULL for 1-level cohort")
  assert_true(is.null(fit$L_cohort),   msg = "L_cohort should be NULL for 1-level cohort")
  assert_true(ncol(fit$EL) == .e$K,    msg = "EL ncol should still be K")
})

# ---------------------------------------------------------------------------
# T5: ELBO recorded — correct length and all finite
# ---------------------------------------------------------------------------
run_test("LBCohort-T5: ELBO is finite and has entries for each iteration", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 10L,
                              verbose = FALSE, cohort_id = .e$cohort2)
  )
  n_iter <- sum(fit$history$elbo_full != 0)
  assert_true(n_iter >= 1L, msg = "at least one ELBO entry recorded")
  assert_finite(fit$history$elbo_full[seq_len(n_iter)],
                msg = "all ELBO entries finite")
})

# ---------------------------------------------------------------------------
# T6: cohort KL contribution is <= 0 (compute_normal_kl sign check on real fit)
# ---------------------------------------------------------------------------
run_test("LBCohort-T6: cohort KL contribution to elbo_full is <= 0", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 8L,
                              verbose = FALSE, cohort_id = .e$cohort2)
  )
  # compute_normal_kl returns a negative scalar; verify on the actual fit output
  kl_val <- compute_normal_kl(fit$EF_cohort, fit$EF2_cohort,
                              sigma_F_cohort = 1.0)
  assert_true(kl_val <= 0,
              msg = "cohort KL contribution must be <= 0 (KL divergence >= 0)")
  assert_finite(kl_val, msg = "cohort KL contribution is finite")
})

# ---------------------------------------------------------------------------
# T7: EF_cohort second moments >= squared means (variance non-negative)
# ---------------------------------------------------------------------------
run_test("LBCohort-T7: EF2_cohort >= EF_cohort^2 everywhere", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L,
                              verbose = FALSE, cohort_id = .e$cohort2)
  )
  assert_true(
    all(fit$EF2_cohort >= fit$EF_cohort^2 - 1e-12),
    msg = "EF2_cohort >= EF_cohort^2 (posterior variance >= 0)"
  )
})

# ---------------------------------------------------------------------------
# T8: Three-cohort fit (C=3) — EF_cohort ncol == 2
# ---------------------------------------------------------------------------
run_test("LBCohort-T8: three-cohort fit returns EF_cohort with 2 columns", {
  fit <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L,
                              verbose = FALSE, cohort_id = .e$cohort3)
  )
  assert_true(ncol(fit$EF_cohort) == 2L,
              msg = "C=3 -> C-1=2 cohort columns")
  assert_true(nrow(fit$EF_cohort) == .e$p,
              msg = "EF_cohort nrow == p for three-cohort fit")
})

# ---------------------------------------------------------------------------
# T9: Small sigma_F_cohort shrinks EF_cohort toward zero
# ---------------------------------------------------------------------------
run_test("LBCohort-T9: small sigma_F_cohort shrinks EF_cohort toward zero", {
  fit_tight <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L, verbose = FALSE,
                              cohort_id = .e$cohort2, sigma_F_cohort = 1e-4)
  )
  fit_loose <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L, verbose = FALSE,
                              cohort_id = .e$cohort2, sigma_F_cohort = 10.0)
  )
  assert_true(
    mean(abs(fit_tight$EF_cohort)) < mean(abs(fit_loose$EF_cohort)),
    msg = "tight prior -> smaller |EF_cohort| on average"
  )
})

# ---------------------------------------------------------------------------
# T10: character vs factor cohort_id gives same EF_cohort
# ---------------------------------------------------------------------------
run_test("LBCohort-T10: character and factor cohort_id give same EF_cohort", {
  set.seed(7)
  fit_chr <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L, verbose = FALSE,
                              cohort_id = .e$cohort2)
  )
  set.seed(7)
  fit_fac <- suppressMessages(
    fit_supervised_mf_modular(.e$Y, .e$time, .e$status,
                              K = .e$K, max_iter = 5L, verbose = FALSE,
                              cohort_id = factor(.e$cohort2))
  )
  assert_near(fit_chr$EF_cohort, fit_fac$EF_cohort, tol = 1e-10,
              msg = "character and factor cohort_id give identical EF_cohort")
})

# ---------------------------------------------------------------------------
# T11: Cohort absorption — adding cohort_id does not crash; ELBO is finite
#       (Sanity: model runs end-to-end with strong platform offset in Y)
# ---------------------------------------------------------------------------
run_test("LBCohort-T11: strong platform offset absorbed without error", {
  set.seed(55)
  n <- 40L; p <- 20L; K <- 2L
  offset <- c(rep(0, n / 2), rep(5, n / 2))   # large cohort offset
  Y_off  <- matrix(rnorm(n * p), n, p) + offset
  time   <- rexp(n, 0.1)
  status <- rbinom(n, 1, 0.7)
  cid    <- rep(c("A", "B"), each = n / 2)

  fit <- suppressMessages(
    fit_supervised_mf_modular(Y_off, time, status,
                              K = K, max_iter = 10L,
                              verbose = FALSE, cohort_id = cid)
  )
  n_iter <- sum(fit$history$elbo_full != 0)
  assert_finite(fit$history$elbo_full[seq_len(n_iter)],
                msg = "ELBO finite under strong platform offset")
  assert_true(!is.null(fit$EF_cohort), msg = "EF_cohort present")
})
