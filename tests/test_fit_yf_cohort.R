# =============================================================================
# tests/test_fit_yf_cohort.R
# Integration tests for cohort_id extension in fit_cox_on_yf().
# Mirrors test_fit_modular_cohort.R but includes YFB-specific invariants:
#   - EF_norms length == K (cohort columns must not enter ZF computation)
#   - EBeta length == K
# =============================================================================

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------
local({
  set.seed(88)
  n <- 50L; p <- 25L; K <- 3L
  .f <<- list(
    n = n, p = p, K = K,
    Y      = matrix(rnorm(n * p), n, p),
    time   = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    cohort2 = rep(c("A", "B"), each = n / 2),
    cohort3 = rep(c("X", "Y", "Z"), length.out = n)
  )
})

# ---------------------------------------------------------------------------
# T1: Backward compatibility — cohort_id=NULL returns K-factor EF_norms
# ---------------------------------------------------------------------------
run_test("YFBCohort-T1: cohort_id=NULL gives valid output; EF_norms length K", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE, cohort_id = NULL)
  )
  assert_true(length(fit$EF_norms) == .f$K, msg = "EF_norms length == K")
  assert_true(length(fit$EBeta)    == .f$K, msg = "EBeta length == K")
  assert_true(is.null(fit$EF_cohort),        msg = "EF_cohort NULL for cohort_id=NULL")
  assert_true(is.null(fit$L_cohort),         msg = "L_cohort NULL for cohort_id=NULL")
})

# ---------------------------------------------------------------------------
# T2: Output structure with cohort — required cohort fields present
# ---------------------------------------------------------------------------
run_test("YFBCohort-T2: cohort_id='A'/'B' returns cohort output fields", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE, cohort_id = .f$cohort2)
  )
  assert_true(!is.null(fit$L_cohort),   msg = "L_cohort missing")
  assert_true(!is.null(fit$EF_cohort),  msg = "EF_cohort missing")
  assert_true(!is.null(fit$EF2_cohort), msg = "EF2_cohort missing")
  assert_true(!is.null(fit$cohort_id),  msg = "cohort_id missing")
})

# ---------------------------------------------------------------------------
# T3: K-factor matrices retain K columns; EF_cohort is p x (C-1)
# ---------------------------------------------------------------------------
run_test("YFBCohort-T3: EF has K cols; EF_cohort is p x (C-1)", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE, cohort_id = .f$cohort2)
  )
  assert_true(ncol(fit$EF) == .f$K,          msg = "EF ncol should be K")
  assert_true(ncol(fit$EL) == .f$K,          msg = "EL ncol should be K")
  assert_true(nrow(fit$EF_cohort) == .f$p,   msg = "EF_cohort nrow should be p")
  assert_true(ncol(fit$EF_cohort) == 1L,     msg = "EF_cohort ncol should be C-1=1")
})

# ---------------------------------------------------------------------------
# T4: EF_norms length == K — cohort columns must NOT contaminate ZF
# ---------------------------------------------------------------------------
run_test("YFBCohort-T4: EF_norms length == K (not K + C_cols)", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE, cohort_id = .f$cohort2)
  )
  assert_true(length(fit$EF_norms) == .f$K,
              msg = sprintf("EF_norms length should be K=%d, got %d",
                            .f$K, length(fit$EF_norms)))
})

# ---------------------------------------------------------------------------
# T5: Single-cohort input degrades gracefully — no cohort fields returned
# ---------------------------------------------------------------------------
run_test("YFBCohort-T5: nlevels(cohort_id)==1 returns no cohort fields", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE,
                  cohort_id = rep("A", .f$n))
  )
  assert_true(is.null(fit$EF_cohort), msg = "EF_cohort should be NULL for 1-level cohort")
  assert_true(is.null(fit$L_cohort),  msg = "L_cohort should be NULL for 1-level cohort")
  assert_true(length(fit$EF_norms) == .f$K, msg = "EF_norms length still == K")
})

# ---------------------------------------------------------------------------
# T6: ELBO finite and properly recorded
# ---------------------------------------------------------------------------
run_test("YFBCohort-T6: ELBO is finite for two-cohort fit", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 8L,
                  verbose = FALSE, cohort_id = .f$cohort2)
  )
  n_iter <- sum(fit$history$elbo_full != 0)
  assert_true(n_iter >= 1L, msg = "at least one ELBO entry")
  assert_finite(fit$history$elbo_full[seq_len(n_iter)],
                msg = "all recorded ELBO values finite")
})

# ---------------------------------------------------------------------------
# T7: Three-cohort fit (C=3) — EF_cohort ncol == 2; EF_norms still length K
# ---------------------------------------------------------------------------
run_test("YFBCohort-T7: three-cohort fit: EF_cohort 2 cols; EF_norms K", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE, cohort_id = .f$cohort3)
  )
  assert_true(ncol(fit$EF_cohort) == 2L,
              msg = "C=3 -> C-1=2 cohort columns")
  assert_true(length(fit$EF_norms) == .f$K,
              msg = "EF_norms length still == K with three cohorts")
})

# ---------------------------------------------------------------------------
# T8: EF_cohort second moments >= squared means
# ---------------------------------------------------------------------------
run_test("YFBCohort-T8: EF2_cohort >= EF_cohort^2 everywhere", {
  fit <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L,
                  verbose = FALSE, cohort_id = .f$cohort2)
  )
  assert_true(
    all(fit$EF2_cohort >= fit$EF_cohort^2 - 1e-12),
    msg = "EF2_cohort >= EF_cohort^2 (posterior variance >= 0)"
  )
})

# ---------------------------------------------------------------------------
# T9: Small sigma_F_cohort shrinks EF_cohort toward zero
# ---------------------------------------------------------------------------
run_test("YFBCohort-T9: small sigma_F_cohort shrinks EF_cohort toward zero", {
  fit_tight <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L, verbose = FALSE,
                  cohort_id = .f$cohort2, sigma_F_cohort = 1e-4)
  )
  fit_loose <- suppressMessages(
    fit_cox_on_yf(.f$Y, .f$time, .f$status,
                  K = .f$K, max_iter = 5L, verbose = FALSE,
                  cohort_id = .f$cohort2, sigma_F_cohort = 10.0)
  )
  assert_true(
    mean(abs(fit_tight$EF_cohort)) < mean(abs(fit_loose$EF_cohort)),
    msg = "tight prior -> smaller |EF_cohort| on average"
  )
})
