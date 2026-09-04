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

# =============================================================================
# Integration tests for beta_cohort_id (Stage 2, DECISIONS.md 2026-09-04) --
# COHORT-SPECIFIC SURVIVAL COEFFICIENTS, a distinct extension from cohort_id
# above (genomics offsets in L/F). Separate fixture, separate section: these
# tests exercise fit_cox_on_yf()'s beta_cohort_id argument end-to-end,
# including predict_cox_on_yf()'s cohort_id_test path. Unit tests for the
# underlying arithmetic (update_beta_cohort_k/_all, compute_pooled_beta) live
# in tests/test_update_beta_cohort.R; these are the integration layer.
# =============================================================================

local({
  set.seed(4242)
  n <- 80L; p <- 40L; K <- 3L
  .fb <<- list(
    n = n, p = p, K = K,
    Y      = matrix(rnorm(n * p), n, p),
    time   = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7),
    cohort2 = rep(c("A", "B"), each = n / 2)
  )
})

# ---------------------------------------------------------------------------
# BC-INT-T1: beta_cohort_id = NULL is the regression gate -- must reproduce
# the current (non-cohort) fit bit-for-bit. Non-negotiable per the 9/4 plan.
# ---------------------------------------------------------------------------
run_test("YFBBetaCohort-T1: beta_cohort_id=NULL reproduces the fit BIT-FOR-BIT", {
  fit_ref <- suppressMessages(
    fit_cox_on_yf(.fb$Y, .fb$time, .fb$status, K = .fb$K, max_iter = 15L,
                  prior_beta = "normal", alpha = 0.5, sign_correction = TRUE,
                  verbose = FALSE)
  )
  fit_null <- suppressMessages(
    fit_cox_on_yf(.fb$Y, .fb$time, .fb$status, K = .fb$K, max_iter = 15L,
                  prior_beta = "normal", alpha = 0.5, sign_correction = TRUE,
                  verbose = FALSE, beta_cohort_id = NULL)
  )
  assert_equal(fit_null$EBeta, fit_ref$EBeta, msg = "EBeta must be bit-for-bit identical")
  assert_equal(fit_null$EF, fit_ref$EF, msg = "EF must be bit-for-bit identical")
  assert_equal(fit_null$EL, fit_ref$EL, msg = "EL must be bit-for-bit identical")
  assert_equal(fit_null$history$elbo_full, fit_ref$history$elbo_full,
               msg = "elbo_full trajectory must be bit-for-bit identical")
})

# ---------------------------------------------------------------------------
# BC-INT-T2: shape and structure of a cohort-beta fit
# ---------------------------------------------------------------------------
run_test("YFBBetaCohort-T2: beta_cohort_id gives a K x C EBeta and the expected extra fields", {
  fit <- suppressMessages(
    fit_cox_on_yf(.fb$Y, .fb$time, .fb$status, K = .fb$K, max_iter = 10L,
                  verbose = FALSE, beta_cohort_id = .fb$cohort2)
  )
  assert_true(is.matrix(fit$EBeta), msg = "EBeta should be a matrix under beta_cohort_id")
  assert_true(all(dim(fit$EBeta) == c(.fb$K, 2L)), msg = "EBeta should be K x C")
  assert_equal(colnames(fit$EBeta), c("A", "B"), msg = "EBeta columns should be named by cohort level")
  assert_length(fit$EBeta_pooled, .fb$K, msg = "EBeta_pooled should be a K-vector")
  assert_finite(fit$EBeta_pooled, msg = "EBeta_pooled must be finite")
  assert_true(all(fit$EBeta2 >= fit$EBeta^2 - 1e-8), msg = "EBeta2 >= EBeta^2 elementwise")
})

# ---------------------------------------------------------------------------
# BC-INT-T3: end-to-end recovery of genuinely cohort-specific survival signal
# ---------------------------------------------------------------------------
run_test("YFBBetaCohort-T3: recovers opposite-signed survival signal across two cohorts", {
  set.seed(77)
  n2 <- 300L; p2 <- 60L; K2 <- 2L
  # Strong single latent factor, opposite prognostic direction per cohort.
  L_true <- matrix(abs(rnorm(n2)), n2, 1)
  F_true <- matrix(abs(rnorm(p2)), p2, 1)
  Y2 <- L_true %*% t(F_true) + matrix(rnorm(n2 * p2, sd = 0.2), n2, p2)
  Y2 <- cbind(Y2, matrix(rnorm(n2 * (p2)), n2, p2))  # decoy noise genes
  cohort2 <- rep(c("A", "B"), each = n2 / 2)
  beta_true <- ifelse(cohort2 == "A", 2.5, -2.5)
  risk <- as.vector(L_true) * beta_true
  time2   <- rexp(n2, rate = exp(risk - mean(risk)))
  status2 <- rbinom(n2, 1, 0.85)

  fit <- suppressMessages(
    fit_cox_on_yf(Y2, time2, status2, K = K2, max_iter = 60L, tol = 1e-4,
                  prior_beta = "normal", alpha = 1.0, sign_correction = TRUE,
                  verbose = FALSE, beta_cohort_id = cohort2)
  )
  # The two cohort columns of the survival-carrying factor should have
  # opposite signs (which factor carries the signal depends on CAVI's own
  # ordering, so check across all K factors for the largest |beta| spread).
  beta_range_per_factor <- apply(fit$EBeta, 1, function(r) diff(range(r)))
  best_k <- which.max(beta_range_per_factor)
  assert_true(sign(fit$EBeta[best_k, 1]) != sign(fit$EBeta[best_k, 2]),
              "the factor with the largest cross-cohort beta spread should have opposite-signed cohort coefficients")
})

# ---------------------------------------------------------------------------
# BC-INT-T4: predict_cox_on_yf() -- unseen cohort errors, NULL falls back to pooled
# ---------------------------------------------------------------------------
run_test("YFBBetaCohort-T4: predict_cox_on_yf errors loudly on an unseen cohort level", {
  fit <- suppressMessages(
    fit_cox_on_yf(.fb$Y, .fb$time, .fb$status, K = .fb$K, max_iter = 10L,
                  verbose = FALSE, beta_cohort_id = .fb$cohort2)
  )
  err <- tryCatch({
    predict_cox_on_yf(.fb$Y[1:5, , drop = FALSE], fit$EF, fit$EBeta, EF_norms = fit$EF_norms,
                       cohort_id_test = rep("UNSEEN_COHORT", 5))
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error on an unseen cohort_id_test level")
})

run_test("YFBBetaCohort-T5: predict_cox_on_yf falls back to EBeta_pooled with cohort_id_test=NULL", {
  fit <- suppressMessages(
    fit_cox_on_yf(.fb$Y, .fb$time, .fb$status, K = .fb$K, max_iter = 10L,
                  verbose = FALSE, beta_cohort_id = .fb$cohort2)
  )
  pred <- predict_cox_on_yf(.fb$Y[1:10, , drop = FALSE], fit$EF, fit$EBeta_pooled,
                             EF_norms = fit$EF_norms)
  assert_length(pred$risk_scores, 10L, "risk_scores should have length n_test")
  assert_finite(pred$risk_scores, "risk_scores must be finite")

  # Passing the raw K x C matrix without cohort_id_test should fail loudly
  # rather than silently misuse it.
  err <- tryCatch({
    predict_cox_on_yf(.fb$Y[1:10, , drop = FALSE], fit$EF, fit$EBeta, EF_norms = fit$EF_norms)
    NULL
  }, error = function(e) e)
  assert_true(!is.null(err), "expected an error passing a K x C EBeta with cohort_id_test=NULL")
})

run_test("YFBBetaCohort-T6: predict_cox_on_yf scores matched cohort levels correctly by name", {
  fit <- suppressMessages(
    fit_cox_on_yf(.fb$Y, .fb$time, .fb$status, K = .fb$K, max_iter = 10L,
                  verbose = FALSE, beta_cohort_id = .fb$cohort2)
  )
  pred <- predict_cox_on_yf(.fb$Y, fit$EF, fit$EBeta, EF_norms = fit$EF_norms,
                             cohort_id_test = .fb$cohort2)
  # Manually recompute risk scores by looking up each patient's own column.
  ZF_manual <- .fb$Y %*% sweep(fit$EF, 2, fit$EF_norms, "/")
  cohort_idx <- match(.fb$cohort2, colnames(fit$EBeta))
  expected <- rowSums(ZF_manual * t(fit$EBeta)[cohort_idx, , drop = FALSE])
  assert_near(pred$risk_scores, expected, tol = 1e-8,
              "predict_cox_on_yf's cohort-matched risk scores must match a manual per-patient lookup")
})

# =============================================================================
# Section: EBeta_pooled pre/post-Phase-C coherence, and EBeta2 sign-flip
# invariance (review finding, Step 3, fixed 2026-09-04 -- see DECISIONS.md).
#
# Fixture below (seed=18, with a short sample()-selected data-generating
# regime) was found by search specifically because it makes fit_cox_on_yf()'s
# Phase C sign-check trigger a flip (a genuine, if rare, event -- most random
# data converges beta to the correct sign on its own since beta itself is
# unconstrained). A flip is REQUIRED to exercise this fix: without one,
# EBeta_pooled would trivially match between sign_correction=TRUE/FALSE under
# both the old buggy code and the fix, since neither path would ever touch a
# post-flip EBeta.
# =============================================================================

local({
  n <- 80L; p <- 20L
  cohort2 <- rep(c("A", "B"), each = n / 2)
  set.seed(18)
  K <- sample(3:8, 1)
  L_true <- matrix(rexp(n * K, 1), n, K)
  F_true <- matrix(rexp(p * K, 1), p, K)
  beta_true <- rnorm(K, 0, sample(c(1, 2, 3), 1))
  Y <- L_true %*% t(F_true) + matrix(rnorm(n * p, sd = sample(c(0.3, 0.5, 0.8), 1)), n, p)
  eta_true <- L_true %*% beta_true
  time   <- rexp(n, rate = exp(sample(c(0.3, 0.5, 0.7, 1), 1) * scale(eta_true)[, 1]))
  status <- rbinom(n, 1L, 0.85)
  mi     <- sample(1:3, 1)
  .pc <<- list(n = n, p = p, K = K, Y = Y, time = time, status = status,
               cohort2 = cohort2, max_iter = mi)
})

run_test("YFBBetaCohort-T7: this fixture's Phase C sign check actually flips (fixture validity check)", {
  fit_true  <- suppressMessages(
    fit_cox_on_yf(.pc$Y, .pc$time, .pc$status, K = .pc$K, max_iter = .pc$max_iter,
                  verbose = FALSE, sign_correction = TRUE, beta_cohort_id = .pc$cohort2)
  )
  fit_false <- suppressMessages(
    fit_cox_on_yf(.pc$Y, .pc$time, .pc$status, K = .pc$K, max_iter = .pc$max_iter,
                  verbose = FALSE, sign_correction = FALSE, beta_cohort_id = .pc$cohort2)
  )
  assert_near(fit_true$EBeta, -fit_false$EBeta, tol = 1e-9,
              "fixture must actually trigger a Phase C flip (sign_correction=TRUE's EBeta must be the exact negation of sign_correction=FALSE's) -- otherwise T8/T9 below are not meaningful")
})

run_test("YFBBetaCohort-T8: EBeta_pooled tracks EBeta's final orientation exactly (coherent-state fix, round 2)", {
  # Before the FIRST fix, EBeta_pooled was computed from rowMeans(EBeta) AFTER
  # Phase C's potential flip, while w/z/ZF (never touched by Phase C) stayed
  # at their pre-flip state -- an incoherent init that fed compute_pooled_beta()
  # a mismatched-sign partial residual, producing a genuinely different
  # (not just sign-flipped) result, not just here but on the cached D4 fit
  # (Program 7: 0.0404 vs. 0.0204, a 2x difference from this alone). Fixed by
  # computing EBeta_pooled's Gauss-Seidel sweep from the PRE-Phase-C EBeta
  # snapshot, so that sweep's inputs are always mutually consistent.
  #
  # That first fix alone left a SECOND bug (caught in code review, round 2):
  # the sweep's OUTPUT was never re-flipped to match Phase C's actual
  # decision, so EBeta_pooled stayed frozen at the pre-flip orientation while
  # EBeta itself got flipped -- an orientation MISMATCH between the two
  # (both are used interchangeably by predict_cox_on_yf() on external
  # cohorts, depending on whether the cohort's beta_cohort_id is known).
  # The correct invariant is therefore NOT "EBeta_pooled is identical
  # regardless of sign_correction" -- it's "EBeta_pooled is negated exactly
  # when EBeta is negated," i.e. it tracks the SAME final global sign.
  fit_true  <- suppressMessages(
    fit_cox_on_yf(.pc$Y, .pc$time, .pc$status, K = .pc$K, max_iter = .pc$max_iter,
                  verbose = FALSE, sign_correction = TRUE, beta_cohort_id = .pc$cohort2)
  )
  fit_false <- suppressMessages(
    fit_cox_on_yf(.pc$Y, .pc$time, .pc$status, K = .pc$K, max_iter = .pc$max_iter,
                  verbose = FALSE, sign_correction = FALSE, beta_cohort_id = .pc$cohort2)
  )
  # T7 already confirms this fixture makes Phase C flip (fit_true$EBeta ==
  # -fit_false$EBeta exactly) -- EBeta_pooled must show the same relationship.
  assert_near(fit_true$EBeta_pooled, -fit_false$EBeta_pooled, tol = 1e-10,
              "EBeta_pooled must be negated exactly when Phase C flips EBeta -- same final orientation, not independent of sign_correction")
})

run_test("YFBBetaCohort-T9: EBeta2 is unchanged by a Phase C sign flip (E[beta^2] invariant to negating beta)", {
  # Before the fix, Phase C set EBeta2 <- EBeta^2 using the NEWLY negated
  # EBeta -- i.e. replaced the posterior second moment (Var + mean^2) with
  # just the squared point estimate, silently discarding the posterior
  # variance. A sign flip must leave EBeta2 exactly as it was.
  fit_true  <- suppressMessages(
    fit_cox_on_yf(.pc$Y, .pc$time, .pc$status, K = .pc$K, max_iter = .pc$max_iter,
                  verbose = FALSE, sign_correction = TRUE, beta_cohort_id = .pc$cohort2)
  )
  fit_false <- suppressMessages(
    fit_cox_on_yf(.pc$Y, .pc$time, .pc$status, K = .pc$K, max_iter = .pc$max_iter,
                  verbose = FALSE, sign_correction = FALSE, beta_cohort_id = .pc$cohort2)
  )
  assert_near(fit_true$EBeta2, fit_false$EBeta2, tol = 1e-12,
              "EBeta2 must be identical whether or not Phase C flips EBeta's sign")
})
