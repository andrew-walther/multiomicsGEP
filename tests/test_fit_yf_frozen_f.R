# =============================================================================
# tests/test_fit_yf_frozen_f.R
# Unit tests for the N_frozen (frozen-F β pre-conditioning) parameter in
# fit_cox_on_yf().
#
# The frozen-F mechanism holds EF fixed at its SVD initialization for the first
# N_frozen CAVI iterations, letting β and L update freely. This breaks the
# β=0 ↔ B_beta=0 fixed-point that traps the YFB model on merged multi-platform
# data where genomics-only EF converges to a platform-contrast direction.
# =============================================================================

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------
local({
  set.seed(77)
  n <- 40L; p <- 20L; K <- 2L
  .ff <<- list(
    n = n, p = p, K = K,
    Y      = matrix(rnorm(n * p), n, p),
    time   = rexp(n, 0.1),
    status = rbinom(n, 1, 0.7)
  )
})

# ---------------------------------------------------------------------------
# T1: N_frozen=0 is backward compatible — runs without error
# ---------------------------------------------------------------------------
run_test("FrozenF-T1: N_frozen=0 (default) runs without error", {
  fit <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 5L,
                  N_frozen = 0L, verbose = FALSE)
  )
  assert_true(length(fit$EBeta)    == .ff$K, msg = "EBeta length K")
  assert_true(nrow(fit$EF)         == .ff$p, msg = "EF rows == p")
  assert_true(ncol(fit$EF)         == .ff$K, msg = "EF cols == K")
  assert_true(is.numeric(fit$EBeta),          msg = "EBeta is numeric")
})

# ---------------------------------------------------------------------------
# T2: N_frozen>0 runs without error and returns valid output
# ---------------------------------------------------------------------------
run_test("FrozenF-T2: N_frozen=3 runs and returns valid output", {
  fit <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 10L,
                  N_frozen = 3L, verbose = FALSE)
  )
  assert_true(length(fit$EBeta)    == .ff$K, msg = "EBeta length K")
  assert_true(nrow(fit$EF)         == .ff$p, msg = "EF rows == p")
  assert_true(ncol(fit$EF)         == .ff$K, msg = "EF cols == K")
  assert_true(length(fit$EF_norms) == .ff$K, msg = "EF_norms length K")
  assert_true(!any(is.na(fit$EBeta)),         msg = "EBeta has no NAs")
  assert_true(!any(is.na(fit$EF)),            msg = "EF has no NAs")
})

# ---------------------------------------------------------------------------
# T3: ELBO is finite after a run with N_frozen>0
# ---------------------------------------------------------------------------
run_test("FrozenF-T3: ELBO is finite with N_frozen=5", {
  fit <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 8L,
                  N_frozen = 5L, verbose = FALSE)
  )
  assert_true(all(is.finite(fit$history$elbo_full)),
              msg = "All ELBO values are finite")
})

# ---------------------------------------------------------------------------
# T4: N_frozen >= max_iter is allowed — EF never unfreezes; result is valid
# ---------------------------------------------------------------------------
run_test("FrozenF-T4: N_frozen >= max_iter runs without error", {
  fit <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 5L,
                  N_frozen = 100L, verbose = FALSE)
  )
  assert_true(length(fit$EBeta) == .ff$K, msg = "EBeta length K")
  assert_true(!any(is.na(fit$EBeta)),      msg = "EBeta has no NAs")
})

# ---------------------------------------------------------------------------
# T5: N_frozen=-1 raises an error
# ---------------------------------------------------------------------------
run_test("FrozenF-T5: N_frozen=-1 raises an error", {
  caught <- tryCatch({
    suppressMessages(
      fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                    K = .ff$K, max_iter = 3L,
                    N_frozen = -1L, verbose = FALSE)
    )
    FALSE
  }, error = function(e) TRUE)
  assert_true(caught, msg = "N_frozen=-1 should raise an error")
})

# ---------------------------------------------------------------------------
# T6: N_frozen="a" (non-numeric) raises an error
# ---------------------------------------------------------------------------
run_test("FrozenF-T6: N_frozen='a' raises an error", {
  caught <- tryCatch({
    suppressMessages(
      fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                    K = .ff$K, max_iter = 3L,
                    N_frozen = "a", verbose = FALSE)
    )
    FALSE
  }, error = function(e) TRUE)
  assert_true(caught, msg = "N_frozen='a' should raise an error")
})

# ---------------------------------------------------------------------------
# T7: N_frozen=0 vs N_frozen=5 produce different EF matrices
#     (frozen-F phase changes the trajectory of EF adaptation post-unfreeze)
# ---------------------------------------------------------------------------
run_test("FrozenF-T7: N_frozen=0 and N_frozen=5 produce different EF", {
  fit0 <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 15L,
                  N_frozen = 0L, verbose = FALSE)
  )
  fit5 <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 15L,
                  N_frozen = 5L, verbose = FALSE)
  )
  ef_diff <- max(abs(fit0$EF - fit5$EF))
  assert_true(ef_diff > 1e-10, msg = "EF should differ between N_frozen=0 and N_frozen=5")
})

# ---------------------------------------------------------------------------
# T9: convergence cannot fire while F is still frozen. Before the fix, beta/L
#     can plateau during the frozen phase (F held fixed, nothing forcing
#     further movement) and the ELBO convergence check -- guarded only by
#     `iter > 5`, not `iter > N_frozen` -- declares convergence before F ever
#     gets a chance to unfreeze, silently defeating the entire mechanism.
#     Reproduced on real D4-preprocessed PDAC data (DECISIONS.md 2026-08-20):
#     N_frozen in {10,20,30} all converged at iter=8, beta stuck at exactly 0.
# ---------------------------------------------------------------------------
run_test("FrozenF-T9: converged run must have n_iter > N_frozen", {
  fit <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 60L,
                  N_frozen = 50L, verbose = FALSE)
  )
  assert_true(fit$history$converged,
              msg = "Fixture no longer converges within max_iter -- update this test's max_iter/tol so it keeps exercising the guard")
  assert_true(fit$history$n_iter > 50L,
              msg = "Convergence fired at n_iter <= N_frozen -- F never unfroze")
})

# ---------------------------------------------------------------------------
# T8: N_frozen works with cohort_id — cohort F columns also unfrozen properly
# ---------------------------------------------------------------------------
run_test("FrozenF-T8: N_frozen=3 with cohort_id runs without error", {
  cohort <- rep(c("A", "B"), each = .ff$n / 2)
  fit <- suppressMessages(
    fit_cox_on_yf(.ff$Y, .ff$time, .ff$status,
                  K = .ff$K, max_iter = 10L,
                  N_frozen = 3L, verbose = FALSE,
                  cohort_id = cohort)
  )
  assert_true(length(fit$EBeta)    == .ff$K, msg = "EBeta length K")
  assert_true(!is.null(fit$EF_cohort),        msg = "EF_cohort present")
  assert_true(nrow(fit$EF_cohort)  == .ff$p,  msg = "EF_cohort rows == p")
})
