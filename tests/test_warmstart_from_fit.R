# ============================================================
# Script:  test_warmstart_from_fit.R
# Purpose: Tests for extract_top_k_by_pve() (code/warmstart_from_fit.R)
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Usage:   Rscript tests/run_tests.R  (included via source)
#          Rscript tests/test_warmstart_from_fit.R  (standalone)
# ============================================================

if (!exists("run_test")) source("tests/test_helpers.R")
if (!exists("extract_top_k_by_pve")) source("code/warmstart_from_fit.R")

cat("\n--- test_warmstart_from_fit.R ---\n")

# Build a fake converged-fit object: EL/EF are n x 5 / p x 5, with columns
# 1..5 given decreasing PVE by construction (so ranking is unambiguous).
make_fake_fit <- function(n = 10, p = 12, K = 5) {
  set.seed(1)
  EL <- matrix(rnorm(n * K), n, K)
  EF <- matrix(rnorm(p * K), p, K)
  # Normalize both EL and EF to unit-norm columns first, so PVE
  # (sum(EL[,k]^2) * sum(EF[,k]^2)) is controlled EXACTLY (not just in
  # expectation) by the scale_k multiplier applied below -- with n, p this
  # small, unnormalized random column norms have enough sampling variance to
  # overturn the intended ranking.
  EL <- sweep(EL, 2, sqrt(colSums(EL^2)), "/")
  EF <- sweep(EF, 2, sqrt(colSums(EF^2)), "/")
  # Scale columns so PVE is strictly decreasing: column 1 largest, column 5 smallest.
  scale_k <- c(10, 8, 6, 4, 2)
  EL <- sweep(EL, 2, scale_k, "*")
  pve <- vapply(seq_len(K), function(k) sum(EL[, k]^2) * sum(EF[, k]^2), numeric(1))
  list(
    EL = EL, EF = EF,
    history = list(factor_pve = matrix(pve, nrow = 1, ncol = K))
  )
}

# T1: K_target == K_source returns all columns (dimensions preserved)
run_test("WS-T1: K_target == K_source preserves dimensions", {
  fit <- make_fake_fit()
  out <- extract_top_k_by_pve(fit, 5)
  assert_equal(dim(out$EL_init), c(10L, 5L), "EL_init should be n x K_source")
  assert_equal(dim(out$EF_init), c(12L, 5L), "EF_init should be p x K_source")
})

# T2: K_target < K_source picks the highest-PVE columns, by construction
#     columns 1,2,3 (in that order, since PVE strictly decreases 1..5)
run_test("WS-T2: K_target < K_source selects top-PVE columns in order", {
  fit <- make_fake_fit()
  out <- extract_top_k_by_pve(fit, 3)
  assert_near(out$EL_init, fit$EL[, 1:3], tol = 1e-10,
              "EL_init should be fit$EL columns 1:3 (highest PVE)")
  assert_near(out$EF_init, fit$EF[, 1:3], tol = 1e-10,
              "EF_init should be fit$EF columns 1:3 (highest PVE)")
})

# T3: ranking is by PVE, not column order -- shuffle PVE and confirm the
#     right columns (not just the first K) are picked
run_test("WS-T3: ranking follows PVE values, not column position", {
  fit <- make_fake_fit()
  # Reverse the PVE order so column 5 is now highest, column 1 lowest.
  fit$history$factor_pve <- matrix(rev(as.vector(fit$history$factor_pve)), nrow = 1)
  out <- extract_top_k_by_pve(fit, 2)
  assert_near(out$EL_init, fit$EL[, c(5, 4)], tol = 1e-10,
              "Should select columns 5,4 (highest PVE after reversal)")
})

# T4: K_target > K_source raises an error
run_test("WS-T4: K_target > K_source raises error", {
  fit <- make_fake_fit()
  err <- tryCatch({ extract_top_k_by_pve(fit, 6); NULL },
                  error = function(e) conditionMessage(e))
  assert_true(is.character(err), "K_target > K_source should raise an error")
})

# T5: non-integer / non-positive K_target raises an error
run_test("WS-T5: invalid K_target raises error", {
  fit <- make_fake_fit()
  err1 <- tryCatch({ extract_top_k_by_pve(fit, 0); NULL },
                   error = function(e) conditionMessage(e))
  err2 <- tryCatch({ extract_top_k_by_pve(fit, 2.5); NULL },
                   error = function(e) conditionMessage(e))
  assert_true(is.character(err1), "K_target=0 should raise an error")
  assert_true(is.character(err2), "K_target=2.5 should raise an error")
})

# T6: uses the FINAL iteration's PVE row when history has multiple rows
#     (e.g. a fit that ran several iterations) -- not the first row.
run_test("WS-T6: uses final-iteration PVE row, not the first", {
  fit <- make_fake_fit()
  # Prepend a bogus first-iteration row with the opposite ranking; only the
  # last row should be used for selection.
  bogus_first_row <- rev(as.vector(fit$history$factor_pve))
  fit$history$factor_pve <- rbind(bogus_first_row, fit$history$factor_pve)
  out <- extract_top_k_by_pve(fit, 3)
  assert_near(out$EL_init, fit$EL[, 1:3], tol = 1e-10,
              "Should rank by the LAST row of factor_pve, ignoring earlier rows")
})
