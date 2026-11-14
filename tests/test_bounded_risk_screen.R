# ============================================================
# Script: test_bounded_risk_screen.R
# Purpose: Check S-17 stability, prior specifications and checkpoint behavior.
# Author: Andrew Walther
# Created: 2026-11-14
# Dependencies: ashr, ebnm
# ============================================================
local({
source("results/benchmark_sim/run_bounded_risk_screen.R", local = TRUE)

#' Create a deterministic two-program trace for gate checks
#' @param n Number of traced sweeps.
#' @return Minimal fit with a stable trace and two coefficients.
#' @examples
#' # test_fit(50)
test_fit <- function(n = 50L) {
  tr <- data.frame(iteration = rep(seq_len(n), each = 2L),
    eta_rel_change = .0001, eta_start_max = 1, eta_end_max = 1, beta_max = .3)
  list(history = list(s16 = tr, converged = TRUE, native_elbo_converged = TRUE), EBeta = c(.2, -.3))
}
x <- test_fit()
stopifnot(s17_stability(x, .1)$pass)
x$history$s16$eta_rel_change[1:2] <- .001
stopifnot(!s17_stability(x, .1)$pass)
x <- test_fit(); x$history$s16$eta_end_max[99:100] <- 1.01
stopifnot(!s17_stability(x, .1)$pass)
x <- test_fit(49); stopifnot(!s17_stability(x, .1)$assessable, !s17_stability(x, .1)$pass)
x <- test_fit(); x$history$converged <- FALSE; stopifnot(!s17_stability(x, .1)$pass)
x <- test_fit(); x$history$native_elbo_converged <- FALSE; stopifnot(!s17_stability(x, .1)$pass)
x <- test_fit(); x$EBeta[1] <- .40001; stopifnot(!s17_stability(x, .1)$beta_pass)
x$EBeta[1] <- .4; stopifnot(s17_stability(x, .1)$beta_pass)
x <- test_fit(); x$history$s16$eta_rel_change[1] <- .5
stopifnot(inherits(try(s17_stability(x, .1), silent = TRUE), "try-error"))
a <- s17_controls("normal", .5); b <- s17_controls("laplace", .5)
stopifnot(a$prior_beta == "normal", a$beta_prior_init$sd == .5,
          b$prior_beta == "point_laplace", b$beta_prior_init$scale == .5 / sqrt(2),
          a$beta_prior_fix, b$beta_prior_fix, a$trace_s16, b$trace_s16,
          a$stable_convergence_window == 50L, b$stable_convergence_window == 50L)
cat("S-17 screen gate and prior checks passed\n")

# Check that all calibration passes, including a nonconverged last pass, are
# retained and that a stale source signature cannot silently reuse the result.
s16_fit <- function(Y, time, status, weight, cfg, initial = NULL, fit_control = list())
  list(fit = test_fit(), warnings = character(), elapsed = 1)
s16_calibrate <- function(pilot, Y, time, status, target, cfg, fit_control,
                          retain_pass_fits) {
  stopifnot(retain_pass_fits, fit_control$beta_prior_fix)
  final <- test_fit(); final$history$converged <- FALSE
  list(fit = final, pass_fits = list(test_fit(), final), no_signal_guard = FALSE,
    history = list(list(weight = 100, capped = FALSE, warnings = character(), elapsed = 2),
                   list(weight = 1e5, capped = TRUE, warnings = "limit", elapsed = 3)))
}
s16_checkpoint <- function(fit, Y, time, status)
  list(fit = fit, cox_weights = rep(1, length(time)), shares = data.frame(program = 1))
s16_atomic_save <- function(object, path) saveRDS(object, path)
tmp <- tempfile(); dir.create(tmp)
d <- list(Y = matrix(0, 2, 2), time = 1:2, status = c(1, 1))
obj <- s17_case("real", "normal", .5, d, list(), list(hash = "a"), tmp)
stopifnot(length(obj$arms$s0.05$fits) == 2L,
          !obj$arms$s0.05$fits[[2]]$history$converged,
          obj$arms$s0.05$capped[2], length(obj$arms$s0.05$checkpoints) == 2L)
out <- s17_export(obj, tmp)
stopifnot(nrow(out$rows) == 7L, sum(out$rows$final) == 4L,
          all(out$rows$capped[out$rows$pass_number == 2]),
          !any(out$rows$pass[out$rows$pass_number == 2]))
stopifnot(inherits(try(s17_case("real", "normal", .5, d, list(),
                               list(hash = "b"), tmp), silent = TRUE), "try-error"))
s16_fit <- function(Y, time, status, weight, cfg, initial = NULL, fit_control = list()) {
  fit <- test_fit(); fit$history$converged <- FALSE
  list(fit = fit, warnings = character(), elapsed = 1)
}
failed <- s17_case("weak", "normal", .5, d, list(), list(hash = "a"), tmp)
failed_rows <- s17_export(failed, tmp)$rows
stopifnot(nrow(failed_rows) == 4L, all(!failed_rows$pass),
          all(failed_rows$error[failed_rows$arm != "pilot"] == "pilot did not converge"))
s16_fit <- function(Y, time, status, weight, cfg, initial = NULL, fit_control = list())
  list(fit = test_fit(), warnings = character(), elapsed = 1)
s16_calibrate <- function(pilot, Y, time, status, target, cfg, fit_control, retain_pass_fits)
  list(fit = pilot, pass_fits = list(), no_signal_guard = TRUE, history = list())
guarded <- s17_case("real", "laplace", .5, d, list(), list(hash = "a"), tmp)
guard_rows <- s17_export(guarded, tmp)$rows
stopifnot(all(guard_rows$pass), all(guard_rows$guard[guard_rows$arm != "pilot"]),
          all(guard_rows$weight == 1))
unlink(tmp, recursive = TRUE)
cat("S-17 pass retention, pilot failure, guard, and checkpoint checks passed\n")

})
