# ============================================================
# Script: test_share_targeted_runner.R
# Purpose: Check false-active classification and risk variance calculations.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: results/benchmark_sim/run_share_targeted_weight.R
# ============================================================
source("results/benchmark_sim/run_share_targeted_weight.R")
stopifnot(identical(s16_real_arm_names(numeric(0)), c("w1", "twostep")),
          identical(s16_real_arm_names(.01), c("w1", "twostep", "s0.01")))
x <- cbind(c(-1, 0, 1), c(1, 0, -1) / 2)
a <- s16_false_variance(x, c(TRUE, FALSE))
stopifnot(isTRUE(all.equal(unname(a["ratio"]), 4)),
          isTRUE(all.equal(unname(a["signed"]), 2)))
stopifnot(identical(s16_false_active(c(TRUE, TRUE), c(.8, .1)), c(FALSE, TRUE)),
          identical(s16_false_active(c(TRUE, TRUE), c(.6, .1)), c(TRUE, TRUE)),
          identical(s16_false_active(c(TRUE, TRUE), c(.8, .1), TRUE), c(TRUE, TRUE)))
stopifnot(all(s16_false_variance(x, c(FALSE, FALSE)) == 0),
          all(is.na(s16_false_variance(cbind(x[, 1], -x[, 1]), c(TRUE, FALSE)))))
cat("S-16 runner metric checks passed\n")
null_rows <- data.frame(cell = rep("weak_null", 10), analysis = "primary", target = .1,
  seed = 1:10, capped_duplicate = FALSE, error = "", converged = TRUE,
  false_active = c(rep(1, 9), 0), false_risk_variance = rep(.009, 10))
stopifnot(s16_null_gate(null_rows, .1)$pass)
null_rows$false_risk_variance[1] <- NA_real_
stopifnot(!s16_null_gate(null_rows, .1)$pass)
null_rows$false_risk_variance[1] <- .02
stopifnot(!s16_null_gate(null_rows, .1)$pass)
# Isolate calibration calls to verify warm-start and cap-transition behavior.
local({
  source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
  calls <- list()
  share_targeted_ws <- function(fit, Y, time, status, target_share, w_current = 1, pilot_loading = NULL) {
    if (is.null(pilot_loading)) return(list(w_S = 1e5, capped = TRUE,
      no_signal_guard = FALSE, pilot_loading = c(1, 2)))
    if (w_current == 1e5) return(list(w_S = 100, capped = FALSE, ratio = .01,
      lost_activity = FALSE, matching_correlation = 1, selected_program = 1,
      diagnostics = data.frame(per_gene_median_share = .001)))
    list(w_S = 100, capped = FALSE, ratio = 1, lost_activity = FALSE,
      matching_correlation = 1, selected_program = 1,
      diagnostics = data.frame(per_gene_median_share = .002))
  }
  s16_fit <- function(Y, time, status, weight, cfg, initial, fit_control) {
    calls[[length(calls) + 1L]] <<- list(weight = weight, initial = initial)
    list(fit = list(history = list(converged = TRUE), posterior_state = paste0("state", length(calls))),
         warnings = character(), elapsed = 1)
  }
  pilot <- list(posterior_state = "pilot")
  x <- s16_calibrate(pilot, matrix(1, 2, 2), 1:2, c(1, 0), .5, list())
  stopifnot(x$passes == 2, x$capped, x$history[[1]]$weight == 1e5,
            x$history[[2]]$weight == 100, !is.null(x$history[[1]]$cap_state),
            identical(calls[[1]]$initial, "pilot"), identical(calls[[2]]$initial, "state1"))
  calls <- list()
  y <- s16_calibrate(pilot, matrix(1, 2, 2), 1:2, c(1, 0), .5, list(), previous_cap = "0.1")
  stopifnot(y$capped_duplicate, y$alias == "0.1", length(calls) == 0)
  s16_dataset <- function(d, cfg, targets, null, fit_control) {
    list(d = d, targets = targets, null = null, controls = fit_control)
  }
  smoke <- s16_smoke_run("seed99", list(), c(.01, .1), list(test = TRUE))
  stopifnot(identical(smoke$d, "seed99"), identical(smoke$targets, c(.01, .1)),
            identical(smoke$null, FALSE), isTRUE(smoke$controls$test))
})
local({
  source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
  s16_fit <- function(Y, time, status, weight, cfg, fit_control) {
    if (weight == 2) stop("deliberate paired failure")
    list(fit = list(history = list(converged = TRUE)), elapsed = 1, warnings = character())
  }
  s16_checkpoint <- function(fit, Y, time, status) list(fit = fit)
  s16_sim_score <- function(fit, d, null) data.frame(false_active = 0)
  sig <- list(targets = list(
    "0.1" = list(fit = list(history = list(converged = TRUE)), no_signal_guard = FALSE,
      history = list(list(weight = 2)), pilot_loading = c(1, 2)),
    "0.5" = list(fit = list(history = list(converged = TRUE)), no_signal_guard = FALSE,
      history = list(list(weight = 3)), pilot_loading = c(2, 1))))
  d <- list(Ytr = matrix(1, 2, 2), time = 1:2, event = c(1, 0))
  x <- s16_paired_null(sig, d, list())
  stopifnot(grepl("deliberate", x[[1]]$error), x[[2]]$weight == 3,
            identical(x[[2]]$pilot_loading, c(2, 1)))
})
# Fold centering, held-out projection and arm failures without fitting a model.
local({
  source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
  seen <- list()
  fake <- function(Y) list(EF = matrix(1, ncol(Y), 1), EBeta = 1,
    EBeta2 = 1.01, history = list(converged = TRUE, n_iter = 4))
  s16_fit <- function(Y, time, status, weight, cfg, initial = NULL, fit_control = list()) {
    seen[["train"]] <<- Y
    list(fit = fake(Y), elapsed = 1, warnings = character())
  }
  s16_twostep <- function(Y, time, status, cfg, fit_control) stop("two-step deliberate error")
  s16_calibrate <- function(...) stop("target deliberate error")
  s16_checkpoint <- function(fit, Y, time, status) list(fit = fit, shares = data.frame(
    active = TRUE, standardized_beta = 1, directional_share = .1, per_gene_share = .01))
  frozen_reverse_cindex <- function(risk, time, status) .6
  calc_cox_taylor_yf <- function(eta, time, status) list(logPL = sum(eta))
  d <- list(Y = matrix(c(1, 3, 5, 7, 11, 13), ncol = 1), time = 1:6,
    status = c(1, 0, 1, 0, 1, 0), outer = c(1, 1, 2, 3, 4, 5))
  x <- s16_real_task(1, d, list(), .1, list())
  stopifnot(abs(mean(seen$train)) < 1e-12, x$center == 9,
    identical(as.numeric(x$arms$w1$risk), c(-8, -6, -4, -2, 2, 4)),
    x$arms$w1$cvpl == -14, x$arms$w1$C == .6,
    grepl("two-step deliberate", x$arms$twostep$error),
    grepl("target deliberate", x$arms$s0.1$error))
  tab <- s16_real_rows(x, 1, c("w1", "twostep", "s0.1"))
  stopifnot(tab$cvpl[1] == -14, tab$achieved[1] == .1, tab$n_iter[1] == 4,
            nzchar(tab$error[2]), nzchar(tab$error[3]))
})
# Eligible weaker arm wins the registered .5 mean-CVPL tie, including w < 1.
local({
  source("results/benchmark_sim/run_share_targeted_weight.R", local = TRUE)
  rows <- expand.grid(fold = 0:5, arm = c("w1", "s0.01", "twostep"), stringsAsFactors = FALSE)
  rows$converged <- TRUE; rows$error <- ""; rows$C <- .6
  rows$cvpl <- ifelse(rows$arm == "w1", -1, ifelse(rows$arm == "s0.01", -1.4, -.9))
  rows$weight <- ifelse(rows$arm == "w1", 1, ifelse(rows$arm == "s0.01", .5, 0))
  stopifnot(s16_real_summary(rows)$selected == "s0.01")
  rows$converged[rows$arm == "s0.01" & rows$fold %in% 1:2] <- FALSE
  stopifnot(s16_real_summary(rows)$selected == "w1")
  rows$converged[rows$arm == "s0.01"] <- TRUE
  rows$cvpl[rows$arm == "s0.01" & rows$fold == 1] <- NA_real_
  stopifnot(s16_real_summary(rows)$selected == "w1")
})
