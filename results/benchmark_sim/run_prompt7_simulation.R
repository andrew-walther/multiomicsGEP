# ============================================================
# Script: run_prompt7_simulation.R
# Purpose: Audit saved classification checkpoints and run the mixed-signal simulation.
# Created: 2026-11-06
# Dependencies: survival, ebnm, yaml, parallel and S-13/S-14 runners
# ============================================================
# Synthetic data only. One fit per fork worker, one BLAS thread per worker.

source("results/benchmark_sim/run_prompt7_real.R")
source("results/benchmark_sim/run_yfb_survival_weight_p2000.R")
s9_load()
source("results/benchmark_sim/run_yfb_supervision_strength.R")
source("results/benchmark_sim/run_yfb_near_threshold.R")

#' Generate the S-14-style mixed-signal expression and outcomes.
#' @param seed Replicate number.
#' @return Training and validation matrices and outcomes, true directions and shares.
#' @examples
#' # prompt7_mixed_data(1)
prompt7_mixed_data <- function(seed) {
  d <- simulate_multimodal_yfb_data("low_variance_prognostic_factor",
    n_train = 400L, n_validation = 120L, p_expression = 2000L, p_methylation = 50L,
    seed = seed, truncate = FALSE, noise_scale = .25)
  F <- d$truth$F$expression
  idx <- c(1L, 3L)
  set.seed(seed + 4000L)
  extra <- simulate_multimodal_yfb_loadings(nrow(F), 1L)[, 1]
  S <- lapply(idx, function(k) outer(d$training$L[, k], F[, k]))
  Sv <- lapply(idx, function(k) outer(d$validation$L[, k], F[, k]))
  base <- d$training$Y$expression - S[[1]] - S[[2]] +
    outer(rexp(nrow(d$training$L)), extra)
  val <- d$validation$Y$expression - Sv[[1]] - Sv[[2]] +
    outer(rexp(nrow(d$validation$L)), extra)
  target <- c(.05, .001)
  centered <- lapply(S, scale, scale = FALSE)
  residual <- scale(base, scale = FALSE)
  a <- c(1, 1)
  # Alternate the S-14 scalar quadratic updates until both shares agree.
  for (pass in seq_len(12L)) for (k in 1:2) {
    other <- residual + a[3L - k] * centered[[3L - k]]
    aa <- (1 - target[k]) * sum(centered[[k]]^2)
    bb <- -2 * target[k] * sum(centered[[k]] * other)
    cc <- -target[k] * sum(other^2)
    a[k] <- (-bb + sqrt(bb^2 - 4 * aa * cc)) / (2 * aa)
  }
  raw <- base + a[1] * S[[1]] + a[2] * S[[2]]
  raw_v <- val + a[1] * Sv[[1]] + a[2] * Sv[[2]]
  share <- a^2 * vapply(centered, function(x) sum(x^2), numeric(1)) /
    sum(scale(raw, scale = FALSE)^2)
  stopifnot(max(abs(share - target)) < 1e-10)
  truth <- sweep(F[, idx, drop = FALSE], 1, apply(raw, 2, sd), "/")
  truth <- sweep(truth, 2, sqrt(colSums(truth^2)), "/")
  Y <- s9_zscore(raw)
  Yv <- s9_zscore(raw_v)
  z <- Y %*% truth
  zv <- Yv %*% truth
  zsd <- apply(z, 2, sd)
  eta <- rowSums(sweep(z, 2, zsd, "/"))
  etav <- rowSums(sweep(zv, 2, zsd, "/"))
  multiplier <- 1.1 / sd(eta)
  set.seed(seed + 5000L)
  draw <- function(risk) {
    event_time <- (-log(runif(length(risk))) / (.01 * exp(risk)))^(1 / 1.5)
    censor <- rexp(length(risk))
    censor <- censor * multimodal_yfb_calibrate_censoring(event_time, censor, .3)
    list(time = pmin(event_time, censor), event = as.integer(event_time <= censor))
  }
  tr <- draw(eta * multiplier)
  va <- draw(etav * multiplier)
  list(Y = Y, Yv = Yv, truth = truth, share = share, tr = tr, va = va, oracle = etav * multiplier)
}

#' Measure fixed top-50 gene precision for one matched direction.
#' @param loading Fitted loading vector. @param truth Planted loading vector.
#' @return Number of top genes in the planted support divided by 50.
#' @examples
#' # prompt7_top50(seq_len(100), c(rep(0, 50), rep(1, 50)))
prompt7_top50 <- function(loading, truth) {
  stopifnot(length(loading) == length(truth), length(truth) >= 50L)
  top <- head(order(-abs(loading), seq_along(loading)), 50)
  sum(truth[top] != 0) / 50
}

#' Fit and score one mixed-signal replicate and arm.
#' @param task Design row. @param cfg Project configuration.
#' @return Two rows per fit, one per planted program, with metrics and fit status.
#' @examples
#' # prompt7_mixed_fit(task, cfg)
prompt7_mixed_fit <- function(task, cfg) {
  d <- prompt7_mixed_data(task$seed)
  notices <- character()
  arm <- task$arm
  ans <- tryCatch(withCallingHandlers({
    set.seed(task$seed)
    w <- c(twostep = 0, A1_w1 = 1, A1_w5 = 5, A2_w1 = 1)[[arm]]
    f <- s13_fit(d$Y, d$tr$time, d$tr$event, w,
                 if (arm == "A2_w1") "joint" else "meanfield", cfg)
    F <- f$EF
    if (arm == "twostep") {
      live <- which(colSums(F^2) > 1e-8)
      z <- d$Y %*% F[, live, drop = FALSE]
      cf <- coxph(Surv(d$tr$time, d$tr$event) ~ z, control = coxph.control(timefix = FALSE))
      stopifnot(all(is.finite(coef(cf))))
      beta <- numeric(ncol(F))
      beta[live] <- coef(cf)
      active <- seq_along(beta) %in% live[abs(coef(cf) / sqrt(diag(vcov(cf)))) >= 1.96]
    } else {
      beta <- f$EBeta
      active <- prompt7_active(f)
    }
    m <- prompt7_matches(d$truth, F)
    m <- m[match(1:2, m$true_factor), ]
    matched <- m$estimated_factor
    risk <- drop(d$Yv %*% F %*% beta)
    qualifies <- !is.na(matched) & m$abs_correlation >= .7
    recovered <- qualifies & active[matched] & beta[matched] * m$correlation > 0
    precision <- vapply(1:2, function(k)
      if (is.na(matched[k])) NA_real_ else prompt7_top50(F[, matched[k]], d$truth[, k]),
      numeric(1))
    data.frame(direction = c("strong", "near_threshold"), planted_share = d$share,
      correlation = m$abs_correlation, recovered = recovered, top50_overlap = 50 * precision,
      top50_precision = precision, false_active = sum(active) - sum(active[matched[qualifies]]),
      C = frozen_reverse_cindex(risk, d$va$time, d$va$event),
      oracle_C = frozen_reverse_cindex(d$oracle, d$va$time, d$va$event),
      converged = f$history$converged, sweeps = f$history$n_iter, error = "")
  }, warning = function(w) notices <<- c(notices, conditionMessage(w))),
  error = function(e) data.frame(direction = c("strong", "near_threshold"), planted_share = d$share,
    correlation = NA_real_, recovered = FALSE, top50_overlap = NA_real_, top50_precision = NA_real_,
    false_active = NA_integer_, C = NA_real_, oracle_C = NA_real_,
    converged = FALSE, sweeps = NA_integer_, error = conditionMessage(e)))
  cbind(task, ans, warnings = paste(unique(notices), collapse = " | "))
}

#' Apply the coefficient-and-share rule to a saved single-signal fit.
#' @param Y Prepared training matrix. @param F,beta Saved estimates.
#' @param truth Planted loading. @param active Activity mask or NULL when unavailable.
#' @return One row with planted retention, overlap, and false-active coverage.
#' @examples
#' # prompt7_saved_rule(Y, F, beta, truth, active)
prompt7_saved_rule <- function(Y, F, beta, truth, active) {
  p <- prompt7_contribution(Y, F, beta)
  partial <- with(p, !is.na(ratio) & ratio >= .25 & variance_share >= .05)
  m <- prompt7_matches(matrix(truth, ncol = 1L), F)
  k <- if (is.null(m)) NA_integer_ else m$estimated_factor[1]
  qualifies <- !is.na(k) && m$abs_correlation[1] >= .7
  known <- !is.null(active)
  false <- if (known) setdiff(which(active), if (qualifies) k else integer()) else integer()
  data.frame(matched = k, matched_ge_70 = qualifies,
    planted_pass_partial = qualifies && partial[k],
    planted_active = if (known) qualifies && active[k] else NA,
    false_active_denominator = if (known) length(false) else NA_integer_,
    false_active_dropped_partial = if (known) sum(!partial[false]) else NA_integer_,
    top50_precision = if (!is.na(k)) prompt7_top50(F[, k], truth) else NA_real_,
    fold_evidence = FALSE)
}

if (sys.nframe() == 0L) {
  cfg <- yaml::read_yaml("config/globals.yml")
  dir.create(OUT7, recursive = TRUE, showWarnings = FALSE)
  if ("mixed" %in% commandArgs(TRUE)) {
    design <- expand.grid(seed = 1:10, arm = c("twostep", "A1_w1", "A1_w5", "A2_w1"),
                          stringsAsFactors = FALSE)
    stopifnot(prompt7_top50(1:100, c(rep(0, 50), rep(1, 50))) == 1)
    work <- "data/prompt7_mixed"
    dir.create(work, showWarnings = FALSE)
    results <- parallel::mclapply(seq_len(nrow(design)), function(i) {
      path <- file.path(work, paste0("task_", i, ".rds"))
      if (file.exists(path)) return(readRDS(path))
      out <- prompt7_mixed_fit(design[i, ], cfg)
      saveRDS(out, path)
      out
    }, mc.cores = 8L, mc.preschedule = FALSE)
    stopifnot(all(vapply(results, is.data.frame, logical(1))))
    write.csv(do.call(rbind, results), file.path(OUT7, "mixed_simulation.csv"), row.names = FALSE)
  } else {
    s13 <- "results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05"
    s14 <- "results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05"
    pick <- read.csv(file.path(s13, "cv_selected.csv"))
    x <- do.call(rbind, lapply(seq_len(nrow(pick)), function(i) {
      t <- pick[i, ]
      f <- readRDS(file.path(s13, "fits", sprintf("sig_%s_s%02d_%s_w%g.rds",
        t$cell, t$seed, t$arm, t$w_S)))
      d <- s13_data(t$seed, t$cell)
      a <- prompt7_active(f$state)
      cbind(t[, c("cell", "seed", "arm", "w_S", "converged")],
        prompt7_saved_rule(d$Ytr, f$state$EF, f$state$EBeta, d$truth, a))
    }))
    metrics <- read.csv(file.path(s14, "metrics.csv"))
    y <- do.call(rbind, lapply(seq_len(nrow(metrics)), function(i) {
      t <- metrics[i, ]
      f <- readRDS(file.path(s14, "fits", sprintf("%s_sh%g_s%02d_%s.rds",
        t$strength, t$share, t$seed, t$arm)))
      d <- s14_data(t$seed, t$share, t$strength)
      a <- NULL
      if (t$arm == "twostep") {
        live <- which(colSums(f$EF^2) > 1e-8)
        z <- d$Ytr %*% f$EF[, live, drop = FALSE]
        cf <- coxph(Surv(d$time, d$event) ~ z, control = coxph.control(timefix = FALSE))
        stopifnot(max(abs(coef(cf) - f$beta[live])) < 1e-8)
        a <- seq_along(f$beta) %in% live[abs(coef(cf) / sqrt(diag(vcov(cf)))) >= 1.96]
      }
      cbind(t[, c("strength", "share", "seed", "arm", "converged", "false_active")],
        prompt7_saved_rule(d$Ytr, f$EF, f$beta, d$truth, a))
    }))
    write.csv(x, file.path(OUT7, "s13_rule.csv"), row.names = FALSE)
    write.csv(y, file.path(OUT7, "s14_rule.csv"), row.names = FALSE)
  }
}
