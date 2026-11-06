# ============================================================
# Script: run_prompt7_ebmf.R
# Purpose: Compare Cox models on unsupervised EBMF programs at three ranks.
# Created: 2026-11-06
# Dependencies: flashier, survival and Prompt 6 helpers
# ============================================================

source("results/benchmark_sim/run_prompt6_comparators.R")
if (sys.nframe() == 0L) {
  d <- prompt6_reference()
  ranks <- c(7L, 15L, 40L)
  fit <- parallel::mclapply(ranks, function(k) {
    set.seed(42L)
    notices <- character()
    fit <- tryCatch(withCallingHandlers({
      f <- flashier::flash(d$Y, var_type = 2, greedy_Kmax = k, backfit = TRUE, verbose = 0)
      stopifnot(f$n_factors > 0)
      F <- flashier::ldf(f, type = "2")$F
      score <- d$Y %*% F
      cox <- coxph(Surv(d$time, d$status) ~ score)
      stopifnot(all(is.finite(coef(cox))))
      list(F = F, beta = coef(cox), rank = ncol(F), converged = f$converged)
    }, warning = function(w) notices <<- c(notices, conditionMessage(w))),
    error = function(e) list(error = conditionMessage(e)))
    fit$warnings <- paste(unique(notices), collapse = " | ")
    fit
  }, mc.cores = length(ranks), mc.set.seed = FALSE)
  stopifnot(all(vapply(fit, is.list, logical(1))))
  names(fit) <- paste0("EBMF_", ranks)
  status <- do.call(rbind, lapply(seq_along(fit), function(i) {
    f <- fit[[i]]
    data.frame(model = names(fit)[i], requested = ranks[i],
      achieved = if (is.null(f$rank)) NA_integer_ else f$rank,
      converged = if (is.null(f$converged)) NA else all(f$converged),
      error = if (is.null(f$error)) "" else f$error, warnings = f$warnings)
  }))
  cos <- cfg$pdac$external_cohorts
  ext <- setNames(lapply(cos, function(co) {
    x <- prompt6_external(co, d$genes)
    x$risks <- c(prompt6_risk(x$Y, d)[1], lapply(fit, function(f)
      if (is.null(f$F)) rep(NA_real_, nrow(x$Y)) else drop(x$Y %*% f$F %*% f$beta)))
    x
  }), cos)
  ans <- c(prompt6_evaluate(ext), list(ebmf_fits = status))
  out <- "results/benchmark_sim/outputs/prompt7_desurv_style"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(ans)) write.csv(ans[[nm]], file.path(out, paste0("rank_", nm, ".csv")), row.names = FALSE)
}
