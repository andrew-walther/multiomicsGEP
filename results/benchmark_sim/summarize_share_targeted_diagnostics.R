# ============================================================
# Script: summarize_share_targeted_diagnostics.R
# Purpose: Rebuild aggregate S-16 diagnostic and tracked-program tables.
# Author: Andrew Walther
# Created: 2026-11-13
# Dependencies: yaml, saved S-16 diagnostics and prepared D4 training data
# ============================================================
# Patient-level fits and input stay under data/. Only aggregate tables are written.
private <- "data/share_targeted_ws_2026-11-13/diagnostics"
out <- "results/benchmark_sim/outputs/share_targeted_ws_2026-11-13"
tol <- yaml::read_yaml("config/globals.yml")$cavi$tol
paths <- list.files(private, pattern = paste0(
  "^(simulation_w[0-9]+|real_pilot_w1|real_target_(0\\.01|0\\.05|0\\.1|0\\.25|0\\.5))",
  "(_ramp50)?\\.rds$"), full.names = TRUE)
if (length(paths) != 16L) stop("Expected 16 complete S-16 diagnostic fit checkpoints")

#' Summarize one saved diagnostic fit and its applicable objective changes.
#' @param path Checkpoint path. @return One aggregate data-frame row.
#' @examples
#' # s16_summary_row("data/share_targeted_ws_2026-11-13/diagnostics/simulation_w10.rds")
s16_summary_row <- function(path) {
  x <- readRDS(path); f <- x$fit
  a <- f$history$s16[f$history$s16$program == 1, ]
  post <- if ("ramp" %in% names(a) && any(a$ramp)) a[!a$ramp, ] else a
  if (nrow(a) != f$history$n_iter || !any(is.finite(post$elbo_rel_signed)))
    stop("Incomplete objective trace: ", path)
  smallest <- min(post$elbo_rel_signed, na.rm = TRUE)
  decrease <- max(0, -smallest)
  data.frame(case = sub("\\.rds$", "", basename(path)), weight = x$weight,
    converged = f$history$converged, sweeps = f$history$n_iter, seconds = x$elapsed,
    warning_count = length(x$warnings), min_signed_elbo_change = smallest,
    max_eta = max(a$eta_end_max), max_beta = max(a$beta_max),
    last_eta_relative_change = tail(a$eta_rel_change, 1),
    max_relative_decrease = decrease,
    meets_convergence_and_decrease = f$history$converged && decrease <= tol)
}

#' Track a real-data pilot loading and recompute its final directional share.
#' @param path Fit checkpoint. @param Y Prepared training matrix.
#' @return One aggregate data-frame row. @examples
#' # s16_tracked_row(path, Y)
s16_tracked_row <- function(path, Y) {
  x <- readRDS(path); f <- x$fit; name <- sub("\\.rds$", "", basename(path))
  target <- sub("_ramp50$", "", sub("^real_target_", "", name))
  calibration <- readRDS(file.path(private, paste0("real_target_", target, "_calibration.rds")))
  corr <- abs(as.vector(cor(calibration$pilot_loading, f$EF)))
  if (all(!is.finite(corr))) stop("No defined pilot loading match: ", path)
  k <- which.max(replace(corr, !is.finite(corr), -Inf))
  beta <- f$EBeta[k]; beta_sd <- sqrt(max(0, f$EBeta2[k] - beta^2))
  direction <- f$EF[, k] / sqrt(sum(f$EF[, k]^2))
  expression <- sum(f$EL2[, k]) * sum(f$Tau * direction^2)
  survival <- x$weight * f$EBeta2[k] * sum(f$CoxW * as.vector(Y %*% direction)^2)
  share <- survival / (expression + survival)
  recorded <- tail(f$history$s16$directional_share[f$history$s16$program == k], 1)
  if (!isTRUE(all.equal(share, recorded, tolerance = 1e-12)))
    stop("Saved Cox weights disagree with final trace: ", path)
  data.frame(case = name, target = as.numeric(target), pilot_program = calibration$selected_program,
    matched_program = k, abs_loading_correlation = corr[k], directional_share = share,
    active = beta_sd > 0 && abs(beta) >= 1.96 * beta_sd, converged = f$history$converged)
}

summary <- do.call(rbind, lapply(paths, s16_summary_row))
z <- readRDS("data/yfb_supervised_real_2026-11-04.rds")
Y <- z$Y[, match(z$sets$D4, z$genes), drop = FALSE]
if (any(!is.finite(Y))) stop("D4 features missing or non-finite")
Y <- scale(Y, center = TRUE, scale = FALSE)
real_paths <- grep("/real_target_.*[0-9](_ramp50)?\\.rds$", paths, value = TRUE)
tracked <- do.call(rbind, lapply(real_paths, s16_tracked_row, Y = Y))
write.csv(summary, file.path(out, "diagnostic_summary.csv"), row.names = FALSE)
write.csv(tracked, file.path(out, "diagnostics_tracked_programs.csv"), row.names = FALSE)
cat(sprintf("Summarized %d fits and %d tracked real-data fits.\n", nrow(summary), nrow(tracked)))
