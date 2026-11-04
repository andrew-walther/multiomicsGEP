# ============================================================
# Script: summarize_yfb_survival_weight_p2000.R
# Purpose: Attempt S-9 tables, Gate S verdict, CV-chosen arms and per-weight curves.
# Author: Andrew Walther
# Created: 2026-11-04
# Dependencies: base R, grDevices
# Run: Rscript results/benchmark_sim/summarize_yfb_survival_weight_p2000.R [grid|smoke]
# ============================================================
# Reads checkpoints only; refits nothing. Gate S as registered in
# docs/plans/2026-11-03-survival-weight-p2000-plan.md: at p = 2,000 a single
# arm qualifies in >= 2 of 3 variance cells, each with recovery >= 3/10 above
# the reference, mean false-active <= 1 and >= 9/10 converged.

args <- commandArgs(TRUE)
sub  <- if (length(args)) args[1] else "grid"
root <- "results/benchmark_sim/outputs/yfb_survival_weight_2026-11-04"
dir  <- file.path(root, sub)
s7   <- "results/benchmark_sim/outputs/yfb_supervision_2026-10-30/attempt7"

# Load ----

#' Read every S-9 checkpoint's metric row, plus the expected design
#' @param dir Output directory.
#' @return Data frame of fits; registered fits that are missing become failure rows.
#' @examples
#' # s9_read(dir)
s9_read <- function(dir) {
  files <- list.files(dir, pattern = "_(w[0-9]+|p_over_50|share0\\.2)\\.rds$", full.names = TRUE)
  m <- do.call(rbind, lapply(files, function(f) readRDS(f)$metrics))
  design <- read.csv(file.path(dir, "design.csv"))
  arms <- c("w0", paste0("w", c(1, 10, 100, 1000, 10000)), "p_over_50", "share0.2")
  exp <- merge(design, data.frame(arm = arms))
  exp <- exp[!(exp$p == 50 & exp$arm == "p_over_50"), ]
  key <- function(x) paste(x$seed, x$p, x$variance, x$arm)
  miss <- exp[!key(exp) %in% key(m), ]
  if (nrow(miss)) {
    warning(nrow(miss), " registered fits are missing; counted as failed fits")
    fill <- m[rep(1, nrow(miss)), ]; fill[] <- NA
    fill[c("seed", "p", "variance", "arm")] <- miss[c("seed", "p", "variance", "arm")]
    fill$recovered <- FALSE; fill$false_active <- 7; fill$converged <- FALSE; fill$error <- "missing checkpoint"
    m <- rbind(m, fill)
  }
  # w_S = p/50 equals w_S = 1 at p = 50: that fit is reused, labelled for the curves
  dup <- subset(m, p == 50 & arm == "w1"); dup$arm <- rep("p_over_50", nrow(dup))
  rbind(m, dup)
}

#' S-7 unsupervised reference metrics for the same datasets
#' @param design Design rows.
#' @return Data frame of seed, p, variance, recovered, false_active, converged.
#' @examples
#' # s9_reference(design)
s9_reference <- function(design) {
  s7m <- read.csv(file.path(s7, "metrics.csv"))
  s7m <- subset(s7m, arm == "unsupervised")[, c("seed", "p", "variance", "recovered", "false_active", "converged", "c_validation")]
  merge(design, s7m)
}

# Gate ----

#' Per-cell Gate S qualification at p = 2,000
#' @param m Fit rows; @param ref Reference rows (recovered per dataset).
#' @return Data frame per arm x variance with counts and qualification.
#' @examples
#' # s9_gate_cells(m, ref)
s9_gate_cells <- function(m, ref) {
  m2 <- subset(m, p == 2000 & arm != "w0")
  out <- do.call(rbind, lapply(split(m2, list(m2$arm, m2$variance), drop = TRUE), function(x) {
    r <- subset(ref, p == 2000 & variance == x$variance[1])
    data.frame(arm = x$arm[1], variance = x$variance[1], n = nrow(x),
               recovered = sum(x$recovered), reference_recovered = sum(r$recovered),
               gain = sum(x$recovered) - sum(r$recovered),
               mean_false_active = mean(x$false_active), converged = sum(x$converged))
  }))
  out$qualifies <- out$gain >= 3 & out$mean_false_active <= 1 & out$converged >= 0.9 * out$n
  out[order(out$arm, out$variance), ]
}

# Main ----
m <- s9_read(dir)
design <- read.csv(file.path(dir, "design.csv"))
ref7 <- s9_reference(design)
ref0 <- subset(m, arm == "w0")[, c("seed", "p", "variance", "recovered")]

cells7 <- s9_gate_cells(m, ref7); cells0 <- s9_gate_cells(m, ref0)
verdict <- function(cells) {
  a <- aggregate(qualifies ~ arm, cells, sum)
  data.frame(arm = a$arm, qualifying_cells = a$qualifies, passes = a$qualifies >= 2)
}
g7 <- verdict(cells7); g0 <- verdict(cells0)
write.csv(m, file.path(dir, "metrics.csv"), row.names = FALSE)
write.csv(cells7, file.path(dir, "gate_cells.csv"), row.names = FALSE)
write.csv(cells0, file.path(dir, "gate_cells_vs_w0.csv"), row.names = FALSE)
write.csv(g7, file.path(dir, "gate.csv"), row.names = FALSE)
write.csv(g0, file.path(dir, "gate_vs_w0.csv"), row.names = FALSE)

## Per-arm summary (curves) ----
num <- function(x) as.numeric(x)
summ <- do.call(rbind, lapply(split(m, list(m$arm, m$p, m$variance), drop = TRUE), function(x)
  data.frame(arm = x$arm[1], p = x$p[1], variance = x$variance[1], n = nrow(x),
             w_S_median = median(num(x$w_S), na.rm = TRUE),
             recovered = sum(x$recovered), mean_false_active = mean(x$false_active),
             converged = sum(x$converged), errors = sum(x$error != "", na.rm = TRUE),
             cor_planted_median = median(x$cor_planted, na.rm = TRUE),
             share_match_median = median(x$share_match, na.rm = TRUE),
             share_live_median = median(x$share_live, na.rm = TRUE),
             cal_slope_median = median(x$cal_slope, na.rm = TRUE),
             beta_covered = sum(x$beta_covered, na.rm = TRUE),
             beta_ratio_median = median(x$beta_match / x$beta_star, na.rm = TRUE),
             c_validation_mean = mean(x$c_validation, na.rm = TRUE),
             cv_c_mean = mean(x$cv_c, na.rm = TRUE),
             sweeps_median = median(x$sweeps, na.rm = TRUE),
             max_beta_max = max(x$max_beta, na.rm = TRUE),
             recon_rel_max = max(x$recon_rel, na.rm = TRUE),
             min_tau_min = min(x$min_tau, na.rm = TRUE),
             lost_variance = sum(x$lost_variance, na.rm = TRUE),
             share_centered = median(x$share_centered), share_uncentered = median(x$share_uncentered))))
summ <- summ[order(summ$p, summ$variance, summ$w_S_median), ]
write.csv(summ, file.path(dir, "summary.csv"), row.names = FALSE)

## CV-chosen arm per dataset ----
sup <- subset(m, arm != "w0" & !(p == 50 & arm == "p_over_50"))
cvpick <- do.call(rbind, lapply(split(sup, list(sup$seed, sup$p, sup$variance), drop = TRUE), function(x) {
  x <- x[is.finite(x$cv_c), ]
  if (!nrow(x)) return(NULL)
  x[which.max(x$cv_c), c("seed", "p", "variance", "arm", "w_S", "cv_c", "recovered", "false_active", "converged")]
}))
write.csv(cvpick, file.path(dir, "cv_chosen.csv"), row.names = FALSE)
cvcell <- do.call(rbind, lapply(split(cvpick, list(cvpick$p, cvpick$variance), drop = TRUE), function(x)
  data.frame(p = x$p[1], variance = x$variance[1], n = nrow(x),
             arms_chosen = paste(names(table(x$arm)), table(x$arm), sep = ":", collapse = " "),
             recovered = sum(x$recovered), mean_false_active = mean(x$false_active),
             converged = sum(x$converged))))
write.csv(cvcell, file.path(dir, "cv_chosen_cells.csv"), row.names = FALSE)

## Curves ----
png(file.path(dir, "weight_curves.png"), width = 1500, height = 1000, res = 120)
op <- par(mfrow = c(2, 4), mar = c(4, 4, 2.5, 1))
metrics <- c(recovered = "Recovered (of 10)", mean_false_active = "Mean false-active",
             share_match_median = "Survival share, matched", cal_slope_median = "Calibration slope",
             beta_covered = "beta* covered (of 10)", converged = "Converged (of 10)",
             c_validation_mean = "Validation C", cv_c_mean = "Training CV C")
fx <- subset(summ, arm %in% paste0("w", c(1, 10, 100, 1000, 10000)))
cols <- c("0.0025" = "#1b9e77", "0.01" = "#d95f02", "0.05" = "#7570b3")
for (mt in names(metrics)) {
  plot(NA, xlim = c(1, 1e4), ylim = range(fx[[mt]], na.rm = TRUE), log = "x",
       xlab = "w_S", ylab = metrics[[mt]], main = metrics[[mt]])
  for (pp in c(50, 2000)) for (v in names(cols)) {
    z <- fx[fx$p == pp & fx$variance == as.numeric(v), ]
    lines(z$w_S_median, z[[mt]], col = cols[[v]], lty = if (pp == 50) 2 else 1, type = "b", pch = 19)
  }
}
legend("bottomleft", c("p=2000", "p=50", paste("var", names(cols))), lty = c(1, 2, 1, 1, 1),
       col = c(1, 1, cols), bty = "n", cex = .8)
par(op); dev.off()

cat("Gate S vs S-7 reference:\n"); print(g7, row.names = FALSE)
cat("\nGate S vs canonical w_S = 0:\n"); print(g0, row.names = FALSE)
cat("\nCV-chosen arms:\n"); print(cvcell, row.names = FALSE)
