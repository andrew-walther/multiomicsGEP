# ============================================================
# Script: summarize_sim_panel.R
# Purpose: Summarize the w = 1 simulation panel against O1-O4 and O7
#          (definitions in results/ssbmf/sim_design.md).
# Author: Andrew Walther
# Created: 2026-10-26
# Dependencies: base R
# Usage: Rscript results/ssbmf/summarize_sim_panel.R   (from the repository root)
# ============================================================

out_dir <- "results/ssbmf/outputs/sim_panel"
rows <- read.csv(file.path(out_dir, "sim_panel_rows.csv"))
rows$fit_id <- paste(rows$cell, rows$seed, rows$K_init, rows$arm)
fits <- rows[!duplicated(rows$fit_id), ]
pct <- function(x) sprintf("%.0f%% (%d/%d)", 100 * mean(x), sum(x), length(x))

# Per cell and arm ----
by <- split(rows, list(rows$cell, rows$arm), drop = TRUE)
summ <- do.call(rbind, lapply(by, function(x) {
  f <- x[!duplicated(x$fit_id), ]
  hv <- x[x$pve_true >= 0.02, ]
  lv <- x[x$beta_true != 0 & x$pve_true < 0.02, ]
  cov_rows <- x[(x$pve_true >= 0.02 & x$best_r >= 0.9) | (x$pve_true < 0.02 & x$best_r >= 0.7), ]
  data.frame(cell = x$cell[1], arm = x$arm[1], fits = nrow(f),
             O1_no_decrease = mean(f$step_decreases == 0),
             worst_rel_step = min(f$worst_relative_step),
             O2_converged = mean(f$converged),
             median_sweeps = median(f$sweeps), max_sweeps = max(f$sweeps),
             O3_highvar_recovered = mean(hv$best_r >= 0.9),
             O4_lowvar_recovered = if (nrow(lv)) mean(lv$best_r >= 0.7 & lv$active & lv$sign_ok) else NA,
             lowvar_matched = if (nrow(lv)) mean(lv$best_r >= 0.7) else NA,
             O7_coverage90 = if (nrow(cov_rows)) mean(cov_rows$covered) else NA,
             O7_n = nrow(cov_rows),
             mean_val_C = mean(f$val_C),
             rejections = sum(f$rejections),
             median_runtime_s = median(f$runtime_s),
             surv_share_prognostic = median(x$surv_share_matched[x$beta_true != 0], na.rm = TRUE),
             surv_share_null = median(x$surv_share_matched[x$beta_true == 0], na.rm = TRUE))
}))
summ <- summ[order(summ$cell, summ$arm), ]
write.csv(summ, file.path(out_dir, "sim_panel_summary.csv"), row.names = FALSE)

# Low-variance recovery by K_init ----
lv <- rows[rows$beta_true != 0 & rows$pve_true < 0.02, ]
lvk <- aggregate(cbind(recovered = best_r >= 0.7 & active & sign_ok, matched = best_r >= 0.7, best_r = best_r) ~ cell + arm + K_init,
                 data = lv, FUN = mean)
write.csv(lvk, file.path(out_dir, "sim_panel_lowvar_by_K.csv"), row.names = FALSE)

# Overall verdicts ----
ss <- rows[grepl("^ssbmf", rows$arm), ]
ssf <- ss[!duplicated(ss$fit_id), ]
hv <- ss[ss$pve_true >= 0.02, ]
cov_rows <- ss[(ss$pve_true >= 0.02 & ss$best_r >= 0.9) | (ss$pve_true < 0.02 & ss$best_r >= 0.7), ]
verdict <- c(
  O1 = sprintf("ssbmf fits with no step lowering the objective by > 1e-8 |obj|: %s; worst relative step %.1e",
               pct(ssf$step_decreases == 0), min(ssf$worst_relative_step)),
  O2 = sprintf("ssbmf fits converged within 1,000 sweeps: %s", pct(ssf$converged)),
  O3 = sprintf("high-variance programs recovered at |r| >= 0.9: %s", pct(hv$best_r >= 0.9)),
  O7 = sprintf("90%% intervals covering the true coefficient: %s", pct(cov_rows$covered))
)
for (a in unique(lv$arm)) {
  for (cl in unique(lv$cell)) {
    x <- lv[lv$arm == a & lv$cell == cl, ]
    verdict[paste("O4", cl, a)] <- sprintf("cell %s, %s: low-variance program recovered (|r| >= 0.7, survival-active, correct sign): %s",
                                          cl, a, pct(x$best_r >= 0.7 & x$active & x$sign_ok))
  }
}
writeLines(c("# Simulation panel verdicts (w_S = 1)", "", paste0("- ", verdict)), file.path(out_dir, "sim_panel_verdicts.md"))
cat(paste0(verdict, collapse = "\n"), "\n\n")
print(format(summ, digits = 3), row.names = FALSE)
cat("\nLow-variance recovery by K_init:\n"); print(format(lvk, digits = 3), row.names = FALSE)
