# ============================================================
# Script:  docs/progress_book/figs/make_percohort_beta_fig.R
# Purpose: Bar chart of the beta_cohort_id fit's per-cohort survival
#          coefficients across all 7 factors, for the 9/4 progress-book
#          chapter (§3, "Interpretability: seeing how cohorts differ").
#
#          Recreates a figure that was originally produced ad hoc (no
#          committed script was found for it) and had gone stale relative
#          to fit_cox_on_yf()'s Phase C orientation fix (review finding,
#          Step 2, 2026-09-04, DECISIONS.md) -- the joint coefficients'
#          sign convention changed, so the original PNG showed Program 3
#          positive / Program 7 negative, the opposite of the corrected
#          fit. Magnitudes and which cohort dominates which program are
#          unaffected; only the sign is.
#
#   Inputs:
#     results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_fits.rds
#   Output:
#     docs/progress_book/figs/2026-09-04_percohort_beta.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Dependencies: ggplot2
# Usage:   Rscript docs/progress_book/figs/make_percohort_beta_fig.R
# ============================================================

suppressPackageStartupMessages(library(ggplot2))

if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../code/fit_modular.R")) {
  ROOT <- "../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

FITS_RDS <- file.path(ROOT, "results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_fits.rds")
OUT_PNG  <- file.path(ROOT, "docs/progress_book/figs/2026-09-04_percohort_beta.png")

if (!file.exists(FITS_RDS)) stop("Missing cached fits: ", FITS_RDS, " -- run run_cohort_beta_comparison.R first.")

fits <- readRDS(FITS_RDS)
fit  <- fits$joint_yfb_beta_c
EBeta <- fit$EBeta  # K x C matrix, colnames = cohort levels
K <- nrow(EBeta)

d <- data.frame(
  program = factor(rep(paste("Program", seq_len(K)), ncol(EBeta)), levels = paste("Program", seq_len(K))),
  cohort  = rep(colnames(EBeta), each = K),
  beta    = as.vector(EBeta)
)

# UNC palette, matching make_k_sweep_6panel_fig.R: Carolina blue / a muted
# red. Cohort order (CPTAC first) matches the original figure's legend.
cohort_levels <- colnames(EBeta)
col_map <- setNames(c("#4B9CD3", "#B31B1B")[seq_along(cohort_levels)], cohort_levels)

fig <- ggplot(d, aes(program, beta, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = col_map, name = "Cohort") +
  labs(title = "Per-cohort survival coefficients, all 7 factors (beta_cohort_id fit)",
       x = NULL, y = expression(beta[k]^{(c)}),
       caption = "Programs 3 and 7 are the two ARD-kept survival-active factors; the rest are near-zero in both cohorts.") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0, size = 10))

ggsave(OUT_PNG, fig, width = 10, height = 6, dpi = 170, bg = "white")
cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
