# ============================================================
# Script:  docs/progress_book/figs/make_3way_percohort_cindex_fig.R
# Purpose: Grouped bar chart of external C-index by cohort for the 9/4
#          progress-book chapter's 3-way performance comparison (§3):
#          joint_yfb (current model), joint_yfb_cohort_L (+ cohort
#          indicator), and the two-step EBMF->Cox baseline.
#
#          Recreates a figure that was originally produced ad hoc (no
#          committed script was found for it) and had gone stale relative
#          to the review-findings fixes (Steps 1-4, 2026-09-04,
#          DECISIONS.md): run_cohort_beta_comparison.R was re-run with
#          fresh fits, and joint_yfb_cohort_L's mean external C moved from
#          the original figure's ~0.605 to the corrected 0.5431.
#
#   Inputs:
#     results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_results.csv
#   Output:
#     docs/progress_book/figs/2026-09-04_3way_percohort_cindex.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Dependencies: ggplot2, tidyr
# Usage:   Rscript docs/progress_book/figs/make_3way_percohort_cindex_fig.R
# ============================================================

suppressPackageStartupMessages({ library(ggplot2); library(tidyr) })

if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../code/fit_modular.R")) {
  ROOT <- "../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

CBC_CSV <- file.path(ROOT, "results/benchmark_sim/outputs/cohort_beta_comparison/cohort_beta_comparison_results.csv")
OUT_PNG <- file.path(ROOT, "docs/progress_book/figs/2026-09-04_3way_percohort_cindex.png")

if (!file.exists(CBC_CSV)) stop("Missing: ", CBC_CSV, " -- run run_cohort_beta_comparison.R first.")

d <- read.csv(CBC_CSV, stringsAsFactors = FALSE)

ARMS <- c(joint_yfb = "Joint YFB (current model)",
          joint_yfb_cohort_L = "Joint YFB + cohort indicator",
          two_step_ebmf_cox = "Two-step EBMF->Cox")
COHORT_COLS <- c("c_Dijk", "c_Moffitt_GEO_array", "c_PACA_AU_array", "c_PACA_AU_seq", "c_Puleo_array")
COHORT_LABELS <- c(c_Dijk = "Dijk", c_Moffitt_GEO_array = "Moffitt_GEO_array",
                    c_PACA_AU_array = "PACA_AU_array", c_PACA_AU_seq = "PACA_AU_seq",
                    c_Puleo_array = "Puleo_array")

d3 <- d[d$arm %in% names(ARMS), c("arm", COHORT_COLS, "mean_external_c")]
long <- pivot_longer(d3, cols = c(all_of(COHORT_COLS), "mean_external_c"),
                      names_to = "cohort", values_to = "c_index")
long$cohort <- ifelse(long$cohort == "mean_external_c", "Mean", COHORT_LABELS[long$cohort])
long$cohort <- factor(long$cohort, levels = c(unname(COHORT_LABELS), "Mean"))
long$arm_label <- factor(ARMS[long$arm], levels = unname(ARMS))

# UNC palette, matching make_k_sweep_6panel_fig.R / make_percohort_beta_fig.R.
col_map <- setNames(c("#4B9CD3", "#13294B", "#B31B1B"), unname(ARMS))

fig <- ggplot(long, aes(cohort, c_index, fill = arm_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_fill_manual(values = col_map, name = NULL) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(title = "External C-index by cohort: 3-way performance comparison",
       x = NULL, y = "External C-index",
       caption = "Dashed line = chance (C=0.5). All three arms trained on TCGA+CPTAC (n=273).") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 20, hjust = 1),
        plot.caption = element_text(hjust = 0, size = 10))

ggsave(OUT_PNG, fig, width = 12, height = 7, dpi = 170, bg = "white")
cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
