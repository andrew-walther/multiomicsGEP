# ============================================================
# Script:  presentation/walther_lab_meeting_08_27_2026/figs/make_k_eff_step_fig.R
# Purpose: Build a compact K_eff_total-vs-K_init figure for the appendix
#          K_init-sweep slide, showing the "step up" pattern: K_eff_total
#          (ARD-kept factor count) climbs from 2 at K_init=2 to a plateau
#          around 5 by K_init~9, rather than staying fixed as K_init grows.
#          Stacked bars break K_eff_total into survival-active vs.
#          genomics-only so it's visible that the step-up is driven by
#          additional genomics-only factors, not more survival-active ones.
#
#   Inputs:
#     results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv
#   Output:
#     presentation/walther_lab_meeting_08_27_2026/assets/k_eff_step.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-27
# Dependencies: ggplot2, dplyr, tidyr
# Usage:   Rscript presentation/walther_lab_meeting_08_27_2026/figs/make_k_eff_step_fig.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

KSWEEP_CSV <- file.path(ROOT, "results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv")
OUT_PNG    <- file.path(ROOT, "presentation/walther_lab_meeting_08_27_2026/assets/k_eff_step.png")

if (!file.exists(KSWEEP_CSV)) stop("Missing K_init sweep CSV: ", KSWEEP_CSV)

d <- read.csv(KSWEEP_CSV, stringsAsFactors = FALSE) |> filter(fit_status == "ok")

K_REC <- 7L

d_long <- d |>
  select(K_init, K_survival_active, K_genomics_only) |>
  pivot_longer(c(K_survival_active, K_genomics_only), names_to = "type", values_to = "count") |>
  mutate(type = recode(type, K_survival_active = "Survival-active",
                       K_genomics_only  = "Genomics-only"),
         type = factor(type, levels = c("Genomics-only", "Survival-active")))

col_surv <- "#B31B1B"   # UNC red -- survival-active
col_geno <- "#4B9CD3"   # Carolina blue -- genomics-only
col_rec  <- "#13294B"   # navy -- recommended K_init marker

p <- ggplot(d_long, aes(K_init, count, fill = type)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = K_REC, linetype = "dotted", color = col_rec, linewidth = 0.6) +
  geom_text(data = d, aes(K_init, K_eff_total, label = K_eff_total),
            inherit.aes = FALSE, vjust = -0.6, size = 3.4, color = "gray25") +
  scale_fill_manual(values = c("Genomics-only" = col_geno, "Survival-active" = col_surv), name = NULL) +
  scale_x_continuous(breaks = seq(2, 20, 1)) +
  scale_y_continuous(breaks = seq(0, 6, 1), limits = c(0, 6.2)) +
  labs(x = expression(K[init]), y = expression(K[eff]~"(ARD-kept factors)"),
       title = expression(K[eff]~"climbs with"~K[init]~"before plateauing, rather than staying fixed"),
       caption = "Dotted line = recommended K_init=7. Bar labels = K_eff_total.") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 15),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

ggsave(OUT_PNG, p, width = 10, height = 5, dpi = 170, bg = "white")

cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
