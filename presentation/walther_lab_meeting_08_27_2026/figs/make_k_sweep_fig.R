# ============================================================
# Script:  presentation/walther_lab_meeting_08_27_2026/figs/make_k_sweep_fig.R
# Purpose: Build the K_init=2..20 consensus-sweep figure for the 8/27 lab
#          meeting deck: four panels (ELBO, BIC, ELBO-style joint
#          log-likelihood, external C-index) over K_init, each marking the
#          criterion-preferred K_init and the recommended K_init=7, with the
#          resulting K_eff_total (ARD-kept factor count) annotated per point.
#
#          This is the figure Stage 1 of the two-stage K-selection framework
#          (methods slide) refers to -- the point of showing all four curves
#          together, not just the recommendation, is to make the ELBO/BIC vs.
#          external-C disagreement visible rather than picking a winner and
#          hiding the rest.
#
#   Inputs:
#     results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv
#   Output:
#     presentation/walther_lab_meeting_08_27_2026/assets/k_init_sweep_4panel.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-08-27
# Dependencies: ggplot2, dplyr, patchwork
# Usage:   Rscript presentation/walther_lab_meeting_08_27_2026/figs/make_k_sweep_fig.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

KSWEEP_CSV <- file.path(ROOT, "results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv")
OUT_PNG    <- file.path(ROOT, "presentation/walther_lab_meeting_08_27_2026/assets/k_init_sweep_4panel.png")

if (!file.exists(KSWEEP_CSV)) stop("Missing K_init sweep CSV: ", KSWEEP_CSV)

d <- read.csv(KSWEEP_CSV, stringsAsFactors = FALSE) |> filter(fit_status == "ok")

K_REC <- 7L
elbo_pref <- d$K_init[which.max(d$elbo_full)]
bic_pref  <- d$K_init[which.min(d$bic)]
c_pref    <- d$K_init[which.max(d$mean_external_c)]

# UNC palette: Carolina blue for the actual curve/points, navy for the
# recommended-K_init marker, a muted red for each panel's own criterion-best.
col_curve <- "#4B9CD3"
col_rec   <- "#13294B"
col_best  <- "#B31B1B"

mark_layer <- function(pref_K) {
  list(
    geom_vline(xintercept = K_REC, linetype = "dotted", color = col_rec, linewidth = 0.5),
    if (pref_K != K_REC)
      geom_vline(xintercept = pref_K, linetype = "dashed", color = col_best, linewidth = 0.5)
  )
}

base_theme <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank())

p_elbo <- ggplot(d, aes(K_init, elbo_full)) +
  mark_layer(elbo_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 2) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "ELBO", x = "K_init", y = "elbo_full") +
  base_theme

p_bic <- ggplot(d, aes(K_init, bic)) +
  mark_layer(bic_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 2) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "BIC (lower = better)", x = "K_init", y = "BIC") +
  base_theme

p_ll <- ggplot(d, aes(K_init, loglik_joint)) +
  mark_layer(elbo_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 2) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "Log-likelihood (genomics + survival)", x = "K_init", y = "loglik_joint") +
  base_theme

p_c <- ggplot(d, aes(K_init, mean_external_c)) +
  mark_layer(c_pref) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray55") +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 2) +
  geom_text(aes(label = K_eff_total), vjust = -1.1, size = 3.1, color = "gray30") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "External C-index (5 held-out cohorts)", x = "K_init", y = "mean C-index",
       caption = "Point labels = K_eff_total (ARD-kept factors). Dotted line = recommended K_init=7; dashed = this panel's own criterion-best.") +
  base_theme

fig <- (p_elbo + p_bic) / (p_ll + p_c) +
  plot_annotation(
    title = sprintf("K_init consensus sweep (K_init = 2..20): ELBO/BIC prefer K_init=%d, external C prefers K_init=%d",
                     elbo_pref, c_pref),
    theme = theme(plot.title = element_text(size = 13, face = "bold"))
  )

ggsave(OUT_PNG, fig, width = 11.5, height = 8, dpi = 170, bg = "white")

cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
cat(sprintf("ELBO-preferred K_init=%d | BIC-preferred K_init=%d | C-preferred K_init=%d | recommended K_init=%d\n",
            elbo_pref, bic_pref, c_pref, K_REC))
