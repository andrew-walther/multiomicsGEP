# ============================================================
# Script:  docs/progress_book/figs/make_k_sweep_6panel_fig.R
# Purpose: Build the K_init=2..20 consensus-sweep figure for the 9/4
#          progress-book chapter: six panels (ELBO, in-sample BIC, in-sample
#          joint log-likelihood, held-out survival log-likelihood, bi-CV
#          genomics log-likelihood, external C-index) over K_init, each
#          marking the criterion-preferred K_init and the recommended
#          K_init=7, with K_eff_total annotated on the external-C panel.
#
#          Extends the 8/27 4-panel figure
#          (presentation/walther_lab_meeting_08_27_2026/figs/make_k_sweep_fig.R)
#          with the two genuinely held-out criteria added in Stage 1 of the
#          9/4 plan (code/compute_cv_loglik.R). The point of showing all six
#          curves together, not just the recommendation, is unchanged: make
#          disagreement between criteria visible rather than picking a
#          winner and hiding the rest.
#
#   Inputs:
#     results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv
#   Output:
#     docs/progress_book/figs/2026-09-04_k_init_sweep_6panel.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Dependencies: ggplot2, dplyr, patchwork
# Usage:   Rscript docs/progress_book/figs/make_k_sweep_6panel_fig.R
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(patchwork)
})

if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../code/fit_modular.R")) {
  ROOT <- "../.."
} else stop("Cannot locate project root (code/fit_modular.R).")

KSWEEP_CSV <- file.path(ROOT, "results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv")
OUT_PNG    <- file.path(ROOT, "docs/progress_book/figs/2026-09-04_k_init_sweep_6panel.png")

if (!file.exists(KSWEEP_CSV)) stop("Missing K_init sweep CSV: ", KSWEEP_CSV)

d <- read.csv(KSWEEP_CSV, stringsAsFactors = FALSE) |> filter(fit_status == "ok")

K_REC <- 7L
elbo_pref     <- d$K_init[which.max(d$elbo_full)]
bic_pref      <- d$K_init[which.min(d$bic)]
c_pref        <- d$K_init[which.max(d$mean_external_c)]
cv_surv_pref  <- if (any(!is.na(d$cv_survival_logPL_per_event)))
                   d$K_init[which.max(d$cv_survival_logPL_per_event)] else NA_integer_
bicv_pref     <- if (any(!is.na(d$bicv_genomics_total_loglik)))
                   d$K_init[which.max(d$bicv_genomics_total_loglik)] else NA_integer_

# UNC palette: Carolina blue for the actual curve/points, navy for the
# recommended-K_init marker, a muted red for each panel's own criterion-best.
col_curve <- "#4B9CD3"
col_rec   <- "#13294B"
col_best  <- "#B31B1B"

mark_layer <- function(pref_K) {
  list(
    geom_vline(xintercept = K_REC, linetype = "dotted", color = col_rec, linewidth = 0.5),
    if (!is.na(pref_K) && pref_K != K_REC)
      geom_vline(xintercept = pref_K, linetype = "dashed", color = col_best, linewidth = 0.5)
  )
}

base_theme <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        panel.grid.minor = element_blank())

p_elbo <- ggplot(d, aes(K_init, elbo_full)) +
  mark_layer(elbo_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 1.8) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "ELBO (higher = better)", x = "K_init", y = "elbo_full") +
  base_theme

p_bic <- ggplot(d, aes(K_init, bic)) +
  mark_layer(bic_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 1.8) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "In-sample BIC (lower = better)", x = "K_init", y = "BIC") +
  base_theme

p_ll <- ggplot(d, aes(K_init, loglik_joint)) +
  mark_layer(elbo_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 1.8) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "In-sample joint log-likelihood (higher = better)", x = "K_init", y = "loglik_joint") +
  base_theme

p_cv_surv <- ggplot(d, aes(K_init, cv_survival_logPL_per_event)) +
  mark_layer(cv_surv_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 1.8) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "Held-out survival log-lik./event (higher = better)",
       x = "K_init", y = "CV logPL / event") +
  base_theme

p_bicv <- ggplot(d, aes(K_init, bicv_genomics_total_loglik)) +
  mark_layer(bicv_pref) +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 1.8) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "Bi-CV genomics log-likelihood (higher = better)",
       x = "K_init", y = "bi-CV log-lik.") +
  base_theme

p_c <- ggplot(d, aes(K_init, mean_external_c)) +
  mark_layer(c_pref) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray55") +
  geom_line(color = col_curve, linewidth = 0.9) +
  geom_point(color = col_curve, size = 1.8) +
  geom_text(aes(label = K_eff_total), vjust = -1.1, size = 2.8, color = "gray30") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "External C-index (higher = better)", x = "K_init", y = "mean C-index",
       caption = "Point labels = K_eff_total (ARD-kept factors). Dotted line = recommended K_init=7; dashed = this panel's own criterion-best.") +
  base_theme

fig <- (p_elbo + p_bic + p_ll) / (p_cv_surv + p_bicv + p_c) +
  plot_annotation(
    title = sprintf("K_init consensus sweep (K_init = 2..20): ELBO/BIC prefer K_init=%d, external C prefers K_init=%d",
                     elbo_pref, c_pref),
    theme = theme(plot.title = element_text(size = 12.5, face = "bold"))
  )

ggsave(OUT_PNG, fig, width = 16, height = 8, dpi = 170, bg = "white")

cat(sprintf("Wrote %s\n", normalizePath(OUT_PNG)))
cat(sprintf("ELBO-pref K_init=%d | BIC-pref K_init=%d | in-sample-LL-pref K_init=%d | held-out-surv-LL-pref K_init=%s | bi-CV-genomics-LL-pref K_init=%s | C-pref K_init=%d | recommended K_init=%d\n",
            elbo_pref, bic_pref, elbo_pref,
            ifelse(is.na(cv_surv_pref), "NA", cv_surv_pref),
            ifelse(is.na(bicv_pref), "NA", bicv_pref),
            c_pref, K_REC))
