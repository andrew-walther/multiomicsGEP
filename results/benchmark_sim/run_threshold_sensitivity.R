# ============================================================
# Script:  results/benchmark_sim/run_threshold_sensitivity.R
# Purpose: Stage 4 of the 9/4 plan. classify_factors()'s beta_thresh
#          (0.001) and pve_thresh (0.01) are not derived from a citable
#          source (ROADMAP.md "Literature grounding..." item, raised
#          2026-08-27). This converts "the thresholds are uncited" into
#          "here is the range over which the K_eff answer does not
#          change" -- no re-fitting, just re-running classify_factors()
#          over a grid on the 19 cached K_init sweep fits.
#
#          Also: (1) sweeps rel_thresh (added August, simulation-
#          calibrated at 0.65) and confirms it must not be applied to the
#          real fit, where Program 3's ratio-to-max is ~0.28 yet it is
#          externally validated (DECISIONS.md 2026-07-15); (2) checks
#          beta comparability across factors -- classify_factors() tests
#          raw |E_q[beta_k]| against a fixed cutoff, which assumes every
#          factor's projection score is on a comparable scale. F_k is
#          unit-L2-normalized (fit$EF_norms), removing the dominant scale
#          gap, but Var(ZF_k) is never standardized anywhere in the
#          codebase -- this computes it directly and reports a
#          variance-standardized |beta_k|*sd(ZF_k) ranking alongside the
#          raw one.
#
#   Inputs:
#     results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_fits.rds
#     (19 cached fits, K_init=2..20 -- no re-fitting)
#   Outputs:
#     results/benchmark_sim/outputs/threshold_sensitivity/
#       threshold_grid_K7.csv        -- K_survival_active over the full
#                                        beta_thresh x pve_thresh grid at K_init=7
#       threshold_vs_kinit.csv       -- K_survival_active vs K_init at 3 beta_thresh values
#       rel_thresh_sweep.csv         -- rel_thresh sensitivity at K_init=7
#       beta_comparability_K7.csv    -- raw vs variance-standardized beta ranking at K_init=7
#     docs/progress_book/figs/2026-09-04_threshold_sensitivity_heatmap.png
#     docs/progress_book/figs/2026-09-04_threshold_vs_kinit.png
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Usage:   Rscript results/benchmark_sim/run_threshold_sensitivity.R
#          (fast -- no fitting, just re-classifying 19 cached fits over a grid)
# ============================================================

if (!file.exists("code/fit_modular.R") && file.exists("../../code/fit_modular.R"))
  setwd("../../")

suppressPackageStartupMessages({ library(yaml); library(survival); library(ggplot2); library(dplyr) })

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/update_beta.R"); source("code/update_beta_cohort.R")
source("code/update_L.R"); source("code/update_F.R"); source("code/update_tau.R")
source("code/compute_elbo.R"); source("code/update_F_cohort.R")
source("code/predict.R"); source("code/predict_cox_on_yf.R")
tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/preprocess_desurv.R")
source("code/select_K.R")

OUT_DIR <- "results/benchmark_sim/outputs/threshold_sensitivity"
FIG_DIR <- "docs/progress_book/figs"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

FITS_RDS <- "results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_fits.rds"
if (!file.exists(FITS_RDS)) stop("Missing cached fits: ", FITS_RDS, " -- run run_k_init_sweep.R first.")
fits <- readRDS(FITS_RDS)
cat(sprintf("Loaded %d cached fits (K_init = %s)\n\n", length(fits),
            paste(sort(as.integer(names(fits))), collapse = ", ")))

# --------------------------------------------------------------------------
# 1. Rebuild the SAME D4-preprocessed Y_train used to produce these fits
#    (deterministic given the code -- needed for compute_pve()/ZF, not for
#    re-fitting).
# --------------------------------------------------------------------------

cat("--- Rebuilding D4 training matrix (for PVE/ZF only, no re-fitting) ---\n")
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
pp <- preprocess_merged_cohorts(
  cohort_raw_list = train_raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank"
)
Y_train <- pp$Y
cat(sprintf("  n=%d, p=%d\n\n", nrow(Y_train), ncol(Y_train)))

BETA_THRESHOLDS <- c(0, 1e-4, 5e-4, 0.001, 0.005, 0.01, 0.05)
PVE_THRESHOLDS  <- c(0.001, 0.005, 0.01, 0.02, 0.05)
K_INIT_VALUES   <- sort(as.integer(names(fits)))

# --------------------------------------------------------------------------
# 2. Full grid at K_init=7 -- heatmap of K_survival_active.
# --------------------------------------------------------------------------

cat("=== Threshold grid at K_init=7 ===\n")
fit7 <- fits[["7"]]
grid_rows <- list()
for (bt in BETA_THRESHOLDS) {
  for (pt in PVE_THRESHOLDS) {
    cls <- classify_factors(fit7, Y_train, beta_thresh = bt, pve_thresh = pt)
    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      beta_thresh = bt, pve_thresh = pt,
      K_survival_active = sum(cls$category == "survival_active"),
      K_genomics_only   = sum(cls$category == "genomics_only"),
      K_dead            = sum(cls$category == "dead"),
      stringsAsFactors = FALSE
    )
  }
}
grid7 <- do.call(rbind, grid_rows)
write.csv(grid7, file.path(OUT_DIR, "threshold_grid_K7.csv"), row.names = FALSE)
print(grid7)

p_heat <- ggplot(grid7, aes(x = factor(beta_thresh), y = factor(pve_thresh), fill = K_survival_active)) +
  geom_tile(color = "white") +
  geom_text(aes(label = K_survival_active), size = 4) +
  scale_fill_gradient(low = "#DCE6F1", high = "#13294B") +
  labs(title = "K_survival_active at K_init=7 across the threshold grid",
       x = "beta_thresh", y = "pve_thresh", fill = "K_survival_active",
       caption = "config/globals.yml defaults: beta_thresh=0.001, pve_thresh=0.01") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "2026-09-04_threshold_sensitivity_heatmap.png"),
       p_heat, width = 7.5, height = 5.5, dpi = 170, bg = "white")

# --------------------------------------------------------------------------
# 3. K_survival_active vs K_init at three beta_thresh values (pve_thresh
#    fixed at its default 0.01).
# --------------------------------------------------------------------------

cat("\n=== K_survival_active vs K_init at 3 beta_thresh values ===\n")
BT_TRIPLE <- c(0.001, 0.01, 0.05)
vs_rows <- list()
for (K_init in K_INIT_VALUES) {
  fit_k <- fits[[as.character(K_init)]]
  for (bt in BT_TRIPLE) {
    cls <- classify_factors(fit_k, Y_train, beta_thresh = bt, pve_thresh = 0.01)
    vs_rows[[length(vs_rows) + 1]] <- data.frame(
      K_init = K_init, beta_thresh = bt,
      K_survival_active = sum(cls$category == "survival_active"),
      stringsAsFactors = FALSE
    )
  }
}
vs_kinit <- do.call(rbind, vs_rows)
write.csv(vs_kinit, file.path(OUT_DIR, "threshold_vs_kinit.csv"), row.names = FALSE)

p_vs <- ggplot(vs_kinit, aes(K_init, K_survival_active, color = factor(beta_thresh))) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(title = "K_survival_active vs K_init at three beta_thresh values (pve_thresh=0.01)",
       x = "K_init", y = "K_survival_active", color = "beta_thresh") +
  theme_minimal(base_size = 12)
ggsave(file.path(FIG_DIR, "2026-09-04_threshold_vs_kinit.png"),
       p_vs, width = 8, height = 5, dpi = 170, bg = "white")

# --------------------------------------------------------------------------
# 4. rel_thresh sweep at K_init=7 -- confirm it must NOT be applied to the
#    real fit (Program 3's ratio-to-max is ~0.28 yet externally validated,
#    DECISIONS.md 2026-07-15).
# --------------------------------------------------------------------------

cat("\n=== rel_thresh sweep at K_init=7 (beta_thresh=0.001, pve_thresh=0.01) ===\n")
REL_THRESHOLDS <- c(NA, seq(0.1, 0.9, by = 0.1))
rel_rows <- list()
ab_beta7 <- abs(fit7$EBeta)
for (rt in REL_THRESHOLDS) {
  cls <- classify_factors(fit7, Y_train, beta_thresh = 0.001, pve_thresh = 0.01,
                           rel_thresh = if (is.na(rt)) NULL else rt)
  rel_rows[[length(rel_rows) + 1]] <- data.frame(
    rel_thresh = rt, K_survival_active = sum(cls$category == "survival_active"),
    stringsAsFactors = FALSE
  )
}
rel_sweep <- do.call(rbind, rel_rows)
write.csv(rel_sweep, file.path(OUT_DIR, "rel_thresh_sweep.csv"), row.names = FALSE)
print(rel_sweep)
cat(sprintf("\nRatio of the smaller survival-active |beta| to max(|beta|) at K_init=7: %.3f\n",
            min(ab_beta7[ab_beta7 > 0.001]) / max(ab_beta7)))
cat("Any rel_thresh >= this ratio would drop the smaller externally-validated factor -- confirms\n")
cat("rel_thresh must not be applied here (DECISIONS.md 2026-07-15 / 2026-08-19).\n")

# --------------------------------------------------------------------------
# 5. Beta comparability: raw |beta_k| vs variance-standardized |beta_k|*sd(ZF_k).
# --------------------------------------------------------------------------

cat("\n=== Beta comparability at K_init=7: raw vs variance-standardized ranking ===\n")
ZF7 <- Y_train %*% sweep(fit7$EF, 2, fit7$EF_norms, "/")
sd_ZF7 <- apply(ZF7, 2, sd)
comparability <- data.frame(
  factor        = seq_len(ncol(fit7$EF)),
  abs_beta_raw  = abs(fit7$EBeta),
  sd_ZF         = sd_ZF7,
  abs_beta_std  = abs(fit7$EBeta) * sd_ZF7,
  stringsAsFactors = FALSE
)
comparability$rank_raw <- rank(-comparability$abs_beta_raw, ties.method = "min")
comparability$rank_std <- rank(-comparability$abs_beta_std, ties.method = "min")
write.csv(comparability, file.path(OUT_DIR, "beta_comparability_K7.csv"), row.names = FALSE)
print(comparability)
cat(sprintf("\nsd(ZF_k) range across factors: [%.4f, %.4f] (ratio %.2fx)\n",
            min(sd_ZF7), max(sd_ZF7), max(sd_ZF7) / min(sd_ZF7)))
cat(sprintf("Raw and variance-standardized rankings %s.\n",
            if (identical(comparability$rank_raw, comparability$rank_std))
              "AGREE exactly" else "DISAGREE -- see beta_comparability_K7.csv"))

cat("\n=== Stage 4 outputs written to ", OUT_DIR, " and ", FIG_DIR, " ===\n")
