# ============================================================
# Script:  presentation/walther_lab_meeting_06_18_2026/figs/make_factor_figs.R
# Purpose: Regenerate the real-data gene-program heatmaps and Kaplan-Meier
#          figures for the 6/18 deck with audience-facing titles — i.e. WITHOUT
#          the internal D4/D5 configuration labels baked into the report PNGs.
#
#          Reuses the exact fits, gene sets, and plotting logic from
#          docs/reports/desurv_factor_diagnostics_05_27_26.qmd, changing only the
#          figure/panel titles to describe each model by its characteristics
#          (projection predictor, DeSurv-aligned, ± cohort indicator).
#
#   Recommended model (report config "D4"): projection (YF)β, DeSurv-aligned, no cohort indicator
#   Cohort-indicator model (report config "D5"): same + cohort membership indicator
#
#   Inputs:
#     results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds
#     config/globals.yml, code/*, PDAC raw data (PDAC_DATA_ROOT)
#   Outputs (into the deck assets/ dir):
#     real_gep_heatmap.png, real_km.png   (recommended model)
#     cohort_heatmap.png,  cohort_km.png  (cohort-indicator contrast, appendix)
#
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-06-15
# Dependencies: yaml, survival, ggplot2, dplyr, tidyr, patchwork
# Usage:   PDAC_DATA_ROOT=<path> Rscript \
#            presentation/walther_lab_meeting_06_18_2026/figs/make_factor_figs.R
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(survival); library(ggplot2)
  library(dplyr); library(tidyr); library(patchwork)
})

# Resolve project root.
if (file.exists("code/fit_modular.R")) {
  ROOT <- "."
} else if (file.exists("../../../../code/fit_modular.R")) {
  ROOT <- "../../../.."
} else stop("Cannot locate project root (code/fit_modular.R).")
setwd(ROOT)

ASSETS <- "presentation/walther_lab_meeting_06_18_2026/assets"
dir.create(ASSETS, recursive = TRUE, showWarnings = FALSE)

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R")

b <- cfg$benchmark; p <- cfg$preprocessing
BETA_THRESH  <- cfg$k_selection$beta_threshold
TOP_N_DESURV <- p$top_n_genes_desurv

# --------------------------------------------------------------------------
# 1. Load training data + DeSurv-aligned preprocessing (matches the report)
# --------------------------------------------------------------------------
TRAIN_COHORTS <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN_COHORTS, TRAIN_COHORTS), function(ds)
  load_pdac_raw(ds, PDAC_DATA_ROOT))
time_train   <- c(train_raw$TCGA_PAAD$time,   train_raw$CPTAC$time)
status_train <- c(train_raw$TCGA_PAAD$status, train_raw$CPTAC$status)

invisible(capture.output(
  pp_desurv <- preprocess_merged_cohorts(
    cohort_raw_list          = train_raw,
    log_transform_flags      = PLATFORM_LOG_TRANSFORM[TRAIN_COHORTS],
    top_n                    = TOP_N_DESURV,
    rank_transform           = FALSE,
    per_platform_standardize = TRUE,
    normalize_method         = "none",
    selection_per_cohort     = TRUE,
    selection_method         = "combined_rank"
  ), type = "output"
))
genes_desurv <- pp_desurv$gene_names

fits <- readRDS("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds")

# Model metadata — audience-facing titles only (no D-labels).
MODELS <- list(
  recommended = list(fit_id = "D4", n_genes = 30,
                     heatmap_title = "Recommended model · projection (YF)β, DeSurv-aligned — top 30 genes",
                     heatmap_out   = file.path(ASSETS, "real_gep_heatmap.png"),
                     km_out        = file.path(ASSETS, "real_km.png")),
  cohort      = list(fit_id = "D5", n_genes = 40,
                     heatmap_title = "Cohort-indicator model · projection (YF)β — top 40 genes",
                     heatmap_out   = file.path(ASSETS, "cohort_heatmap.png"),
                     km_out        = file.path(ASSETS, "cohort_km.png"))
)

# --------------------------------------------------------------------------
# 2. Plotting helpers (copied verbatim from the diagnostics report, except the
#    titles, which are passed in / stripped of the internal config id)
# --------------------------------------------------------------------------
logrank_p <- function(EL_k, time, status) {
  grp <- ifelse(EL_k > median(EL_k), "High", "Low")
  if (length(unique(grp)) < 2) return(NA_real_)
  lr <- tryCatch(survdiff(Surv(time, status) ~ grp), error = function(e) NULL)
  if (is.null(lr)) return(NA_real_)
  1 - pchisq(lr$chisq, df = 1)
}

make_km_plot <- function(EL_k, time, status, beta_k, lrp, factor_label) {
  grp <- ifelse(EL_k > median(EL_k), "High loading", "Low loading")
  step_data <- do.call(rbind, lapply(c("High loading", "Low loading"), function(g) {
    idx    <- grp == g
    km_sub <- survfit(Surv(time[idx], status[idx]) ~ 1)
    data.frame(time = c(0, km_sub$time), surv = c(1, km_sub$surv),
               group = g, stringsAsFactors = FALSE)
  }))
  direction <- if (beta_k > 0) "Adverse" else "Protective"
  p_lab <- if (!is.na(lrp)) sprintf("Log-rank p = %.4f", lrp) else "Log-rank p = NA"
  ggplot(step_data, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 0.9) +
    scale_color_manual(values = c("High loading" = "#d62728", "Low loading" = "#1f77b4"),
                       name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, label = p_lab, size = 4) +
    labs(title = sprintf("%s  (β̂ = %+.4f, %s)", factor_label, beta_k, direction),
         x = "Follow-up time", y = "Survival probability") +
    ylim(0, 1) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 12, face = "bold"))
}

make_heatmap <- function(fit_id, title, n_genes = 40) {
  EF       <- fits[[fit_id]]$EF
  EBeta    <- fits[[fit_id]]$EBeta
  K        <- ncol(EF)
  active_k <- which(abs(EBeta) > BETA_THRESH)
  top_idx  <- order(apply(abs(EF), 1, max), decreasing = TRUE)[1:min(n_genes, length(genes_desurv))]
  top_g    <- genes_desurv[top_idx]

  df_wide <- cbind(
    data.frame(Gene = top_g, stringsAsFactors = FALSE),
    as.data.frame(EF[top_idx, , drop = FALSE]) |> setNames(paste0("F", seq_len(K)))
  )
  df_long <- pivot_longer(df_wide, -Gene, names_to = "Factor", values_to = "Weight")
  df_long$Factor <- factor(df_long$Factor, levels = paste0("F", seq_len(K)))
  df_long$Gene   <- factor(df_long$Gene,   levels = rev(top_g))

  flabels <- setNames(paste0("F", seq_len(K)), paste0("F", seq_len(K)))
  flabels[paste0("F", active_k)] <- paste0("F", active_k, " *")
  face_vec <- ifelse(paste0("F", seq_len(K)) %in% paste0("F", active_k), "bold", "plain")
  wmax <- max(abs(df_long$Weight), na.rm = TRUE)

  ggplot(df_long, aes(x = Factor, y = Gene, fill = Weight)) +
    geom_tile(color = "white", linewidth = 0.1) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#d6604d",
                         midpoint = 0, limits = c(-wmax, wmax), name = "Weight") +
    scale_x_discrete(labels = flabels) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(size = 11, face = face_vec),
          axis.text.y = element_text(size = 9.5),
          plot.title  = element_text(size = 12, face = "bold"),
          legend.key.height = unit(0.5, "cm"))
}

# --------------------------------------------------------------------------
# 3. Build and save the four figures
# --------------------------------------------------------------------------
for (nm in names(MODELS)) {
  m       <- MODELS[[nm]]
  fit     <- fits[[m$fit_id]]
  EBeta   <- fit$EBeta
  active  <- which(abs(EBeta) > BETA_THRESH)
  if (length(active) == 0) stop(sprintf("Model %s has no survival-active factors.", nm))

  # KM panels — one per active factor, clean "Factor k" labels (no D-label).
  km_panels <- lapply(active, function(k) {
    lrp <- logrank_p(fit$EL[, k], time_train, status_train)
    make_km_plot(EL_k = fit$EL[, k], time = time_train, status = status_train,
                 beta_k = EBeta[k], lrp = lrp, factor_label = sprintf("Factor %d", k))
  })
  km_combined <- patchwork::wrap_plots(km_panels, ncol = length(km_panels))
  ggsave(m$km_out, km_combined, width = 14, height = 8, dpi = 110)
  cat(sprintf("Wrote %s  (active factors: %s)\n", m$km_out, paste(active, collapse = ", ")))

  # Gene-program heatmap — characteristic-based title.
  hm <- make_heatmap(m$fit_id, title = m$heatmap_title, n_genes = m$n_genes)
  ggsave(m$heatmap_out, hm, width = 10, height = 10, dpi = 115)
  cat(sprintf("Wrote %s\n", m$heatmap_out))
}

# --------------------------------------------------------------------------
# 4. Appendix figures: convergence traces + β barplots, with CLEAN labels
#    (replace the report's "LB orig (M4)" / "YFB orig (M5)" with characteristic
#    descriptions so no internal M-codes are visible to the audience).
# --------------------------------------------------------------------------

# Map each report config to an audience-facing, characteristic-based label.
CLEAN_LABELS <- c(
  D1 = "LB · variance top-2000",
  D2 = "YFB · variance top-2000",
  D3 = "LB · DeSurv-aligned",
  D4 = "YFB · DeSurv-aligned",
  D5 = "YFB · DeSurv + cohort"
)
label_levels <- unname(CLEAN_LABELS)

# Convergence (RMSE) trace data frame.
rmse_df <- do.call(rbind, lapply(names(CLEAN_LABELS), function(id) {
  h <- fits[[id]]$history$rmse
  data.frame(Label = CLEAN_LABELS[[id]], Iter = seq_along(h), RMSE = h,
             stringsAsFactors = FALSE)
}))
rmse_df$Label <- factor(rmse_df$Label, levels = label_levels)

rmse_colors <- setNames(
  c("#4c78a8", "#f58518", "#72b7b2", "#e45756", "#b279a2"), label_levels)

rmse_plot <- ggplot(rmse_df, aes(x = Iter, y = RMSE, color = Label)) +
  geom_line(linewidth = 1.0) +
  scale_color_manual(values = rmse_colors, name = NULL) +
  labs(x = "Iteration", y = "Training RMSE") +
  theme_bw(base_size = 13) +
  theme(legend.position = "right", legend.text = element_text(size = 11))
ggsave(file.path(ASSETS, "rmse_traces.png"), rmse_plot, width = 9, height = 5.2, dpi = 120)
cat(sprintf("Wrote %s\n", file.path(ASSETS, "rmse_traces.png")))

# β coefficient data frame (Active = |β| above the survival threshold).
beta_df <- do.call(rbind, lapply(names(CLEAN_LABELS), function(id) {
  EBeta  <- fits[[id]]$EBeta
  active <- abs(EBeta) > BETA_THRESH
  data.frame(Label = CLEAN_LABELS[[id]], Factor = seq_along(EBeta),
             Beta = round(EBeta, 4), Active = active, stringsAsFactors = FALSE)
}))
beta_df$Label <- factor(beta_df$Label, levels = label_levels)
beta_df$Fill  <- ifelse(beta_df$Active,
                        ifelse(beta_df$Beta > 0, "Adverse", "Protective"), "Inactive")

beta_plot <- ggplot(beta_df, aes(x = factor(Factor), y = Beta, fill = Fill)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Adverse" = "#d62728", "Protective" = "#1f77b4",
                               "Inactive" = "#cccccc"), name = NULL) +
  facet_wrap(~ Label, scales = "free_x", ncol = 3) +
  labs(x = "Factor", y = expression(hat(beta)[k])) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(size = 11))
ggsave(file.path(ASSETS, "beta_plots.png"), beta_plot, width = 11, height = 6.4, dpi = 120)
cat(sprintf("Wrote %s\n", file.path(ASSETS, "beta_plots.png")))

cat("\nDone — figures use characteristic-based titles (no D/M-labels).\n")
