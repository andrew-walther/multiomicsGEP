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
                     heatmap_title = "Recommended model · projection (YF)β, survival-ranked gene selection — top 30 genes",
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
logrank_p <- function(score_k, time, status) {
  grp <- ifelse(score_k > median(score_k), "High", "Low")
  if (length(unique(grp)) < 2) return(NA_real_)
  lr <- tryCatch(survdiff(Surv(time, status) ~ grp), error = function(e) NULL)
  if (is.null(lr)) return(NA_real_)
  1 - pchisq(lr$chisq, df = 1)
}

# Empirical survival direction of a program's activation score.
#   We stratify by the (YF) PROJECTION SCORE — the quantity the risk model
#   actually scores on (eta = (YF)beta), NOT the genomics loading E[L]. The
#   Adverse/Protective label is read off the DATA (sign of the Cox coefficient
#   on the projection score): coef > 0 => higher activation = worse survival
#   = Adverse. This is the marginal direction the KM curve displays; we report
#   it rather than the joint-model beta sign, which can differ under suppression
#   when programs are correlated.
direction_label <- function(score_k, time, status) {
  cc <- tryCatch(coef(survival::coxph(Surv(time, status) ~ score_k)),
                 error = function(e) NA_real_)
  if (is.na(cc)) return("Survival-active")
  if (cc > 0) "Adverse" else "Protective"
}

make_km_plot <- function(score_k, time, status, lrp, factor_label, direction) {
  grp <- ifelse(score_k > median(score_k), "High activation", "Low activation")
  step_data <- do.call(rbind, lapply(c("High activation", "Low activation"), function(g) {
    idx    <- grp == g
    km_sub <- survfit(Surv(time[idx], status[idx]) ~ 1)
    data.frame(time = c(0, km_sub$time), surv = c(1, km_sub$surv),
               group = g, stringsAsFactors = FALSE)
  }))
  p_lab <- if (!is.na(lrp)) sprintf("Log-rank  p = %s",
                                    formatC(lrp, format = "g", digits = 2)) else "Log-rank p = NA"
  ggplot(step_data, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 1.2) +
    scale_color_manual(values = c("High activation" = "#d62728", "Low activation" = "#1f77b4"),
                       name = NULL) +
    annotate("text", x = Inf, y = Inf, hjust = 1.08, vjust = 1.5, label = p_lab,
             size = 7, fontface = "bold") +
    labs(title = sprintf("%s  —  %s program", factor_label, direction),
         x = "Follow-up time", y = "Survival probability") +
    ylim(0, 1) +
    theme_bw(base_size = 17) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 16),
          plot.title  = element_text(size = 18, face = "bold"))
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
# DeSurv-aligned preprocessed training matrix → (YF) projection scores.
# ZF = Y · EF is the n×K matrix of OBSERVED projection scores the YFB risk model
# scores on (eta = ZF · beta). We stratify KM curves and read survival direction
# from ZF[,k], NOT from the genomics loading E[L].
Y_train <- pp_desurv$Y

active_gene_rows <- list()   # accumulate top genes per active factor for the deck

for (nm in names(MODELS)) {
  m       <- MODELS[[nm]]
  fit     <- fits[[m$fit_id]]
  EBeta   <- fit$EBeta
  EF      <- fit$EF
  active  <- which(abs(EBeta) > BETA_THRESH)   # survival-active: |beta_hat| > threshold
  if (length(active) == 0) stop(sprintf("Model %s has no survival-active factors.", nm))
  ZF <- Y_train %*% EF                          # (YF) projection scores, n×K

  # KM panels — one per survival-active factor, stratified by the (YF) projection
  # score and labeled by the empirical survival direction.
  km_panels <- lapply(active, function(k) {
    score_k <- ZF[, k]
    lrp <- logrank_p(score_k, time_train, status_train)
    dir <- direction_label(score_k, time_train, status_train)
    make_km_plot(score_k = score_k, time = time_train, status = status_train,
                 lrp = lrp, factor_label = sprintf("Program %d", k), direction = dir)
  })
  km_combined <- patchwork::wrap_plots(km_panels, ncol = length(km_panels))
  ggsave(m$km_out, km_combined, width = 15, height = 8.2, dpi = 120)
  dirs <- vapply(active, function(k) direction_label(ZF[, k], time_train, status_train), character(1))
  cat(sprintf("Wrote %s  (active factors: %s; directions: %s)\n",
              m$km_out, paste(active, collapse = ", "), paste(dirs, collapse = ", ")))

  # Export top-weighted genes per active factor (recommended model only) for the deck.
  if (nm == "recommended") {
    for (k in active) {
      ord  <- order(abs(EF[, k]), decreasing = TRUE)[1:8]
      active_gene_rows[[length(active_gene_rows) + 1]] <- data.frame(
        factor    = k,
        direction = direction_label(ZF[, k], time_train, status_train),
        abs_beta  = round(abs(EBeta[k]), 4),
        logrank_p = signif(logrank_p(ZF[, k], time_train, status_train), 2),
        top_genes = paste(genes_desurv[ord], collapse = ", "),
        stringsAsFactors = FALSE)
    }
  }

  # Gene-program heatmap — characteristic-based title.
  hm <- make_heatmap(m$fit_id, title = m$heatmap_title, n_genes = m$n_genes)
  ggsave(m$heatmap_out, hm, width = 10, height = 10, dpi = 115)
  cat(sprintf("Wrote %s\n", m$heatmap_out))
}

# Persist the active-factor gene table so the .qmd can render it without needing
# PDAC raw data at render time.
genes_csv <- file.path(ASSETS, "active_factor_genes.csv")
write.csv(do.call(rbind, active_gene_rows), genes_csv, row.names = FALSE)
cat(sprintf("Wrote %s\n", genes_csv))

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

# β coefficient data frame. We plot |β̂| (activity MAGNITUDE) and flag which
# factors clear the survival-activity threshold. We deliberately do NOT colour by
# the sign of β̂ here: in the (YF)β parameterization the joint coefficient sign can
# differ from a program's marginal survival direction (suppression among
# correlated programs), so the Adverse/Protective call is read off the KM curves
# (sign of the Cox coefficient on the projection score), not off β̂'s sign.
beta_df <- do.call(rbind, lapply(names(CLEAN_LABELS), function(id) {
  EBeta  <- fits[[id]]$EBeta
  active <- abs(EBeta) > BETA_THRESH
  data.frame(Label = CLEAN_LABELS[[id]], Factor = seq_along(EBeta),
             AbsBeta = round(abs(EBeta), 4),
             Status = ifelse(active, "Survival-active", "Expression-only"),
             stringsAsFactors = FALSE)
}))
beta_df$Label <- factor(beta_df$Label, levels = label_levels)

beta_plot <- ggplot(beta_df, aes(x = factor(Factor), y = AbsBeta, fill = Status)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = BETA_THRESH, color = "black", linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c("Survival-active" = "#d62728",
                               "Expression-only" = "#cccccc"), name = NULL) +
  facet_wrap(~ Label, scales = "free_x", ncol = 3) +
  labs(x = "Factor", y = expression("|" * hat(beta)[k] * "|  (survival-activity magnitude)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.text = element_text(size = 11))
ggsave(file.path(ASSETS, "beta_plots.png"), beta_plot, width = 11, height = 6.4, dpi = 120)
cat(sprintf("Wrote %s\n", file.path(ASSETS, "beta_plots.png")))

cat("\nDone — figures use characteristic-based titles (no D/M-labels).\n")
