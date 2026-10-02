# ============================================================
# Script: make_km_figures.R
# Purpose: Kaplan-Meier figures for the recommended single-modality YFB fit
#          (D4) in the five held-out PDAC cohorts: (1) risk-score tertiles per
#          cohort; (2) Program 7 (adverse) and Program 3 (protective)
#          projection tertiles, pooled across cohorts.
# Author: Andrew Walther
# Created: 2026-10-02
# Dependencies: survival, ggplot2, yaml; results/benchmark_sim/benchmark_helpers.R,
#               code/preprocess_desurv.R, code/predict_cox_on_yf.R
# Run: PDAC_DATA_ROOT=data/PDAC_data Rscript results/benchmark_sim/make_km_figures.R
# ============================================================
#
# Orientation. The saved D4 fit predates the 2026-09-04 sign-correction fix,
# so its eta = Z beta is a good-prognosis score (beta_7 = -0.040 on the
# adverse program, beta_3 = +0.012 on the protective one). Raw concordance
# with reverse = TRUE is below 0.5 in all five cohorts. The risk score used
# here is therefore -eta: one global orientation, fixed before looking at
# any validation outcome. It reproduces the reported external C-indices.
#
# Program direction. Adverse/protective labels follow the marginal survival
# association of each program's projection (DECISIONS.md 2026-06-16), so
# "high Program 7" and "high Program 3" are plotted directly.

suppressPackageStartupMessages({ library(survival); library(ggplot2) })
cfg <- yaml::read_yaml("config/globals.yml")   # benchmark_helpers.R reads cfg
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R"); source("code/predict_cox_on_yf.R")

# load_pdac_raw() symlinks the root into a temp dir, so it must be absolute
pdac_root <- normalizePath(Sys.getenv("PDAC_DATA_ROOT", "data/PDAC_data"))
out_dir <- "results/benchmark_sim/outputs/km_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dc <- "results/benchmark_sim/outputs/desurv_comparison"
fit <- readRDS(file.path(dc, "desurv_comparison_fits.rds"))[["D4"]]
train_genes <- readRDS(file.path(dc, "d4_gene_names.rds"))
saved <- readRDS(file.path(dc, "desurv_comparison_riskscores.rds"))[["D4"]]
cohort_labels <- c(Dijk = "Dijk", Moffitt_GEO_array = "Moffitt (array)",
                   PACA_AU_array = "PACA-AU (array)", PACA_AU_seq = "PACA-AU (RNA-seq)",
                   Puleo_array = "Puleo (array)")

# Rebuild projections exactly as run_desurv_comparison.R did ----
rows <- list()
for (cohort in cfg$pdac$external_cohorts) {
  raw <- load_pdac_raw(cohort, pdac_root)
  pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
                                  log_transform = cfg$pdac$platform_log_transform[[cohort]],
                                  cohort_name = cohort, rank_transform = FALSE,
                                  per_platform_standardize = TRUE)
  common <- intersect(train_genes, pre$gene_names)
  Y <- pre$Y[, match(common, pre$gene_names), drop = FALSE]
  EF <- fit$EF[match(common, train_genes), , drop = FALSE]
  eta <- predict_cox_on_yf(Y, EF, fit$EBeta, EF_norms = fit$EF_norms)$risk_scores
  # Fail loud if the rebuild does not reproduce the saved validation scores
  if (length(eta) != length(saved[[cohort]]$risk) ||
      max(abs(eta - saved[[cohort]]$risk)) > 1e-8) {
    stop(cohort, ": rebuilt risk scores differ from the saved D4 scores (max diff ",
         signif(max(abs(eta - saved[[cohort]]$risk)), 3), ").")
  }
  Z <- sweep(Y %*% EF, 2, fit$EF_norms, "/")       # normalized projections, n x K
  rows[[cohort]] <- data.frame(cohort = cohort_labels[[cohort]], time = raw$time,
                               status = raw$status, risk = -eta,
                               program7 = Z[, 7], program3 = Z[, 3])
}
km_data <- do.call(rbind, rows)

# Tertiles within each cohort (cohorts are on different platforms)
tertile <- function(x, labels) {
  cut(x, stats::quantile(x, c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = labels)
}
km_data <- do.call(rbind, lapply(split(km_data, km_data$cohort), function(d) {
  d$risk_group <- tertile(d$risk, c("Low", "Middle", "High"))
  d$p7_group <- tertile(d$program7, c("Low", "Middle", "High"))
  d$p3_group <- tertile(d$program3, c("Low", "Middle", "High"))
  d
}))

#' Kaplan-Meier step data for ggplot
km_steps <- function(d, group, facet = NULL) {
  f <- survfit(stats::as.formula(paste("Surv(time, status) ~", group)), data = d)
  s <- summary(f, censored = TRUE)
  out <- data.frame(time = c(0, s$time), surv = c(1, s$surv),
                    strata = c(as.character(s$strata[1]), as.character(s$strata)))
  # restart each stratum at (0, 1)
  out <- do.call(rbind, lapply(split(data.frame(time = s$time, surv = s$surv,
                                                n.censor = s$n.censor, strata = s$strata),
                                     s$strata), function(x) {
    rbind(data.frame(time = 0, surv = 1, n.censor = 0, strata = x$strata[1]), x)
  }))
  out$group <- sub(".*=", "", as.character(out$strata))
  if (!is.null(facet)) out$facet <- facet
  out
}
cols3 <- c(Low = "#2a78d6", Middle = "#eda100", High = "#e34948")
# Program activity is a magnitude: one hue, light -> dark (validated blue ramp)
seq3 <- c(Low = "#9ec5f4", Middle = "#3987e5", High = "#104281")
fmt_p <- function(p) if (p < 1e-4) "p < 0.0001" else sprintf("p = %s", format(signif(p, 2), scientific = FALSE))
theme_set(theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                                legend.position = "bottom"))

# Figure 1: risk tertiles per held-out cohort ----
lab <- do.call(rbind, lapply(split(km_data, km_data$cohort), function(d) {
  lr <- survdiff(Surv(time, status) ~ risk_group, data = d)
  p <- stats::pchisq(lr$chisq, df = 2, lower.tail = FALSE)
  C <- concordance(Surv(time, status) ~ risk, data = d, reverse = TRUE)$concordance
  data.frame(facet = d$cohort[1], n = nrow(d),
             label = sprintf("n = %d   C = %.3f\nlog-rank %s", nrow(d), C, fmt_p(p)))
}))
steps <- do.call(rbind, lapply(split(km_data, km_data$cohort), function(d)
  km_steps(d, "risk_group", d$cohort[1])))
steps$group <- factor(steps$group, levels = c("Low", "Middle", "High"))
p1 <- ggplot(steps, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.7) +
  geom_point(data = steps[steps$n.censor > 0, ], shape = 3, size = 1.2, show.legend = FALSE) +
  geom_text(data = lab, aes(x = Inf, y = 1, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.1, size = 3) +
  facet_wrap(~facet, nrow = 2) +
  scale_colour_manual(values = cols3, name = "Risk tertile (within cohort)") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Months", y = "Overall survival",
       title = "Frozen single-modality model in five held-out cohorts")
ggsave(file.path(out_dir, "km_risk_tertiles_external.png"), p1, width = 9, height = 6, dpi = 150)

# Figure 2: Program 7 and Program 3 tertiles, pooled across cohorts ----
prog <- rbind(
  transform(km_steps(km_data, "p7_group", "Program 7 (adverse; basal-like)")),
  transform(km_steps(km_data, "p3_group", "Program 3 (protective; classical)")))
prog$group <- factor(prog$group, levels = c("Low", "Middle", "High"))
plab <- do.call(rbind, lapply(c("p7_group", "p3_group"), function(g) {
  # log-rank stratified by cohort: tertiles are within cohort
  lr <- survdiff(stats::as.formula(paste("Surv(time, status) ~", g, "+ strata(cohort)")),
                 data = km_data)
  p <- stats::pchisq(lr$chisq, df = 2, lower.tail = FALSE)
  data.frame(facet = if (g == "p7_group") "Program 7 (adverse; basal-like)" else
               "Program 3 (protective; classical)",
             label = sprintf("n = %d (5 cohorts)\ncohort-stratified log-rank %s",
                             nrow(km_data), fmt_p(p)))
}))
p2 <- ggplot(prog, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.7) +
  geom_text(data = plab, aes(x = Inf, y = 1, label = label), inherit.aes = FALSE,
            hjust = 1.05, vjust = 1.1, size = 3) +
  facet_wrap(~facet, nrow = 1) +
  scale_colour_manual(values = seq3, name = "Program projection tertile (within cohort)") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Months", y = "Overall survival",
       title = "Survival by program activity, held-out cohorts pooled")
ggsave(file.path(out_dir, "km_programs_pooled_external.png"), p2, width = 9, height = 4.2, dpi = 150)

write.csv(lab[, c("facet", "n", "label")], file.path(out_dir, "km_risk_tertile_stats.csv"), row.names = FALSE)
print(lab[, c("facet", "label")]); print(plab)
