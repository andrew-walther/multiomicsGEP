# ============================================================
# Script: score_yfb_supervised_external.R
# Purpose: Attempt S-12: score the fixed S-11 fits (w_S = 1 YFB, and w_S = 0 loadings then Cox)
#          once on the five external PDAC cohorts.
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: survival, yaml; code/predict_ssbmf_yfb.R, code/preprocess_desurv.R,
#               code/concordance_ci.R, results/benchmark_sim/benchmark_helpers.R
# Run: Rscript results/benchmark_sim/score_yfb_supervised_external.R
# ============================================================
# Registration: docs/plans/2026-11-05-supervised-external-scoring-plan.md. No fit, no tuning:
# both models come from saved S-11 checkpoints and the local training matrix (git-ignored).
# Only aggregate rows are written to the output folder.

suppressPackageStartupMessages({ library(survival); library(yaml) })
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R"); source("code/predict_ssbmf_yfb.R"); source("code/concordance_ci.R")

S11  <- "results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04"
OUT  <- "results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05"
z    <- readRDS("data/yfb_supervised_real_2026-11-04.rds")   # training matrix (local only)
w1   <- readRDS(file.path(S11, "fits", "size_D4_N2218_w1_f0.rds"))
w0   <- readRDS(file.path(S11, "fits", "size_D4_N2218_w0_f0.rds"))
stopifnot(identical(w1$genes, w0$genes), length(w1$genes) == 2218L)
genes <- w1$genes

# Two-step coefficients: Cox of training survival on the w_S = 0 live program scores ----
# (same construction as twostep() in summarize_yfb_supervised_real.R)
EF0  <- w0$state$EF; live <- which(colSums(EF0^2) > 1e-8)
Ytr  <- z$Y[, match(genes, z$genes), drop = FALSE]
cf   <- coxph(Surv(z$time, z$status) ~ I(Ytr %*% EF0[, live, drop = FALSE]))
b0   <- rep(0, ncol(EF0)); b0[live] <- ifelse(is.finite(coef(cf)), coef(cf), 0)

#' Risk scores of both fixed models on one prepared cohort matrix
#' @param Y Cohort patients x D4 genes, z-scored within the cohort (missing genes 0).
#' @return List with `w1` (supervised) and `twostep` risk vectors; higher = higher hazard.
#' @examples
#' # s12_risk(Y_ext)
s12_risk <- function(Y) {
  list(w1 = as.numeric(predict_ssbmf_yfb(Y, w1$state$EF, w1$state$EBeta,
                                         EF_norms = w1$state$EF_norms)$risk_scores),
       twostep = drop(Y %*% EF0 %*% b0))
}

#' Calibration slope: Cox coefficient of survival on a frozen risk score
#' @param risk Risk score. @param time,status Outcome. @return Slope and its SE.
#' @examples
#' # cal_slope(risk, time, status)
cal_slope <- function(risk, time, status) {
  f <- coxph(Surv(time, status) ~ risk)
  c(slope = unname(coef(f)), se = sqrt(vcov(f)[1, 1]))
}

# Score each cohort once ----
rows <- list(); keep <- list()
for (co in cfg$pdac$external_cohorts) {
  raw <- load_pdac_raw(co, PDAC_DATA_ROOT)
  pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
                                  log_transform = PLATFORM_LOG_TRANSFORM[[co]], cohort_name = co,
                                  rank_transform = FALSE, per_platform_standardize = TRUE)
  Y <- matrix(0, raw$n, length(genes))                 # absent gene = cohort mean (0)
  hit <- match(genes, pre$gene_names)
  Y[, !is.na(hit)] <- pre$Y[, hit[!is.na(hit)], drop = FALSE]
  r <- s12_risk(Y)
  keep[[co]] <- list(w1 = r$w1, twostep = r$twostep, time = raw$time, status = raw$status)
  for (arm in names(r)) {
    ci <- bootstrap_concordance_ci(r[[arm]], raw$time, raw$status, B = 2000, seed = 1, flip = FALSE)
    cs <- cal_slope(r[[arm]], raw$time, raw$status)
    rows[[length(rows) + 1]] <- data.frame(cohort = co, model = arm, n = raw$n,
      events = sum(raw$status), missing_genes = sum(is.na(hit)),
      C = frozen_reverse_cindex(r[[arm]], raw$time, raw$status), lo = ci$lower, hi = ci$upper,
      cal_slope = cs[["slope"]], cal_se = cs[["se"]], risk_sd = sd(r[[arm]]))
  }
}
res <- do.call(rbind, rows)

# Paired differences (supervised - two-step), per cohort and pooled ----
diff <- do.call(rbind, lapply(names(keep), function(co) {
  k <- keep[[co]]
  d <- bootstrap_concordance_diff_ci(k$w1, k$twostep, k$time, k$status, B = 2000, seed = 1,
                                     flip_a = FALSE, flip_b = FALSE)
  data.frame(cohort = co, diff = d$estimate, lo = d$lower, hi = d$upper, se = d$se)
}))
pooled <- function(cos, label) {
  d <- bootstrap_concordance_diff_ci_stratified(lapply(keep[cos], `[[`, "w1"),
         lapply(keep[cos], `[[`, "twostep"), lapply(keep[cos], `[[`, "time"),
         lapply(keep[cos], `[[`, "status"), B = 2000, seed = 1)
  data.frame(cohort = label, diff = d$estimate, lo = d$lower, hi = d$upper, se = d$se)
}
diff <- rbind(diff, pooled(names(keep), "pooled (stratified, 5 cohorts)"),
              pooled(setdiff(names(keep), "PACA_AU_seq"), "pooled without PACA_AU_seq"))

# Pooled HR per SD of the risk score (added after registration, 2026-11-05) ----
# DeSurv's headline external metric: stratified Cox (strata = cohort) of survival on the risk
# score divided by its pooled SD across all external patients. Descriptive only.
pooled_hr <- do.call(rbind, lapply(c("w1", "twostep"), function(arm) {
  df <- do.call(rbind, lapply(names(keep), function(co)
    data.frame(time = keep[[co]]$time, event = keep[[co]]$status, cohort = co, risk = keep[[co]][[arm]])))
  df$risk_z <- df$risk / sd(df$risk)
  f <- summary(coxph(Surv(time, event) ~ risk_z + strata(cohort), data = df))
  data.frame(model = arm, n = nrow(df), hr_per_sd = f$conf.int[1, "exp(coef)"],
             lo = f$conf.int[1, "lower .95"], hi = f$conf.int[1, "upper .95"],
             p = f$coefficients[1, "Pr(>|z|)"])
}))

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(pooled_hr, file.path(OUT, "external_pooled_hr.csv"), row.names = FALSE)
print(pooled_hr, digits = 3)
write.csv(res, file.path(OUT, "external_cindex.csv"), row.names = FALSE)
write.csv(diff, file.path(OUT, "external_paired_diff.csv"), row.names = FALSE)
writeLines(c(system("git rev-parse HEAD", intern = TRUE), format(Sys.time()),
             sprintf("two-step Cox on %d live programs", length(live))),
           file.path(OUT, "provenance.txt"))
print(res, digits = 3); print(diff, digits = 3)
