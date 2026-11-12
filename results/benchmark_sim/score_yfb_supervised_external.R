# ============================================================
# Script: score_yfb_supervised_external.R
# Purpose: Attempt S-12: score the fixed S-11 fits (w_S = 1 YFB, and w_S = 0 loadings then Cox)
#          once on the five external PDAC cohorts.
# Author: Andrew Walther
# Created: 2026-11-05
# Dependencies: survival, yaml; code/predict_ssbmf_yfb.R, code/preprocess_desurv.R,
#               code/concordance_ci.R, results/benchmark_sim/benchmark_helpers.R
# Run: Rscript results/benchmark_sim/score_yfb_supervised_external.R [--donor-only]
# ============================================================
# Registration: docs/plans/2026-11-05-supervised-external-scoring-plan.md. No fit, no tuning:
# both models come from saved S-11 checkpoints and the local training matrix (git-ignored).
# Only aggregate rows are written to the output folder.

suppressPackageStartupMessages({ library(survival); library(yaml) })
cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R"); source("code/predict_ssbmf_yfb.R"); source("code/concordance_ci.R")

donor_only <- "--donor-only" %in% commandArgs(trailingOnly = TRUE)
if (!dir.exists(PDAC_DATA_ROOT)) stop("PDAC_DATA_ROOT must point to local PDAC data.")

source("results/benchmark_sim/exhibit_rank.R")   # EXHIBIT_K: score the fits of another rank
S11  <- exhibit_out("results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04", "s11")
OUT  <- exhibit_out("results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05", "s12")
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
  ids <- as.character(raw$sampID)
  if (length(ids) != raw$n || anyNA(ids) || any(!nzchar(ids)))
    stop("Missing or unaligned sample identifiers in ", co)
  # Resolve sample IDs to biological donors before combining the two assays.
  # The loader appends _seq to RNA-seq sample IDs after reading this metadata.
  if (co %in% c("PACA_AU_array", "PACA_AU_seq")) {
    meta <- readRDS(file.path(PDAC_DATA_ROOT, "original", paste0(co, ".rds")))$sampInfo
    stopifnot(all(c("icgc_sample_id", "icgc_donor_id") %in% names(meta)),
              !anyDuplicated(meta$icgc_sample_id))
    sample_row <- match(sub("_seq$", "", ids), meta$icgc_sample_id)
    donor_id <- as.character(meta$icgc_donor_id[sample_row])
    if (length(donor_id) != raw$n || anyNA(donor_id) || any(!nzchar(donor_id)))
      stop("Could not match all PACA-AU samples to donor metadata in ", co)
    donor <- paste0("PACA_AU:", donor_id)
  } else {
    donor <- paste0(co, ":", ids)
  }
  keep[[co]] <- list(w1 = r$w1, twostep = r$twostep, time = raw$time,
                     status = raw$status, donor = donor)
  if (donor_only) next
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

#' Summarize the equal-cohort-weighted paired difference.
#' @param cos Cohort keys in the in-memory score list.
#' @param label Display label for the output row.
#' @param donors Optional donor ID vectors in the order of cos.
#' @return One aggregate row with difference, percentile limits and SE.
#' @examples
#' # pooled(names(keep), "pooled", lapply(keep, `[[`, "donor"))
pooled <- function(cos, label, donors = NULL) {
  d <- bootstrap_concordance_diff_ci_stratified(lapply(keep[cos], `[[`, "w1"),
         lapply(keep[cos], `[[`, "twostep"), lapply(keep[cos], `[[`, "time"),
         lapply(keep[cos], `[[`, "status"), B = 2000, seed = 1, donor_list = donors)
  data.frame(cohort = label, diff = d$estimate, lo = d$lower, hi = d$upper, se = d$se)
}
# Compare the same equal-cohort-weighted estimand using a donor bootstrap.
donor_ids <- lapply(keep, `[[`, "donor")
shared <- length(intersect(donor_ids$PACA_AU_array, donor_ids$PACA_AU_seq))
if (shared != 47L) stop("PACA-AU donor overlap differs from the preprocessing audit: ", shared)
donor_diff <- pooled(names(keep), "pooled (donors, 5 cohorts)", donor_ids)
donor_diff$n_donors <- length(unique(unlist(donor_ids, use.names = FALSE)))
donor_diff$n_shared_paca_au <- shared
if (donor_only) {
  prior <- read.csv(file.path(OUT, "external_paired_diff.csv"))
  prior_diff <- prior$diff[prior$cohort == "pooled (stratified, 5 cohorts)"]
  stopifnot(length(prior_diff) == 1L, abs(prior_diff - donor_diff$diff) < 1e-10)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
write.csv(donor_diff, file.path(OUT, "external_paired_diff_donor.csv"), row.names = FALSE)
print(donor_diff, digits = 4)
if (donor_only) quit(save = "no")

# Paired differences (supervised - two-step), per cohort and pooled ----
diff <- do.call(rbind, lapply(names(keep), function(co) {
  k <- keep[[co]]
  d <- bootstrap_concordance_diff_ci(k$w1, k$twostep, k$time, k$status, B = 2000, seed = 1,
                                     flip_a = FALSE, flip_b = FALSE)
  data.frame(cohort = co, diff = d$estimate, lo = d$lower, hi = d$upper, se = d$se)
}))
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
