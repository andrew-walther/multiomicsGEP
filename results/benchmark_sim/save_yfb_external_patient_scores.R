# ============================================================
# Script: save_yfb_external_patient_scores.R
# Purpose: Save per-patient external risk scores and program scores from the fixed S-11 fits
#          (same procedure as S-12) under the git-ignored data/ folder, for KM figures.
# Author: Andrew Walther
# Created: 2026-11-06
# Dependencies: survival, yaml; code/predict_ssbmf_yfb.R, code/preprocess_desurv.R,
#               code/concordance_ci.R, results/benchmark_sim/benchmark_helpers.R
# Run: Rscript results/benchmark_sim/save_yfb_external_patient_scores.R
# ============================================================
# Project the saved S-11 loadings using the S-12 preprocessing. The two-step
# Cox coefficients are reconstructed from training data, as in S-12.
# Patient-level output excludes identifiers and stays under ignored data/.

suppressPackageStartupMessages({ library(survival); library(yaml) })

# Parameters ----
source("results/benchmark_sim/exhibit_rank.R")   # EXHIBIT_K: cache the scores of another rank
S11 <- exhibit_out("results/benchmark_sim/outputs/yfb_supervised_real_2026-11-04", "s11")
FIT_W1 <- file.path(S11, "fits", "size_D4_N2218_w1_f0.rds")
FIT_W0 <- file.path(S11, "fits", "size_D4_N2218_w0_f0.rds")
TRAIN <- "data/yfb_supervised_real_2026-11-04.rds"
CINDEX_REF <- file.path(exhibit_out("results/benchmark_sim/outputs/yfb_supervised_external_2026-11-05",
                                    "s12"), "external_cindex.csv")
OUT_RDS <- exhibit_data("data/yfb_external_patient_scores_ws1.rds", "external_patient_scores_ws1.rds")

cfg <- yaml::read_yaml("config/globals.yml")
source("results/benchmark_sim/benchmark_helpers.R")
source("code/preprocess_desurv.R"); source("code/predict_ssbmf_yfb.R"); source("code/concordance_ci.R")

z  <- readRDS(TRAIN)
w1 <- readRDS(FIT_W1)
w0 <- readRDS(FIT_W0)
stopifnot(identical(w1$genes, w0$genes))
genes <- w1$genes

# Two-step coefficients (as in S-12): Cox of training survival on w_S = 0 live program scores ----
EF0  <- w0$state$EF; live <- which(colSums(EF0^2) > 1e-8)
Ytr  <- z$Y[, match(genes, z$genes), drop = FALSE]
cf   <- coxph(Surv(z$time, z$status) ~ I(Ytr %*% EF0[, live, drop = FALSE]))
b0   <- rep(0, ncol(EF0)); b0[live] <- ifelse(is.finite(coef(cf)), coef(cf), 0)

# Score each external cohort, keep patient-level rows ----
rows <- list()
for (co in cfg$pdac$external_cohorts) {
  raw <- load_pdac_raw(co, PDAC_DATA_ROOT)
  pre <- preprocess_desurv_cohort(raw$Y, raw$gene_names, top_n = NULL,
                                  log_transform = PLATFORM_LOG_TRANSFORM[[co]], cohort_name = co,
                                  rank_transform = FALSE, per_platform_standardize = TRUE)
  Y <- matrix(0, raw$n, length(genes))
  hit <- match(genes, pre$gene_names)
  Y[, !is.na(hit)] <- pre$Y[, hit[!is.na(hit)], drop = FALSE]
  p <- predict_ssbmf_yfb(Y, w1$state$EF, w1$state$EBeta, EF_norms = w1$state$EF_norms)
  Zm <- p$ZF_test; colnames(Zm) <- paste0("Z", seq_len(ncol(Zm)))
  rows[[co]] <- data.frame(cohort = co, time = raw$time, event = raw$status,
                           risk = as.numeric(p$risk_scores), risk_twostep = drop(Y %*% EF0 %*% b0),
                           Zm, stringsAsFactors = FALSE)
}
sc <- do.call(rbind, rows); rownames(sc) <- NULL

# Verify frozen scores against the registered counts and concordance ----
ref <- read.csv(CINDEX_REF)
for (i in seq_len(nrow(ref))) {
  d <- rows[[ref$cohort[i]]]
  risk <- if (ref$model[i] == "w1") d$risk else d$risk_twostep
  stopifnot(nrow(d) == ref$n[i], sum(d$event) == ref$events[i],
            abs(frozen_reverse_cindex(risk, d$time, d$event) - ref$C[i]) < 1e-10)
}

if (file.exists(OUT_RDS)) {
  previous <- readRDS(OUT_RDS)
  if (!isTRUE(all.equal(previous[names(sc)], sc, tolerance = 1e-10)))
    stop("The existing ignored cache differs from the fixed S-12 projections.")
  cat("Verified existing ignored score cache:", OUT_RDS, "\n")
} else {
  saveRDS(sc, OUT_RDS)
  cat("Saved", nrow(sc), "external records to", OUT_RDS, "\n")
}
