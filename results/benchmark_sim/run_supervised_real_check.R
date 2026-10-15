# ============================================================
# Script: run_supervised_real_check.R
# Purpose: Real-data criterion for the survival-supervised single-modality
#          configuration (Workstream C, C7): lower 95% limit of the external
#          Delta C against the reported alpha_F = 0 model (D4) above -0.02.
# Author: Andrew Walther
# Created: 2026-10-15
# Dependencies: survival, ebnm, yaml; code/fit_cox_on_yf.R,
#               code/preprocess_desurv.R, code/predict_cox_on_yf.R,
#               code/concordance_ci.R
# Run: Rscript results/benchmark_sim/run_supervised_real_check.R
#      (PDAC_DATA_ROOT or the globals.yml default must point at the data)
# ============================================================
#
# Training data and preprocessing are D4's (run_desurv_comparison.R): TCGA_PAAD
# + CPTAC, per-platform z-standardization, combined-rank selection of the top
# 3000 genes per cohort, no further normalization; K = 7.
#
# Arms:
#   supervised  fit_cox_on_yf(canonicalize, w_S = 1, init_supervised,
#               shared point-Laplace beta prior, 10-iteration warm-up,
#               point-Laplace L/F)  -- DECISIONS.md 2026-10-15
#   alpha0_pl   alpha_F = 0 with the point-Laplace L/F prior (separates the
#               effect of supervision from the change of prior)
#   D4          the reported alpha_F = 0 fit (cached risk scores; frozen
#               orientation risk = -eta, DECISIONS.md 2026-10-02)
#
# Primary estimand (fixed before running): mean over the 5 external cohorts
# of the per-cohort Delta C (supervised - D4), 95% percentile interval from a
# cohort-stratified paired bootstrap (B = 2000). Criterion: lower limit > -0.02.
# Per-cohort and concatenated differences are reported alongside.
#
# Memory/runtime: n ~ 280, p = 3000, K = 7; a few minutes per fit on a laptop.

suppressPackageStartupMessages({ library(survival); library(yaml) })
source("code/predict_cox_on_yf.R"); source("code/concordance_ci.R")
suppressMessages(tryCatch(source("code/fit_cox_on_yf.R"), error = function(e) NULL))
source("code/preprocess_desurv.R")
cfg <- yaml::read_yaml("config/globals.yml")
# load_pdac_raw(), PDAC_DATA_ROOT, PLATFORM_LOG_TRANSFORM (needs cfg)
source("results/benchmark_sim/benchmark_helpers.R")
OUT_DIR <- "results/benchmark_sim/outputs/supervised_real_check"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
K <- cfg$benchmark$k_merged_yfb_desurv   # 7

# 1. Training data (D4 preprocessing) ----
TRAIN <- cfg$pdac$training_cohorts
train_raw <- lapply(setNames(TRAIN, TRAIN), function(ds) load_pdac_raw(ds, PDAC_DATA_ROOT))
time_train   <- unlist(lapply(TRAIN, function(ds) train_raw[[ds]]$time))
status_train <- unlist(lapply(TRAIN, function(ds) train_raw[[ds]]$status))
pp <- preprocess_merged_cohorts(
  cohort_raw_list = train_raw, log_transform_flags = PLATFORM_LOG_TRANSFORM[TRAIN],
  top_n = cfg$preprocessing$top_n_genes_desurv, rank_transform = FALSE,
  per_platform_standardize = TRUE, normalize_method = "none",
  selection_per_cohort = TRUE, selection_method = "combined_rank")
d4_genes <- readRDS("results/benchmark_sim/outputs/desurv_comparison/d4_gene_names.rds")
stopifnot(identical(pp$gene_names, d4_genes))   # same gene set as the reported fit
Y_train <- pp$Y
cat(sprintf("Training: n = %d, p = %d, events = %d\n", nrow(Y_train), ncol(Y_train),
            sum(status_train)))

# 2. Fits ----
common <- list(Y_train, time_train, status_train, K = K, max_iter = cfg$cavi$max_iter,
               prior_LF = "point_laplace", verbose = FALSE, trace_supervision = TRUE)
arms <- list(
  supervised = c(common, list(canonicalize = TRUE, w_S = 1, init_supervised = TRUE,
                              beta_prior_mode = "shared", prior_beta = "point_laplace",
                              beta_warmup = 10L)),
  alpha0_pl  = c(common, list(alpha_F = 0, prior_beta = "normal", alpha = cfg$benchmark$alpha))
)
fits <- lapply(names(arms), function(a) {
  set.seed(42L); t0 <- proc.time()[["elapsed"]]
  f <- do.call(fit_cox_on_yf, arms[[a]])
  cat(sprintf("%-10s %.0f s | converged %s in %d | beta %s\n", a,
              proc.time()[["elapsed"]] - t0, f$history$converged, f$history$n_iter,
              paste(sprintf("%+.3f", f$EBeta), collapse = " ")))
  f
})
names(fits) <- names(arms)
saveRDS(fits, file.path(OUT_DIR, "fits.rds"))

# 3. External scoring ----
d4_risk <- readRDS("results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_riskscores.rds")[["D4"]]
risk <- list()
for (co in cfg$pdac$external_cohorts) {
  raw <- load_pdac_raw(co, PDAC_DATA_ROOT)
  pre <- preprocess_desurv_cohort(Y = raw$Y, gene_names = raw$gene_names, top_n = NULL,
                                  log_transform = PLATFORM_LOG_TRANSFORM[[co]], cohort_name = co,
                                  rank_transform = FALSE, per_platform_standardize = TRUE)
  g  <- intersect(d4_genes, pre$gene_names)
  Ye <- pre$Y[, match(g, pre$gene_names), drop = FALSE]
  idx <- match(g, d4_genes)
  stopifnot(identical(raw$time, d4_risk[[co]]$time))
  risk[[co]] <- c(lapply(fits, function(f)
    predict_cox_on_yf(Ye, f$EF[idx, , drop = FALSE], f$EBeta, EF_norms = f$EF_norms)$risk_scores),
    list(D4 = -d4_risk[[co]]$risk, time = raw$time, status = raw$status))
}
saveRDS(risk, file.path(OUT_DIR, "riskscores.rds"))

# 4. Per-cohort C and paired differences (orientation frozen at fit time) ----
cindex <- function(r, co) frozen_reverse_cindex(r, risk[[co]]$time, risk[[co]]$status)
rows <- list()
for (co in names(risk)) for (a in c("supervised", "alpha0_pl")) {
  ci <- bootstrap_concordance_diff_ci(risk[[co]][[a]], risk[[co]]$D4, risk[[co]]$time,
                                      risk[[co]]$status, B = 2000, seed = 1,
                                      flip_a = FALSE, flip_b = FALSE)
  rows[[length(rows) + 1]] <- data.frame(cohort = co, arm = a, n = length(risk[[co]]$time),
    C_arm = cindex(risk[[co]][[a]], co), C_D4 = cindex(risk[[co]]$D4, co),
    dC = ci$estimate, dC_lower = ci$lower, dC_upper = ci$upper)
}
per <- do.call(rbind, rows)

# Primary: mean over cohorts of per-cohort Delta C, cohort-stratified bootstrap
strat_mean_diff <- function(a, B = 2000, seed = 1) {
  set.seed(seed)
  point <- mean(per$dC[per$arm == a])
  boot <- replicate(B, mean(vapply(names(risk), function(co) {
    i <- sample.int(length(risk[[co]]$time), replace = TRUE)
    t <- risk[[co]]$time[i]; s <- risk[[co]]$status[i]
    frozen_reverse_cindex(risk[[co]][[a]][i], t, s) - frozen_reverse_cindex(risk[[co]]$D4[i], t, s)
  }, numeric(1))))
  q <- stats::quantile(boot, c(0.025, 0.975), names = FALSE)
  data.frame(cohort = "MEAN_OF_COHORTS", arm = a, n = sum(sapply(risk, function(x) length(x$time))),
             C_arm = mean(per$C_arm[per$arm == a]), C_D4 = mean(per$C_D4[per$arm == a]),
             dC = point, dC_lower = q[1], dC_upper = q[2])
}
concat <- function(a) {
  r1 <- unlist(lapply(risk, `[[`, a)); r0 <- unlist(lapply(risk, `[[`, "D4"))
  t <- unlist(lapply(risk, `[[`, "time")); s <- unlist(lapply(risk, `[[`, "status"))
  ci <- bootstrap_concordance_diff_ci(r1, r0, t, s, B = 2000, seed = 1, flip_a = FALSE, flip_b = FALSE)
  data.frame(cohort = "CONCATENATED", arm = a, n = length(t), C_arm = frozen_reverse_cindex(r1, t, s),
             C_D4 = frozen_reverse_cindex(r0, t, s), dC = ci$estimate, dC_lower = ci$lower,
             dC_upper = ci$upper)
}
res <- rbind(per, do.call(rbind, lapply(c("supervised", "alpha0_pl"), strat_mean_diff)),
             do.call(rbind, lapply(c("supervised", "alpha0_pl"), concat)))
res$criterion_pass <- ifelse(res$cohort == "MEAN_OF_COHORTS", res$dC_lower > -0.02, NA)
write.csv(res, file.path(OUT_DIR, "delta_c.csv"), row.names = FALSE)
print(res, digits = 3, row.names = FALSE)

# 5. Programs: survival-active (|E beta| >= 1.96 sd) and variance shares ----
prog <- do.call(rbind, lapply(names(fits), function(a) {
  f <- fits[[a]]; sdb <- sqrt(pmax(f$EBeta2 - f$EBeta^2, 0))
  pve <- vapply(seq_len(K), function(k) sum(f$EL[, k]^2) * sum(f$EF[, k]^2), numeric(1)) / sum(Y_train^2)
  data.frame(arm = a, program = seq_len(K), EBeta = f$EBeta, beta_sd = sdb,
             survival_active = abs(f$EBeta) >= 1.96 * sdb & sdb > 0,
             alive = colSums(f$EF^2) > 1e-8, share_rank1 = pve)
}))
write.csv(prog, file.path(OUT_DIR, "programs.csv"), row.names = FALSE)
print(prog, digits = 3, row.names = FALSE)
