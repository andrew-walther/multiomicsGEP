# ============================================================
# Script:  run_K_cv.R
# Purpose: Run select_K_cv() on TCGA_PAAD (tcga_only) for both
#          LB (eta = L*beta) and YFB (eta = YF*beta) model structures,
#          each with point_normal and normal priors.
#          K_grid = {2,...,10,15,20}, 5-fold stratified CV, 1-SE rule.
#          Alpha fixed at 0.50 (CV-selected in full LB benchmark).
#          Results are printed side-by-side and saved to CSV.
#
# Usage (from project root):
#   Rscript results/benchmark_sim/run_K_cv.R
#   PDAC_DATA_ROOT=/path/to/data Rscript results/benchmark_sim/run_K_cv.R
#
# Output: results/benchmark_sim/outputs/K_cv/
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-05-06
# Dependencies: code/fit_modular.R, code/predict.R,
#               code/fit_cox_on_yf.R, code/predict_cox_on_yf.R,
#               code/train_test_split.R, code/select_K.R,
#               code/preprocess_desurv.R,
#               results/benchmark_sim/benchmark_helpers.R
# ============================================================

suppressPackageStartupMessages({
  library(yaml); library(ebnm); library(survival)
})

# cfg must be loaded before benchmark_helpers.R
cfg <- yaml::read_yaml("config/globals.yml")

source("code/update_beta.R"); source("code/update_L.R")
source("code/update_F.R");    source("code/update_tau.R")
source("code/compute_elbo.R")
suppressMessages(tryCatch(
  source("code/fit_modular.R"),
  error = function(e) invisible(NULL)
))
source("code/predict.R")
tryCatch(source("code/fit_cox_on_yf.R"),    error = function(e) invisible(NULL))
tryCatch(source("code/predict_cox_on_yf.R"), error = function(e) invisible(NULL))
source("code/train_test_split.R")
source("code/select_K.R")
source("code/preprocess_desurv.R")
source("results/benchmark_sim/benchmark_helpers.R")

PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT",
                              unset = file.path(path.expand("~"),
                                "OneDrive - University of North Carolina at Chapel Hill",
                                "UNC Dissertation (Liu)", "PDAC_data"))

# --- global settings ---
K_GRID   <- c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 15L, 20L)
N_FOLDS  <- 5L
ALPHA    <- 0.50
MAX_ITER <- cfg$cavi$max_iter
TOL      <- cfg$cavi$tol
PRIOR_LF <- "point_exponential"
PRIORS   <- c("point_normal", "normal")

out_dir <- "results/benchmark_sim/outputs/K_cv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --- load and preprocess TCGA_PAAD (shared across LB and YFB K-CV) ---
if (!dir.exists(PDAC_DATA_ROOT))
  stop(sprintf("PDAC_DATA_ROOT not found: %s", PDAC_DATA_ROOT))

cat("  Loading TCGA_PAAD ...\n")
raw <- load_pdac_raw("TCGA_PAAD", PDAC_DATA_ROOT)
pre <- preprocess_desurv_cohort(
  Y             = raw$Y,
  gene_names    = raw$gene_names,
  top_n         = cfg$preprocessing$top_n_genes,
  log_transform = TRUE,
  cohort_name   = "TCGA_PAAD"
)
Y_train      <- pre$Y
time_train   <- raw$time
status_train <- raw$status
cat(sprintf("  Training matrix: n=%d, p=%d | events=%d (%.0f%%)\n\n",
            nrow(Y_train), ncol(Y_train),
            sum(status_train), 100 * mean(status_train)))

# ============================================================
# Helper: run K-CV for one model × prior, print, save, return
# ============================================================
run_one_kcv <- function(model_name, prior_b, label) {
  cat(sprintf("============================================================\n"))
  cat(sprintf("  Model: %s | Prior: %s\n", model_name, prior_b))
  cat(sprintf("  Total fits: %d K × %d folds = %d\n\n",
              length(K_GRID), N_FOLDS, length(K_GRID) * N_FOLDS))

  t0 <- proc.time()

  if (model_name == "LB") {
    kcv <- select_K_cv(
      Y_train, time_train, status_train,
      K_grid     = K_GRID,
      n_folds    = N_FOLDS,
      use_1se    = TRUE,
      seed       = 42L,
      verbose    = TRUE,
      model      = "LB",
      alpha      = ALPHA,
      max_iter   = MAX_ITER,
      tol        = TOL,
      prior_LF   = PRIOR_LF,
      prior_beta = prior_b
    )
  } else {
    kcv <- select_K_cv(
      Y_train, time_train, status_train,
      K_grid     = K_GRID,
      n_folds    = N_FOLDS,
      use_1se    = TRUE,
      seed       = 42L,
      verbose    = TRUE,
      model      = "YFB",
      alpha      = ALPHA,
      max_iter   = MAX_ITER,
      tol        = TOL,
      prior_LF   = PRIOR_LF,
      prior_beta = prior_b
    )
  }
  elapsed <- (proc.time() - t0)["elapsed"]

  cat(sprintf("\n  Elapsed: %.1f min\n", elapsed / 60))
  cat(sprintf("  Selected K: %d  (rule: %s)\n\n",
              kcv$K_opt, kcv$selection_rule))

  # annotated table
  tbl <- kcv$cv_table
  best_K <- tbl$K[which.max(tbl$mean_cindex)]
  tbl$note <- ifelse(
    tbl$K == kcv$K_opt & tbl$K == best_K, "<-- selected (= best)",
    ifelse(tbl$K == kcv$K_opt,            "<-- selected (1-SE)",
    ifelse(tbl$K == best_K,               "    (best mean C)", "")))
  tbl$mean_cindex <- round(tbl$mean_cindex, 4)
  tbl$se_cindex   <- round(tbl$se_cindex,   4)
  print(tbl, row.names = FALSE)
  cat("\n")

  # save per-model-prior CSVs
  cv_out   <- kcv$cv_table; cv_out$model   <- model_name; cv_out$prior_beta   <- prior_b
  fold_out <- kcv$fold_results; fold_out$model <- model_name; fold_out$prior_beta <- prior_b
  write.csv(cv_out,
            file.path(out_dir, sprintf("K_cv_table_%s_%s.csv", model_name, prior_b)),
            row.names = FALSE)
  write.csv(fold_out,
            file.path(out_dir, sprintf("K_cv_folds_%s_%s.csv", model_name, prior_b)),
            row.names = FALSE)

  list(kcv = kcv, elapsed = elapsed)
}

# ============================================================
# LB K-CV
# ============================================================
cat("============================================================\n")
cat(" K-CV Sweep — LB model, TCGA_PAAD (tcga_only)\n")
cat(sprintf("  K_grid : %s\n", paste(K_GRID, collapse = ", ")))
cat(sprintf("  n_folds: %d | alpha=%.2f | max_iter=%d\n\n",
            N_FOLDS, ALPHA, MAX_ITER))

lb_results <- list()
for (prior_b in PRIORS)
  lb_results[[prior_b]] <- run_one_kcv("LB", prior_b, prior_b)

# LB side-by-side
cat("============================================================\n")
cat(" LB — side-by-side prior comparison\n")
cat("============================================================\n\n")
pn_lb <- lb_results[["point_normal"]]$kcv$cv_table
nm_lb <- lb_results[["normal"]]$kcv$cv_table
cmp_lb <- data.frame(
  K         = pn_lb$K,
  PN_mean_C = round(pn_lb$mean_cindex, 4),  PN_se = round(pn_lb$se_cindex, 4),
  NM_mean_C = round(nm_lb$mean_cindex, 4),  NM_se = round(nm_lb$se_cindex, 4),
  selected  = paste0(
    ifelse(pn_lb$K == lb_results[["point_normal"]]$kcv$K_opt, "PN ", ""),
    ifelse(nm_lb$K == lb_results[["normal"]]$kcv$K_opt,       "NM",  "")),
  stringsAsFactors = FALSE
)
print(cmp_lb, row.names = FALSE)
cat(sprintf("\n  LB point_normal K_opt = %d\n", lb_results[["point_normal"]]$kcv$K_opt))
cat(sprintf("  LB normal       K_opt = %d\n\n", lb_results[["normal"]]$kcv$K_opt))
write.csv(cmp_lb, file.path(out_dir, "K_cv_comparison_LB.csv"), row.names = FALSE)

# ============================================================
# YFB K-CV
# ============================================================
cat("============================================================\n")
cat(" K-CV Sweep — YFB model, TCGA_PAAD (tcga_only)\n")
cat(sprintf("  K_grid : %s\n", paste(K_GRID, collapse = ", ")))
cat(sprintf("  n_folds: %d | alpha=%.2f | max_iter=%d\n\n",
            N_FOLDS, ALPHA, MAX_ITER))

yfb_results <- list()
for (prior_b in PRIORS)
  yfb_results[[prior_b]] <- run_one_kcv("YFB", prior_b, prior_b)

# YFB side-by-side
cat("============================================================\n")
cat(" YFB — side-by-side prior comparison\n")
cat("============================================================\n\n")
pn_yfb <- yfb_results[["point_normal"]]$kcv$cv_table
nm_yfb <- yfb_results[["normal"]]$kcv$cv_table
cmp_yfb <- data.frame(
  K         = pn_yfb$K,
  PN_mean_C = round(pn_yfb$mean_cindex, 4), PN_se = round(pn_yfb$se_cindex, 4),
  NM_mean_C = round(nm_yfb$mean_cindex, 4), NM_se = round(nm_yfb$se_cindex, 4),
  selected  = paste0(
    ifelse(pn_yfb$K == yfb_results[["point_normal"]]$kcv$K_opt, "PN ", ""),
    ifelse(nm_yfb$K == yfb_results[["normal"]]$kcv$K_opt,        "NM",  "")),
  stringsAsFactors = FALSE
)
print(cmp_yfb, row.names = FALSE)
cat(sprintf("\n  YFB point_normal K_opt = %d\n", yfb_results[["point_normal"]]$kcv$K_opt))
cat(sprintf("  YFB normal       K_opt = %d\n\n", yfb_results[["normal"]]$kcv$K_opt))
write.csv(cmp_yfb, file.path(out_dir, "K_cv_comparison_YFB.csv"), row.names = FALSE)

# ============================================================
# LB vs YFB head-to-head (normal prior — the working prior for both)
# ============================================================
cat("============================================================\n")
cat(" LB vs YFB — head-to-head (normal prior)\n")
cat("============================================================\n\n")
cmp_h2h <- data.frame(
  K            = nm_lb$K,
  LB_mean_C    = round(nm_lb$mean_cindex, 4),  LB_se  = round(nm_lb$se_cindex, 4),
  YFB_mean_C   = round(nm_yfb$mean_cindex, 4), YFB_se = round(nm_yfb$se_cindex, 4),
  LB_selected  = nm_lb$K == lb_results[["normal"]]$kcv$K_opt,
  YFB_selected = nm_yfb$K == yfb_results[["normal"]]$kcv$K_opt,
  stringsAsFactors = FALSE
)
print(cmp_h2h, row.names = FALSE)
write.csv(cmp_h2h, file.path(out_dir, "K_cv_comparison_h2h.csv"), row.names = FALSE)

# ============================================================
# Meta / combined CSVs (backwards-compatible with old K_cv_comparison.csv)
# ============================================================
# Overwrite old LB-only comparison with new LB comparison for backwards compat
write.csv(cmp_lb, file.path(out_dir, "K_cv_comparison.csv"), row.names = FALSE)

meta_rows <- lapply(list(
  list(model="LB",  prior="point_normal", res=lb_results[["point_normal"]]),
  list(model="LB",  prior="normal",       res=lb_results[["normal"]]),
  list(model="YFB", prior="point_normal", res=yfb_results[["point_normal"]]),
  list(model="YFB", prior="normal",       res=yfb_results[["normal"]])
), function(x) {
  data.frame(
    model          = x$model,
    prior_beta     = x$prior,
    K_opt          = x$res$kcv$K_opt,
    selection_rule = x$res$kcv$selection_rule,
    elapsed_min    = round(x$res$elapsed / 60, 2),
    K_grid         = paste(K_GRID, collapse = ","),
    n_folds        = N_FOLDS,
    alpha          = ALPHA,
    run_date       = as.character(Sys.Date()),
    stringsAsFactors = FALSE
  )
})
meta <- do.call(rbind, meta_rows)
write.csv(meta, file.path(out_dir, "K_cv_meta.csv"), row.names = FALSE)

cat(sprintf("  Results saved to %s/\n", out_dir))
cat("============================================================\n")
