# ============================================================
# Script:  run_K_cv.R
# Purpose: Run select_K_cv() on TCGA_PAAD (tcga_only) for both
#          point_normal and normal priors to identify the optimal K
#          via 5-fold CV C-index. K_grid = {2,...,10,15,20}.
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
source("code/train_test_split.R")
source("code/select_K.R")
source("code/preprocess_desurv.R")
source("results/benchmark_sim/benchmark_helpers.R")

PDAC_DATA_ROOT <- Sys.getenv("PDAC_DATA_ROOT",
                              unset = file.path(path.expand("~"),
                                "OneDrive - University of North Carolina at Chapel Hill",
                                "UNC Dissertation (Liu)", "PDAC_data"))

cat("============================================================\n")
cat(" K-CV Sweep — LB model, TCGA_PAAD (tcga_only)\n")
cat(" Priors: point_normal vs normal\n")
cat("============================================================\n\n")

# --- settings ---
K_GRID   <- c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 15L, 20L)
N_FOLDS  <- 5L
ALPHA    <- 0.50
LAMBDA   <- cfg$cavi$lambda
MAX_ITER <- cfg$cavi$max_iter
TOL      <- cfg$cavi$tol
PRIOR_LF <- "point_exponential"
PRIORS   <- c("point_normal", "normal")

cat(sprintf("  K_grid : %s\n", paste(K_GRID, collapse = ", ")))
cat(sprintf("  n_folds: %d | alpha=%.2f | lambda=%.2f | max_iter=%d\n\n",
            N_FOLDS, ALPHA, LAMBDA, MAX_ITER))

# --- load and preprocess TCGA_PAAD ---
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

# --- run K-CV for each prior ---
out_dir <- "results/benchmark_sim/outputs/K_cv"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

results_list <- list()

for (prior_b in PRIORS) {
  cat(sprintf("============================================================\n"))
  cat(sprintf("  Prior: %s\n", prior_b))
  cat(sprintf("  Total fits: %d K × %d folds = %d\n\n",
              length(K_GRID), N_FOLDS, length(K_GRID) * N_FOLDS))

  t0 <- proc.time()
  kcv <- select_K_cv(
    Y_train, time_train, status_train,
    K_grid     = K_GRID,
    n_folds    = N_FOLDS,
    use_1se    = TRUE,
    seed       = 42L,
    verbose    = TRUE,
    alpha      = ALPHA,
    lambda     = LAMBDA,
    max_iter   = MAX_ITER,
    tol        = TOL,
    prior_LF   = PRIOR_LF,
    prior_beta = prior_b
  )
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
    ifelse(tbl$K == best_K,               "    (best mean C)",
                                           "")))
  tbl$mean_cindex <- round(tbl$mean_cindex, 4)
  tbl$se_cindex   <- round(tbl$se_cindex,   4)
  print(tbl, row.names = FALSE)
  cat("\n")

  # save per-prior CSVs
  cv_tbl_out <- kcv$cv_table
  cv_tbl_out$prior_beta <- prior_b
  fold_out <- kcv$fold_results
  fold_out$prior_beta <- prior_b

  write.csv(cv_tbl_out,
            file.path(out_dir, sprintf("K_cv_table_%s.csv", prior_b)),
            row.names = FALSE)
  write.csv(fold_out,
            file.path(out_dir, sprintf("K_cv_folds_%s.csv", prior_b)),
            row.names = FALSE)

  results_list[[prior_b]] <- list(kcv = kcv, elapsed = elapsed)
}

# --- side-by-side comparison ---
cat("============================================================\n")
cat(" Side-by-side comparison\n")
cat("============================================================\n\n")

pn <- results_list[["point_normal"]]$kcv$cv_table
nm <- results_list[["normal"]]$kcv$cv_table

cmp <- data.frame(
  K              = pn$K,
  PN_mean_C      = round(pn$mean_cindex, 4),
  PN_se          = round(pn$se_cindex,   4),
  NM_mean_C      = round(nm$mean_cindex, 4),
  NM_se          = round(nm$se_cindex,   4),
  stringsAsFactors = FALSE
)
# flag selected K for each prior
cmp$selected <- paste0(
  ifelse(cmp$K == results_list[["point_normal"]]$kcv$K_opt, "PN ", ""),
  ifelse(cmp$K == results_list[["normal"]]$kcv$K_opt,       "NM",  "")
)
print(cmp, row.names = FALSE)

cat(sprintf("\n  point_normal K_opt = %d\n",
            results_list[["point_normal"]]$kcv$K_opt))
cat(sprintf("  normal       K_opt = %d\n\n",
            results_list[["normal"]]$kcv$K_opt))

# save combined table
write.csv(cmp, file.path(out_dir, "K_cv_comparison.csv"), row.names = FALSE)

meta <- data.frame(
  prior_beta     = PRIORS,
  K_opt          = sapply(PRIORS, function(p) results_list[[p]]$kcv$K_opt),
  selection_rule = sapply(PRIORS, function(p) results_list[[p]]$kcv$selection_rule),
  elapsed_min    = sapply(PRIORS, function(p) round(results_list[[p]]$elapsed / 60, 2)),
  K_grid         = paste(K_GRID, collapse = ","),
  n_folds        = N_FOLDS,
  alpha          = ALPHA,
  run_date       = as.character(Sys.Date()),
  stringsAsFactors = FALSE
)
write.csv(meta, file.path(out_dir, "K_cv_meta.csv"), row.names = FALSE)

cat(sprintf("  Results saved to %s/\n", out_dir))
cat("============================================================\n")
