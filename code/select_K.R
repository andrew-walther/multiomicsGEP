# ============================================================
# Script: select_K.R
# Purpose: Data-driven selection of the number of latent factors K
#          for the SBMF model.
#          Option A: fit large K and count active (non-pruned) factors.
#          Option B: cross-validated C-index (stub for Longleaf HPC).
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-04-01
# Dependencies: code/fit_modular.R (must be sourced first)
# ============================================================

# ============================================================
# compute_pve() ----
# ============================================================

#' Compute per-factor proportion of variance explained (PVE).
#'
#' PVE_k = ||EL[,k] %*% t(EF[,k])||_F^2 / ||Y||_F^2
#'
#' Note: this uses the posterior means EL and EF, not the full
#' posterior (which would require second moments and a larger calculation).
#' It is a useful diagnostic but not identical to the model ELBO.
#'
#' @param res     list with $EL (n × K), $EF (p × K) — from fit_supervised_mf_modular()
#' @param Y       numeric matrix (n × p) — the original data
#'
#' @return numeric vector of length K: PVE per factor (values in [0, 1])
compute_pve <- function(res, Y) {
  K         <- ncol(res$EL)
  total_var <- sum(Y^2)
  if (total_var == 0) return(rep(0, K))
  sapply(seq_len(K), function(k) {
    sum(outer(res$EL[, k], res$EF[, k])^2) / total_var
  })
}

# ============================================================
# auto_prune_K() ----
# ============================================================

#' Select K by fitting a large model and counting active factors.
#'
#' Fits the SBMF model with K = K_max latent factors.  The point-normal
#' prior drives uninformative factors toward exactly zero, producing
#' sparse β and F columns.  An "active" factor is one where either:
#'   - |β_k| > beta_thresh (the factor has prognostic signal), OR
#'   - PVE_k > pve_thresh (the factor explains >1% of genomic variance)
#'
#' The count of active factors is the recommended effective K.  This is
#' analogous to the EBMF/flashr approach to rank selection.
#'
#' **Computational cost:** One fit at K_max instead of one fit per K
#' candidate — much cheaper than cross-validation.
#'
#' @param Y          numeric matrix (n × p)
#' @param time       numeric vector (n)
#' @param status     integer vector (n)
#' @param K_max      integer: maximum K to fit (default 10)
#' @param beta_thresh  numeric: |β_k| threshold for "active" (default 0.05)
#' @param pve_thresh   numeric: PVE_k threshold for "active" (default 0.01 = 1%)
#' @param ...        additional arguments passed to fit_supervised_mf_modular()
#'
#' @return Named list:
#'   $K_effective   integer: number of active factors
#'   $pve           numeric vector (K_max): PVE per factor
#'   $beta          numeric vector (K_max): |EBeta| per factor
#'   $active        logical vector (K_max): TRUE if factor is active
#'   $fit           full fit object from fit_supervised_mf_modular()
auto_prune_K <- function(Y, time, status, K_max = 10,
                          beta_thresh = 0.05, pve_thresh = 0.01, ...) {

  cat(sprintf("  [auto_prune_K] Fitting K=%d to identify active factors...\n", K_max))

  fit <- fit_supervised_mf_modular(Y, time, status, K = K_max, ...)

  pve     <- compute_pve(fit, Y)
  ab_beta <- abs(fit$EBeta)

  # A factor is "active" if it either has prognostic signal (|β| > thresh)
  # or explains non-trivial genomic variance (PVE > 1%)
  active <- (ab_beta > beta_thresh) | (pve > pve_thresh)

  K_eff <- sum(active)
  cat(sprintf("  [auto_prune_K] K_effective = %d / %d (beta_thresh=%.2f, pve_thresh=%.3f)\n",
              K_eff, K_max, beta_thresh, pve_thresh))
  cat(sprintf("  [auto_prune_K] |beta|: [%s]\n",
              paste(sprintf("%.3f", ab_beta), collapse = ", ")))
  cat(sprintf("  [auto_prune_K] PVE%%:   [%s]\n",
              paste(sprintf("%.1f%%", pve * 100), collapse = ", ")))

  list(
    K_effective = K_eff,
    pve         = pve,
    beta        = ab_beta,
    active      = active,
    fit         = fit
  )
}

# ============================================================
# select_K_cv() ----
# ============================================================

#' K selection via cross-validated held-out C-index with the 1-SE rule.
#'
#' For each candidate K in K_grid, runs n_folds-fold stratified CV:
#'   1. Fit the model on the (n_folds-1)/n_folds training portion.
#'   2. Project the held-out fold (LB: predict_supervised_mf;
#'      YFB: predict_cox_on_yf).
#'   3. Compute held-out C-index: concordance(Surv(time, status) ~ I(-risk)).
#'
#' After collecting all fold-level C-indices, summarises by mean ± SE per K.
#' The 1-SE rule (use_1se = TRUE) selects the **smallest K** whose mean C-index
#' is within one SE of the maximum — preferring parsimony when evidence for a
#' larger K is weak.  Set use_1se = FALSE to select the maximising K directly.
#'
#' **Computational cost:** n_folds × |K_grid| full fits.
#' For K_grid = 2:10 + 15 + 20 (11 values) and n_folds = 5, that is 55 fits.
#' Each fit runs up to max_iter iterations on (n_folds-1)/n_folds of n subjects.
#' On Longleaf HPC, parallelise over K values with SLURM array jobs; locally,
#' expect ~15–30 min for n ~ 150, p = 2000, K_max = 20.
#'
#' sign_correction is always disabled inside CV folds (same rationale as in
#' select_alpha_cv): fold-level sign correction produces inconsistent signs
#' across folds and inflates apparent C-index variance.
#'
#' @param Y          numeric matrix (n × p)
#' @param time       numeric vector length n: survival/censoring times
#' @param status     integer vector length n: event indicators (1=event, 0=censored)
#' @param K_grid     integer vector: candidate K values
#'                   (default c(2,3,4,5,6,7,8,9,10,15,20))
#' @param n_folds    integer >= 2: number of CV folds (default 5)
#' @param use_1se    logical: apply 1-SE rule (default TRUE)
#' @param seed       integer: RNG seed for stratified fold creation (default 42)
#' @param verbose    logical: print per-fold progress (default FALSE)
#' @param model      character: "LB" (default, calls fit_supervised_mf_modular +
#'                   predict_supervised_mf) or "YFB" (calls fit_cox_on_yf +
#'                   predict_cox_on_yf).  Both models must be sourced before
#'                   calling; see dependency checks below.
#' @param ...        additional arguments passed to the fitting function.
#'                   Do NOT pass K (overridden per candidate) or sign_correction
#'                   (disabled inside folds for both models).
#'
#' @return Named list:
#'   $K_opt         integer: selected K (1-SE rule or max)
#'   $cv_table      data.frame with K, mean_cindex, se_cindex, n_folds
#'   $fold_results  data.frame with K, fold, n_train, n_test, n_event_test,
#'                  n_censored_test, cindex — one row per (K, fold) combination
#'   $selection_rule character: "1se" or "max"
#'   $model         character: "LB" or "YFB" (echoed for caller bookkeeping)
#'
#' @examples
#' \dontrun{
#' # LB model
#' source("code/fit_modular.R"); source("code/predict.R")
#' source("code/train_test_split.R"); source("code/select_K.R")
#' res_lb  <- select_K_cv(Y, time, status, model = "LB",
#'                         K_grid = c(2,3,5,8,10), n_folds = 3)
#' # YFB model
#' source("code/fit_cox_on_yf.R"); source("code/predict_cox_on_yf.R")
#' res_yfb <- select_K_cv(Y, time, status, model = "YFB",
#'                         K_grid = c(2,3,5,8,10), n_folds = 3)
#' }
#'
#' @seealso \code{\link{auto_prune_K}} for a cheaper single-fit alternative.
select_K_cv <- function(Y, time, status,
                         K_grid   = c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 15L, 20L),
                         n_folds  = 5L,
                         use_1se  = TRUE,
                         seed     = 42L,
                         verbose  = FALSE,
                         model    = "LB",
                         ...) {

  # --- model validation ---
  model <- match.arg(model, c("LB", "YFB"))

  # --- dependency checks (branch on model) ---
  if (model == "LB") {
    if (!exists("fit_supervised_mf_modular", mode = "function"))
      stop("fit_supervised_mf_modular() must be sourced before select_K_cv(model='LB').")
    if (!exists("predict_supervised_mf", mode = "function"))
      stop("predict_supervised_mf() must be sourced before select_K_cv(model='LB').")
  } else {
    if (!exists("fit_cox_on_yf", mode = "function"))
      stop("fit_cox_on_yf() must be sourced before select_K_cv(model='YFB').")
    if (!exists("predict_cox_on_yf", mode = "function"))
      stop("predict_cox_on_yf() must be sourced before select_K_cv(model='YFB').")
  }
  if (!exists("create_stratified_folds", mode = "function"))
    stop("create_stratified_folds() must be sourced before select_K_cv().")

  # --- input validation ---
  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix.")
  n <- nrow(Y)
  if (length(time) != n)
    stop(sprintf("time must have length %d (nrow(Y)).", n))
  if (length(status) != n)
    stop(sprintf("status must have length %d (nrow(Y)).", n))
  if (!all(status %in% c(0, 1)))
    stop("status must contain only 0 and 1.")
  K_grid <- as.integer(sort(unique(K_grid)))
  if (any(K_grid < 1))
    stop("K_grid values must be positive integers.")
  n_folds <- as.integer(n_folds)
  if (n_folds < 2)
    stop("n_folds must be >= 2.")

  # Strip K (overridden per candidate) and sign_correction (disabled in folds
  # for LB; not a formal arg in fit_cox_on_yf, so must be stripped for YFB too).
  extra <- list(...)
  extra[["K"]]               <- NULL
  extra[["sign_correction"]] <- NULL

  # --- create stratified folds once (shared across all K) ---
  fold_obj <- create_stratified_folds(status, n_folds = n_folds, seed = seed)

  n_combos  <- length(K_grid) * n_folds
  fold_rows <- vector("list", n_combos)
  row_idx   <- 1L

  for (K in K_grid) {
    if (verbose)
      cat(sprintf("  [select_K_cv] K = %d ...\n", K))

    for (fold_id in seq_len(n_folds)) {
      test_idx  <- fold_obj$folds[[fold_id]]
      train_idx <- setdiff(seq_len(n), test_idx)

      # --- fit on training fold ---
      if (model == "LB") {
        # LB: sign_correction=FALSE hardcoded — consistent fold signs
        fit_args <- c(
          list(Y      = Y[train_idx, , drop = FALSE],
               time   = time[train_idx],
               status = status[train_idx],
               K      = K,
               sign_correction = FALSE,
               verbose = FALSE),
          extra
        )
        fit <- do.call(fit_supervised_mf_modular, fit_args)

        # Project held-out fold
        pred <- predict_supervised_mf(Y[test_idx, , drop = FALSE],
                                      fit$EF, fit$EBeta)
      } else {
        # YFB: fit_cox_on_yf does not have sign_correction; verbose=FALSE
        fit_args <- c(
          list(Y      = Y[train_idx, , drop = FALSE],
               time   = time[train_idx],
               status = status[train_idx],
               K      = K,
               verbose = FALSE),
          extra
        )
        fit <- do.call(fit_cox_on_yf, fit_args)

        # Project held-out fold; pass EF_norms so column scaling is consistent
        pred <- predict_cox_on_yf(Y[test_idx, , drop = FALSE],
                                  fit$EF, fit$EBeta,
                                  EF_norms = fit$EF_norms)
      }

      cindex <- tryCatch(
        as.numeric(survival::concordance(
          survival::Surv(time[test_idx], status[test_idx]) ~ I(-pred$risk_scores)
        )$concordance),
        error = function(e) NA_real_
      )

      if (verbose)
        cat(sprintf("    fold %d/%d: n_train=%d, n_test=%d, C=%.3f\n",
                    fold_id, n_folds,
                    length(train_idx), length(test_idx),
                    ifelse(is.na(cindex), -1, cindex)))

      fold_rows[[row_idx]] <- data.frame(
        K               = K,
        fold            = fold_id,
        n_train         = length(train_idx),
        n_test          = length(test_idx),
        n_event_test    = sum(status[test_idx] == 1),
        n_censored_test = sum(status[test_idx] == 0),
        cindex          = cindex,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
    }
  }

  fold_results <- do.call(rbind, fold_rows)
  rownames(fold_results) <- NULL

  # --- summarise per K ---
  cv_table <- do.call(rbind, lapply(K_grid, function(K) {
    cidx <- fold_results$cindex[fold_results$K == K]
    cidx_obs <- cidx[!is.na(cidx)]
    data.frame(
      K           = K,
      mean_cindex = if (length(cidx_obs) > 0) mean(cidx_obs) else NA_real_,
      se_cindex   = if (length(cidx_obs) > 1)
                      stats::sd(cidx_obs) / sqrt(length(cidx_obs))
                    else NA_real_,
      n_folds     = length(cidx_obs),
      stringsAsFactors = FALSE
    )
  }))
  rownames(cv_table) <- NULL

  # --- K selection ---
  # Only consider K values with complete fold results
  complete <- !is.na(cv_table$mean_cindex)
  if (!any(complete))
    stop("All K values produced NA C-indices — check data and model arguments.")

  best_idx <- which.max(cv_table$mean_cindex[complete])
  best_K   <- cv_table$K[complete][best_idx]
  best_c   <- cv_table$mean_cindex[complete][best_idx]

  if (use_1se) {
    # SE at the best K (may be NA if only one non-NA fold)
    se_best   <- cv_table$se_cindex[cv_table$K == best_K]
    se_best   <- ifelse(is.na(se_best), 0, se_best)
    threshold <- best_c - se_best
    # Smallest K whose mean C-index >= threshold
    eligible  <- cv_table$K[complete & cv_table$mean_cindex >= threshold]
    K_opt     <- min(eligible)
    selection_rule <- "1se"
  } else {
    K_opt <- best_K
    selection_rule <- "max"
  }

  list(
    K_opt          = K_opt,
    cv_table       = cv_table,
    fold_results   = fold_results,
    selection_rule = selection_rule,
    model          = model
  )
}
