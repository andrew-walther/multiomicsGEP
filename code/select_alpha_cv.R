# ============================================================
# Script: select_alpha_cv.R
# Purpose: Cross-validated alpha selection for the SBMF model.
#          Evaluates alpha values by held-out C-index and supports
#          either max-C-index or 1-SE rule selection.
# Author: Codex
# Created: 2026-04-23
# Dependencies: code/fit_modular.R, code/predict.R, code/train_test_split.R
#               must be sourced first.
# ============================================================

# ============================================================
# select_alpha_cv() ----
# ============================================================

#' Select alpha via stratified K-fold cross-validation.
#'
#' For each alpha in `alpha_grid`, fits the SBMF model on `n_folds - 1` folds,
#' projects the held-out fold with `predict_supervised_mf()`, and evaluates
#' held-out C-index using the project convention:
#'
#'   concordance(Surv(time, status) ~ I(-risk_scores))
#'
#' When `use_1se = TRUE`, applies the 1-SE rule by selecting the largest alpha
#' whose mean held-out C-index is at least `max_mean - se_at_best`.
#'
#' @param Y          numeric matrix (n x p)
#' @param time       numeric vector length n
#' @param status     integer/numeric vector length n with entries in {0, 1}
#' @param alpha_grid numeric vector of candidate alpha values in [0, 1]
#' @param n_folds    integer >= 2, default 5
#' @param K_max      integer >= 1: number of latent factors to fit in each fold
#' @param use_1se    logical: use 1-SE rule (default TRUE)
#' @param seed       integer random seed used for fold creation
#' @param ...        additional arguments passed to fit_supervised_mf_modular()
#'
#' @return Named list:
#'   $alpha_opt      selected alpha
#'   $cv_table       data.frame with alpha, mean_cindex, se_cindex, n_folds
#'   $fold_results   data.frame with per-alpha, per-fold C-index
#'   $selection_rule character: "1se" or "max"
#'
#' @examples
#' # Requires fit_modular.R, predict.R, and train_test_split.R to be sourced first.
#' # res <- select_alpha_cv(Y, time, status, alpha_grid = c(0.1, 0.5, 0.9), n_folds = 3, K_max = 3)
select_alpha_cv <- function(Y, time, status,
                            alpha_grid = seq(0.1, 0.9, by = 0.1),
                            n_folds    = 5,
                            K_max      = 10,
                            use_1se    = TRUE,
                            seed       = 42,
                            ...) {

  if (!exists("fit_supervised_mf_modular", mode = "function"))
    stop("fit_supervised_mf_modular() must be sourced before select_alpha_cv().")
  if (!exists("predict_supervised_mf", mode = "function"))
    stop("predict_supervised_mf() must be sourced before select_alpha_cv().")
  if (!exists("create_stratified_folds", mode = "function"))
    stop("create_stratified_folds() must be sourced before select_alpha_cv().")

  if (!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix.")
  n <- nrow(Y)
  if (length(time) != n)
    stop(sprintf("time must have length %d to match nrow(Y).", n))
  if (length(status) != n)
    stop(sprintf("status must have length %d to match nrow(Y).", n))
  if (!all(status %in% c(0, 1)))
    stop("status must contain only 0 and 1.")
  if (length(K_max) != 1 || !is.finite(K_max) || K_max < 1 || K_max != as.integer(K_max))
    stop("K_max must be an integer >= 1.")
  if (length(n_folds) != 1 || !is.finite(n_folds) || n_folds < 2 || n_folds != as.integer(n_folds))
    stop("n_folds must be an integer >= 2.")
  if (!is.numeric(alpha_grid) || length(alpha_grid) < 1)
    stop("alpha_grid must be a non-empty numeric vector.")
  if (any(!is.finite(alpha_grid)) || any(alpha_grid < 0) || any(alpha_grid > 1))
    stop("alpha_grid must contain only finite values in [0, 1].")

  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  n_folds    <- as.integer(n_folds)
  K_max      <- as.integer(K_max)

  fold_obj <- create_stratified_folds(status, n_folds = n_folds, seed = seed)

  fold_rows <- vector("list", length(alpha_grid) * n_folds)
  row_idx <- 1L

  for (alpha in alpha_grid) {
    for (fold_id in seq_len(n_folds)) {
      test_idx  <- fold_obj$folds[[fold_id]]
      train_idx <- setdiff(seq_len(n), test_idx)

      fit <- fit_supervised_mf_modular(
        Y[train_idx, , drop = FALSE],
        time[train_idx],
        status[train_idx],
        K = K_max,
        alpha = alpha,
        verbose = FALSE,
        ...
      )

      pred <- predict_supervised_mf(Y[test_idx, , drop = FALSE], fit$EF, fit$EBeta)
      cindex <- survival::concordance(
        survival::Surv(time[test_idx], status[test_idx]) ~ I(-pred$risk_scores)
      )$concordance

      fold_rows[[row_idx]] <- data.frame(
        alpha = alpha,
        fold = fold_id,
        n_train = length(train_idx),
        n_test = length(test_idx),
        n_event_test = sum(status[test_idx] == 1),
        n_censored_test = sum(status[test_idx] == 0),
        cindex = as.numeric(cindex)
      )
      row_idx <- row_idx + 1L
    }
  }

  fold_results <- do.call(rbind, fold_rows)

  cv_table <- do.call(
    rbind,
    lapply(alpha_grid, function(alpha) {
      cidx <- fold_results$cindex[fold_results$alpha == alpha]
      data.frame(
        alpha = alpha,
        mean_cindex = mean(cidx),
        se_cindex = stats::sd(cidx) / sqrt(length(cidx)),
        n_folds = length(cidx)
      )
    })
  )
  rownames(cv_table) <- NULL

  best_idx <- which.max(cv_table$mean_cindex)
  if (use_1se) {
    threshold <- cv_table$mean_cindex[best_idx] - cv_table$se_cindex[best_idx]
    eligible  <- cv_table$alpha[cv_table$mean_cindex >= threshold]
    alpha_opt <- max(eligible)
    selection_rule <- "1se"
  } else {
    alpha_opt <- max(cv_table$alpha[cv_table$mean_cindex == cv_table$mean_cindex[best_idx]])
    selection_rule <- "max"
  }

  list(
    alpha_opt = alpha_opt,
    cv_table = cv_table,
    fold_results = fold_results,
    selection_rule = selection_rule
  )
}
