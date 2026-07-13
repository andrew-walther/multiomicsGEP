# ============================================================
# Script:  select_k_alpha_bo.R
# Purpose: Joint (K, alpha) tuning via Bayesian optimization -- Step 3 of
#          docs/plans/ssbmf_k_parsimony_followup_plan_07_13_2026.md. Tests
#          whether jointly tuning K and alpha (rather than fixing alpha=0.5
#          and tuning K alone) can reach a smaller K with comparable external
#          performance -- the DeSurv-aligned comparison sketched in
#          docs/plans/joint_k_alpha_bayesopt_plan_07_12_2026.md.
#
#          Reuses code/select_K.R's select_K_cv() as the per-point CV
#          objective (called with a single-element K_grid), rather than
#          duplicating its fold-fitting logic -- select_K_cv() already
#          threads arbitrary extra arguments (including `alpha`) through to
#          the underlying fit function via `...`.
# Author:  Claude Code (reviewed by Andrew Walther)
# Created: 2026-07-13
# Dependencies: code/select_K.R (select_K_cv), code/train_test_split.R
#               (create_stratified_folds), rBayesianOptimization package
# ============================================================

#' Joint (K, alpha) CV Objective for Bayesian Optimization
#'
#' Scores a single (K, alpha) point by mean cross-validated external
#' concordance, via a single-K call to select_K_cv(). K is rounded to the
#' nearest integer before fitting, since the Gaussian-process optimizer
#' proposes continuous values for all dimensions including K.
#'
#' @param Y, time, status  training data (see select_K_cv())
#' @param K       numeric (rounded to integer before use)
#' @param alpha   numeric in [0, 1]
#' @param n_folds,seed,model,...  passed through to select_K_cv()
#'
#' @return list(Score = mean CV concordance, Pred = 0) -- the return
#'   contract rBayesianOptimization::BayesianOptimization() requires.
select_k_alpha_bo_objective <- function(Y, time, status, K, alpha,
                                         n_folds = 5, seed = 42,
                                         model = "YFB", ...) {
  if (!exists("select_K_cv", mode = "function"))
    stop("select_K_cv() must be sourced before select_k_alpha_bo_objective().")
  K_int <- round(K)
  res <- select_K_cv(Y, time, status,
                     K_grid = K_int, n_folds = n_folds, use_1se = FALSE,
                     seed = seed, verbose = FALSE, model = model,
                     alpha = alpha, ...)
  list(Score = res$cv_table$mean_cindex[1], Pred = 0)
}

#' Joint (K, alpha) Selection via Bayesian Optimization
#'
#' Runs rBayesianOptimization::BayesianOptimization() over the joint (K,
#' alpha) space, using select_k_alpha_bo_objective() (mean CV concordance
#' via select_K_cv()) as the objective. K is treated as continuous by the
#' underlying Gaussian process and rounded to the nearest integer inside the
#' objective; init_points random draws seed the GP before n_iter acquisition-
#' guided evaluations follow.
#'
#' @param Y, time, status  training data
#' @param K_bounds      numeric length-2 vector: c(lower, upper) for K
#'                      (integers; lower < upper).
#' @param alpha_bounds  numeric length-2 vector: c(lower, upper) for alpha,
#'                      both within [0, 1] (lower < upper).
#' @param n_folds       integer: CV folds per evaluation (default 5).
#' @param seed          integer: fold-assignment seed, held fixed across all
#'                      evaluations so every (K, alpha) point is scored on
#'                      the same folds (default 42).
#' @param model         "LB" or "YFB" (default "YFB").
#' @param init_points   integer: random evaluations before GP fitting
#'                      (default 10).
#' @param n_iter        integer: acquisition-guided evaluations after the
#'                      random init (default 20).
#' @param bo_seed       integer: seed set immediately before calling
#'                      BayesianOptimization() (default 42).
#'                      rBayesianOptimization::BayesianOptimization() has no
#'                      seed argument of its own -- its init_points draws
#'                      come from the ambient RNG state, so without fixing
#'                      this the search is not reproducible run-to-run (and,
#'                      with few init_points in a small search space, can
#'                      occasionally draw near-duplicate points that make
#'                      GPfit::GP_fit()'s covariance matrix singular).
#' @param verbose       logical: print BO progress? (default TRUE)
#' @param ...           additional arguments passed through to
#'                      select_K_cv() / the underlying fit function
#'                      (e.g. prior_beta, max_iter).
#'
#' @return The list returned by rBayesianOptimization::BayesianOptimization():
#'   $Best_Par (named vector with K, alpha), $Best_Value, $History, $Pred.
#'
#' @examples
#' \dontrun{
#' res <- select_k_alpha_bayesopt(Y_train, time_train, status_train,
#'          K_bounds = c(2L, 10L), alpha_bounds = c(0, 1),
#'          init_points = 10, n_iter = 20)
#' res$Best_Par   # c(K = ..., alpha = ...)
#' }
select_k_alpha_bayesopt <- function(Y, time, status,
                                     K_bounds     = c(2L, 10L),
                                     alpha_bounds = c(0, 1),
                                     n_folds      = 5,
                                     seed         = 42,
                                     model        = "YFB",
                                     init_points  = 10,
                                     n_iter       = 20,
                                     bo_seed      = 42,
                                     verbose      = TRUE,
                                     ...) {
  if (length(K_bounds) != 2 || K_bounds[1] >= K_bounds[2])
    stop("K_bounds must be c(lower, upper) with lower < upper.")
  if (length(alpha_bounds) != 2 || alpha_bounds[1] >= alpha_bounds[2] ||
      alpha_bounds[1] < 0 || alpha_bounds[2] > 1)
    stop("alpha_bounds must be c(lower, upper) with 0 <= lower < upper <= 1.")

  extra_args <- list(...)

  objective <- function(K, alpha) {
    do.call(select_k_alpha_bo_objective,
            c(list(Y = Y, time = time, status = status, K = K, alpha = alpha,
                   n_folds = n_folds, seed = seed, model = model),
              extra_args))
  }

  set.seed(bo_seed)
  rBayesianOptimization::BayesianOptimization(
    FUN         = objective,
    bounds      = list(K = K_bounds, alpha = alpha_bounds),
    init_points = init_points,
    n_iter      = n_iter,
    acq         = "ucb",
    verbose     = verbose
  )
}

#' Select a Trustworthy Winner From a BO History Table
#'
#' rBayesianOptimization's raw Best_Par can land on a degenerate alpha
#' extreme -- alpha near 0 (model ignores survival entirely) or alpha near 1
#' (ignores genomics entirely) -- that scores well on the CV objective by
#' incidental unsupervised-reconstruction alignment with the outcome, not
#' genuine survival modeling. This is a documented real failure mode in this
#' project (DECISIONS.md 2026-05-05: alpha=1.0 degenerate K-CV selection with
#' K_eff=0; the "lucky PCA direction alignment" archived-baseline finding),
#' and the BO objective itself has no defense against it since raw CV
#' concordance doesn't know or care whether beta is actually non-zero.
#'
#' This function re-fits the top `n_candidates` (K, alpha) points from a BO
#' History table on the full training data and returns the best-scoring
#' candidate that has at least one active survival factor (K_eff > 0 at
#' `beta_threshold`), instead of blindly trusting BO's own Best_Par.
#'
#' @param history       data.frame with columns K, alpha, Value (the shape
#'                      of BayesianOptimization()$History).
#' @param Y, time, status  training data for the refit.
#' @param n_candidates  integer: how many top-scoring History rows to
#'                      re-check (default 5).
#' @param beta_threshold  scalar: |beta| threshold for counting K_eff
#'                      (default 0.001, matching this project's convention
#'                      elsewhere -- e.g. config/globals.yml's
#'                      k_selection$beta_threshold).
#' @param fit_fn        function used to refit each candidate; must accept
#'                      (Y, time, status, K, alpha, verbose, ...) and return
#'                      a list with an $EBeta element (default
#'                      fit_cox_on_yf). Overridable for testing the
#'                      selection logic in isolation from CAVI's numerical
#'                      behavior on any particular dataset.
#' @param ...           passed through to fit_fn() (e.g. max_iter, prior_beta).
#'
#' @return list(K, alpha, cv_value, fit, k_eff, candidates_checked) where
#'   candidates_checked is a data.frame (K, alpha, cv_value, k_eff) for every
#'   candidate examined, ordered by cv_value descending. Throws an
#'   informative error if none of the top n_candidates has K_eff > 0 -- a
#'   genuine finding to report, not something to silently paper over.
pick_trustworthy_bo_winner <- function(history, Y, time, status,
                                        n_candidates   = 5,
                                        beta_threshold = 0.001,
                                        fit_fn         = fit_cox_on_yf,
                                        ...) {
  if (!is.function(fit_fn))
    stop("fit_fn must be a function (default fit_cox_on_yf; override only for testing).")
  extra_args <- list(...)
  extra_args[["verbose"]] <- NULL   # this function always fits quietly

  ord <- order(history$Value, decreasing = TRUE)
  top <- history[ord[seq_len(min(n_candidates, nrow(history)))], , drop = FALSE]

  fits    <- vector("list", nrow(top))
  checked <- vector("list", nrow(top))
  for (i in seq_len(nrow(top))) {
    K_i   <- round(top$K[i])
    set.seed(42L)
    fit_i <- suppressMessages(do.call(fit_fn, c(
      list(Y = Y, time = time, status = status, K = K_i, alpha = top$alpha[i], verbose = FALSE),
      extra_args
    )))
    fits[[i]]    <- fit_i
    checked[[i]] <- data.frame(
      K = K_i, alpha = top$alpha[i], cv_value = top$Value[i],
      k_eff = sum(abs(fit_i$EBeta) > beta_threshold),
      stringsAsFactors = FALSE
    )
  }
  candidates_checked <- do.call(rbind, checked)

  trustworthy_idx <- which(candidates_checked$k_eff > 0)
  if (length(trustworthy_idx) == 0) {
    stop(sprintf(paste(
      "None of the top %d BO candidates have any active survival factor",
      "(K_eff > 0) -- the search's top region appears to be dominated by",
      "degenerate alpha rather than genuine survival signal. Inspect the",
      "returned candidates_checked before proceeding."
    ), nrow(top)))
  }
  winner_i <- trustworthy_idx[1]   # candidates_checked is already cv_value-descending

  list(
    K        = candidates_checked$K[winner_i],
    alpha    = candidates_checked$alpha[winner_i],
    cv_value = candidates_checked$cv_value[winner_i],
    fit      = fits[[winner_i]],
    k_eff    = candidates_checked$k_eff[winner_i],
    candidates_checked = candidates_checked
  )
}
