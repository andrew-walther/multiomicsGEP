# ============================================================
# Script: code/compute_cv_loglik.R
# Purpose: Genuine held-out (cross-validated) log-likelihood criteria for
#          K_init selection, alongside the in-sample criteria in
#          code/compute_bic.R. Two separate quantities, reported
#          separately and never summed:
#            - cv_survival_loglik():  held-out Cox partial log-likelihood
#              (leakage-free: eta_new = (Y_test EF) beta is the exact
#              YFB prediction formula, not an approximation).
#            - bicv_genomics_loglik(): bi-cross-validated genomics
#              log-likelihood (Owen & Perry 2009). Ordinary row-wise
#              held-out genomics likelihood is not well-defined here --
#              predict.R's L_test = Y_test F (F'F)+ reuses the held-out
#              patient's own expression to score that same patient, so the
#              score falls monotonically in K by construction. There is
#              also no missing-data (NA/mask) support anywhere in
#              update_L.R / update_F.R / update_tau.R or the SVD
#              initializers, so Wold-style held-out-cell masking is not
#              available. Bi-cross-validation avoids both problems by
#              holding out rows and columns jointly.
# Author: Claude Code (reviewed by Andrew Walther)
# Created: 2026-09-04
# Dependencies: code/fit_cox_on_yf.R (fit_cox_on_yf, calc_cox_taylor_yf),
#               code/predict_cox_on_yf.R (predict_cox_on_yf),
#               code/train_test_split.R (create_stratified_folds),
#               code/update_L.R (compute_R_k),
#               code/update_L_surv_YFB.R (update_L_surv_YFB_k),
#               code/update_F_surv_YFB.R (update_F_surv_YFB_k),
#               code/update_tau.R (update_tau).
#               All must be sourced first (fit_cox_on_yf.R sources the
#               update_* files itself; sourcing it is sufficient).
# ============================================================

# ============================================================
# gaussian_matrix_loglik() ----
# ============================================================

#' Gaussian log-likelihood of a residual matrix under column-specific precision
#'
#' No reusable Gaussian log-density evaluator existed anywhere in the repo
#' before this file (compute_bic.R inlines its own Gaussian constant rather
#' than exposing a general one). This is that evaluator, used by
#' `bicv_genomics_loglik()` below and available for future scoring code.
#'
#'   loglik = sum_j [ (n_row/2) * log(tau_j / (2*pi)) ] - (1/2) * sum_ij [ tau_j * resid_ij^2 ]
#'
#' @param resid n x p numeric matrix of residuals (observed - predicted)
#' @param tau   p-vector of per-column noise precisions (> 0)
#'
#' @return scalar: total Gaussian log-likelihood of `resid` under N(0, 1/tau_j)
#'   per column j, independent across rows and columns.
gaussian_matrix_loglik <- function(resid, tau) {
  if (!is.matrix(resid)) stop("resid must be a matrix.")
  if (length(tau) != ncol(resid))
    stop(sprintf("tau has length %d, expected ncol(resid) = %d.", length(tau), ncol(resid)))
  if (any(tau <= 0) || any(!is.finite(tau)))
    stop("tau must be strictly positive and finite.")

  n_row <- nrow(resid)
  sum(0.5 * log(tau / (2 * pi))) * n_row - 0.5 * sum(sweep(resid^2, 2, tau, "*"))
}

# ============================================================
# .fit_genomics_only() ----
# ============================================================

#' Fit L, F, Tau via CAVI with no survival term at all (internal helper)
#'
#' `bicv_genomics_loglik()` below needs a genomics reconstruction fit on a
#' row x column submatrix, and wants nothing from the survival side of the
#' model. Under the Cox-on-YF (YFB) model at `alpha_F = 0` (the recommended
#' default), `fit_cox_on_yf()`'s q(L) and q(F) updates are already pure
#' genomics -- L never appears in the Cox likelihood under YFB, and F's
#' survival term is zero-weighted at alpha_F=0 (see update_L_surv_YFB.R,
#' update_F_surv_YFB.R). This helper reuses exactly those two update
#' functions plus update_tau(), omitting the beta/Cox machinery entirely
#' rather than calling fit_cox_on_yf() with a fabricated time/status vector.
#'
#' SVD initialization follows fit_cox_on_yf.R's convention: `abs()` rather
#' than `pmax(..., 0)`, since pmax can zero out an entire column when the
#' SVD vector points all-negative, producing a degenerate ZF column later.
#' Convergence is checked on the genomics elbo_proxy only (no survival term
#' exists to include).
#'
#' @param Y        n x p numeric matrix
#' @param K        integer: number of latent factors
#' @param max_iter integer: maximum CAVI iterations (default 100)
#' @param tol      numeric: relative elbo_proxy convergence threshold (default 1e-5)
#' @param prior_LF character: EBNM prior family for L and F (default "point_exponential")
#'
#' @return Named list: $EL, $EL2, $EF, $EF2, $Tau, $n_iter
.fit_genomics_only <- function(Y, K, max_iter = 100, tol = 1e-5,
                                prior_LF = "point_exponential") {

  required_fns <- c("compute_R_k", "update_L_surv_YFB_k", "update_F_surv_YFB_k", "update_tau")
  missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function")]
  if (length(missing_fns) > 0) {
    stop("Missing required function(s): ", paste(missing_fns, collapse = ", "),
         ". Source code/fit_cox_on_yf.R (which sources all of them) first.")
  }

  svd_init <- svd(Y, nu = K, nv = K)
  d_k <- sqrt(pmax(svd_init$d[1:K], 0))
  EL  <- abs(svd_init$u %*% diag(d_k, K, K))
  EF  <- abs(svd_init$v %*% diag(d_k, K, K))
  EL2 <- EL^2
  EF2 <- EF^2
  Tau <- 1.0 / pmax(apply(Y, 2, var), 1e-8)

  elbo_prev <- NA_real_
  n_iter <- max_iter
  for (iter in seq_len(max_iter)) {
    for (k in seq_len(K)) {
      R_k   <- compute_R_k(Y, EL, EF, k)
      res_L <- update_L_surv_YFB_k(Tau, EF[, k], EF2[, k], R_k, prior_family = prior_LF)
      EL[, k]  <- res_L$mean
      EL2[, k] <- res_L$second

      R_k   <- compute_R_k(Y, EL, EF, k)
      res_F <- update_F_surv_YFB_k(Tau, EL[, k], EL2[, k], R_k,
                                    prior_family = prior_LF, alpha = 0)
      EF[, k]  <- res_F$mean
      EF2[, k] <- res_F$second
    }

    res_tau <- update_tau(Y, EL, EL2, EF, EF2)
    Tau <- res_tau$Tau
    elbo <- res_tau$elbo_proxy

    if (iter > 5 && is.finite(elbo) && is.finite(elbo_prev) && elbo_prev != 0) {
      if (abs((elbo - elbo_prev) / abs(elbo_prev)) < tol) {
        n_iter <- iter
        break
      }
    }
    elbo_prev <- elbo
  }

  list(EL = EL, EL2 = EL2, EF = EF, EF2 = EF2, Tau = Tau, n_iter = n_iter)
}

# ============================================================
# cv_survival_loglik() ----
# ============================================================

#' Held-out Cox partial log-likelihood for a YFB fit, by stratified K-fold CV
#'
#' Answers "how well does this K generalize to held-out patients' survival
#' outcomes," as a genuine cross-validated companion to compute_bic.R's
#' in-sample survival term. Under YFB, eta_new = (Y_test EF) beta is the
#' SAME formula used at training time (predict_cox_on_yf()), so held-out
#' scoring introduces no train/test formula mismatch -- unlike the
#' genomics side (see bicv_genomics_loglik() below).
#'
#' \strong{Leakage guard (orientation).} Each fold is fit with
#' `sign_correction = FALSE`, then oriented from a correct-direction
#' (`reverse = TRUE`) training concordance check computed ONLY on that
#' fold's training rows: if training concordance is below 0.5, EBeta is
#' flipped before scoring the held-out fold. `reverse = TRUE` matters here:
#' `fit_cox_on_yf()`'s own post-loop sign rule (fit_cox_on_yf.R:705-720)
#' omits it, which is inverted for a Cox risk score (see
#' `survival::concordance()`'s own documentation and code/compute_bic.R's
#' matching orientation-correction comment, DECISIONS.md 2026-09-04) --
#' copying that rule verbatim here would silently re-introduce the same
#' bug into a held-out metric, where it is not sign-invariant and would
#' corrupt every score. Unlike C-index, partial log-likelihood is NOT
#' sign-invariant, so orienting against held-out outcomes (rather than
#' training outcomes) would leak held-out label information into the
#' score -- hence training-only, not "no orientation at all."
#'
#' \strong{Risk-set size.} Each fold's partial log-likelihood is computed
#' within that fold's own risk sets, so `logPL` here is not on the same
#' scale as a full-data logPL from compute_bic.R -- comparisons are only
#' valid across K at the same fold assignment (which is exactly this
#' function's use case).
#'
#' @param Y        n x p genomics matrix
#' @param time     n-vector: survival/censoring times
#' @param status   n-vector: event indicators (1 = event, 0 = censored)
#' @param K        integer: number of latent factors
#' @param n_folds  integer >= 2: stratified CV folds (default 5, matches
#'                 config/globals.yml cv_loglik$n_folds)
#' @param seed     integer: fold-assignment seed (default 42, matches
#'                 config/globals.yml cv_loglik$seed)
#' @param max_iter,tol,prior_LF,prior_beta,alpha: forwarded to
#'                 `fit_cox_on_yf()` for each fold's training fit.
#' @param cohort_id,strata_id,beta_cohort_id NULL (default), or n-vectors
#'                 forwarded to `fit_cox_on_yf()` -- row-subsetted to each
#'                 fold's training indices internally (the same pattern
#'                 `select_K.R::select_K_cv()` uses), since passing a
#'                 full-length vector straight through would trip
#'                 `fit_cox_on_yf()`'s own length check against the
#'                 fold-subsetted `Y`. `strata_id` is ALSO row-subsetted to
#'                 the test fold and passed to the held-out
#'                 `calc_cox_taylor_yf()` call (fixed 2026-09-04, DECISIONS.md
#'                 -- previously a model fit under stratified risk sets was
#'                 scored under one pooled risk set on the held-out side).
#' @param cv_scoring "within_cohort" (default) or "unseen_cohort" -- only
#'                 meaningful when `beta_cohort_id` is supplied; ignored
#'                 otherwise. These are two DIFFERENT questions and are not
#'                 interchangeable (fixed 2026-09-04, DECISIONS.md -- every
#'                 fold previously scored with `EBeta_pooled` regardless):
#'                 "within_cohort" answers "how well does this model predict
#'                 a held-out PATIENT from a cohort it has already seen,
#'                 using that patient's own (held-out but KNOWN) cohort
#'                 label's beta column" -- the natural reading of an
#'                 ordinary K-fold CV, since folds split patients within the
#'                 same fixed set of training cohorts, not cohorts
#'                 themselves. "unseen_cohort" answers "how well does the
#'                 pooled/shared beta generalize to a cohort with NO
#'                 beta^(c) of its own" -- scores every fold with
#'                 `EBeta_pooled`, ignoring the test fold's own (known)
#'                 cohort labels; this is the right question for comparing
#'                 against genuinely external cohorts (see
#'                 `predict_cox_on_yf()`'s `cohort_id_test = NULL` path
#'                 elsewhere in this project), not for an ordinary CV fold.
#' @param ...      additional arguments forwarded to `fit_cox_on_yf()`
#'                 (e.g. init_method). Do not pass K or sign_correction.
#'
#' @return Named list:
#'   $K                     echoed input K
#'   $cv_scoring            echoed input cv_scoring, or NA_character_ when
#'                          beta_cohort_id is NULL (the argument is unused)
#'   $fold_results          data.frame: fold, n_train, n_test, n_event_test,
#'                          c_train, logPL, logPL_per_event
#'   $total_logPL           sum of held-out logPL across folds
#'   $total_events          sum of held-out event counts across folds
#'   $mean_logPL_per_event  total_logPL / total_events
#'   $sd_logPL, $se_logPL   spread of per-fold logPL across folds
cv_survival_loglik <- function(Y, time, status, K,
                                n_folds  = 5L,
                                seed     = 42L,
                                max_iter = 100,
                                tol      = 1e-5,
                                prior_LF   = "point_exponential",
                                prior_beta = "normal",
                                alpha      = 0.5,
                                cohort_id       = NULL,
                                strata_id       = NULL,
                                beta_cohort_id  = NULL,
                                cv_scoring      = c("within_cohort", "unseen_cohort"),
                                ...) {

  cv_scoring <- match.arg(cv_scoring)
  required_fns <- c("fit_cox_on_yf", "predict_cox_on_yf", "calc_cox_taylor_yf",
                     "create_stratified_folds")
  missing_fns <- required_fns[!vapply(required_fns, exists, logical(1), mode = "function")]
  if (length(missing_fns) > 0) {
    stop("Missing required function(s): ", paste(missing_fns, collapse = ", "),
         ". Source code/fit_cox_on_yf.R and code/train_test_split.R first.")
  }
  if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix.")
  n <- nrow(Y)
  if (length(time) != n)   stop(sprintf("time must have length %d (nrow(Y)).", n))
  if (length(status) != n) stop(sprintf("status must have length %d (nrow(Y)).", n))

  extra <- list(...)
  extra[["sign_correction"]] <- NULL  # this function controls orientation itself

  fold_obj <- create_stratified_folds(status, n_folds = n_folds, seed = seed)

  fold_rows <- vector("list", n_folds)
  for (fold_id in seq_len(n_folds)) {
    test_idx  <- fold_obj$folds[[fold_id]]
    train_idx <- setdiff(seq_len(n), test_idx)

    fit_args <- c(
      list(Y = Y[train_idx, , drop = FALSE], time = time[train_idx], status = status[train_idx],
           K = K, max_iter = max_iter, tol = tol,
           prior_LF = prior_LF, prior_beta = prior_beta, alpha = alpha,
           sign_correction = FALSE, verbose = FALSE,
           cohort_id      = if (is.null(cohort_id)) NULL else cohort_id[train_idx],
           strata_id      = if (is.null(strata_id)) NULL else strata_id[train_idx],
           beta_cohort_id = if (is.null(beta_cohort_id)) NULL else beta_cohort_id[train_idx]),
      extra
    )
    fit <- do.call(fit_cox_on_yf, fit_args)
    use_cohort <- !is.null(beta_cohort_id)

    # Orientation from TRAINING concordance only -- see roxygen leakage guard
    # above. reverse = TRUE: eta is a Cox risk score (larger eta -> shorter
    # survival), the opposite of concordance()'s formula-method default
    # assumption. EBeta_for_orient/eta_train/the resulting flip decision
    # depend on cv_scoring -- see the @param cv_scoring roxygen above for why
    # these are genuinely different questions, not interchangeable:
    #   - not use_cohort:      unchanged (fit$EBeta, a K-vector).
    #   - "unseen_cohort":     fit$EBeta_pooled, ignoring the fold's own
    #                          (known) cohort labels -- the pre-2026-09-04
    #                          behavior, unconditional on cv_scoring.
    #   - "within_cohort":     fit$EBeta (the K x C matrix), oriented via a
    #                          single GLOBAL flip decided from the pooled
    #                          training eta using each training patient's own
    #                          known cohort column -- mirrors fit_cox_on_yf()'s
    #                          own Phase C convention (one global flip, applied
    #                          uniformly to every cohort column).
    ZF_train <- Y[train_idx, , drop = FALSE] %*% sweep(fit$EF, 2, fit$EF_norms, "/")
    if (!use_cohort) {
      EBeta_for_orient <- fit$EBeta
      eta_train <- as.vector(ZF_train %*% fit$EBeta)
    } else if (cv_scoring == "unseen_cohort") {
      EBeta_for_orient <- fit$EBeta_pooled
      eta_train <- as.vector(ZF_train %*% fit$EBeta_pooled)
    } else {  # within_cohort
      EBeta_for_orient <- fit$EBeta
      train_cohort_idx <- match(beta_cohort_id[train_idx], colnames(fit$EBeta))
      eta_train <- rowSums(ZF_train * t(fit$EBeta)[train_cohort_idx, , drop = FALSE])
    }
    c_train   <- tryCatch(
      as.numeric(survival::concordance(
        survival::Surv(time[train_idx], status[train_idx]) ~ eta_train, reverse = TRUE
      )$concordance),
      error = function(e) NA_real_
    )
    EBeta_oriented <- EBeta_for_orient
    if (is.finite(c_train) && c_train < 0.5) EBeta_oriented <- -EBeta_oriented

    # Held-out scoring: eta_new = (Y_test EF) beta -- exact YFB prediction
    # formula. "within_cohort" scores each test patient with their own
    # (held-out but known) cohort's column via cohort_id_test; the other two
    # cases pass a K-vector, so cohort_id_test stays NULL (its only valid
    # value for a non-matrix EBeta).
    pred <- if (use_cohort && cv_scoring == "within_cohort") {
      predict_cox_on_yf(Y[test_idx, , drop = FALSE], fit$EF, EBeta_oriented,
                         EF_norms = fit$EF_norms, cohort_id_test = beta_cohort_id[test_idx])
    } else {
      predict_cox_on_yf(Y[test_idx, , drop = FALSE], fit$EF, EBeta_oriented,
                         EF_norms = fit$EF_norms)
    }
    # strata_id row-subsetted to the test fold (fixed 2026-09-04, DECISIONS.md):
    # a model fit under stratified risk sets must be scored under stratified
    # risk sets, not one pooled risk set.
    surv <- calc_cox_taylor_yf(pred$risk_scores, time[test_idx], status[test_idx],
                                strata = if (is.null(strata_id)) NULL else strata_id[test_idx])

    n_event_test <- sum(status[test_idx] == 1)
    fold_rows[[fold_id]] <- data.frame(
      fold            = fold_id,
      n_train         = length(train_idx),
      n_test          = length(test_idx),
      n_event_test    = n_event_test,
      c_train         = c_train,
      logPL           = surv$logPL,
      logPL_per_event = if (n_event_test > 0) surv$logPL / n_event_test else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  fold_results <- do.call(rbind, fold_rows)
  rownames(fold_results) <- NULL

  list(
    K                    = K,
    cv_scoring           = if (is.null(beta_cohort_id)) NA_character_ else cv_scoring,
    fold_results         = fold_results,
    total_logPL          = sum(fold_results$logPL),
    total_events         = sum(fold_results$n_event_test),
    mean_logPL_per_event = sum(fold_results$logPL) / sum(fold_results$n_event_test),
    sd_logPL             = stats::sd(fold_results$logPL),
    se_logPL             = stats::sd(fold_results$logPL) / sqrt(n_folds)
  )
}

# ============================================================
# bicv_genomics_loglik() ----
# ============================================================

#' Bi-cross-validated genomics log-likelihood (Owen & Perry 2009)
#'
#' Row-wise held-out genomics likelihood is not well-defined for this model:
#' predict.R's `L_test = Y_test F (F'F)+` reconstructs a held-out patient's
#' loadings FROM that same patient's expression, so scoring Y_test against
#' L_test F' reuses the data being scored and the apparent fit improves
#' monotonically in K regardless of true structure. There is also no
#' missing-data (NA/mask) support anywhere in the CAVI loop, so Wold-style
#' single-cell masking is not available either. Bi-cross-validation holds
#' out row AND column blocks jointly, so genes used to project the held-out
#' patients are never the genes being scored, and patients used to project
#' the held-out genes are never the patients being scored.
#'
#' For each (row-fold r, col-fold c):
#'   1. Fit genomics-only CAVI on the train x train block `Y[-r, -c]`
#'      (`.fit_genomics_only()` -- no survival term at all).
#'   2. Project held-out rows using TRAIN genes:
#'      `L_r = Y[r, -c] EF_train (EF_train'EF_train)+` (same ridge
#'      pseudoinverse as predict.R's predict_supervised_mf()).
#'   3. Project held-out genes using TRAIN rows:
#'      `F_c = (EL_train'EL_train)+ EL_train' Y[-r, c]`.
#'   4. Estimate "training" noise precision for the held-out genes from
#'      TRAIN ROWS only: `Tau_c = n_train / colSums((Y[-r,c] - EL_train F_c)^2)`,
#'      floored the same way as update_tau.R's tau_floor. This uses F_c
#'      (fit in step 3, on train rows) but never touches held-out rows, so
#'      it is a genuine training-side quantity despite scoring held-out genes.
#'   5. Score the held-out block `Y[r, c]` against `L_r F_c'` under `Tau_c`
#'      via `gaussian_matrix_loglik()`.
#'
#' Every (row-fold, col-fold) pair partitions Y into disjoint blocks whose
#' union is the full n x p matrix, so the sum across all blocks is a
#' held-out log-likelihood for every entry of Y, each scored exactly once.
#'
#' `update_tau.R`'s `n / pmax(col_sums, n * tau_floor)` hardcodes the
#' observed-entry count per column as `n` (no ragged columns). That is
#' exactly why Wold-style single-cell masking was not chosen here instead
#' of bi-cross-validation: bi-CV always fits and scores complete
#' submatrices, so this hardcoded `n` is never violated.
#'
#' @param Y           n x p genomics matrix
#' @param status      n-vector: event indicators, used only to stratify row
#'                    folds (`create_stratified_folds()`); no survival term
#'                    enters this function's fitting or scoring at all.
#' @param K           integer: number of latent factors
#' @param n_row_folds integer >= 2: patient folds (default 5, matches
#'                    config/globals.yml cv_loglik$n_row_folds)
#' @param n_col_folds integer >= 2: gene folds (default 5, matches
#'                    config/globals.yml cv_loglik$n_col_folds)
#' @param seed        integer: fold-assignment seed for both row and column
#'                    folds (default 42, matches config/globals.yml cv_loglik$seed)
#' @param max_iter,tol,prior_LF: forwarded to `.fit_genomics_only()` for
#'                    each train x train block fit.
#' @param lambda      numeric: ridge threshold for the pseudoinverse
#'                    projections (default 1e-8, matches predict.R's default).
#'
#' @return Named list:
#'   $K                total (summed) held-out genomics log-likelihood
#'   $block_results    data.frame: row_fold, col_fold, n_test_rows,
#'                     n_test_cols, block_loglik
#'   $total_loglik     sum of block_loglik across all blocks (all n*p cells)
#'   $n_cells_scored   sanity check: should equal n*p
bicv_genomics_loglik <- function(Y, status, K,
                                  n_row_folds = 5L,
                                  n_col_folds = 5L,
                                  seed        = 42L,
                                  max_iter    = 100,
                                  tol         = 1e-5,
                                  prior_LF    = "point_exponential",
                                  lambda      = 1e-8) {

  if (!exists("create_stratified_folds", mode = "function"))
    stop("create_stratified_folds() must be sourced (code/train_test_split.R) first.")
  if (!is.matrix(Y) || !is.numeric(Y)) stop("Y must be a numeric matrix.")

  n <- nrow(Y); p <- ncol(Y)
  tau_floor <- 1e-8  # matches update_tau.R's default

  row_fold_obj <- create_stratified_folds(status, n_folds = n_row_folds, seed = seed)
  set.seed(seed)
  col_fold_id <- sample(rep_len(seq_len(n_col_folds), p))

  .ridge_pinv <- function(M, K) {
    sv    <- svd(M)
    d_inv <- ifelse(sv$d > lambda * max(sv$d), 1 / sv$d, 0)
    sv$v %*% diag(d_inv, nrow = K, ncol = K) %*% t(sv$u)
  }

  block_rows <- vector("list", n_row_folds * n_col_folds)
  idx <- 1L
  for (r in seq_len(n_row_folds)) {
    row_test  <- row_fold_obj$folds[[r]]
    row_train <- setdiff(seq_len(n), row_test)

    for (c in seq_len(n_col_folds)) {
      col_test  <- which(col_fold_id == c)
      col_train <- which(col_fold_id != c)

      # ---- 1. Fit genomics-only CAVI on the train x train block ----------
      fit_tt <- .fit_genomics_only(Y[row_train, col_train, drop = FALSE], K = K,
                                    max_iter = max_iter, tol = tol, prior_LF = prior_LF)
      EF_train <- fit_tt$EF   # length(col_train) x K
      EL_train <- fit_tt$EL   # length(row_train) x K

      # ---- 2. Held-out rows, train genes -> L_r ---------------------------
      FtF_pinv <- .ridge_pinv(crossprod(EF_train), K)
      L_r <- Y[row_test, col_train, drop = FALSE] %*% EF_train %*% FtF_pinv

      # ---- 3. Held-out genes, train rows -> F_c ---------------------------
      LtL_pinv <- .ridge_pinv(crossprod(EL_train), K)
      F_c <- LtL_pinv %*% t(EL_train) %*% Y[row_train, col_test, drop = FALSE]  # K x length(col_test)

      # ---- 4. "Training" Tau for held-out genes, from TRAIN ROWS only ----
      resid_train_heldoutgenes <- Y[row_train, col_test, drop = FALSE] - EL_train %*% F_c
      n_train <- length(row_train)
      col_sums_c <- colSums(resid_train_heldoutgenes^2)
      Tau_c <- n_train / pmax(col_sums_c, n_train * tau_floor)

      # ---- 5. Score the held-out block -------------------------------------
      pred_block  <- L_r %*% F_c
      resid_block <- Y[row_test, col_test, drop = FALSE] - pred_block
      block_ll    <- gaussian_matrix_loglik(resid_block, Tau_c)

      block_rows[[idx]] <- data.frame(
        row_fold     = r,
        col_fold     = c,
        n_test_rows  = length(row_test),
        n_test_cols  = length(col_test),
        block_loglik = block_ll,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  block_results <- do.call(rbind, block_rows)
  rownames(block_results) <- NULL

  list(
    K              = K,
    block_results  = block_results,
    total_loglik   = sum(block_results$block_loglik),
    n_cells_scored = sum(block_results$n_test_rows * block_results$n_test_cols)
  )
}
