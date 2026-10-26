# Cross-validation over K ----

#' Subset an ssbmf_data object by sample index
#' @keywords internal
#' @noRd
.subset_data <- function(data, idx) {
  data$Y <- lapply(data$Y, function(A) A[idx, , drop = FALSE])
  data$time <- data$time[idx]; data$status <- data$status[idx]
  if (!is.null(data$strata)) data$strata <- droplevels(data$strata[idx])
  data$n <- length(idx); data$samples <- data$samples[idx]
  data
}

#' Harrell's C of a risk score (higher = worse), stratified if strata given
#' @keywords internal
#' @noRd
.cindex <- function(risk, time, status, strata = NULL) {
  if (is.null(strata)) {
    cc <- survival::concordance(survival::Surv(time, status) ~ risk, reverse = TRUE)
  } else {
    cc <- survival::concordance(survival::Surv(time, status) ~ risk + strata(strata), reverse = TRUE)
  }
  unname(cc$concordance)
}

#' Fold labels balanced on events (and strata)
#' @keywords internal
#' @noRd
.make_folds <- function(status, strata, nfolds, seed) {
  set.seed(seed)
  grp <- interaction(status, if (is.null(strata)) 1 else strata, drop = TRUE)
  folds <- integer(length(status))
  for (g in levels(grp)) {
    idx <- which(grp == g)
    folds[idx] <- sample(rep_len(seq_len(nfolds), length(idx)))
  }
  folds
}

#' Choose K by training cross-validated concordance
#'
#' For every K in `K_grid`, fits [ssbmf_fit()] on each training fold and
#' scores the held-out fold by projection (no refitting). The selected K is
#' the smallest whose mean C is within one standard error of the best
#' (one-standard-error rule). The model is then refitted on all of `data`.
#'
#' @param data an [ssbmf_data()] object (training data only).
#' @param K_grid integer vector of K values.
#' @param nfolds number of folds (default 5).
#' @param control settings from [ssbmf_control()].
#' @param seed seed for the fold assignment.
#' @param cores number of cores for the fold fits (`parallel::mclapply`;
#'   1 = sequential).
#' @return list with `cv` (data.frame: `K`, `mean_C`, `se_C`,
#'   `converged_folds`, `max_step_decreases`), `fold_C` (per fold), `folds`,
#'   `best` (selected K) and `fit` (refit on all data).
#' @family fit
#' @export
ssbmf_cv <- function(data, K_grid, nfolds = 5, control = ssbmf_control(), seed = 1, cores = 1) {
  if (!inherits(data, "ssbmf_data")) stop("`data` must be created with ssbmf_data().")
  folds <- .make_folds(data$status, data$strata, nfolds, seed)
  jobs <- expand.grid(K = K_grid, fold = seq_len(nfolds))
  run <- function(r) {
    K <- jobs$K[r]; f <- jobs$fold[r]
    tr <- which(folds != f); te <- which(folds == f)
    fit <- ssbmf_fit(.subset_data(data, tr), K, control)
    dte <- .subset_data(data, te)
    risk <- predict(fit, newdata = dte)
    data.frame(K = K, fold = f,
               C = .cindex(risk, dte$time, dte$status, dte$strata),
               converged = fit$converged, step_decreases = fit$step_decreases, sweeps = fit$sweeps)
  }
  res <- if (cores > 1) parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = cores) else lapply(seq_len(nrow(jobs)), run)
  bad <- vapply(res, function(x) inherits(x, "try-error") || !is.data.frame(x), TRUE)
  if (any(bad)) stop("Some cross-validation fits failed: ", paste(vapply(res[bad], as.character, ""), collapse = "; "))
  fold_C <- do.call(rbind, res)
  cv <- do.call(rbind, lapply(split(fold_C, fold_C$K), function(x) {
    data.frame(K = x$K[1], mean_C = mean(x$C), se_C = stats::sd(x$C) / sqrt(nrow(x)),
               converged_folds = sum(x$converged), max_step_decreases = max(x$step_decreases))
  }))
  i_best <- which.max(cv$mean_C)
  best <- min(cv$K[cv$mean_C >= cv$mean_C[i_best] - cv$se_C[i_best]])
  list(cv = cv, fold_C = fold_C, folds = folds, best = best, fit = ssbmf_fit(data, best, control))
}
