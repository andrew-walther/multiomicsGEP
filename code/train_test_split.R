# ============================================================
# Script: train_test_split.R
# Purpose: Stratified train/test splitting for hold-out evaluation
#          of the SBMF model on survival data.
# Author: Andrew Walther
# Created: 2026-03-31
# Dependencies: none (base R only)
# ============================================================

# ============================================================
# stratified_split() ----
# ============================================================

#' Create a stratified train/test split preserving event rate.
#'
#' Splits n patients into training and test sets such that the
#' proportion of events (status == 1) is approximately equal in both.
#' This is critical for survival data because a test set with very
#' few events yields an unreliable C-index.
#'
#' **Strategy:**
#' 1. Separate patient indices into events (status == 1) and censored (status == 0)
#' 2. Sample test_frac from each group independently
#' 3. Guarantee at least 1 event and 1 censored patient in the test set
#'    (required for C-index computation)
#'
#' @param status   integer vector of length n: event indicators (1 = event, 0 = censored)
#' @param test_frac  numeric in (0, 1): fraction of data to hold out (default 0.2)
#' @param seed     integer: random seed for reproducibility (default 42)
#'
#' @return Named list:
#'   $train_idx   integer vector of training sample indices (sorted)
#'   $test_idx    integer vector of test sample indices (sorted)
#'   $n_train     number of training samples
#'   $n_test      number of test samples
#'   $event_rate_train  proportion of events in training set
#'   $event_rate_test   proportion of events in test set
#'
#' @examples
#' status <- c(1,1,1,0,0,0,0,1,0,1)
#' sp <- stratified_split(status, test_frac = 0.3, seed = 1)
#' # sp$train_idx, sp$test_idx are disjoint and cover 1:10
stratified_split <- function(status, test_frac = 0.2, seed = 42) {

  # --- Input validation ---
  n <- length(status)
  if (n < 4)
    stop(sprintf("Need at least 4 samples for a meaningful split (got %d).", n))
  if (!all(status %in% c(0L, 1L)))
    stop("status must contain only 0 and 1.")
  if (test_frac <= 0 || test_frac >= 1)
    stop("test_frac must be strictly between 0 and 1.")

  events   <- which(status == 1)
  censored <- which(status == 0)

  if (length(events) < 2)
    stop(sprintf("Need at least 2 events for stratified split (got %d).", length(events)))
  if (length(censored) < 2)
    stop(sprintf("Need at least 2 censored for stratified split (got %d).", length(censored)))

  # --- Stratified sampling ---
  set.seed(seed)

  # Number to hold out per group — at least 1 from each

  n_test_ev <- max(1L, round(length(events) * test_frac))
  n_test_cn <- max(1L, round(length(censored) * test_frac))

  # Ensure we don't exhaust either group (leave at least 1 for training)
  n_test_ev <- min(n_test_ev, length(events) - 1L)
  n_test_cn <- min(n_test_cn, length(censored) - 1L)

  test_ev <- sample(events, n_test_ev)
  test_cn <- sample(censored, n_test_cn)

  test_idx  <- sort(c(test_ev, test_cn))
  train_idx <- sort(setdiff(seq_len(n), test_idx))

  list(
    train_idx        = train_idx,
    test_idx         = test_idx,
    n_train          = length(train_idx),
    n_test           = length(test_idx),
    event_rate_train = mean(status[train_idx] == 1),
    event_rate_test  = mean(status[test_idx] == 1)
  )
}

# ============================================================
# create_stratified_folds() ----
# ============================================================

#' Create stratified K-fold splits preserving event rate.
#'
#' Assigns each sample to one of `n_folds` folds while approximately
#' preserving the event/censoring mix in every fold. This is intended for
#' cross-validation workflows where each held-out fold should contain both
#' events and censored samples so survival metrics remain well-defined.
#'
#' **Strategy:**
#' 1. Separate indices into events and censored groups
#' 2. Shuffle each group independently with the same seed
#' 3. Distribute each shuffled group round-robin across folds
#' With cohort_id, shuffle within joint cohort/event cells instead and continue
#' allocation across cohorts within each status. Sparse cohorts can be rejected
#' even if a different allocation could work. No automatic search is attempted.
#'
#' @param status   integer vector of length n: event indicators (1 = event, 0 = censored)
#' @param n_folds  integer >= 2: number of folds (default 5)
#' @param seed     integer: random seed for reproducibility (default 42)
#' @param cohort_id Optional vector of cohort labels, length n with no missing
#'   values. NULL preserves event-only folds. Otherwise balance each cohort/event
#'   cell across folds and require every training complement to retain all cohorts.
#'   Sparse cells are allowed, but an allocation losing a training cohort errors.
#'
#' @return Named list:
#'   $folds              list of length n_folds; each element is a sorted integer vector
#'                       of held-out sample indices for that fold
#'   $fold_id            integer vector of length n assigning each sample to a fold
#'   $event_rate_by_fold numeric vector of length n_folds with per-fold event rates
#'
#' @examples
#' status <- c(rep(1, 10), rep(0, 10))
#' fd <- create_stratified_folds(status, n_folds = 5, seed = 1)
#' lengths(fd$folds)
#' fd <- create_stratified_folds(status, 5, 1, cohort_id = rep(c("A", "B"), 10))
#' @family cross-validation
#' @seealso select_K_cv, select_alpha_cv, cv_survival_loglik, bicv_genomics_loglik
create_stratified_folds <- function(status, n_folds = 5, seed = 42,
                                    cohort_id = NULL) {

  n <- length(status)
  if (n < 4)
    stop(sprintf("Need at least 4 samples for stratified folds (got %d).", n))
  if (!all(status %in% c(0L, 1L)))
    stop("status must contain only 0 and 1.")
  if (length(n_folds) != 1 || !is.finite(n_folds) || n_folds < 2 || n_folds != as.integer(n_folds))
    stop("n_folds must be an integer >= 2.")

  n_folds <- as.integer(n_folds)
  events   <- which(status == 1)
  censored <- which(status == 0)

  if (length(events) < n_folds) {
    stop(sprintf("Need at least %d events to place >=1 event in every fold (got %d).",
                 n_folds, length(events)))
  }
  if (length(censored) < n_folds) {
    stop(sprintf("Need at least %d censored samples to place >=1 censored in every fold (got %d).",
                 n_folds, length(censored)))
  }

  folds <- vector("list", n_folds)
  for (k in seq_len(n_folds)) folds[[k]] <- integer(0)

  if (!is.null(cohort_id) &&
      (!is.atomic(cohort_id) || !is.null(dim(cohort_id)) ||
       length(cohort_id) != n || anyNA(cohort_id)))
    stop("cohort_id must be a vector of length status with no missing values.")

  set.seed(seed)
  if (is.null(cohort_id)) {
    # Keep the historical shuffle and assignment exactly for existing callers.
    events   <- sample(events, length(events))
    censored <- sample(censored, length(censored))
    for (i in seq_along(events)) {
      fold_k <- ((i - 1L) %% n_folds) + 1L
      folds[[fold_k]] <- c(folds[[fold_k]], events[i])
    }
    for (i in seq_along(censored)) {
      fold_k <- ((i - 1L) %% n_folds) + 1L
      folds[[fold_k]] <- c(folds[[fold_k]], censored[i])
    }
  } else {
    cohort_id <- as.character(cohort_id)
    # Continue round-robin allocation across cohorts within each event status.
    # Each joint cell and each pooled status differ by at most one across folds.
    # sample.int avoids sample(single_index)'s special 1:index behavior in R.
    for (st in c(1L, 0L)) {
      offset <- 0L
      for (co in unique(cohort_id)) {
        idx <- which(cohort_id == co & status == st)
        if (!length(idx)) next
        idx <- idx[sample.int(length(idx))]
        assignment <- ((seq_along(idx) - 1L + offset) %% n_folds) + 1L
        for (k in seq_len(n_folds))
          folds[[k]] <- c(folds[[k]], idx[assignment == k])
        offset <- (offset + length(idx)) %% n_folds
      }
    }
    # Cohort/event balance alone cannot ensure training coverage for tiny cohorts.
    # Do not return a split that would lose a Cox stratum or coefficient cohort.
    for (k in seq_len(n_folds)) {
      missing <- setdiff(unique(cohort_id), cohort_id[setdiff(seq_len(n), folds[[k]])])
      if (length(missing))
        stop(sprintf("Cohort(s) %s absent from training fold %d; this sparse dataset needs a different fold design.",
                     paste(missing, collapse = ", "), k))
    }
  }

  fold_id <- integer(n)
  for (k in seq_len(n_folds)) {
    folds[[k]] <- sort(folds[[k]])
    fold_id[folds[[k]]] <- k
  }

  list(
    folds = folds,
    fold_id = fold_id,
    event_rate_by_fold = vapply(folds, function(idx) mean(status[idx] == 1), numeric(1))
  )
}
