# ============================================================
# Script: train_test_split.R
# Purpose: Stratified train/test splitting for hold-out evaluation
#          of the SBMF model on survival data.
# Author: Claude Code (reviewed by Andrew Walther)
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
#'
#' @param status   integer vector of length n: event indicators (1 = event, 0 = censored)
#' @param n_folds  integer >= 2: number of folds (default 5)
#' @param seed     integer: random seed for reproducibility (default 42)
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
create_stratified_folds <- function(status, n_folds = 5, seed = 42) {

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

  set.seed(seed)
  events   <- sample(events, length(events))
  censored <- sample(censored, length(censored))

  folds <- vector("list", n_folds)
  for (k in seq_len(n_folds)) folds[[k]] <- integer(0)

  for (i in seq_along(events)) {
    fold_k <- ((i - 1L) %% n_folds) + 1L
    folds[[fold_k]] <- c(folds[[fold_k]], events[i])
  }
  for (i in seq_along(censored)) {
    fold_k <- ((i - 1L) %% n_folds) + 1L
    folds[[fold_k]] <- c(folds[[fold_k]], censored[i])
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
