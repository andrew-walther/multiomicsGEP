# ============================================================
# Script: preprocess_multimodal_yfb.R
# Purpose: Validate and freeze matched multimodal YFB data ordering.
# Author: Andrew Walther
# Created: 2026-09-17
# Dependencies: base R
# ============================================================

.multimodal_yfb_modalities <- c("expression", "methylation")

.validate_multimodal_yfb_blocks <- function(Y, label, nonnegative = NULL,
                                            modalities = NULL) {
  # Y may hold both modalities or one of them (e.g. expression alone, to fit
  # the same model to a single modality); blocks keep the canonical order.
  if (!is.list(Y) || is.null(names(Y)) || length(Y) == 0L ||
      !all(names(Y) %in% .multimodal_yfb_modalities) || anyDuplicated(names(Y))) {
    stop(label, " Y must be a named list of expression and/or methylation blocks.")
  }
  modalities <- modalities %||% .multimodal_yfb_modalities[.multimodal_yfb_modalities %in% names(Y)]
  if (!identical(names(Y), modalities)) {
    stop(label, " Y must contain exactly the blocks ", paste(modalities, collapse = ", "),
         ", in that order.")
  }
  # nonnegative: named logical by modality. Blocks fit with a point-exponential
  # loading prior must be nonnegative; signed priors (point-Laplace, Normal)
  # allow centered or otherwise signed data. Default: all blocks nonnegative.
  if (is.null(nonnegative)) {
    nonnegative <- stats::setNames(rep(TRUE, length(modalities)), modalities)
  }
  for (modality in modalities) {
    block <- Y[[modality]]
    if (!is.matrix(block) || !is.numeric(block) || is.null(rownames(block)) ||
        is.null(colnames(block)) || anyDuplicated(rownames(block)) ||
        anyDuplicated(colnames(block)) || any(!is.finite(block)) ||
        (isTRUE(nonnegative[[modality]]) && any(block < 0))) {
      stop(label, " ", modality, " must be a finite ",
           if (isTRUE(nonnegative[[modality]])) "nonnegative " else "",
           "numeric matrix with unique row and column names.")
    }
  }
  subject_ids <- rownames(Y[[1]])
  for (modality in modalities[-1]) {
    if (!setequal(subject_ids, rownames(Y[[modality]]))) {
      stop(label, " all modalities must contain exactly the same subject IDs.")
    }
    Y[[modality]] <- Y[[modality]][subject_ids, , drop = FALSE]
  }
  Y
}

#' Validate training data and store frozen multimodal feature ordering
#'
#' @param Y Named expression and methylation subject-by-feature matrices.
#' @param time Named positive follow-up-time vector.
#' @param event Named 0/1 event-indicator vector.
#' @param nonnegative Named logical by modality: must the block be
#'   nonnegative? (Default: all TRUE, as required by point-exponential loadings.)
#' @return List with aligned `Y`, outcome vectors, and `training_spec` metadata.
#' @examples
#' Y <- list(expression = matrix(1, 1, 1, dimnames = list("s1", "g1")),
#'           methylation = matrix(1, 1, 1, dimnames = list("s1", "cg1")))
#' preprocess_multimodal_yfb_training(Y, c(s1 = 1), c(s1 = 1))
#' @family multimodal_yfb_preprocessing
preprocess_multimodal_yfb_training <- function(Y, time, event, nonnegative = NULL) {
  Y <- .validate_multimodal_yfb_blocks(Y, "Training", nonnegative)
  subject_ids <- rownames(Y[[1]])
  if (!is.numeric(time) || !is.numeric(event) || is.null(names(time)) ||
      is.null(names(event)) || anyDuplicated(names(time)) || anyDuplicated(names(event)) ||
      !setequal(names(time), subject_ids) || !setequal(names(event), subject_ids) ||
      any(!is.finite(time)) || any(time <= 0) || any(!(event %in% c(0, 1)))) {
    stop("Training time and event must be named, finite, and exactly match subject IDs; time must be positive and event 0/1.")
  }
  list(
    Y = Y,
    time = unname(time[subject_ids]),
    event = unname(event[subject_ids]),
    training_spec = list(
      subject_ids = subject_ids,
      feature_names = lapply(Y, colnames),
      modalities = names(Y),
      nonnegative = nonnegative
    )
  )
}

#' Align validation blocks to frozen multimodal training features
#'
#' @param Y Named expression and methylation validation matrices.
#' @param training_spec Metadata returned by `preprocess_multimodal_yfb_training()`.
#' @return List with aligned `Y` and ignored extra feature names by modality.
#' @examples
#' spec <- list(feature_names = list(expression = "g1", methylation = "cg1"))
#' Y <- list(expression = matrix(1, 1, 1, dimnames = list("v1", "g1")),
#'           methylation = matrix(1, 1, 1, dimnames = list("v1", "cg1")))
#' align_multimodal_yfb_prediction(Y, spec)
#' @family multimodal_yfb_preprocessing
align_multimodal_yfb_prediction <- function(Y, training_spec) {
  modalities <- training_spec$modalities %||% .multimodal_yfb_modalities
  if (!is.list(training_spec$feature_names) ||
      !identical(names(training_spec$feature_names), modalities)) {
    stop("training_spec must contain feature names for every training modality.")
  }
  # Validation data may carry extra modalities; keep the training ones only
  Y <- Y[intersect(modalities, names(Y))]
  Y <- .validate_multimodal_yfb_blocks(Y, "Validation", training_spec$nonnegative, modalities)
  ignored_features <- vector("list", length(modalities))
  names(ignored_features) <- modalities
  for (modality in modalities) {
    required <- training_spec$feature_names[[modality]]
    missing <- setdiff(required, colnames(Y[[modality]]))
    if (length(missing) > 0L) {
      stop("Validation ", modality, " is missing training feature(s): ",
           paste(missing, collapse = ", "), ".")
    }
    ignored_features[[modality]] <- setdiff(colnames(Y[[modality]]), required)
    Y[[modality]] <- Y[[modality]][, required, drop = FALSE]
  }
  list(Y = Y, ignored_features = ignored_features)
}
