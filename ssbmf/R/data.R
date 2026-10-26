# Data constructor ----

#' Package data blocks and survival outcomes for ssbmf
#'
#' Validates one or more data blocks measured on the same patients, together
#' with survival times, event indicators and optional Cox strata. Blocks are
#' stored internally as samples x features matrices.
#'
#' @param X A numeric matrix (one block) or a named list of numeric matrices
#'   (one per block). Each matrix must have sample names (column names when
#'   `orientation = "features_by_samples"`, row names otherwise) and feature
#'   names. A single matrix becomes block `"block1"`.
#' @param time Positive survival or censoring times, one per sample.
#' @param status Event indicator (1 = event, 0 = censored), one per sample.
#' @param strata Optional stratum label per sample (e.g. cohort); used as Cox
#'   strata, so each stratum has its own baseline hazard.
#' @param orientation Whether `X` is `"features_by_samples"` (default, as in
#'   Bioconductor and DeSurv) or `"samples_by_features"`.
#' @return An object of class `ssbmf_data`: a list with `Y` (named list of
#'   samples x features matrices), `time`, `status`, `strata` (factor or
#'   `NULL`), `n`, `p` (named integer), `samples`, `features` (named list).
#' @family data
#' @examples
#' X <- matrix(rnorm(200), 20, 10,
#'             dimnames = list(paste0("g", 1:20), paste0("s", 1:10)))
#' dd <- ssbmf_data(X, time = rexp(10), status = rbinom(10, 1, 0.7))
#' dd
#' @export
ssbmf_data <- function(X, time, status, strata = NULL,
                       orientation = c("features_by_samples", "samples_by_features")) {
  orientation <- match.arg(orientation)
  if (is.matrix(X) || is.data.frame(X)) X <- list(block1 = X)
  if (!is.list(X) || length(X) == 0) stop("`X` must be a matrix or a non-empty named list of matrices.")
  if (is.null(names(X)) || any(names(X) == "")) {
    if (length(X) == 1) names(X) <- "block1" else stop("A list `X` must have a name for every block.")
  }
  Y <- lapply(names(X), function(m) {
    A <- as.matrix(X[[m]])
    if (!is.numeric(A)) stop(sprintf("Block '%s' is not numeric.", m))
    if (orientation == "features_by_samples") A <- t(A)
    if (is.null(rownames(A))) stop(sprintf("Block '%s' has no sample names.", m))
    if (is.null(colnames(A))) colnames(A) <- paste0(m, "_f", seq_len(ncol(A)))
    if (anyDuplicated(colnames(A))) stop(sprintf("Block '%s' has duplicated feature names.", m))
    if (anyDuplicated(rownames(A))) stop(sprintf("Block '%s' has duplicated sample names.", m))
    if (anyNA(A)) {
      bad <- colnames(A)[colSums(is.na(A)) > 0]
      stop(sprintf("Block '%s' has missing values in %d feature(s) (e.g. %s); impute or remove them first.",
                   m, length(bad), paste(utils::head(bad, 3), collapse = ", ")))
    }
    sds <- apply(A, 2, stats::sd)
    if (any(!is.finite(sds) | sds == 0)) {
      bad <- colnames(A)[!is.finite(sds) | sds == 0]
      stop(sprintf("Block '%s' has %d constant feature(s) (e.g. %s); remove them first.",
                   m, length(bad), paste(utils::head(bad, 3), collapse = ", ")))
    }
    A
  })
  names(Y) <- names(X)
  samples <- rownames(Y[[1]])
  for (m in names(Y)[-1]) {
    if (!setequal(rownames(Y[[m]]), samples) || nrow(Y[[m]]) != length(samples)) {
      stop(sprintf("Blocks '%s' and '%s' do not have the same sample names.", names(Y)[1], m))
    }
    Y[[m]] <- Y[[m]][samples, , drop = FALSE]
  }
  n <- length(samples)
  if (length(time) != n || length(status) != n) stop("`time` and `status` must have one value per sample.")
  time <- as.numeric(time); status <- as.numeric(status)
  if (anyNA(time) || any(time <= 0)) stop("`time` must be positive and not missing.")
  if (anyNA(status) || !all(status %in% c(0, 1))) stop("`status` must be 0 (censored) or 1 (event).")
  if (sum(status) == 0) stop("There are no events.")
  if (!is.null(strata)) {
    if (length(strata) != n) stop("`strata` must have one value per sample.")
    if (anyNA(strata)) stop("`strata` has missing values.")
    strata <- factor(strata)
  }
  structure(list(Y = Y, time = time, status = status, strata = strata, n = n,
                 p = vapply(Y, ncol, 1L), samples = samples,
                 features = lapply(Y, colnames)),
            class = "ssbmf_data")
}

#' @export
print.ssbmf_data <- function(x, ...) {
  cat(sprintf("ssbmf_data: %d samples, %d events, %d block(s)\n", x$n, sum(x$status), length(x$Y)))
  for (m in names(x$Y)) cat(sprintf("  %-14s %d features\n", m, x$p[[m]]))
  if (!is.null(x$strata)) cat(sprintf("  strata: %s\n", paste(levels(x$strata), collapse = ", ")))
  invisible(x)
}

#' Align new data to the training features of a fit
#'
#' @param newdata An `ssbmf_data` object, a matrix (one block) or a named list
#'   of matrices.
#' @param features Named list of training feature names by block.
#' @param orientation Orientation of matrix/list `newdata`.
#' @return Named list of samples x features matrices with training feature order.
#' @keywords internal
#' @noRd
.align_newdata <- function(newdata, features, orientation) {
  if (inherits(newdata, "ssbmf_data")) {
    Y <- newdata$Y
  } else {
    if (is.matrix(newdata) || is.data.frame(newdata)) {
      if (length(features) != 1) stop("The fit has several blocks; pass `newdata` as a named list.")
      newdata <- stats::setNames(list(newdata), names(features))
    }
    Y <- lapply(newdata, function(A) {
      A <- as.matrix(A)
      if (orientation == "features_by_samples") t(A) else A
    })
  }
  missing_blocks <- setdiff(names(features), names(Y))
  if (length(missing_blocks)) stop(sprintf("`newdata` lacks block(s): %s.", paste(missing_blocks, collapse = ", ")))
  out <- lapply(names(features), function(m) {
    A <- Y[[m]]
    miss <- setdiff(features[[m]], colnames(A))
    if (length(miss)) {
      stop(sprintf("`newdata` block '%s' lacks %d training feature(s), e.g. %s.",
                   m, length(miss), paste(utils::head(miss, 5), collapse = ", ")))
    }
    A[, features[[m]], drop = FALSE]
  })
  names(out) <- names(features)
  n <- vapply(out, nrow, 1L)
  if (length(unique(n)) != 1) stop("Blocks of `newdata` have different numbers of samples.")
  out
}
