# Methods for fitted models ----

#' Risk scores or program scores for new or training samples
#'
#' Scores new patients in closed form from their molecular data alone:
#' \eqn{\eta = \sum_m (y_m - \bar y_m)^\top F_m \beta}, with the training means
#' \eqn{\bar y_m} and posterior means of \eqn{F_m} and \eqn{\beta}. Higher risk
#' means shorter expected survival (the scale is the log hazard). Features are
#' matched by name: extra features are ignored, missing ones are an error.
#'
#' @param object an `ssbmf_fit`.
#' @param newdata `NULL` (training samples), an [ssbmf_data()] object, a matrix
#'   (one-block fits) or a named list of matrices.
#' @param type `"risk"` (eta, the log hazard) or `"scores"` (per-program
#'   projection scores Z, samples x K, on the fitted scale).
#' @param orientation orientation of matrix or list `newdata`.
#' @param ... unused.
#' @return A named numeric vector (`"risk"`) or a matrix (`"scores"`).
#' @family methods
#' @examples
#' sim <- ssbmf_simulate(n = 80, p = c(expression = 60), K = 2, seed = 1)
#' fit <- ssbmf_fit(ssbmf_data(sim$X, sim$time, sim$status), K = 3,
#'                  control = ssbmf_control(max_sweeps = 30))
#' head(predict(fit, newdata = sim$validation$X))
#' @export
predict.ssbmf_fit <- function(object, newdata = NULL, type = c("risk", "scores"),
                              orientation = c("features_by_samples", "samples_by_features"), ...) {
  type <- match.arg(type); orientation <- match.arg(orientation)
  if (is.null(newdata)) {
    Z <- object$scores
  } else {
    Y <- .align_newdata(newdata, object$features, orientation)
    Z <- Reduce(`+`, lapply(names(Y), function(m) {
      sweep(Y[[m]], 2, object$center[[m]], "-") %*% object$F[[m]]$mean
    }))
    rownames(Z) <- rownames(Y[[1]])
  }
  if (type == "scores") return(Z)
  stats::setNames(drop(Z %*% object$beta$mean), rownames(Z))
}

#' Canonical per-program scale: unit standard deviation of the training score
#' @keywords internal
#' @noRd
.canon_scale <- function(object) {
  s <- object$score_sd
  s[!is.finite(s) | s == 0] <- 1
  s
}

#' Survival coefficients with posterior intervals
#'
#' Coefficients are reported per standard deviation of each program's
#' training projection score (the canonical scale), so they are comparable
#' across programs and fits. Under this scale
#' `ssbmf_loadings()` and `coef()` reproduce `predict()` exactly.
#'
#' @param object an `ssbmf_fit`.
#' @param level coverage of the posterior interval (default 0.95).
#' @param ... unused.
#' @return data.frame with `program`, `mean`, `sd`, `lower`, `upper`, `lfsr`
#'   (local false sign rate, the posterior probability of the sign opposite to
#'   the mean).
#' @family methods
#' @export
coef.ssbmf_fit <- function(object, level = 0.95, ...) {
  s <- .canon_scale(object)
  b <- object$beta$mean * s; sdv <- object$beta$sd * s
  z <- stats::qnorm(1 - (1 - level) / 2)
  lfsr <- ifelse(sdv > 0, stats::pnorm(-abs(b) / sdv), ifelse(b == 0, 1, 0))
  data.frame(program = names(object$beta$mean), mean = unname(b), sd = unname(sdv),
             lower = unname(b - z * sdv), upper = unname(b + z * sdv), lfsr = unname(lfsr))
}

#' Posterior-mean loadings and activities on the canonical scale
#'
#' `ssbmf_loadings()` returns F_m divided by each program's score standard
#' deviation; `ssbmf_activities()` returns L multiplied by it, so
#' `L F_m^T` is unchanged.
#'
#' @param fit an `ssbmf_fit`.
#' @param block block name; `NULL` returns a list over blocks.
#' @return a features x K matrix (or a list of them) / a samples x K matrix.
#' @family methods
#' @export
ssbmf_loadings <- function(fit, block = NULL) {
  s <- .canon_scale(fit)
  out <- lapply(fit$F, function(Fm) sweep(Fm$mean, 2, s, "/"))
  if (is.null(block)) out else out[[block]]
}

#' @rdname ssbmf_loadings
#' @export
ssbmf_activities <- function(fit) sweep(fit$L$mean, 2, .canon_scale(fit), "*")

#' Top-loading features of each program
#'
#' @param fit an `ssbmf_fit`.
#' @param n number of features per program and block.
#' @param block block name or `NULL` for all blocks.
#' @return data.frame with `program`, `block`, `rank`, `feature`, `loading`
#'   (canonical scale), ordered by absolute loading.
#' @family methods
#' @export
ssbmf_top_features <- function(fit, n = 50, block = NULL) {
  Ls <- ssbmf_loadings(fit)
  if (!is.null(block)) Ls <- Ls[block]
  do.call(rbind, lapply(names(Ls), function(m) {
    A <- Ls[[m]]
    do.call(rbind, lapply(colnames(A), function(k) {
      o <- utils::head(order(-abs(A[, k])), n)
      data.frame(program = k, block = m, rank = seq_along(o), feature = rownames(A)[o],
                 loading = unname(A[o, k]))
    }))
  }))
}

#' @export
print.ssbmf_fit <- function(x, ...) {
  cat(sprintf("ssbmf_fit: %d samples (%d events), K = %d, blocks: %s\n", x$n, x$events,
              length(x$beta$mean), paste(sprintf("%s (%d)", names(x$features), lengths(x$features)), collapse = ", ")))
  cat(sprintf("  %s after %d sweeps; objective %.4f; steps lowering it: %d; rejected steps: %d\n",
              if (x$converged) "converged" else "NOT converged", x$sweeps,
              x$objective$final[["total"]], x$step_decreases, sum(x$rejections)))
  r <- table(factor(x$roles$role, c("survival_active", "reconstruction_only", "inactive")))
  cat(sprintf("  programs: %d survival-active, %d reconstruction-only, %d inactive\n", r[1], r[2], r[3]))
  invisible(x)
}

#' @export
summary.ssbmf_fit <- function(object, ...) {
  print(object)
  cat("\nCoefficients (per SD of program score, 95% intervals):\n")
  print(coef(object), digits = 3, row.names = FALSE)
  cat("\nRoles:\n")
  print(object$roles, digits = 3, row.names = FALSE)
  invisible(object)
}
