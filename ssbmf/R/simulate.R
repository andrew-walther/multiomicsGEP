# Simulated data ----

#' Simulate data blocks and survival from the ssbmf model
#'
#' Generates `n` training and `n_validation` validation samples:
#' \eqn{Y_m = 1\mu_m^T + L F_m^T + E_m} with nonnegative activities
#' (exponential, 20% zeros) and sparse nonnegative loadings on disjoint
#' feature sets (each program loads on 10% of each block's features).
#' Survival is Weibull (shape 1.5) with log hazard
#' \eqn{\eta = \sum_k \beta_k Z_k / \mathrm{sd}(Z_k)}, where
#' \eqn{Z_k = \tilde Y F_k} is the projection of the training-centered data on
#' the true loadings, so `beta` is the log hazard ratio per standard deviation
#' of the true program score. Censoring is uniform, with its upper limit set
#' to give the requested censored fraction in the training data.
#'
#' @param n number of training samples.
#' @param p named integer vector of features per block.
#' @param K number of true programs.
#' @param beta log hazard ratios per SD of each program score; default
#'   `c(0.7, 0.35, 0, ...)`, or, with `low_variance_prognostic = TRUE`,
#'   0 for every program except the last (1.0).
#' @param low_variance_prognostic if `TRUE`, the last program's loadings are
#'   scaled down so that it explains about 1% of each block's variance.
#' @param feature_correlation `"none"` or `"block"`: with `"block"`, features
#'   come in groups of 20 that share a nuisance factor not related to the
#'   programs (within-group correlation about 0.3 in the noise), as real
#'   expression and methylation data do.
#' @param noise `"lower"` (noise SD 1) or `"higher"` (noise SD 2); loadings
#'   have mean 1, so a loaded feature's signal variance is about 2.
#' @param n_validation number of validation samples.
#' @param censoring target censored fraction (default 0.3).
#' @param seed random seed.
#' @return list with `X` (features x samples matrices by block), `time`,
#'   `status`, `truth` (`L`, `F`, `beta`, `pve` (program x block), `score_sd`),
#'   and `validation` (`X`, `time`, `status`).
#' @family simulate
#' @examples
#' sim <- ssbmf_simulate(n = 100, p = c(expression = 200, methylation = 300), K = 3, seed = 1)
#' str(sim$X)
#' @export
ssbmf_simulate <- function(n, p, K, beta = NULL, low_variance_prognostic = FALSE,
                           feature_correlation = c("none", "block"), noise = c("lower", "higher"),
                           n_validation = n, censoring = 0.3, seed = 1) {
  feature_correlation <- match.arg(feature_correlation); noise <- match.arg(noise)
  if (is.null(names(p))) names(p) <- paste0("block", seq_along(p))
  if (is.null(beta)) {
    beta <- if (low_variance_prognostic) c(rep(0, K - 1), 1.0) else c(0.7, 0.35, rep(0, max(0, K - 2)))[seq_len(K)]
  }
  stopifnot(length(beta) == K, K <= 10)
  set.seed(seed)
  N <- n + n_validation
  sd_e <- if (noise == "lower") 1 else 2
  L <- matrix(stats::rexp(N * K) * (stats::runif(N * K) > 0.2), N, K)
  Fs <- list(); Ys <- list()
  for (m in names(p)) {
    pm <- p[[m]]
    Fm <- matrix(0, pm, K)
    sets <- split(sample(pm), rep(seq_len(10), length.out = pm))   # 10 disjoint sets of 10%
    for (k in seq_len(K)) Fm[sets[[k]], k] <- stats::rexp(length(sets[[k]]))
    E <- matrix(stats::rnorm(N * pm, 0, sd_e), N, pm)
    if (feature_correlation == "block") {
      grp <- rep(seq_len(ceiling(pm / 20)), each = 20, length.out = pm)
      G <- matrix(stats::rnorm(N * max(grp)), N, max(grp))
      E <- sqrt(0.7) * E + sqrt(0.3) * sd_e * G[, grp] * rep(sample(c(-1, 1), pm, TRUE), each = N)
    }
    if (low_variance_prognostic) {
      # scale the last program so that it explains about 1% of the block variance
      signal_other <- L[, -K, drop = FALSE] %*% t(Fm[, -K, drop = FALSE])
      tot_other <- sum(scale(signal_other, scale = FALSE)^2) + sum(scale(E, scale = FALSE)^2)
      lk <- scale(L[, K], scale = FALSE)
      v_k <- sum(lk^2) * sum(Fm[, K]^2)
      Fm[, K] <- Fm[, K] * sqrt(0.01 / 0.99 * tot_other / v_k)
    }
    mu <- stats::rnorm(pm, 5, 1)
    Ym <- outer(rep(1, N), mu) + L %*% t(Fm) + E
    dimnames(Ym) <- list(paste0("s", seq_len(N)), paste0(substr(m, 1, 3), "_", seq_len(pm)))
    Fs[[m]] <- `dimnames<-`(Fm, list(colnames(Ym), paste0("T", seq_len(K))))
    Ys[[m]] <- Ym
  }
  tr <- seq_len(n); va <- n + seq_len(n_validation)
  # true projection scores, centered at training means
  Z <- Reduce(`+`, lapply(names(p), function(m) {
    sweep(Ys[[m]], 2, colMeans(Ys[[m]][tr, , drop = FALSE]), "-") %*% Fs[[m]]
  }))
  zsd <- apply(Z[tr, , drop = FALSE], 2, stats::sd)
  eta <- drop(sweep(Z, 2, zsd, "/") %*% beta)
  Tev <- (-log(stats::runif(N)) / (0.1 * exp(eta)))^(1 / 1.5)
  # censoring upper limit chosen so the training censored fraction is `censoring`
  U <- stats::runif(N)
  frac <- function(cmax) mean(U[tr] * cmax < Tev[tr])
  cmax <- if (censoring <= 0) Inf else stats::uniroot(function(c) frac(c) - censoring,
                                                      c(1e-3, 1e4 * max(Tev)))$root
  Cc <- U * cmax
  time <- pmin(Tev, Cc); status <- as.numeric(Tev <= Cc)
  pve <- vapply(names(p), function(m) {
    Yc <- scale(Ys[[m]][tr, ], scale = FALSE); tot <- sum(Yc^2)
    vapply(seq_len(K), function(k) sum(scale(L[tr, k], scale = FALSE)^2) * sum(Fs[[m]][, k]^2) / tot, 0)
  }, numeric(K))
  pve <- matrix(pve, K, length(p), dimnames = list(paste0("T", seq_len(K)), names(p)))
  X_of <- function(idx) lapply(Ys, function(A) t(A[idx, , drop = FALSE]))
  list(X = X_of(tr), time = time[tr], status = status[tr],
       truth = list(L = L[tr, , drop = FALSE], F = Fs, beta = beta, pve = pve, score_sd = zsd,
                    L_validation = L[va, , drop = FALSE], eta = eta[tr]),
       validation = list(X = X_of(va), time = time[va], status = status[va], eta = eta[va]))
}
