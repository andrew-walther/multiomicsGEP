# Data preparation and initialization (derivation, Section 6) ----

#' Prepare the internal data object used by the updates
#' @param data an `ssbmf_data` object
#' @param control an `ssbmf_control` object
#' @return list with blocks `Y`, centered blocks `Yc` (training means), their
#'   squares `Yc2`, `center`, survival quantities and settings
#' @keywords internal
#' @noRd
.prep_data <- function(data, control) {
  M <- length(data$Y); blocks <- names(data$Y)
  center <- lapply(data$Y, colMeans)
  Yc <- lapply(blocks, function(m) sweep(data$Y[[m]], 2, center[[m]], "-")); names(Yc) <- blocks
  w_block <- rep(1, M); names(w_block) <- blocks
  if (!is.null(control$w_block)) w_block[names(control$w_block)] <- control$w_block
  prior_F <- if (length(control$prior_F) == 1) rep(control$prior_F, M) else control$prior_F[blocks]
  if (anyNA(prior_F)) stop("`prior_F` must name every block.")
  names(prior_F) <- blocks
  list(Y = data$Y, Yc = Yc, Yc2 = lapply(Yc, function(A) A^2), center = center,
       n = data$n, M = M, p = data$p, blocks = blocks,
       status = data$status, time = data$time, strata = data$strata,
       ri = .risk_index(data$time, data$strata),
       tau_max = lapply(Yc, function(A) control$tau_cap / apply(A, 2, stats::var)),
       w_block = w_block, w_surv = control$w_surv,
       prior_L = control$prior_L, prior_F = prior_F,
       max_halvings = control$max_halvings, beta_inner = control$beta_inner)
}

#' Initial state: scaled concatenated SVD, nonnegative activities, least-squares loadings
#' @keywords internal
#' @noRd
.init_state <- function(d, K, init = "svd") {
  # 1. scale each centered block by its root mean square and concatenate
  Xs <- do.call(cbind, lapply(d$Yc, function(A) A / sqrt(mean(A^2))))
  sv <- svd(Xs, nu = K, nv = 0)
  L <- vapply(seq_len(K), function(k) {
    u <- sv$u[, k] * sqrt(sv$d[k])
    if (sum(pmax(u, 0)^2) < sum(pmax(-u, 0)^2)) u <- -u
    pmax(u, 0) + 1e-3
  }, numeric(d$n))
  L <- matrix(L, d$n, K)
  # 4. optional survival-score start for the last program
  if (init == "survival") {
    u0 <- .cox_moment(numeric(d$n), numeric(d$n), d$status, d$ri)$rho   # null martingale residuals
    z <- Reduce(`+`, lapply(d$Yc, function(A) A %*% crossprod(A, u0)))
    z <- drop(z) / stats::sd(drop(z))
    L[, K] <- pmax(z - min(z), 0) + 1e-3
  }
  Lc <- sweep(L, 2, colMeans(L), "-")
  G <- crossprod(Lc) + diag(1e-8, K)
  st <- list(mu = list(), tau = list(), EF = list(), EF2 = list(), R = list(),
             gF = list(), klF = list(), msgF = list())
  for (m in seq_len(d$M)) {
    Fm <- t(solve(G, crossprod(Lc, d$Yc[[m]])))            # p x K least squares
    if (d$prior_F[m] == "point_exponential") Fm <- pmax(Fm, 0)
    st$EF[[m]] <- Fm; st$EF2[[m]] <- Fm^2
    st$mu[[m]] <- d$center[[m]] - drop(Fm %*% colMeans(L))
    st$R[[m]] <- d$Y[[m]] - outer(rep(1, d$n), st$mu[[m]]) - L %*% t(Fm)
    st$tau[[m]] <- pmin(d$n / pmax(colSums(st$R[[m]]^2), 1e-12), d$tau_max[[m]])
    st$gF[[m]] <- vector("list", K); st$klF[[m]] <- numeric(K)
    st$msgF[[m]] <- rep(list(list(c = numeric(d$p[m]), A = numeric(d$p[m]))), K)
  }
  st$EL <- L; st$EL2 <- L^2; st$gL <- vector("list", K); st$klL <- numeric(K)
  st$b <- numeric(K); st$s2 <- numeric(K); st$gbeta <- NULL; st$klB <- 0
  st$msgB <- list(c = numeric(K), A = numeric(K))
  st$Zbar <- Reduce(`+`, lapply(seq_len(d$M), function(m) d$Yc[[m]] %*% st$EF[[m]]))
  st$V <- matrix(0, d$n, K)
  st$counters <- list(halvings = 0, rejections = c(L = 0, F = 0, beta = 0))
  st$obj <- list(rec = numeric(d$M), S = 0, total = -Inf)
  st
}

#' Effective number of independent features of a block
#'
#' From residuals R (n x p, centered), S = R^T R / n_df. The participation
#' ratio p_eff = (tr S)^2 / tr(S^2) equals p for independent features of equal
#' variance; with n_df < p the sample tr(S^2) is inflated by about
#' (tr S)^2 / n_df, which is subtracted (computed through the n x n Gram matrix).
#' @param R residual matrix
#' @param n_df residual degrees of freedom
#' @return p_eff (between 1 and p)
#' @keywords internal
#' @noRd
.effective_features <- function(R, n_df) {
  G <- tcrossprod(R) / n_df
  trS <- sum(diag(G)); trS2 <- sum(G^2) - trS^2 / n_df
  p <- ncol(R)
  if (trS2 <= 0) return(p)
  min(p, max(1, trS^2 / trS2))
}
