# Coordinate updates (derivation, Sections 4-5) ----
#
# Every update returns the state with `obj` refreshed. Updates that call
# ebnm (L, F, beta) are checked against the exact objective: a step that
# lowers it is shortened (F, beta: the survival message is damped) or
# rejected, so every accepted step increases the objective.

#' Accept a candidate state if the objective did not decrease
#' @return list(st, accepted)
#' @keywords internal
#' @noRd
.accept <- function(old, new) {
  new$obj$total >= old$obj$total - 1e-12 * abs(old$obj$total)
}

## Intercept and precisions ----

#' Update the feature intercepts mu_m (exact maximizer)
#' @keywords internal
#' @noRd
.update_mu <- function(st, d) {
  for (m in seq_len(d$M)) {
    Rplus <- sweep(st$R[[m]], 2, st$mu[[m]], "+")       # Y - E[L] E[F]^T
    st$mu[[m]] <- colMeans(Rplus)
    st$R[[m]] <- sweep(Rplus, 2, st$mu[[m]], "-")
  }
  .objective(st, d)
}

#' Update the precisions tau_m (exact maximizer, capped)
#' @keywords internal
#' @noRd
.update_tau <- function(st, d) {
  for (m in seq_len(d$M)) {
    ss <- pmax(.ss_block(st, m), 0)
    st$tau[[m]] <- pmin(d$n / ss, d$tau_max[[m]])
  }
  .objective(st, d)
}

## Activities L ----

#' Update q(L_.k) and g_L,k by EBNM on block-pooled pseudo-observations
#' @keywords internal
#' @noRd
.update_L <- function(st, d, k, safeguard = TRUE) {
  h <- numeric(d$n); a <- 0
  for (m in seq_len(d$M)) {
    tf <- st$tau[[m]] * st$EF[[m]][, k]
    h <- h + d$w_block[m] * (drop(st$R[[m]] %*% tf) + st$EL[, k] * sum(tf * st$EF[[m]][, k]))
    a <- a + d$w_block[m] * sum(st$tau[[m]] * st$EF2[[m]][, k])
  }
  fit <- .ebnm_fit(h / a, rep(1 / sqrt(a), d$n), d$prior_L, g_init = st$gL[[k]])
  new <- st
  dL <- fit$mean - st$EL[, k]
  for (m in seq_len(d$M)) new$R[[m]] <- st$R[[m]] - outer(dL, st$EF[[m]][, k])
  new$EL[, k] <- fit$mean; new$EL2[, k] <- fit$mean2
  new$gL[k] <- list(fit$g); new$klL[k] <- fit$kl
  new <- .objective(new, d)
  if (!safeguard || .accept(st, new)) return(new)
  st$counters$rejections[["L"]] <- st$counters$rejections[["L"]] + 1
  st
}

## Loadings F ----

#' Survival message for column k of block m (derivation, Section 4.2)
#'
#' A_j = (b_k^2 + s_k^2) sum_i P_i ytilde_ij^2,
#' c_j = b_k sum_i rho_i ytilde_ij - s_k^2 sum_i P_i Zbar_ik ytilde_ij + A_j fbar_jk.
#' @keywords internal
#' @noRd
.msg_F <- function(st, d, m, k) {
  e2 <- st$b[k]^2 + st$s2[k]
  if (d$w_surv == 0 || e2 == 0) return(list(c = numeric(d$p[m]), A = numeric(d$p[m])))
  cm <- .surv_state(st, d)
  A <- pmax(e2 * drop(crossprod(d$Yc2[[m]], cm$P)), 0)
  cc <- st$b[k] * drop(crossprod(d$Yc[[m]], cm$rho)) -
    st$s2[k] * drop(crossprod(d$Yc[[m]], cm$P * st$Zbar[, k])) + A * st$EF[[m]][, k]
  list(c = cc, A = A)
}

#' Update q(F_m,.k) and g_m,k: reconstruction (exact) plus survival message
#' @keywords internal
#' @noRd
.update_F <- function(st, d, m, k, safeguard = TRUE) {
  tau <- st$tau[[m]]; wm <- d$w_block[m]
  aR <- wm * tau * sum(st$EL2[, k])
  hR <- wm * tau * (drop(crossprod(st$R[[m]], st$EL[, k])) + st$EF[[m]][, k] * sum(st$EL[, k]^2))
  msg_new <- .msg_F(st, d, m, k)
  msg_old <- st$msgF[[m]][[k]]
  surv_on <- any(msg_new$A > 0) || any(msg_old$A > 0)
  step <- 1
  for (attempt in 0:d$max_halvings) {
    cc <- (1 - step) * msg_old$c + step * msg_new$c
    A <- (1 - step) * msg_old$A + step * msg_new$A
    prec <- aR + d$w_surv * A
    fit <- .ebnm_fit((hR + d$w_surv * cc) / prec, 1 / sqrt(prec), d$prior_F[m], g_init = st$gF[[m]][[k]])
    new <- .set_F(st, d, m, k, fit)
    new$msgF[[m]][[k]] <- list(c = cc, A = A)
    new$obj$rec[m] <- .recon_block(new, d, m)
    new$obj$S <- .surv_state(new, d)$S
    new$obj$total <- .total(new, d)
    if (!safeguard || .accept(st, new)) {
      new$counters$halvings <- new$counters$halvings + attempt
      return(new)
    }
    if (!surv_on) break
    step <- step / 2
  }
  st$counters$rejections[["F"]] <- st$counters$rejections[["F"]] + 1
  st
}

#' Write a new q(F_m,.k) into the state, updating residuals and predictor moments
#' @keywords internal
#' @noRd
.set_F <- function(st, d, m, k, fit) {
  dF <- fit$mean - st$EF[[m]][, k]
  dVar <- (fit$mean2 - fit$mean^2) - (st$EF2[[m]][, k] - st$EF[[m]][, k]^2)
  st$R[[m]] <- st$R[[m]] - outer(st$EL[, k], dF)
  st$Zbar[, k] <- st$Zbar[, k] + drop(d$Yc[[m]] %*% dF)
  st$V[, k] <- pmax(st$V[, k] + drop(d$Yc2[[m]] %*% dVar), 0)
  st$EF[[m]][, k] <- fit$mean; st$EF2[[m]][, k] <- fit$mean2
  st$gF[[m]][k] <- list(fit$g); st$klF[[m]][k] <- fit$kl
  st
}

## Survival coefficients beta ----

#' Survival message for beta (derivation, Section 4.4)
#'
#' A_k = sum_i P_i (Zbar_ik^2 + V_ik); c_k = sum_i rho_i Zbar_ik + b_k sum_i P_i Zbar_ik^2.
#' @keywords internal
#' @noRd
.msg_beta <- function(st, d) {
  cm <- .surv_state(st, d)
  A <- pmax(drop(crossprod(st$Zbar^2 + pmax(st$V, 0), cm$P)), 0)
  cc <- drop(crossprod(st$Zbar, cm$rho)) + st$b * drop(crossprod(st$Zbar^2, cm$P))
  list(c = cc, A = A)
}

#' Update q(beta) and sigma_beta^2 (shared Normal prior fitted by marginal likelihood)
#' @keywords internal
#' @noRd
.update_beta <- function(st, d, safeguard = TRUE) {
  if (d$w_surv == 0) return(st)
  for (inner in seq_len(d$beta_inner)) {
    msg_new <- .msg_beta(st, d)
    msg_old <- st$msgB
    step <- 1; accepted <- FALSE
    for (attempt in 0:d$max_halvings) {
      cc <- (1 - step) * msg_old$c + step * msg_new$c
      A <- (1 - step) * msg_old$A + step * msg_new$A
      fit <- .ebnm_fit(cc / A, 1 / sqrt(d$w_surv * A), "normal", g_init = st$gbeta)
      new <- st
      new$b <- fit$mean; new$s2 <- fit$mean2 - fit$mean^2
      new$gbeta <- fit$g; new$klB <- fit$kl
      new$msgB <- list(c = cc, A = A)
      new$obj$S <- .surv_state(new, d)$S
      new$obj$total <- .total(new, d)
      if (!safeguard || .accept(st, new)) { accepted <- TRUE; break }
      step <- step / 2
    }
    if (!accepted) { st$counters$rejections[["beta"]] <- st$counters$rejections[["beta"]] + 1; break }
    new$counters$halvings <- new$counters$halvings + attempt
    change <- max(abs(new$b - st$b)) / max(1e-8, max(abs(new$b)))
    st <- new
    if (change < 1e-6) break
  }
  st
}
