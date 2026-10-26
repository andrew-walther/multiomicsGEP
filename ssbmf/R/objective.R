# The objective and its components ----
#
#   L(q) = sum_m w_m R_m + w_S S - sum_k KL_L,k - sum_{m,k} KL_m,k - KL_beta
#
# R_m (exact expected Gaussian log-likelihood of block m):
#   R_m = sum_j [ n/2 log(tau_mj / 2 pi) - tau_mj/2 * SS_mj ],
#   SS_mj = sum_i r_ij^2 + sum_k ( sum_i E L_ik^2 E F_jk^2 - sum_i (E L_ik)^2 (E F_jk)^2 )
# S: the Gaussian-moment (Jensen) bound on E_q log PL (derivation, Section 3).
# The KL terms are stored when their column is updated, always against the
# prior fitted in that same update.

#' Expected residual sum of squares by feature for block m
#' @keywords internal
#' @noRd
.ss_block <- function(st, m) {
  colSums(st$R[[m]]^2) +
    drop(st$EF2[[m]] %*% colSums(st$EL2)) - drop(st$EF[[m]]^2 %*% colSums(st$EL^2))
}

#' Reconstruction term R_m of block m
#' @keywords internal
#' @noRd
.recon_block <- function(st, d, m) {
  ss <- .ss_block(st, m)
  sum(0.5 * d$n * log(st$tau[[m]] / (2 * pi)) - 0.5 * st$tau[[m]] * ss)
}

#' Predictor moments under q: eta_bar = Zbar b, v = Zbar^2 s2 + V (b^2 + s2)
#' @keywords internal
#' @noRd
.eta_moments <- function(st) {
  list(eta_bar = drop(st$Zbar %*% st$b),
       v = drop(st$Zbar^2 %*% st$s2 + st$V %*% (st$b^2 + st$s2)))
}

#' Moment-bound Cox term S at the current state
#' @keywords internal
#' @noRd
.surv_state <- function(st, d) {
  em <- .eta_moments(st)
  .cox_moment(em$eta_bar, em$v, d$status, d$ri)
}

#' Total objective from cached components
#' @keywords internal
#' @noRd
.total <- function(st, d) {
  sum(d$w_block * st$obj$rec) + d$w_surv * st$obj$S -
    sum(st$klL) - sum(unlist(st$klF)) - st$klB
}

#' Recompute every component of the objective from scratch
#'
#' @param st state list
#' @param d prepared data (from `.prep_data()`)
#' @return the state with `obj` refreshed; `obj$total` is the objective
#' @keywords internal
.objective <- function(st, d) {
  st$obj$rec <- vapply(seq_len(d$M), function(m) .recon_block(st, d, m), 0)
  st$obj$S <- .surv_state(st, d)$S
  st$obj$total <- .total(st, d)
  st
}
