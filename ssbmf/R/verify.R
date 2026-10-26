# Numerical stationarity check of the updates (used by verify_derivation.R) ----

#' Objective after replacing q(F_m,jk) for selected j by the posterior for perturbed pseudo-data
#'
#' q(F_m,.k) is the ebnm posterior for pseudo-data (x_j, s_j) under prior g;
#' the variational parameters are therefore (x_j, log s_j) with g fixed. This
#' helper recomputes the posterior of the selected entries under the fixed g
#' and returns the exact objective.
#' @keywords internal
#' @noRd
.obj_perturb_F <- function(st, d, m, k, j, x, s) {
  xs <- st$pseudoF; xx <- xs$x; ss <- xs$s
  xx[j] <- x; ss[j] <- s
  fit <- .ebnm_fit(xx, ss, d$prior_F[m], g_init = st$gF[[m]][[k]], fix_g = TRUE)
  new <- .set_F(st, d, m, k, fit)
  .objective(new, d)$obj$total
}

#' Stationarity of the F, L and beta updates by central finite differences
#'
#' Fits a small two-block simulated data set, then for one column of F (each
#' block), one column of L and the beta vector, iterates the update to its
#' fixed point and evaluates the finite-difference gradient of the exact
#' objective in the variational parameters of a few entries. Gradients are
#' reported per unit of the entry's posterior standard deviation. For
#' contrast, the same gradients are reported at the state before the update.
#'
#' @param seed seed for the simulated data.
#' @param verbose print a table.
#' @return data.frame with `update`, `parameter`, `grad_before`, `grad_after`.
#' @keywords internal
.verify_stationarity <- function(seed = 1, verbose = TRUE) {
  sim <- ssbmf_simulate(n = 60, p = c(a = 30, b = 60), K = 3, seed = seed)
  data <- ssbmf_data(sim$X, sim$time, sim$status)
  ctl <- ssbmf_control(max_sweeps = 15)
  d <- .prep_data(data, ctl)
  st <- .init_state(d, 3)
  st <- .update_mu(st, d); st <- .update_tau(st, d)
  for (k in 1:3) st <- .update_L(st, d, k, safeguard = FALSE)
  for (m in 1:2) for (k in 1:3) st <- .update_F(st, d, m, k, safeguard = FALSE)
  st <- .objective(st, d)
  for (it in 1:15) {
    st <- .update_mu(st, d); st <- .update_tau(st, d)
    for (k in 1:3) st <- .update_L(st, d, k)
    for (m in 1:2) for (k in 1:3) st <- .update_F(st, d, m, k)
    st <- .update_beta(st, d)
  }
  out <- list()
  h_rel <- 1e-3

  # F columns: iterate the message update to its fixed point
  pseudo_F <- function(st, m, k) {
    tau <- st$tau[[m]]
    aR <- tau * sum(st$EL2[, k])
    hR <- tau * (drop(crossprod(st$R[[m]], st$EL[, k])) + st$EF[[m]][, k] * sum(st$EL[, k]^2))
    msg <- .msg_F(st, d, m, k)
    prec <- aR + d$w_surv * msg$A
    list(x = (hR + d$w_surv * msg$c) / prec, s = 1 / sqrt(prec))
  }
  fd_F <- function(st, m, k, js) {
    st$pseudoF <- pseudo_F(st, m, k)
    # re-express the current q as the posterior of its own pseudo-data under g
    fit <- .ebnm_fit(st$pseudoF$x, st$pseudoF$s, d$prior_F[m], g_init = st$gF[[m]][[k]], fix_g = TRUE)
    st <- .objective(.set_F(st, d, m, k, fit), d)
    vapply(js, function(j) {
      x0 <- st$pseudoF$x[j]; s0 <- st$pseudoF$s[j]; h <- h_rel * s0
      (.obj_perturb_F(st, d, m, k, j, x0 + h, s0) - .obj_perturb_F(st, d, m, k, j, x0 - h, s0)) / (2 * h) * s0
    }, 0)
  }
  for (m in 1:2) {
    k <- 1
    js <- order(-abs(st$EF[[m]][, k]))[1:4]
    before <- NA_real_
    st2 <- st
    for (it in 1:50) {
      prev <- st2$EF[[m]][, k]
      st2 <- .update_F(st2, d, m, k)
      if (max(abs(st2$EF[[m]][, k] - prev)) < 1e-12 * max(1, max(abs(prev)))) break
    }
    after <- fd_F(st2, m, k, js)
    out[[length(out) + 1]] <- data.frame(update = sprintf("F[%s][,%d]", d$blocks[m], k),
                                         parameter = sprintf("x_%d", js), grad_before = before, grad_after = after)
  }

  # beta: Gaussian q, parameters (b_k, log s_k) with sigma_beta^2 fixed
  kl_beta <- function(b, s2, sig2) sum(0.5 * (s2 / sig2 + b^2 / sig2 - 1 - log(s2 / sig2)))
  obj_beta <- function(st, b, s2) {
    st$b <- b; st$s2 <- s2
    st$klB <- kl_beta(b, s2, st$gbeta$sd[1]^2)
    st$obj$S <- .surv_state(st, d)$S
    .total(st, d)
  }
  fd_beta <- function(st) {
    unlist(lapply(seq_along(st$b), function(k) {
      h <- h_rel * sqrt(st$s2[k]); bp <- bm <- st$b; bp[k] <- bp[k] + h; bm[k] <- bm[k] - h
      gb <- (obj_beta(st, bp, st$s2) - obj_beta(st, bm, st$s2)) / (2 * h) * sqrt(st$s2[k])
      lp <- lm <- st$s2; lp[k] <- lp[k] * exp(2 * h_rel); lm[k] <- lm[k] * exp(-2 * h_rel)
      gs <- (obj_beta(st, st$b, lp) - obj_beta(st, st$b, lm)) / (2 * h_rel)
      c(gb, gs)
    }))
  }
  stb <- st; stb$b <- stb$b * 0.8                        # a non-stationary state for contrast
  before <- fd_beta(stb)
  # with sigma_beta^2 fixed, the stationary q(beta) is reached by the message iteration
  for (it in 1:200) {
    msg <- .msg_beta(st, d)
    sig2 <- st$gbeta$sd[1]^2
    prec <- 1 / sig2 + d$w_surv * msg$A
    b_new <- d$w_surv * msg$c / prec
    done <- max(abs(b_new - st$b)) < 1e-13
    st$b <- b_new; st$s2 <- 1 / prec
    if (done) break
  }
  after <- fd_beta(st)
  out[[length(out) + 1]] <- data.frame(update = "beta", parameter = paste0(rep(c("b_", "log s_"), 3), rep(1:3, each = 2)),
                                       grad_before = before, grad_after = after)

  # L column: exact update; parameters are pseudo-data with g fixed
  k <- 1
  L_pseudo <- function(st) {
    h <- numeric(d$n); a <- 0
    for (m in seq_len(d$M)) {
      tf <- st$tau[[m]] * st$EF[[m]][, k]
      h <- h + drop(st$R[[m]] %*% tf) + st$EL[, k] * sum(tf * st$EF[[m]][, k])
      a <- a + sum(st$tau[[m]] * st$EF2[[m]][, k])
    }
    list(x = h / a, s = rep(1 / sqrt(a), d$n))
  }
  obj_L <- function(st, ps) {
    fit <- .ebnm_fit(ps$x, ps$s, d$prior_L, g_init = st$gL[[k]], fix_g = TRUE)
    dL <- fit$mean - st$EL[, k]
    for (m in seq_len(d$M)) st$R[[m]] <- st$R[[m]] - outer(dL, st$EF[[m]][, k])
    st$EL[, k] <- fit$mean; st$EL2[, k] <- fit$mean2; st$klL[k] <- fit$kl
    .objective(st, d)$obj$total
  }
  stL <- .update_L(st, d, k)
  ps <- L_pseudo(stL)
  is_ <- order(-stL$EL[, k])[1:4]
  after <- vapply(is_, function(i) {
    h <- h_rel * ps$s[i]; p1 <- p2 <- ps; p1$x[i] <- p1$x[i] + h; p2$x[i] <- p2$x[i] - h
    (obj_L(stL, p1) - obj_L(stL, p2)) / (2 * h) * ps$s[i]
  }, 0)
  out[[length(out) + 1]] <- data.frame(update = "L[,1]", parameter = sprintf("x_%d", is_),
                                       grad_before = NA_real_, grad_after = after)
  res <- do.call(rbind, out)
  res$max_grad <- abs(res$grad_after)
  if (verbose) print(format(res, digits = 3), row.names = FALSE)
  res
}
