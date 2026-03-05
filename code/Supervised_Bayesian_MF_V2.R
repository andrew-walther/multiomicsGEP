# ------------------------------------------------------------------------------
# TITLE:       Supervised Bayesian Matrix Factorization — V2 (Corrected)
# AUTHOR:      Andrew Walther
# DATE:        March 2026  (Rebuilt from February 19, 2026)
# DERIVATIONS: derivations/MF_UpdateDerivations/MF_V2_Companion.tex
#
# DESCRIPTION:
#   Coordinate Ascent Variational Inference (CAVI) for the joint model:
#     Genomics:  Y_{n x p} = L_{n x K} F^T_{K x p} + E,  E_{ij} ~ N(0, tau_j^{-1})
#     Survival:  h(t_i | l_i) = h_0(t_i) exp( sum_k l_{ik} beta_k )   [Cox PH]
#
#   Each parameter update is derived as an EBNM problem (see Companion.tex):
#     EBNM(x = B/A, s = 1/sqrt(A), prior_family = "point_normal")
#   where A is the EBNM precision and B is the signal pseudo-numerator.
#
# VARIABLE CONVENTIONS (used throughout this script):
#   EL[i,k]   = E_q[l_{ik}]          posterior mean of loading i on factor k
#   EL2[i,k]  = E_q[l_{ik}^2]        posterior 2nd moment = Var_q(l_{ik}) + mean^2
#   EF[j,k]   = E_q[f_{jk}]          posterior mean of feature weight j on factor k
#   EF2[j,k]  = E_q[f_{jk}^2]        posterior 2nd moment
#   EBeta[k]  = E_q[beta_k]          posterior mean of survival coefficient k
#   EBeta2[k] = E_q[beta_k^2]        posterior 2nd moment
#   Tau[j]    = tau_j                 noise precision (feature-specific)
#   w[i]      = W_{ii}               diagonal Cox Hessian (per-sample weight)
#   z[i]      = working response      z_i = eta_i + u_i / W_{ii}
#
# KEY CHANGES FROM V1 (Supervised_Bayesian_MF.R):
#   [A1] z_no_k computed from CURRENT EL, EBeta inside the k-loop (true
#        Gauss-Seidel CAVI); same z_no_k reused for both L and beta updates.
#   [A2] Orthogonalisation is OFF by default (orthogonalize=FALSE).
#        When ON: EL2/EF2 reset to squared means (known approximation;
#        variance recovered next EBNM call).  See Companion.tex Sec. 10.
#   [A3] pmax() floors on all EBNM precision inputs (prevents degenerate solves).
#   [A4] Dual convergence: both mean|delta_L| and mean|delta_Beta| < tol.
#   [A5] ELBO proxy (genomics log-likelihood) tracked per iteration.
#   [A6] refresh_taylor flag: re-expand Cox Taylor at each factor k
#        (default FALSE = standard IRLS-within-VI).
# ------------------------------------------------------------------------------

DATA_MODE <- "simulated"

real_genomics_mat <- NULL
real_clinical_df  <- NULL
real_time_col     <- "time"
real_status_col   <- "status"

# ==============================================================================
# PART 1: LIBRARIES & HELPER FUNCTIONS
# ==============================================================================

library(survival)
library(ebnm)

# ------------------------------------------------------------------------------
#' Calculate Cox Score and Diagonal Hessian (Taylor Expansion)
#'
#' Transforms the non-conjugate Cox partial likelihood into a locally Gaussian
#' weighted-least-squares form centred at eta_hat = L_bar %*% beta_bar.
#'
#' DERIVATION: Companion.tex Sec. 3 (Eq. working-response).
#' The 2nd-order Taylor expansion gives:
#'   ell(eta) ~ sum_i [-1/2 * W_{ii} * (eta_i - z_i)^2] + C
#' where z_i = eta_hat_i + u_i / W_{ii} is the working response.
#'
#' @param eta    n-vector: current linear predictor eta_i = sum_k l_bar_{ik} beta_bar_k
#' @param time   n-vector: observed survival / censoring times
#' @param status n-vector: event indicator (1 = event, 0 = censored)
#' @return list(u = n-vector score, w = n-vector neg-diagonal Hessian)
# ------------------------------------------------------------------------------
calc_cox_taylor <- function(eta, time, status) {
  n   <- length(time)
  ord <- order(time)
  time_s   <- time[ord]
  status_s <- status[ord]
  eta_s    <- eta[ord]

  theta    <- exp(eta_s)                           # exp(eta_i) for each sample

  risk_sum <- rev(cumsum(rev(theta)))               # sum_{m: t_m >= t_i} exp(eta_m)

  h <- status_s / risk_sum                          # hazard increment delta_i / risk_sum_i
  H <- cumsum(h)                                    # cumulative hazard H_i

  # Score:       u_i = delta_i - exp(eta_i) * H_i
  # Neg-Hessian: W_{ii} = exp(eta_i) * H_i  (positive)
  u_s <- status_s - theta * H
  w_s <- theta * H
  w_s[w_s < 1e-6] <- 1e-6                          # floor to prevent 0-division

  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_s;   w[ord] <- w_s
  list(u = u, w = w)
}

# ------------------------------------------------------------------------------
#' C-Index Comparison: Supervised Latent L vs. PCA Top Components
# ------------------------------------------------------------------------------
get_cindex_comparison <- function(EL, data) {
  pca_y  <- prcomp(data$Y, rank. = min(5, ncol(EL)))
  fit_pc <- coxph(Surv(data$time, data$status) ~ pca_y$x)
  fit_l  <- coxph(Surv(data$time, data$status) ~ EL)
  list(
    c_original = round(summary(fit_pc)$concordance[1], 3),
    c_latent   = round(summary(fit_l)$concordance[1], 3)
  )
}

# ------------------------------------------------------------------------------
#' Extract Top Influential Features per Factor
# ------------------------------------------------------------------------------
get_top_features <- function(EF, n_top = 10) {
  lapply(1:ncol(EF), function(k) {
    weights   <- EF[, k]
    order_idx <- order(abs(weights), decreasing = TRUE)
    data.frame(
      FeatureID = order_idx[1:n_top],
      Weight    = round(weights[order_idx[1:n_top]], 4)
    )
  })
}

# ------------------------------------------------------------------------------
#' Generate Factor Summary Table (Log-rank p, Sparsity, PVE)
# ------------------------------------------------------------------------------
get_factor_summary_table <- function(res, data) {
  K <- ncol(res$L)
  p_vals <- sapply(1:K, function(k) {
    grp     <- ifelse(res$L[, k] > median(res$L[, k]), "High", "Low")
    sd_test <- survdiff(Surv(data$time, data$status) ~ grp)
    1 - pchisq(sd_test$chisq, 1)
  })
  total_var <- sum(data$Y^2)
  pve <- sapply(1:K, function(k) {
    sum((outer(res$L[, k], res$F[, k]))^2) / total_var * 100
  })
  sparsity <- colMeans(res$F != 0) * 100
  data.frame(
    Factor       = 1:K,
    Beta         = round(res$Beta, 3),
    LogRank_P    = round(p_vals, 4),
    Sparsity_Pct = round(sparsity, 2),
    PVE_Pct      = round(pve, 2)
  )
}

# ==============================================================================
# PART 2: THE FITTING ALGORITHM
# ==============================================================================

# ------------------------------------------------------------------------------
#' Fit Supervised Matrix Factorization via CAVI
#'
#' Maximises the Evidence Lower Bound (ELBO) over variational posteriors
#' q(L), q(F), q(beta) and priors g(L), g(F), g(beta) using Coordinate
#' Ascent Variational Inference.  Each coordinate update reduces to an
#' EBNM problem (Companion.tex Secs. 4-6).
#'
#' @param Y              n x p genomics data matrix
#' @param time           n-vector of survival / censoring times
#' @param status         n-vector of event indicators (1=event, 0=censored)
#' @param K              Number of latent factors
#' @param max_iter       Maximum CAVI outer iterations
#' @param tol            Convergence threshold: both mean|dL| and mean|dBeta| < tol
#' @param orthogonalize  Logical (default FALSE).  SVD-rotate L every 10 iters.
#'                       APPROXIMATION: resets EL2/EF2 (Companion.tex Sec. 10).
#' @param refresh_taylor Logical (default FALSE).  If TRUE, recompute (z, w) at
#'                       every factor k (true CAVI; slower).
#' @param verbose        Logical: print iteration logs?
#' @return list(L, F, Beta, Beta2, Tau, history)
# ------------------------------------------------------------------------------
fit_supervised_mf <- function(Y, time, status,
                               K              = 5,
                               max_iter       = 100,
                               tol            = 1e-5,
                               orthogonalize  = FALSE,
                               refresh_taylor = FALSE,
                               verbose        = TRUE) {

  n <- nrow(Y); p <- ncol(Y)

  # ---------------------------------------------------------------------------
  # Initialization  (Companion.tex Sec. 8 — Algorithm Box)
  # ---------------------------------------------------------------------------

  # SVD of Y: deterministic high-variance starting subspace.
  # EL = U sqrt(D),  EF = V sqrt(D)  so  EL %*% t(EF) = Y_rank-K
  svd_init <- svd(Y, nu = K, nv = K)
  d_k <- sqrt(pmax(svd_init$d[1:K], 0))
  EL  <- svd_init$u %*% diag(d_k, K, K)      # n x K
  EF  <- svd_init$v %*% diag(d_k, K, K)      # p x K

  # Second moments: initialised to squared means (zero posterior variance).
  # Posterior variance is populated after the first EBNM call.
  EL2 <- EL^2
  EF2 <- EF^2

  # Warm-start beta via Cox regression on initial loadings.
  df_cox <- as.data.frame(EL)
  colnames(df_cox) <- paste0("L", 1:K)
  df_cox$time   <- time
  df_cox$status <- status
  cox_init <- tryCatch(
    coxph(as.formula("Surv(time, status) ~ ."), data = df_cox, x = FALSE),
    error = function(e) NULL
  )
  if (!is.null(cox_init)) {
    cx_coef <- coef(cox_init)
    cx_coef[is.na(cx_coef)] <- 0
    EBeta <- cx_coef
  } else {
    EBeta <- rep(0, K)
  }
  EBeta2 <- EBeta^2      # zero variance initially; updated at first EBNM call

  # Column-specific precision tau_j (Companion.tex Eq. tau-colspec):
  # initialised from sample variance of each column of Y.
  col_var <- pmax(apply(Y, 2, var), 1e-8)
  Tau     <- 1.0 / col_var                     # p-vector

  # Tracking history
  history <- list(
    rmse       = numeric(max_iter),
    elbo_proxy = numeric(max_iter),
    converged  = FALSE,
    n_iter     = max_iter
  )

  if (verbose) {
    cat("=== Supervised Bayesian MF (V2) — CAVI ===\n")
    cat(sprintf("    n=%d, p=%d, K=%d | max_iter=%d | tol=%.1e\n\n",
                n, p, K, max_iter, tol))
  }

  # ===========================================================================
  # Main CAVI Loop
  # ===========================================================================
  for (iter in 1:max_iter) {

    EL_old    <- EL
    EBeta_old <- EBeta

    # -------------------------------------------------------------------------
    # STEP 1: Cox Taylor Expansion (Companion.tex Sec. 3)
    #
    # Linearise the Cox partial likelihood around the current linear predictor
    # eta_hat_i = sum_k l_bar_{ik} * beta_bar_k.
    #
    # Working response:  z_i  = eta_hat_i + u_i / W_{ii}
    # Weight:            W_{ii}  (negative diagonal Hessian, positive)
    #
    # (z, w) are held FIXED for this outer iteration unless refresh_taylor=TRUE.
    # Standard IRLS-within-VI; exactness restored at next outer iteration.
    # -------------------------------------------------------------------------
    eta    <- as.vector(EL %*% EBeta)              # n-vector: eta_hat_i
    taylor <- calc_cox_taylor(eta, time, status)
    z      <- eta + taylor$u / taylor$w            # n-vector: working response z_i
    w      <- taylor$w                             # n-vector: W_{ii}

    # Reconstruction RMSE at posterior means (monitoring only).
    history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))

    # =========================================================================
    # STEP 2: Factor-wise Coordinate Ascent  k = 1, ..., K
    # =========================================================================
    for (k in 1:K) {

      # [A6] Optional: re-expand Taylor at each k (true CAVI; slower).
      if (refresh_taylor) {
        eta    <- as.vector(EL %*% EBeta)
        taylor <- calc_cox_taylor(eta, time, status)
        z      <- eta + taylor$u / taylor$w
        w      <- taylor$w
      }

      # -----------------------------------------------------------------------
      # Partial Residuals (Companion.tex Sec. 4.1)
      #
      # R^{-k}_{ij} = Y_{ij} - sum_{k' != k} l_bar_{ik'} f_bar_{jk'}
      # z^{-k}_i    = z_i    - sum_{k' != k} l_bar_{ik'} beta_bar_{k'}
      #
      # [A1] Both use CURRENT EL, EF, EBeta (incorporates updates from k' < k
      #      in this iteration — true Gauss-Seidel CAVI).
      #
      # KEY POINT:  z_no_k does NOT depend on l_{ik} or beta_k — only on k'!=k.
      # Therefore the SAME z_no_k is reused for BOTH the L and beta updates
      # of factor k.  (Companion.tex Sec. 6, "z_no_k reuse rationale")
      # -----------------------------------------------------------------------

      # R^{-k}: remove factor k from the full reconstruction
      Y_hat <- EL %*% t(EF)                                 # n x p
      R_k   <- Y - Y_hat + outer(EL[, k], EF[, k])         # n x p: partial residual

      # z^{-k}: partial linear predictor excluding factor k  (Eq. z-residual)
      eta_no_k <- as.vector(EL %*% EBeta) - EL[, k] * EBeta[k]
      z_no_k   <- z - eta_no_k                               # n-vector

      # =======================================================================
      # (a) Update q_{l_k}: Patient Loadings (Companion.tex Sec. 4)
      #
      #  A_{ik} = sum_j tau_j * E_q[f^2_{jk}]  +  W_{ii} * E_q[beta_k^2]
      #         = [scalar genomics term]          + [n-vector survival term]
      #         -> n-vector (varies across samples i via W_{ii})
      #
      #  B_{ik} = sum_j tau_j * R^{-k}_{ij} * f_bar_{jk}  (genomics)
      #         + W_{ii} * z^{-k}_i * beta_bar_k            (survival)
      #         -> n-vector
      #
      #  EBNM: x_i = B_{ik}/A_{ik},  s_i = 1/sqrt(A_{ik})
      # =======================================================================
      A_L <- sum(Tau * EF2[, k]) + w * EBeta2[k]           # n-vector
      A_L <- pmax(A_L, 1e-10)                               # [A3] floor

      # B^{gen}_{ik} = sum_j tau_j R^{-k}_{ij} f_bar_{jk}
      #   Efficient: R_k %*% (Tau * EF[,k]) avoids n x p sweep
      B_L_gen  <- as.vector(R_k %*% (Tau * EF[, k]))        # n-vector
      B_L_surv <- w * z_no_k * EBeta[k]                     # n-vector
      B_L      <- B_L_gen + B_L_surv

      res_L    <- ebnm(x = B_L / A_L, s = 1 / sqrt(A_L), prior_family = "point_normal")
      EL[, k]  <- res_L$posterior$mean                       # l_bar_{ik}
      EL2[, k] <- res_L$posterior$sd^2 + res_L$posterior$mean^2  # E[l^2_{ik}]

      # =======================================================================
      # (b) Update q_{f_k}: Biological Factors (Companion.tex Sec. 5)
      #
      #  A_{jk} = tau_j * sum_i E_q[l^2_{ik}]     -> p-vector
      #  B_{jk} = tau_j * sum_i R^{-k}_{ij} * l_bar_{ik}  -> p-vector
      #
      #  Uses updated EL[,k] and EL2[,k] from step (a).
      #  R_k is still valid (it excludes factor k, which is unchanged for k'!=k).
      #
      #  EBNM: x_j = B_{jk}/A_{jk},  s_j = 1/sqrt(A_{jk})
      # =======================================================================
      sum_EL2_k <- sum(EL2[, k])                              # scalar
      A_F <- pmax(Tau * sum_EL2_k, 1e-10)                     # p-vector [A3]
      B_F <- Tau * as.vector(t(R_k) %*% EL[, k])              # p-vector

      res_F    <- ebnm(x = B_F / A_F, s = 1 / sqrt(A_F), prior_family = "point_normal")
      EF[, k]  <- res_F$posterior$mean
      EF2[, k] <- res_F$posterior$sd^2 + res_F$posterior$mean^2

      # =======================================================================
      # (c) Update q_{beta_k}: Survival Coefficients (Companion.tex Sec. 6)
      #
      #  A_k = sum_i W_{ii} * E_q[l^2_{ik}]        -> scalar
      #        Error-in-variables: uses E[l^2] not l_bar^2.
      #        Prevents beta from overfitting to uncertain loadings.
      #
      #  B_k = sum_i W_{ii} * z^{-k}_i * l_bar_{ik} -> scalar
      #        Uses updated EL[,k] from step (a).
      #
      #  z_no_k is REUSED from above (does not depend on l_{ik} or beta_k).
      #
      #  EBNM: x_k = B_k/A_k,  s_k = 1/sqrt(A_k)
      # =======================================================================
      A_Beta <- max(sum(w * EL2[, k]), 1e-10)                 # scalar [A3]
      B_Beta <- sum(w * z_no_k * EL[, k])                     # scalar

      res_Beta  <- ebnm(x = B_Beta / A_Beta, s = 1 / sqrt(A_Beta), prior_family = "point_normal")
      EBeta[k]  <- res_Beta$posterior$mean
      EBeta2[k] <- res_Beta$posterior$sd^2 + res_Beta$posterior$mean^2

    }  # end k-loop

    # =========================================================================
    # STEP 3: Precision Update tau (Companion.tex Sec. 7)
    #
    #  Expected squared residual:
    #   R_bar^2_{ij} = (Y_{ij} - sum_k l_bar f_bar)^2
    #                + sum_k [ E[l^2_{ik}]*E[f^2_{jk}] - l_bar^2_{ik}*f_bar^2_{jk} ]
    #
    #  Column-specific MLE:  tau_hat_j = n / sum_i R_bar^2_{ij}
    # =========================================================================
    Var_Term <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))      # n x p: variance correction
    R2_bar   <- (Y - EL %*% t(EF))^2 + Var_Term              # n x p: E[residual^2]
    Tau      <- n / pmax(colSums(R2_bar), n * 1e-8)          # p-vector [A3]

    # [A5] ELBO proxy: genomics log-likelihood term  E_q[log P(Y | L, F, tau)]
    history$elbo_proxy[iter] <- sum(n / 2 * log(Tau) - Tau / 2 * colSums(R2_bar))

    # =========================================================================
    # STEP 4 (optional): Orthogonalization  [A2]
    #
    # SVD rotation of EL every 10 iterations for factor identifiability.
    # APPROXIMATION: EL2/EF2 reset to EL^2/EF^2 (drops posterior variance;
    # recovered at next EBNM call).  See Companion.tex Sec. 10.
    # Default OFF: Point-Normal EBNM already promotes sparsity/distinctness.
    # =========================================================================
    if (orthogonalize && iter %% 10 == 0) {
      svd_rot <- svd(EL, nu = K, nv = K)
      EL  <- svd_rot$u %*% diag(svd_rot$d, K, K)
      EF  <- EF %*% svd_rot$v
      EL2 <- EL^2
      EF2 <- EF^2
    }

    # =========================================================================
    # Convergence Check  [A4]
    # =========================================================================
    delta_L    <- mean(abs(EL - EL_old))
    delta_Beta <- mean(abs(EBeta - EBeta_old))

    if (verbose && iter %% 10 == 0) {
      cat(sprintf("  iter %3d | RMSE: %.4f | ELBO: %+.1f | dL: %.2e | dB: %.2e | beta: [%s]\n",
                  iter, history$rmse[iter], history$elbo_proxy[iter],
                  delta_L, delta_Beta,
                  paste(sprintf("%+.2f", EBeta), collapse = ", ")))
    }

    if (iter > 5 && delta_L < tol && delta_Beta < tol) {
      if (verbose) {
        cat(sprintf("\n  Converged at iteration %d  (dL=%.2e, dBeta=%.2e)\n",
                    iter, delta_L, delta_Beta))
      }
      history$converged  <- TRUE
      history$n_iter     <- iter
      history$rmse       <- history$rmse[1:iter]
      history$elbo_proxy <- history$elbo_proxy[1:iter]
      break
    }
  }  # end CAVI loop

  list(L = EL, F = EF, Beta = EBeta, Beta2 = EBeta2, Tau = Tau, history = history)
}

# ==============================================================================
# PART 3: VISUALIZATION & ANALYTICS
# ==============================================================================

#' GEP Heatmap with Colour Legend
plot_gep_heatmap <- function(res, n_features = 50) {
  top_var_genes <- order(rowSums(abs(res$F)), decreasing = TRUE)[1:n_features]
  F_sub   <- res$F[top_var_genes, ]
  palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_val <- max(abs(F_sub))

  layout(matrix(1:2, ncol = 2), widths = c(5, 1))

  par(mar = c(6, 4, 4, 1))
  image(1:nrow(F_sub), 1:ncol(F_sub), F_sub,
        main = paste("GEP Feature Weights (Top", n_features, "IDs)"),
        xlab = "Feature ID (Row Index)", ylab = "Latent Factors",
        col = palette, axes = FALSE, zlim = c(-max_val, max_val))
  axis(1, at = 1:nrow(F_sub), labels = top_var_genes, las = 2, cex.axis = 0.6)
  axis(2, at = 1:ncol(F_sub), labels = paste0("F", 1:ncol(F_sub)), las = 1)
  box()

  par(mar = c(6, 1, 4, 3))
  legend_image <- as.matrix(seq(-max_val, max_val, length.out = 100))
  image(1, seq(-max_val, max_val, length.out = 100), t(legend_image),
        col = palette, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1, cex.axis = 0.8)
  mtext("Weight", side = 4, line = 2, cex = 0.8)
  layout(1)
}

#' Diagnostic Dashboard (renders sequential plots into RStudio Plots history)
visualize_dashboard <- function(res, data) {
  std_mar <- c(5, 5, 4, 2)

  # --- RMSE Trace ---
  par(mfrow = c(1, 1), mar = std_mar)
  rmse_vals <- res$history$rmse[res$history$rmse > 0]
  plot(rmse_vals, type = "l", lwd = 2, col = "#1f77b4",
       main = "Reconstruction RMSE", xlab = "Iteration", ylab = "RMSE", bty = "n")
  abline(h = 1.0, col = "red", lty = 2)
  grid(col = "lightgray", lty = "dotted")

  # --- ELBO Proxy Trace ---
  par(mfrow = c(1, 1), mar = std_mar)
  elbo_vals <- res$history$elbo_proxy[res$history$elbo_proxy != 0]
  if (length(elbo_vals) > 1) {
    plot(elbo_vals, type = "l", lwd = 2, col = "#2ca02c",
         main = "Genomics ELBO Proxy (should increase)", xlab = "Iteration",
         ylab = "E[log P(Y | L,F,tau)]", bty = "n")
    grid(col = "lightgray", lty = "dotted")
  }

  # --- GEP Heatmap ---
  plot_gep_heatmap(res)

  # --- Kaplan-Meier per Factor ---
  K <- ncol(res$L)
  summary_tab <- get_factor_summary_table(res, data)
  par(mfrow = c(2, ceiling(K / 2)), mar = c(4, 4, 3, 1))
  for (k in 1:K) {
    groups  <- ifelse(res$L[, k] > median(res$L[, k]), "High Score", "Low Score")
    km      <- survfit(Surv(data$time, data$status) ~ groups)
    p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
    plot(km, col = c("#d62728", "#1f77b4"), lwd = 2, bty = "n",
         main = paste("Factor", k, "\n", p_label),
         xlab = "Time", ylab = "Survival Probability")
  }

  # --- Signal Recovery (simulated data only) ---
  if (!is.null(data$L_true)) {
    par(mfrow = c(1, 1), mar = std_mar)
    cors        <- cor(data$L_true, res$L)
    best_match  <- apply(abs(cors), 2, which.max)
    target_est  <- which.max(abs(res$Beta))
    target_true <- best_match[target_est]
    sign_corr   <- sign(cors[target_true, target_est])
    plot(data$L_true[, target_true], res$L[, target_est] * sign_corr,
         main  = paste("Signal Recovery  (Est", target_est, "vs True", target_true, ")"),
         xlab  = "Ground Truth Loading", ylab = "Estimated (sign corrected)",
         col   = rgb(0, 0, 0, 0.4), pch = 16, bty = "n")
    abline(0, 1, col = "#d62728", lwd = 2, lty = 2)
  }
}

# ==============================================================================
# PART 4: SIMULATION & EXECUTION
# ==============================================================================

if (DATA_MODE == "simulated") {

  # Synthetic benchmark (Companion.tex Sec. 9)
  # n=250 patients, p=1000 features, K=5 factors (4 prognostic + 1 structural)
  sim_data_fn <- function(n = 250, p = 1000, k = 5) {
    set.seed(123)

    L     <- matrix(rnorm(n * k), n, k)
    F_mat <- matrix(0, p, k)
    for (i in 1:k) {
      active <- sample(1:p, round(p * 0.05))
      F_mat[active, i] <- rnorm(length(active), 0, 5)
    }

    Y <- L %*% t(F_mat) + matrix(rnorm(n * p), n, p)   # Y = LF' + E

    Beta      <- c(1.5, -1.2, 0.8, -0.5, 0.0)
    eta_true  <- as.vector(L %*% Beta)
    raw_times <- (-log(runif(n)) / (0.01 * exp(eta_true)))^(1 / 1.5)
    cens_times <- rexp(n, rate = 1 / 50)

    list(
      Y      = Y,
      time   = pmin(raw_times, cens_times),
      status = as.integer(raw_times <= cens_times),
      L_true = L,
      F_true = F_mat,
      B_true = Beta
    )
  }

  data <- sim_data_fn()

  cat("=== Running Supervised MF V2 on Simulated Data ===\n")
  cat(sprintf("Censoring rate: %.1f%%\n\n", 100 * mean(data$status == 0)))

  res <- fit_supervised_mf(
    Y = data$Y, time = data$time, status = data$status, K = 5,
    max_iter = 100, tol = 1e-5,
    orthogonalize = FALSE, refresh_taylor = FALSE, verbose = TRUE
  )

  cat("\n=== FACTOR SUMMARY TABLE ===\n")
  print(get_factor_summary_table(res, data))

  cat("\n=== PROPORTIONAL HAZARDS (PH) TEST ===\n")
  print(cox.zph(coxph(Surv(data$time, data$status) ~ res$L)))

  cat("\n=== MODEL PERFORMANCE (C-INDEX) ===\n")
  perf <- get_cindex_comparison(res$L, data)
  cat(sprintf("  Top-5 PCA  C-index: %.3f\n", perf$c_original))
  cat(sprintf("  Supervised C-index: %.3f\n", perf$c_latent))

  cat("\n=== ESTIMATED vs TRUE BETA ===\n")
  print(data.frame(
    Factor    = 1:5,
    Beta_true = data$B_true,
    Beta_est  = round(res$Beta, 3),
    Beta2_est = round(res$Beta2, 4)
  ))

  cat("\n=== TOP 5 FEATURES PER GEP ===\n")
  top_feats <- get_top_features(res$F, 5)
  for (k in seq_along(top_feats)) {
    cat(sprintf("GEP %d:\n", k)); print(top_feats[[k]])
  }

  visualize_dashboard(res, data)
}
