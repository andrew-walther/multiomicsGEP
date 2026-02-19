# ------------------------------------------------------------------------------
# TITLE:       Supervised Bayesian Matrix Factorization (Consolidated Suite)
# AUTHOR:      Andrew Walther
# DATE:        February 19, 2026
# DOCUMENTATION: https://gemini.google.com/share/a29cae54a148
#
# DESCRIPTION:
#   This consolidated script implements a unified probabilistic framework to 
#   decompose high-dimensional genomics data (Y) while simultaneously modeling 
#   survival outcomes (t, delta). The code uses Variational Inference to find 
#   a shared latent space (L) that minimizes reconstruction error and hazard 
#   residuals.It includes a full discovery suite with 
#   PVE calculation, PH testing, Log-rank statistics, and Diverging Heatmaps.
# ------------------------------------------------------------------------------

# --- [USER TOGGLE] ---
# Sets the data source. "simulated" generates a synthetic benchmark; 
# "real" is used for applying the model to actual datasets.
DATA_MODE <- "simulated"  

# --- [REAL DATA INPUTS] ---
# Placeholders for user-supplied data objects when DATA_MODE is set to "real".
real_genomics_mat <- NULL   
real_clinical_df  <- NULL   
real_time_col     <- "time" 
real_status_col   <- "status" 

# ==============================================================================
# PART 1: LIBRARIES & HELPER FUNCTIONS
# ==============================================================================

# 'survival' provides the Cox Proportional Hazards math and PH testing (cox.zph).
# It is necessary for modeling time-to-event outcomes.
library(survival)

# 'ebnm' provides Empirical Bayes Normal Means solvers.
# It is the engine that estimates sparse priors (Point-Normal) for the latent variables.
library(ebnm) 

#' Calculate Cox Gradients and Hessian (Diagonal Approximation)
#' 
#' This function transforms the non-conjugate Cox survival likelihood into 
#' a locally linear "Weighted Least Squares" form via Taylor expansion.
#' This "bridge" allows survival residuals to exert a pull on the latent loadings.
calc_cox_taylor <- function(eta, time, status) {
  n <- length(time)
  # Sorting is mandatory for calculating risk sets in Cox partial likelihoods.
  ord <- order(time) 
  time_sorted <- time[ord]; status_sorted <- status[ord]; eta_sorted <- eta[ord]
  
  theta <- exp(eta_sorted) # Relative risk
  risk_sum <- rev(cumsum(rev(theta))) # Risk set denominator
  
  h <- status_sorted / risk_sum # Hazard increments
  H <- cumsum(h) # Cumulative hazard
  
  # Gradient (u) and Diagonal Hessian (w) for the quadratic approximation.
  u_sorted <- status_sorted - theta * H
  w_sorted <- theta * H 
  
  # A numerical floor prevents division-by-zero for censored or low-information patients.
  w_sorted[w_sorted < 1e-6] <- 1e-6 
  
  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_sorted; w[ord] <- w_sorted
  return(list(u = u, w = w))
}

#' Get C-Index Performance Comparison
#' 
#' Calculates Harrell's Concordance Index (C-index) to evaluate the predictive 
#' power of the latent space compared to standard Principal Components (PCA).
get_cindex_comparison <- function(EL, data) {
  # Baseline: Top 5 Principal Components of the original genomics matrix.
  pca_y <- prcomp(data$Y, rank. = min(5, ncol(EL)))
  fit_pc <- coxph(Surv(data$time, data$status) ~ pca_y$x)
  
  # Comparison: Latent Loadings from our Supervised Model.
  fit_l  <- coxph(Surv(data$time, data$status) ~ EL)
  
  return(list(
    c_original = round(summary(fit_pc)$concordance[1], 3),
    c_latent   = round(summary(fit_l)$concordance[1], 3)
  ))
}

#' Extract Top Influential Features per Factor
#' 
#' Ranks features by their weight in the Factor matrix (F). This identifies 
#' the genomic loci (genes/CpGs) that form the regulatory backbone of each program.
get_top_features <- function(EF, n_top = 10) {
  lapply(1:ncol(EF), function(k) {
    weights <- EF[,k]
    order_idx <- order(abs(weights), decreasing = TRUE)
    data.frame(
      FeatureID = order_idx[1:n_top], 
      Weight = round(weights[order_idx[1:n_top]], 4)
    )
  })
}

#' Generate a Comprehensive Factor Summary Table
#' 
#' Aggregates Log-Rank P-values, biological sparsity, and Proportion 
#' of Variance Explained (PVE) into a single summary view for the user.
get_factor_summary_table <- function(res, data) {
  K <- ncol(res$L)
  # Calculate per-factor survival significance via median stratification.
  p_vals <- lapply(1:K, function(k) {
    group <- ifelse(res$L[,k] > median(res$L[,k]), "High", "Low")
    sd_test <- survdiff(Surv(data$time, data$status) ~ group)
    1 - pchisq(sd_test$chisq, length(sd_test$n) - 1)
  })
  # Calculate sparsity: percentage of genomic features active in each factor.
  sparsity <- colMeans(res$F != 0) * 100
  # Calculate genomic PVE: how much variance in Y each factor explains.
  total_var <- sum(data$Y^2)
  pve <- sapply(1:K, function(k) {
    sum((res$L[,k] %*% t(res$F[,k]))^2) / total_var * 100
  })
  
  data.frame(
    Factor = 1:K, 
    Beta = round(res$Beta, 3), 
    LogRank_P = round(unlist(p_vals), 4),
    Sparsity_Pct = round(sparsity, 2), 
    PVE_Pct = round(pve, 2)
  )
}

# ==============================================================================
# PART 2: THE FITTING ALGORITHM (THE ENGINE)
# ==============================================================================

#' Fit Supervised Matrix Factorization
#' 
#' The master loop using Coordinate Ascent Variational Inference (CAVI).
#' It refines L, F, Beta, and Tau to maximize the Evidence Lower Bound (ELBO).
fit_supervised_mf <- function(Y, time, status, K = 5, max_iter = 100, tol = 1e-5) {
  n <- nrow(Y); p <- ncol(Y)
  
  # Initialization via SVD provides a deterministic, high-variance starting subspace.
  # This is necessary because matrix factorization is non-convex and sensitive to start values.
  svd_init <- svd(Y, nu = K, nv = K)
  EL <- svd_init$u %*% diag(sqrt(svd_init$d[1:K]), K, K) 
  EF <- svd_init$v %*% diag(sqrt(svd_init$d[1:K]), K, K)
  EL2 <- EL^2; EF2 <- EF^2
  
  # Warm-start Beta via Cox on initial loadings to ensure early clinical signal recognition.
  df_surv <- data.frame(time = time, status = status, EL)
  colnames(df_surv)[3:(2+K)] <- paste0("L", 1:K)
  cox_fit <- coxph(as.formula(paste("Surv(time, status) ~ .")), data = df_surv)
  EBeta <- coef(cox_fit); EBeta[is.na(EBeta)] <- 0
  EBeta2 <- EBeta^2
  
  # Tau represents feature-specific noise precision, initialized to unit variance.
  Tau <- rep(1, p) 
  history <- list(rmse = numeric(max_iter))
  
  cat("Starting Supervised Inference Engine...\n")
  
  for(iter in 1:max_iter) {
    EL_old <- EL
    
    # 1. Survival Linearization: Translating survival residuals into the L update.
    # This Taylor expansion creates the "bridge" between the two likelihoods.
    eta <- EL %*% EBeta
    taylor <- calc_cox_taylor(as.vector(eta), time, status)
    z <- eta + taylor$u/taylor$w 
    w <- taylor$w
    
    # Track the noise floor. RMSE ~ 1.0 is expected in this simulation (Noise SD=1).
    history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))
    
    # 2. Sequential Factor Updates (Coordinate Ascent)
    for(k in 1:K) {
      
      # --- Update Patient Loadings ***(L)***: Data Fusion ---
      # Precision A_L combines information from genomics (F^2) and survival relevance (Beta^2).
      A_L <- sum(Tau * EF2[,k]) + w * EBeta2[k]
      Y_hat <- EL %*% t(EF)
      # Subtract all factors except k to isolate the current signal.
      R_k <- Y - Y_hat + outer(EL[,k], EF[,k]) 
      B_L_gen <- as.vector(sweep(R_k, 2, Tau, `*`) %*% EF[,k])
      
      z_no_k <- z - (eta - (EL[,k] * EBeta[k]))
      B_L_surv <- w * z_no_k * EBeta[k]
      
      # ebnm finds the sparse posterior mean using Point-Normal priors.
      res_L <- ebnm(x = (B_L_gen + B_L_surv)/A_L, s = 1/sqrt(A_L), prior_family = "point_normal")
      EL[,k] <- res_L$posterior$mean
      EL2[,k] <- res_L$posterior$sd^2 + res_L$posterior$mean^2 # Incorporate variance for uncertainty.
      
      # --- Update Biological Signatures ***(F)*** ---
      # Factors define the identity of the program and are updated via genomics signal.
      A_F <- Tau * sum(EL2[,k])
      B_F <- Tau * as.vector(t(R_k) %*% EL[,k])
      res_F <- ebnm(x = B_F/A_F, s = 1/sqrt(A_F), prior_family = "point_normal")
      EF[,k] <- res_F$posterior$mean
      EF2[,k] <- res_F$posterior$sd^2 + res_F$posterior$mean^2
      
      # --- Update Survival Coefficients ***(Beta)*** ---
      # Using EL2 implements the "Uncertainty Penalty" (Error-in-Variables adjustment).
      A_Beta <- sum(w * EL2[,k]) 
      z_no_k_beta <- z - (EL %*% EBeta - (EL[,k] * EBeta[k]))
      B_Beta <- sum(w * z_no_k_beta * EL[,k])
      res_Beta <- ebnm(x = B_Beta/A_Beta, s = 1/sqrt(A_Beta), prior_family = "point_normal")
      EBeta[k] <- res_Beta$posterior$mean
      EBeta2[k] <- res_Beta$posterior$sd^2 + res_Beta$posterior$mean^2
    }
    
    # 3. Periodic Orthogonalization: SVD rotation to maintain distinct, unique factors.
    if(iter %% 10 == 0) {
      svd_rot <- svd(EL, nu=K, nv=K)
      EL <- svd_rot$u %*% diag(svd_rot$d, K, K)
      EF <- EF %*% svd_rot$v
      EL2 <- EL^2; EF2 <- EF^2
    }
    
    # 4. Precision Update ***(Tau)***: Feature-specific measurement noise estimation.
    Var_Term <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
    Tau <- n / colSums((Y - EL %*% t(EF))^2 + Var_Term)
    
    if(iter %% 10 == 0) cat(sprintf("  Iteration %d | RMSE: %.4f\n", iter, history$rmse[iter]))
    if(mean(abs(EL - EL_old)) < tol && iter > 10) break
  }
  return(list(L=EL, F=EF, Beta=EBeta, Tau=Tau, history=history))
}

# ==============================================================================
# PART 3: VISUALIZATION & ANALYTICS
# ==============================================================================

#' High-Resolution GEP Heatmap with Legend
#' 
#' Diverging Red-Blue scale identifying influential features for each pathway.
#' This allows the user to see which specific genes/CpGs drive a prognostic program.
plot_gep_heatmap <- function(res, n_features = 50) {
  # Select top 50 features by absolute contribution across all factors.
  top_var_genes <- order(rowSums(abs(res$F)), decreasing = TRUE)[1:n_features]
  F_sub <- res$F[top_var_genes, ]
  
  # Blue = Protective/Negative; Red = Risk/Positive influence.
  palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_val <- max(abs(F_sub))
  
  # layout(matrix(1:2, ncol=2), widths=c(5,1)) splits the pane for the legend.
  layout(matrix(1:2, ncol=2), widths=c(5,1))
  
  # Panel 1: Main Heatmap
  par(mar=c(6, 4, 4, 1))
  image(1:nrow(F_sub), 1:ncol(F_sub), F_sub, 
        main=paste("GEP Feature Weights (Top", n_features, "IDs)"),
        xlab="Feature ID (Row Index)", ylab="Latent Factors",
        col=palette, axes=FALSE, zlim=c(-max_val, max_val))
  
  # Use the original feature row indices for the x-axis labels.
  axis(1, at=1:nrow(F_sub), labels=top_var_genes, las=2, cex.axis=0.6)
  axis(2, at=1:ncol(F_sub), labels=paste0("F", 1:ncol(F_sub)), las=1)
  box()
  
  # Panel 2: Vertical Color Legend
  par(mar=c(6, 1, 4, 3))
  legend_image <- as.matrix(seq(-max_val, max_val, length.out=100))
  image(1, seq(-max_val, max_val, length.out=100), t(legend_image), 
        col=palette, axes=FALSE, xlab="", ylab="")
  axis(4, las=1, cex.axis=0.8)
  mtext("Weight", side=4, line=2, cex=0.8)
  
  layout(1) # Reset layout to default.
}

#' Diagnostic Dashboard
#' 
#' Renders independent panes into RStudio Plots history.
#' Use the back/forward arrows in the RStudio Plots tab to navigate results.
visualize_dashboard <- function(res, data) {
  std_mar <- c(5, 5, 4, 2)
  
  # FIGURE 1: Optimization Trace. RMSE plateauing at 1.0 indicates perfect signal recovery.
  par(mfrow=c(1,1), mar=std_mar)
  plot(res$history$rmse[res$history$rmse > 0], type="l", lwd=2, col="#1f77b4", 
       main="Optimization Trace (Reconstruction Error)", xlab="Iteration", ylab="RMSE", bty="n")
  abline(h = 1.0, col = "red", lty = 2); grid(col="lightgray", lty="dotted")
  
  # FIGURE 2: Diverging Heatmap of Program Signatures.
  plot_gep_heatmap(res)
  
  # FIGURE 3: Kaplan-Meier Clinical Stratification Grid.
  K <- ncol(res$L); summary_tab <- get_factor_summary_table(res, data)
  par(mfrow=c(2, ceiling(K/2)), mar=c(4, 4, 3, 1))
  for(k in 1:K) {
    groups <- ifelse(res$L[,k] > median(res$L[,k]), "High Score", "Low Score")
    km <- survfit(Surv(data$time, data$status) ~ groups)
    p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
    plot(km, col=c("#d62728", "#1f77b4"), lwd=2, bty="n", 
         main=paste("Factor", k, "\n", p_label), xlab="Time", ylab="Prob")
  }
  
  # FIGURE 4: Simulated Signal Recovery (Truth vs. Estimated).
  if(!is.null(data$L_true)) {
    par(mfrow=c(1,1), mar=std_mar)
    # Automatically match the factor permutation and sign-flip.
    cors <- cor(data$L_true, res$L); best_match <- apply(abs(cors), 2, which.max)
    target_est <- which.max(abs(res$Beta)); target_true <- best_match[target_est]
    sign_correction <- sign(cors[target_true, target_est])
    
    plot(data$L_true[,target_true], res$L[,target_est] * sign_correction, 
         main=paste("Signal Recovery (Matched: Est", target_est, "vs True", target_true, ")"), 
         xlab="Ground Truth Signal", ylab="Model Estimated (Sign Corrected)", 
         col=rgb(0,0,0,0.4), pch=16, bty="n")
    abline(0,1, col="#d62728", lwd=2, lty=2)
  }
}

# ==============================================================================
# PART 4: EXECUTION & JUSTIFICATION
# ==============================================================================
if(DATA_MODE == "simulated") {
  # Synthetic Data Justification:
  # - N=250: Simulates a standard clinical cancer cohort.
  # - P=1000: Represents curated high-variance genomic features.
  # - K=5: Complex enough for overlapping pathways, simple enough for interpretation.
  sim_data_fn <- function(n=250, p=1000, k=5) {
    set.seed(123) # Reproducibility.
    L <- matrix(rnorm(n*k), n, k)
    F_mat <- matrix(0, p, k)
    # Sparse biological programs (5% density).
    for(i in 1:k) F_mat[sample(1:p, p*0.05), i] <- rnorm(p*0.05, 0, 5) 
    # Y = LF' + E.
    Y <- L %*% t(F_mat) + matrix(rnorm(n*p), n, p)
    # Factors 1-4 prognostic; Factor 5 structural.
    Beta <- c(1.5, -1.2, 0.8, -0.5, 0) 
    eta <- L %*% Beta
    real_times <- (-log(runif(n)) / (0.01 * exp(eta)))^(1/1.5)
    list(Y=Y, time=pmin(real_times, rexp(n, 1/50)), status=as.numeric(real_times<=rexp(n, 1/50)), 
         L_true=L, F_true=F_mat, B_true=Beta)
  }
  
  data <- sim_data_fn() # data
  res <- fit_supervised_mf(data$Y, data$time, data$status, K=5) # residuals
  
  cat("\n=== FACTOR SUMMARY TABLE (Log-rank & PVE) ===\n")
  print(get_factor_summary_table(res, data))
  
  cat("\n=== PROPORTIONAL HAZARDS (PH) TEST ===\n")
  # Validates the Cox assumption for our discovered factors.
  print(cox.zph(coxph(Surv(data$time, data$status) ~ res$L)))
  
  cat("\n=== MODEL PERFORMANCE (C-INDEX) ===\n")
  perf <- get_cindex_comparison(res$L, data)
  cat(sprintf("Original Top PCs C-index: %.3f\nSupervised Latent L C-index: %.3f\n", 
              perf$c_original, perf$c_latent))
  
  cat("\n=== TOP 5 FEATURES PER GEP ===\n")
  top_feats <- get_top_features(res$F, 5)
  for(k in 1:length(top_feats)) { 
    cat(sprintf("GEP %d:\n", k)); print(top_feats[[k]]) 
  }
  
  # Launch Visualization Dashboard Suite.
  visualize_dashboard(res, data)
}