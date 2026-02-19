# ------------------------------------------------------------------------------
# TITLE:       Supervised Bayesian Matrix Factorization (Enhanced Analytics)
# AUTHOR:      Andrew Walther
# DATE:        February 19, 2026
#
# DESCRIPTION:
#   This script implements a unified probabilistic framework to decompose 
#   high-dimensional genomics data (Y) while simultaneously modeling survival 
#   outcomes (t, delta). By sharing a latent space (L), the model identifies 
#   "Gene Expression Programs" (GEPs) that are explicitly prognostic.
# ------------------------------------------------------------------------------

# --- [USER TOGGLE] ---
# Sets the data source. "simulated" generates a synthetic benchmark with 
# known ground truth; "real" is used for applying the model to actual datasets.
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

# 'survival' handles the Cox Proportional Hazards math and C-index calculations.
library(survival)
# 'ebnm' provides Empirical Bayes Normal Means solvers, which estimate the 
# sparse priors (Spike-and-Slab) during the variational update steps.
library(ebnm) 

#' Calculate Cox Gradients and Hessian (Diagonal Approximation)
#' 
#' This function transforms the complex, non-conjugate Cox survival likelihood 
#' into a locally linear "Weighted Least Squares" form using a 2nd-order Taylor 
#' expansion. This is the mathematical "bridge" that allows us to treat survival 
#' outcomes as pseudo-genomics observations.
calc_cox_taylor <- function(eta, time, status) {
  n <- length(time)
  # Sorting by time is required to calculate the risk sets used in Cox models.
  ord <- order(time) 
  time_sorted <- time[ord]; status_sorted <- status[ord]; eta_sorted <- eta[ord]
  
  # Exponentiated linear predictor represents relative risk.
  theta <- exp(eta_sorted)
  # Cumulative hazard denominator for the partial likelihood.
  risk_sum <- rev(cumsum(rev(theta)))
  
  # Hazard increments used for derivative calculations.
  h <- status_sorted / risk_sum
  H <- cumsum(h)
  
  # Gradient (u) and Diagonal Hessian (w) for the quadratic approximation.
  u_sorted <- status_sorted - theta * H
  w_sorted <- theta * H 
  
  # Numerical stability floor to prevent division-by-zero for censored patients.
  w_sorted[w_sorted < 1e-6] <- 1e-6 
  
  # Re-ordering back to original patient indices for data alignment.
  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_sorted; w[ord] <- w_sorted
  return(list(u = u, w = w))
}

#' Get C-Index Performance Comparison
#' 
#' Calculates Harrell's Concordance Index (C-index) to evaluate how well the 
#' model predicts the order of survival events. It compares our Supervised 
#' Latent space against standard Unsupervised Principal Components (PCA) 
#' to prove that the supervision is successfully denoising the clinical signal.
get_cindex_comparison <- function(EL, data) {
  # Baseline: Top 5 Principal Components of the raw genomics data.
  pca_y <- prcomp(data$Y, rank. = min(5, ncol(EL)))
  fit_pc <- coxph(Surv(data$time, data$status) ~ pca_y$x)
  
  # Comparison: The Latent Loadings (L) produced by our Supervised model.
  fit_l  <- coxph(Surv(data$time, data$status) ~ EL)
  
  return(list(
    c_original = round(summary(fit_pc)$concordance[1], 3),
    c_latent   = round(summary(fit_l)$concordance[1], 3)
  ))
}

#' Extract Top Influential Features per Factor
#' 
#' Ranks genomics features (genes/CpGs) by their absolute magnitude in the 
#' Factor matrix (F). This identifies the biological "identity" of each 
#' program and provides a list of regulatory sites associated with prognosis.
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
#' Aggregates clinical significance (Beta/P-value), biological sparsity, 
#' and genomics variance explained (PVE) into a single diagnostic table.
get_factor_summary_table <- function(res, data) {
  K <- ncol(res$L)
  
  # Calculate Log-Rank p-values for each factor via median stratification.
  p_vals <- lapply(1:K, function(k) {
    group <- ifelse(res$L[,k] > median(res$L[,k]), "High", "Low")
    sd_test <- survdiff(Surv(data$time, data$status) ~ group)
    1 - pchisq(sd_test$chisq, length(sd_test$n) - 1)
  })
  
  # Calculate sparsity (percentage of non-zero gene weights).
  sparsity <- colMeans(res$F != 0) * 100
  
  # Calculate Proportion of Variance Explained (PVE) for the genomics matrix.
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
#' The core engine that uses Coordinate Ascent Variational Inference (CAVI) 
#' to iteratively refine the model parameters. It balances the need to 
#' reconstruct the genomics matrix with the need to predict survival.
fit_supervised_mf <- function(Y, time, status, K = 5, max_iter = 100, tol = 1e-5) {
  n <- nrow(Y); p <- ncol(Y)
  
  # Initialization via SVD provides a deterministic, high-variance starting 
  # subspace, which helps the algorithm avoid poor local optima.
  svd_init <- svd(Y, nu = K, nv = K)
  EL <- svd_init$u %*% diag(sqrt(svd_init$d[1:K]), K, K) # Expected Loadings
  EF <- svd_init$v %*% diag(sqrt(svd_init$d[1:K]), K, K) # Expected Factors
  EL2 <- EL^2; EF2 <- EF^2 # Second moments assume initially zero variance.
  
  # Warm-start Beta using a standard Cox model on the initial SVD loadings 
  # to ensure the clinical signal is recognized early in the process.
  df_surv <- data.frame(time = time, status = status, EL)
  colnames(df_surv)[3:(2+K)] <- paste0("L", 1:K)
  cox_fit <- coxph(as.formula(paste("Surv(time, status) ~ .")), data = df_surv)
  EBeta <- coef(cox_fit); EBeta[is.na(EBeta)] <- 0
  EBeta2 <- EBeta^2
  
  # Tau represents inverse noise variance per feature (heteroskedasticity).
  Tau <- rep(1, p) 
  history <- list(rmse = numeric(max_iter))
  
  cat("Starting Supervised Inference Engine...\n")
  
  for(iter in 1:max_iter) {
    EL_old <- EL
    
    # 1. Survival Linearization: Bridging the Cox likelihood to Gaussian priors.
    eta <- EL %*% EBeta
    taylor <- calc_cox_taylor(as.vector(eta), time, status)
    z <- eta + taylor$u/taylor$w 
    w <- taylor$w
    
    # Monitoring RMSE. In simulations with noise SD=1.0, this should floor at 1.0.
    history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))
    
    # 2. Sequential Coordinate Ascent Updates (One Factor at a Time)
    for(k in 1:K) {
      
      # --- Update Loadings (L): Data Fusion Step ---
      # Precision A_L is the sum of genomics precision and survival precision.
      # This forces L to account for both modalities.
      A_L <- sum(Tau * EF2[,k]) + w * EBeta2[k]
      
      # Use partial residuals to isolate the contribution of the k-th factor.
      Y_hat <- EL %*% t(EF)
      R_k <- Y - Y_hat + outer(EL[,k], EF[,k]) 
      B_L_gen <- as.vector(sweep(R_k, 2, Tau, `*`) %*% EF[,k])
      
      z_no_k <- z - (eta - (EL[,k] * EBeta[k]))
      B_L_surv <- w * z_no_k * EBeta[k]
      
      # ebnm solves the sub-problem by finding the sparse posterior mean.
      res_L <- ebnm(x = (B_L_gen + B_L_surv)/A_L, s = 1/sqrt(A_L), prior_family = "point_normal")
      EL[,k] <- res_L$posterior$mean
      EL2[,k] <- res_L$posterior$sd^2 + res_L$posterior$mean^2 # Captures uncertainty.
      
      # --- Update Factors (F): GEP Extraction ---
      # Factors are updated primarily through the genomics matrix signal.
      A_F <- Tau * sum(EL2[,k])
      B_F <- Tau * as.vector(t(R_k) %*% EL[,k])
      res_F <- ebnm(x = B_F/A_F, s = 1/sqrt(A_F), prior_family = "point_normal")
      EF[,k] <- res_F$posterior$mean
      EF2[,k] <- res_F$posterior$sd^2 + res_F$posterior$mean^2
      
      # --- Update Beta: Clinical Supervision Step ---
      # Accounts for "Error-in-Variables" by incorporating the loading variance.
      # This penalizes Beta if the factor loading is uncertain or noisy.
      A_Beta <- sum(w * EL2[,k]) 
      z_no_k_beta <- z - (EL %*% EBeta - (EL[,k] * EBeta[k]))
      B_Beta <- sum(w * z_no_k_beta * EL[,k])
      res_Beta <- ebnm(x = B_Beta/A_Beta, s = 1/sqrt(A_Beta), prior_family = "point_normal")
      EBeta[k] <- res_Beta$posterior$mean
      EBeta2[k] <- res_Beta$posterior$sd^2 + res_Beta$posterior$mean^2
    }
    
    # 3. Periodic Orthogonalization (Rotation)
    # Rotating the basis ensures that the GEPs remain distinct and non-redundant.
    if(iter %% 10 == 0) {
      svd_rot <- svd(EL, nu=K, nv=K)
      EL <- svd_rot$u %*% diag(svd_rot$d, K, K)
      EF <- EF %*% svd_rot$v
      EL2 <- EL^2; EF2 <- EF^2
    }
    
    # 4. Precision Update (Tau): Estimating feature-specific noise floors.
    Var_Term <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
    Tau <- n / colSums((Y - EL %*% t(EF))^2 + Var_Term)
    
    if(iter %% 10 == 0) cat(sprintf("  Iteration %d | RMSE: %.4f\n", iter, history$rmse[iter]))
    
    # Check for parameter stability to determine convergence.
    if(mean(abs(EL - EL_old)) < tol && iter > 10) break
  }
  return(list(L=EL, F=EF, Beta=EBeta, Tau=Tau, history=history))
}

# ==============================================================================
# PART 3: VISUALIZATION & ANALYTICS
# ==============================================================================

#' Detailed GEP Heatmap with Legend
#' 
#' Visualizes the top 50 genomic features for each factor using a 
#' diverging red-blue scale to show regulatory influence gradients.
plot_gep_heatmap <- function(res, n_features = 50) {
  # Select the features with the highest absolute weight across the factors.
  top_var_genes <- order(rowSums(abs(res$F)), decreasing = TRUE)[1:n_features]
  F_sub <- res$F[top_var_genes, ]
  
  # Blue = Protective/Negative weights; Red = Risk/Positive weights.
  palette <- colorRampPalette(c("blue", "white", "red"))(100)
  max_val <- max(abs(F_sub))
  
  # Split the pane to include a color legend on the right.
  layout(matrix(1:2, ncol=2), widths=c(5,1))
  
  par(mar=c(6, 4, 4, 1))
  image(1:nrow(F_sub), 1:ncol(F_sub), F_sub, 
        main=paste("GEP Feature Weights (Top", n_features, "IDs)"),
        xlab="Feature ID (Index)", ylab="Latent Factors",
        col=palette, axes=FALSE, zlim=c(-max_val, max_val))
  
  # Add axis labels using the actual row indices of the features.
  axis(1, at=1:nrow(F_sub), labels=top_var_genes, las=2, cex.axis=0.6)
  axis(2, at=1:ncol(F_sub), labels=paste0("F", 1:ncol(F_sub)), las=1)
  box()
  
  # Draw the legend color bar.
  par(mar=c(6, 1, 4, 3))
  legend_image <- as.matrix(seq(-max_val, max_val, length.out=100))
  image(1, seq(-max_val, max_val, length.out=100), t(legend_image), 
        col=palette, axes=FALSE, xlab="", ylab="")
  axis(4, las=1, cex.axis=0.8)
  mtext("Weight", side=4, line=2, cex=0.8)
  
  layout(1) # Reset layout.
}

#' Generate Diagnostic Dashboard
#' 
#' Sequentially renders independent panes into the RStudio Plots tab history.
visualize_dashboard <- function(res, data) {
  std_mar <- c(5, 5, 4, 2)
  
  # PANE 1: Convergence diagnostics.
  par(mfrow=c(1,1), mar=std_mar)
  plot(res$history$rmse[res$history$rmse > 0], type="l", lwd=2, col="#1f77b4",
       main="Optimization Trace (Reconstruction Error)", 
       xlab="Iteration", ylab="RMSE", bty="n")
  abline(h = 1.0, col = "red", lty = 2) # Theoretical floor for SD=1 noise.
  grid(col="lightgray", lty="dotted")
  
  # PANE 2: Detailed biological pathway identities.
  plot_gep_heatmap(res)
  
  # PANE 3: Clinical stratification grid for every factor.
  K <- ncol(res$L)
  summary_tab <- get_factor_summary_table(res, data)
  par(mfrow=c(2, ceiling(K/2)), mar=c(4, 4, 3, 1))
  for(k in 1:K) {
    groups <- ifelse(res$L[,k] > median(res$L[,k]), "High", "Low")
    km <- survfit(Surv(data$time, data$status) ~ groups)
    p_label <- sprintf("p = %.4f", summary_tab$LogRank_P[k])
    plot(km, col=c("#d62728", "#1f77b4"), lwd=2, bty="n",
         main=paste("Factor", k, "\n", p_label), xlab="Time", ylab="Prob")
  }
  
  # PANE 4: Truth recovery check (Simulation Only).
  if(!is.null(data$L_true)) {
    par(mfrow=c(1,1), mar=std_mar)
    # Matching algorithm finds the best pair of estimated vs true factors.
    cors <- cor(data$L_true, res$L)
    best_match <- apply(abs(cors), 2, which.max)
    target_est <- which.max(abs(res$Beta))
    target_true <- best_match[target_est]
    sign_correction <- sign(cors[target_true, target_est])
    
    plot(data$L_true[,target_true], res$L[,target_est] * sign_correction, 
         main=paste("Matched Signal Recovery (Est", target_est, "vs True", target_true, ")"), 
         xlab="Ground Truth Signal", ylab="Model Estimated (Sign Corrected)", 
         col=rgb(0,0,0,0.4), pch=16, bty="n")
    abline(0,1, col="#d62728", lwd=2, lty=2)
  }
}

# ==============================================================================
# PART 4: EXECUTION
# ==============================================================================
if(DATA_MODE == "simulated") {
  # Simulation parameters: 250 patients, 1000 features, 5 latent biological programs.
  sim_data_fn <- function(n=250, p=1000, k=5) {
    set.seed(123) # Ensures reproducibility across sessions.
    L <- matrix(rnorm(n*k), n, k)
    F_mat <- matrix(0, p, k)
    # Generate sparse biological programs (5% activity).
    for(i in 1:k) F_mat[sample(1:p, p*0.05), i] <- rnorm(p*0.05, 0, 5) 
    # Add random unit noise.
    Y <- L %*% t(F_mat) + matrix(rnorm(n*p), n, p)
    # Define clinical outcomes driven by factors 1-4.
    Beta <- c(1.5, -1.2, 0.8, -0.5, 0) 
    eta <- L %*% Beta
    real_times <- (-log(runif(n)) / (0.01 * exp(eta)))^(1/1.5)
    list(Y=Y, time=pmin(real_times, rexp(n, 1/50)), status=as.numeric(real_times<=rexp(n, 1/50)), 
         L_true=L, F_true=F_mat, B_true=Beta)
  }
  
  data <- sim_data_fn()
  res <- fit_supervised_mf(data$Y, data$time, data$status, K=5)
  
  # --- CONSOLE DISCOVERY LOGS ---
  cat("\n=== FACTOR SUMMARY TABLE ===\n")
  print(get_factor_summary_table(res, data))
  
  cat("\n=== MODEL PERFORMANCE (C-INDEX) ===\n")
  perf <- get_cindex_comparison(res$L, data)
  cat(sprintf("Original Top PCs C-index: %.3f\nSupervised Latent L C-index: %.3f\n", 
              perf$c_original, perf$c_latent))
  
  cat("\n=== TOP 5 FEATURES PER GEP ===\n")
  top_feats <- get_top_features(res$F, 5)
  for(k in 1:length(top_feats)) {
    cat(sprintf("GEP %d:\n", k)); print(top_feats[[k]])
  }
  
  # Launch Visualization Suite - Navigate using Plot Pane back/forward arrows.
  visualize_dashboard(res, data)
}