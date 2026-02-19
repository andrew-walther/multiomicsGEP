# ------------------------------------------------------------------------------
# TITLE:       Supervised Bayesian Matrix Factorization (Enhanced Analytics)
# AUTHOR:      Andrew Walther
# DATE:        February 19, 2026
#
# DESCRIPTION:
#   Simultaneously decomposes genomics data (Y) and survival outcomes (t, delta).
#   Includes convergence tracking, reconstruction error, and per-factor 
#   survival hypothesis testing (Log-rank tests).
# ------------------------------------------------------------------------------

# --- [USER TOGGLE] ---
DATA_MODE <- "simulated"  

# --- [REAL DATA INPUTS] ---
real_genomics_mat <- NULL   
real_clinical_df  <- NULL   
real_time_col     <- "time" 
real_status_col   <- "status" 

# ==============================================================================
# PART 1: LIBRARIES & HELPER FUNCTIONS
# ==============================================================================
library(survival)
library(ebnm) 

#' Calculate Cox Gradients and Hessian (Diagonal Approximation)
calc_cox_taylor <- function(eta, time, status) {
  n <- length(time)
  ord <- order(time) 
  time_sorted <- time[ord]; status_sorted <- status[ord]; eta_sorted <- eta[ord]
  theta <- exp(eta_sorted)
  risk_sum <- rev(cumsum(rev(theta)))
  h <- status_sorted / risk_sum
  H <- cumsum(h)
  u_sorted <- status_sorted - theta * H
  w_sorted <- theta * H 
  w_sorted[w_sorted < 1e-6] <- 1e-6 
  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_sorted; w[ord] <- w_sorted
  return(list(u = u, w = w))
}

#' Get Individual Factor Survival Statistics
#' 
#' WHAT: Performs a Log-rank test for every factor in the model.
#' WHY: To statistically determine which biological programs are prognostic.
get_factor_survival_stats <- function(L, data) {
  K <- ncol(L)
  stats <- lapply(1:K, function(k) {
    # Stratify by median loading score
    group <- ifelse(L[,k] > median(L[,k]), "High", "Low")
    sd_test <- survdiff(Surv(data$time, data$status) ~ group)
    # Calculate p-value from chi-square
    p_val <- 1 - pchisq(sd_test$chisq, length(sd_test$n) - 1)
    return(p_val)
  })
  return(unlist(stats))
}

#' Get C-Index Performance
get_cindex_comparison <- function(EL, data) {
  pca_y <- prcomp(data$Y, rank. = min(5, ncol(EL)))
  fit_pc <- coxph(Surv(data$time, data$status) ~ pca_y$x)
  fit_l  <- coxph(Surv(data$time, data$status) ~ EL)
  return(list(
    c_original = summary(fit_pc)$concordance[1],
    c_latent   = summary(fit_l)$concordance[1]
  ))
}

# ==============================================================================
# PART 2: THE FITTING ALGORITHM
# ==============================================================================
fit_supervised_mf <- function(Y, time, status, K = 5, max_iter = 100, tol = 1e-5) {
  n <- nrow(Y); p <- ncol(Y)
  svd_init <- svd(Y, nu = K, nv = K)
  EL <- svd_init$u %*% diag(sqrt(svd_init$d[1:K]), K, K) 
  EF <- svd_init$v %*% diag(sqrt(svd_init$d[1:K]), K, K)
  EL2 <- EL^2; EF2 <- EF^2
  
  df_surv <- data.frame(time = time, status = status, EL)
  colnames(df_surv)[3:(2+K)] <- paste0("L", 1:K)
  cox_fit <- coxph(as.formula(paste("Surv(time, status) ~ .")), data = df_surv)
  EBeta <- coef(cox_fit); EBeta[is.na(EBeta)] <- 0
  EBeta2 <- EBeta^2
  
  Tau <- rep(1, p) 
  history <- list(rmse = numeric(max_iter))
  
  cat("Starting Supervised Inference Engine...\n")
  
  for(iter in 1:max_iter) {
    EL_old <- EL
    eta <- EL %*% EBeta
    taylor <- calc_cox_taylor(as.vector(eta), time, status)
    z <- eta + taylor$u/taylor$w 
    w <- taylor$w
    history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))
    
    for(k in 1:K) {
      A_L <- sum(Tau * EF2[,k]) + w * EBeta2[k]
      Y_hat <- EL %*% t(EF)
      R_k <- Y - Y_hat + outer(EL[,k], EF[,k]) 
      B_L_gen <- as.vector(sweep(R_k, 2, Tau, `*`) %*% EF[,k])
      z_no_k <- z - (eta - (EL[,k] * EBeta[k]))
      B_L_surv <- w * z_no_k * EBeta[k]
      res_L <- ebnm(x = (B_L_gen + B_L_surv)/A_L, s = 1/sqrt(A_L), prior_family = "point_normal")
      EL[,k] <- res_L$posterior$mean
      EL2[,k] <- res_L$posterior$sd^2 + res_L$posterior$mean^2
      
      A_F <- Tau * sum(EL2[,k])
      B_F <- Tau * as.vector(t(R_k) %*% EL[,k])
      res_F <- ebnm(x = B_F/A_F, s = 1/sqrt(A_F), prior_family = "point_normal")
      EF[,k] <- res_F$posterior$mean
      EF2[,k] <- res_F$posterior$sd^2 + res_F$posterior$mean^2
      
      A_Beta <- sum(w * EL2[,k]) 
      z_no_k_beta <- z - (EL %*% EBeta - (EL[,k] * EBeta[k]))
      B_Beta <- sum(w * z_no_k_beta * EL[,k])
      res_Beta <- ebnm(x = B_Beta/A_Beta, s = 1/sqrt(A_Beta), prior_family = "point_normal")
      EBeta[k] <- res_Beta$posterior$mean
      EBeta2[k] <- res_Beta$posterior$sd^2 + res_Beta$posterior$mean^2
    }
    
    if(iter %% 10 == 0) {
      svd_rot <- svd(EL, nu=K, nv=K)
      EL <- svd_rot$u %*% diag(svd_rot$d, K, K)
      EF <- EF %*% svd_rot$v
      EL2 <- EL^2; EF2 <- EF^2
    }
    
    Var_Term <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
    Tau <- n / colSums((Y - EL %*% t(EF))^2 + Var_Term)
    
    if(iter %% 10 == 0) cat(sprintf("  Iteration %d | RMSE: %.4f\n", iter, history$rmse[iter]))
    if(mean(abs(EL - EL_old)) < tol && iter > 10) break
  }
  return(list(L=EL, F=EF, Beta=EBeta, Tau=Tau, history=history))
}

# ==============================================================================
# PART 3: SIMULATION & ANALYTICS
# ==============================================================================
sim_data_fn <- function(n=200, p=1000, k=5) {
  set.seed(123) 
  L <- matrix(rnorm(n*k), n, k)
  F_mat <- matrix(0, p, k)
  for(i in 1:k) F_mat[sample(1:p, p*0.05), i] <- rnorm(p*0.05, 0, 5) 
  Y <- L %*% t(F_mat) + matrix(rnorm(n*p), n, p)
  Beta <- c(1.5, -1.2, 0.8, -0.5, 0) 
  eta <- L %*% Beta
  real_times <- (-log(runif(n)) / (0.01 * exp(eta)))^(1/1.5)
  list(Y=Y, time=pmin(real_times, rexp(n, 1/50)), status=as.numeric(real_times<=rexp(n, 1/50)), 
       L_true=L, F_true=F_mat, B_true=Beta)
}

visualize_analytics <- function(res, data) {
  # 1. Main Diagnostics (RMSE and Sparsity)
  par(mfrow=c(2,2), mar=c(4.5, 4.5, 3, 2))
  plot(res$history$rmse[res$history$rmse > 0], type="l", lwd=2, col="#1f77b4",
       main="Optimization Trace", xlab="Iteration", ylab="RMSE", bty="n")
  abline(h = 1.0, col = "red", lty = 2) 
  
  threshold <- quantile(abs(res$F), 0.95)
  image(t(abs(res$F) > threshold), main="GEP Top Signatures", 
        xlab="Features", ylab="GEPs", col=c("#fcfcfc", "#d62728"), axes=FALSE)
  axis(2, at=seq(0, 1, length.out=ncol(res$F)), labels=1:ncol(res$F), las=1)
  
  # 2. MATCHED SIGNAL RECOVERY
  if(!is.null(data$L_true)) {
    cors <- cor(data$L_true, res$L)
    best_match <- apply(abs(cors), 2, which.max)
    target_est <- which.max(abs(res$Beta))
    target_true <- best_match[target_est]
    sign_correction <- sign(cors[target_true, target_est])
    
    plot(data$L_true[,target_true], res$L[,target_est] * sign_correction, 
         main="Latent Signal Recovery", xlab="Truth", ylab="Estimated", 
         col=rgb(0,0,0,0.3), pch=16, bty="n")
    abline(0,1, col="#d62728", lwd=2, lty=2)
  }
  
  # 3. PER-FACTOR SURVIVAL GRID
  # Open a new device or just clear the current one to show KM plots
  K <- ncol(res$L)
  p_vals <- get_factor_survival_stats(res$L, data)
  
  dev.new(width=10, height=7)
  par(mfrow=c(2, ceiling(K/2)), mar=c(4, 4, 3, 1))
  for(k in 1:K) {
    groups <- ifelse(res$L[,k] > median(res$L[,k]), "High", "Low")
    km <- survfit(Surv(data$time, data$status) ~ groups)
    p_label <- sprintf("p = %.4f", p_vals[k])
    plot(km, col=c("#d62728", "#1f77b4"), lwd=2, bty="n",
         main=paste("Factor", k, "\n", p_label), xlab="Time", ylab="Prob")
    if(k == 1) legend("bottomleft", legend=c("High", "Low"), col=c("#d62728", "#1f77b4"), lty=1, bty="n", cex=0.7)
  }
}

# ==============================================================================
# PART 4: EXECUTION
# ==============================================================================
if(DATA_MODE == "simulated") {
  data <- sim_data_fn(n=250, p=1000, k=5)
  res <- fit_supervised_mf(data$Y, data$time, data$status, K=5)
  
  cat("\n=== SURVIVAL HYPOTHESIS TESTS (Log-rank) ===\n")
  p_vals <- get_factor_survival_stats(res$L, data)
  for(k in 1:length(p_vals)) {
    cat(sprintf("Factor %d: p-value = %.4f %s\n", 
                k, p_vals[k], ifelse(p_vals[k] < 0.05, "*SIGNIFICANT*", "")))
  }
  
  cat("\n=== MODEL PERFORMANCE ===\n")
  c_test <- get_cindex_comparison(res$L, data)
  cat(sprintf("Original Top PCs C-index: %.3f\nSupervised Latent L C-index: %.3f\n", 
              c_test$c_original, c_test$c_latent))
  
  visualize_analytics(res, data)
}