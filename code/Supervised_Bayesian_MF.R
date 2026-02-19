# ------------------------------------------------------------------------------
# TITLE:       Supervised Bayesian Matrix Factorization
# AUTHOR:      Andrew Walther
# DATE:        February 18, 2026
#
# DESCRIPTION:
#   This script implements a probabilistic matrix factorization model that
#   simultaneously reduces dimensionality of high-throughput genomics data
#   (Y) and models time-to-event survival outcomes (t, delta).
#
# METHODOLOGY:
#   1. Genomics:  Modeled via Gaussian Matrix Factorization (Y ~ L * F').
#   2. Survival:  Modeled via Cox Proportional Hazards, coupled to L via Beta.
#   3. Approx:    2nd-order Taylor expansion (Newton-Raphson) for survival likelihood.
#   4. Priors:    Sparsity induced via Empirical Bayes Normal Means (EBNM).
#
# DEPENDENCIES:
#   - survival (CRAN)
#   - ebnm     (CRAN)
#
# INPUTS:
#   - Genomics Matrix (n x p)
#   - Survival Time vector (n)
#   - Status indicator vector (n)
# ------------------------------------------------------------------------------

# --- [USER TOGGLE] ---
DATA_MODE <- "simulated"  # Options: "simulated" or "real"

# --- [REAL DATA INPUTS] ---
# Only required if DATA_MODE is "real"
real_genomics_mat <- NULL   # Matrix (n x p)
real_clinical_df  <- NULL   # Dataframe (n rows)
real_time_col     <- "time" # Name of time column
real_status_col   <- "status" # Name of event column

# ==============================================================================
# PART 1: LIBRARIES & HELPER FUNCTIONS
# ==============================================================================
library(survival)
library(ebnm) 

#' Calculate Cox Gradients and Hessian (Diagonal Approximation)
#' 
#' WHY THIS IS NEEDED:
#' Approximates the complex Cox partial likelihood as a Weighted Least Squares
#' problem (Pseudo-response 'z' and weights 'w') to enable closed-form updates.
calc_cox_taylor <- function(eta, time, status) {
  n <- length(time)
  ord <- order(time)
  time_sorted <- time[ord]; status_sorted <- status[ord]; eta_sorted <- eta[ord]
  
  theta <- exp(eta_sorted)
  risk_sum <- rev(cumsum(rev(theta)))
  
  # Gradient (u) and Diagonal Hessian (w)
  h <- status_sorted / risk_sum
  H <- cumsum(h)
  u_sorted <- status_sorted - theta * H
  w_sorted <- theta * H 
  w_sorted[w_sorted < 1e-6] <- 1e-6 
  
  u <- numeric(n); w <- numeric(n)
  u[ord] <- u_sorted; w[ord] <- w_sorted
  return(list(u = u, w = w))
}

#' Inspect Data Structure
#' 
#' PURPOSE:
#' diagnostics tool to view dimensions and content summaries of input data.
#' Essential for ensuring REAL data matches the format of SIMULATED data.
inspect_data_structure <- function(data) {
  cat("\n=========================================\n")
  cat("       DATA INSPECTION REPORT            \n")
  cat("=========================================\n")
  
  # 1. Genomics Check
  cat("[1] GENOMICS MATRIX (Y):\n")
  cat(sprintf("    Dimensions: %d Patients x %d Features\n", nrow(data$Y), ncol(data$Y)))
  cat("    Value Preview (Top-Left 5x5):\n")
  print(data$Y[1:5, 1:5])
  
  # 2. Survival Check
  cat("\n[2] SURVIVAL OUTCOMES:\n")
  cat(sprintf("    Time Vector Length: %d\n", length(data$time)))
  cat(sprintf("    Status Vector Length: %d\n", length(data$status)))
  
  # Censoring Statistics
  events <- sum(data$status == 1)
  censored <- sum(data$status == 0)
  cat(sprintf("    Events: %d (%.1f%%) | Censored: %d (%.1f%%)\n", 
              events, 100*events/length(data$status), 
              censored, 100*censored/length(data$status)))
  
  cat("    Time Summary:\n")
  print(summary(data$time))
  
  # 3. Integrity Check
  if(nrow(data$Y) != length(data$time)) {
    warning(">> CRITICAL WARNING: Row mismatch between Genomics and Survival data!")
  } else {
    cat("\n[3] INTEGRITY CHECK: PASS (Dimensions match)\n")
  }
  cat("=========================================\n\n")
}

# ==============================================================================
# PART 2: THE FITTING ALGORITHM (The Engine)
# ==============================================================================
fit_supervised_mf <- function(Y, time, status, K = 5, max_iter = 100, tol = 1e-4) {
  
  n <- nrow(Y); p <- ncol(Y)
  
  # --- INITIALIZATION ---
  # Use SVD for "Warm Start" of Genomics components
  svd_init <- svd(Y, nu = K, nv = K)
  EL <- svd_init$u %*% diag(sqrt(svd_init$d[1:K]), K, K) # Initialize 1st moment for L
  EF <- svd_init$v %*% diag(sqrt(svd_init$d[1:K]), K, K)
  EL2 <- EL^2; EF2 <- EF^2
  
  # Init Beta using standard Cox regression on SVD factors
  df_surv <- data.frame(time = time, status = status, EL)
  colnames(df_surv)[3:(2+K)] <- paste0("L", 1:K)
  cox_fit <- coxph(as.formula(paste("Surv(time, status) ~ .")), data = df_surv)
  EBeta <- coef(cox_fit); EBeta[is.na(EBeta)] <- 0
  EBeta2 <- EBeta^2
  
  Tau <- rep(1, p) # Noise precision
  
  cat("Starting Variational Inference Fit...\n")
  
  for(iter in 1:max_iter) {
    EL_old <- EL
    
    # 1. SURVIVAL TAYLOR APPROXIMATION
    eta <- EL %*% EBeta
    taylor <- calc_cox_taylor(as.vector(eta), time, status)
    u <- taylor$u; w <- taylor$w
    z <- eta + u/w # Pseudo-response
    
    # 2. COORDINATE ASCENT
    for(k in 1:K) {
      
      # === UPDATE L (Loadings) - FUSION STEP ===
      # Confidence from Genomics + Confidence from Survival
      A_L_gen <- sum(Tau * EF2[,k])
      A_L_surv <- w * EBeta2[k]
      A_L <- A_L_gen + A_L_surv
      
      # Pull from Genomics + Pull from Survival
      Y_hat <- EL %*% t(EF)
      R_k <- Y - Y_hat + outer(EL[,k], EF[,k]) 
      R_k_weighted <- sweep(R_k, 2, Tau, `*`)
      B_L_gen <- as.vector(R_k_weighted %*% EF[,k])
      
      eta_no_k <- eta - (EL[,k] * EBeta[k])
      z_no_k <- z - eta_no_k
      B_L_surv <- w * z_no_k * EBeta[k]
      
      res_L <- ebnm(x = (B_L_gen + B_L_surv)/A_L, s = 1/sqrt(A_L), prior_family = "point_normal")
      EL[,k] <- res_L$posterior$mean
      EL2[,k] <- res_L$posterior$sd^2 + res_L$posterior$mean^2
      
      # === UPDATE F (Factors) - GENOMICS ONLY ===
      A_F <- Tau * sum(EL2[,k])
      B_F <- Tau * as.vector(t(R_k) %*% EL[,k])
      res_F <- ebnm(x = B_F/A_F, s = 1/sqrt(A_F), prior_family = "point_normal")
      EF[,k] <- res_F$posterior$mean
      EF2[,k] <- res_F$posterior$sd^2 + res_F$posterior$mean^2
      
      # === UPDATE BETA (Coefficients) - SURVIVAL ONLY ===
      # Accounts for uncertainty in L (E[L^2])
      A_Beta <- sum(w * EL2[,k]) 
      eta_current <- EL %*% EBeta
      z_no_k_beta <- z - (eta_current - (EL[,k] * EBeta[k]))
      B_Beta <- sum(w * z_no_k_beta * EL[,k])
      
      res_Beta <- ebnm(x = B_Beta/A_Beta, s = 1/sqrt(A_Beta), prior_family = "point_normal")
      EBeta[k] <- res_Beta$posterior$mean
      EBeta2[k] <- res_Beta$posterior$sd^2 + res_Beta$posterior$mean^2
    }
    
    # 3. UPDATE TAU & CHECK CONVERGENCE
    Y_hat <- EL %*% t(EF)
    Resid_Sq <- (Y - Y_hat)^2
    Var_Term <- (EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
    Tau <- n / colSums(Resid_Sq + Var_Term)
    
    diff <- mean(abs(EL - EL_old))
    if(diff < tol) { cat("Converged at iter", iter, "\n"); break }
  }
  return(list(EL=EL, EF=EF, EBeta=EBeta, Tau=Tau))
}

# ==============================================================================
# PART 3: SIMULATION FUNCTION
# ==============================================================================
sim_data_fn <- function(n=200, p=1000, k=5) {
  L <- matrix(rnorm(n*k), n, k)
  F_mat <- matrix(0, p, k) # Sparse Factors
  for(i in 1:k) F_mat[sample(1:p, p*0.1), i] <- rnorm(p*0.1) 
  Y <- L %*% t(F_mat) + matrix(rnorm(n*p), n, p)
  
  Beta <- c(0.8, -0.6, rep(0, k-2)) # Only 2 prognostic factors
  eta <- L %*% Beta
  real_times <- (-log(runif(n)) / (0.01 * exp(eta)))^(1/1.5)
  censor_times <- rexp(n, rate = 1/mean(real_times)*0.4)
  list(Y=Y, time=pmin(real_times, censor_times), status=as.numeric(real_times<=censor_times))
}

# ==============================================================================
# PART 4: EXECUTION LOGIC
# ==============================================================================

if(DATA_MODE == "simulated") {
  message(">>> Generating Simulated Data...")
  data <- sim_data_fn(n=200, p=500, k=5)
  
  # --- NEW: INSPECT DATA ---
  inspect_data_structure(data)
  
  res <- fit_supervised_mf(data$Y, data$time, data$status, K=5)
  print("Estimated Betas (First 2 should be non-zero):")
  print(round(res$EBeta, 3))
  
} else if(DATA_MODE == "real") {
  message(">>> Loading Real Data...")
  if(is.null(real_genomics_mat) || is.null(real_clinical_df)) stop("Missing real data objects.")
  
  data <- list(Y = real_genomics_mat, 
               time = real_clinical_df[[real_time_col]], 
               status = real_clinical_df[[real_status_col]])
  
  # --- NEW: INSPECT DATA ---
  inspect_data_structure(data)
  
  res <- fit_supervised_mf(data$Y, data$time, data$status, K=10)
  print("Top Survival Associated Factors:")
  print(round(res$EBeta, 3))
}