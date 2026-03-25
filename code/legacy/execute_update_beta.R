# 1. Load the function
source("code/update_beta.R")   # loads ebnm too

# 2. Make up plausible inputs (no simulation required)
set.seed(42)
n <- 100

# Cox weights: positive, O(1) — represents neg-diagonal Hessian
w <- abs(rnorm(n, mean = 2, sd = 0.5))

# Loadings for factor k: random normal
EL_k  <- rnorm(n)

# Second moment = mean^2 + posterior variance.
# Setting variance > 0 is more realistic than EL_k^2 alone.
EL2_k <- EL_k^2 + 0.1        # 0.1 = posterior variance per patient

# Partial working response: simulate as if true beta_k = 1.5
# z_no_k ~= beta_true * l_ik + noise
z_no_k <- 1.5 * EL_k + rnorm(n, sd = 0.5)

# 3. Call the function
res <- update_beta_k(w, z_no_k, EL_k, EL2_k)

# 4. Inspect
cat("A_k (precision):  ", round(res$A, 3), "\n")
cat("B_k (signal):     ", round(res$B, 3), "\n")
cat("x_k = B/A:        ", round(res$x, 3), "\n")   # ~1.5 if signal is clean
cat("s_k = 1/sqrt(A):  ", round(res$s, 5), "\n")
cat("Posterior mean:   ", round(res$mean, 4), "\n")   # should be ~1.5
cat("Posterior SD:     ", round(res$sd, 4), "\n")
cat("Second moment:    ", round(res$second, 4), "\n")  # sd^2 + mean^2

# --- K=5 full-loop ---
K <- 5
EL    <- matrix(rnorm(n * K), n, K)
EL2   <- EL^2 + 0.05
EBeta <- rnorm(K) * 0.1          # warm start near zero
z     <- EL %*% c(1.5, -1.2, 0.8, -0.5, 0) + rnorm(n, sd = 0.5)
w     <- abs(rnorm(n, mean = 2, sd = 0.5))

res_all <- update_beta_all(w, z, EL, EL2, EBeta)

cat("Estimated beta: ", round(res_all$EBeta, 3), "\n")
# Should roughly recover signs: (+, -, +, -, ~0)
