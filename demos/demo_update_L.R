# =============================================================================
# demos/demo_update_L.R
#
# Interactive demonstration of code/update_L.R
#
# Run from project root:
#   Rscript demos/demo_update_L.R
#
# Five scenarios illustrating the q(L) CAVI update — a VECTOR EBNM
# problem driven by both genomics and survival likelihoods.
#
# Related files:
#   code/update_L.R                           -- the functions being demonstrated
#   derivations/qL/qL_update_derivation.pdf  -- the math behind them
#   tests/test_update_L.R                    -- correctness tests
# =============================================================================

source("code/update_L.R")    # loads compute_R_k, update_L_k, update_L_all
source("code/update_beta.R") # for compute_z_no_k

sep  <- function() cat(strrep("-", 60), "\n")
sec  <- function(title) { cat("\n"); sep(); cat(title, "\n"); sep() }
val  <- function(label, x, digits = 4) {
  cat(sprintf("  %-36s %s\n", paste0(label, ":"),
              paste(round(x, digits), collapse = "  ")))
}

cat("============================================================\n")
cat(" demo_update_L.R\n")
cat(" Demonstrations of the q(l_{ik}) variational update\n")
cat("============================================================\n")

set.seed(42)
n <- 150   # patients
p <- 100   # genomic features


# ============================================================
# DEMO 1: Anatomy of A_L and B_L — two sources, sample-varying
# ============================================================
sec("DEMO 1: Anatomy of A_L and B_L — two sources of information")

cat("
Unlike the beta update (scalar EBNM), the L update is a VECTOR
EBNM because A_L and B_L vary across patients i:

  A_L[i] = sum_j(tau_j * E[f^2_jk])     <- genomics (scalar, same for all i)
          + W_ii * E[beta_k^2]            <- survival (n-vector, patient-specific)

  B_L[i] = sum_j(tau_j * R_k[i,j] * f_bar_jk)   <- genomics
          + W_ii * z_no_k[i] * beta_bar_k          <- survival
\n")

Tau    <- rep(2.0, p)
EF_k   <- rnorm(p);  EF2_k  <- EF_k^2 + 0.1
w      <- abs(rnorm(n, mean = 1.5, sd = 0.3))   # Cox weights
EBeta_k  <- 0.8;  EBeta2_k <- EBeta_k^2 + 0.05
R_k      <- matrix(rnorm(n * p), n, p)
z_no_k   <- rnorm(n)

res <- update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k)

A_gen  <- sum(Tau * EF2_k)      # the scalar genomics contribution
val("A_gen = sum(Tau*EF2_k)", A_gen)
cat(sprintf("  %-36s min=%.3f  max=%.3f  mean=%.3f\n",
            "A_surv[i] = w[i]*EBeta2_k :",
            min(w * EBeta2_k), max(w * EBeta2_k), mean(w * EBeta2_k)))
cat(sprintf("  %-36s min=%.3f  max=%.3f  mean=%.3f\n",
            "A_L[i] = A_gen + A_surv[i] :",
            min(res$A), max(res$A), mean(res$A)))

cat("\n  A_L varies across patients because the Cox weights W_ii differ.\n")
cat("  The genomics part is constant (scalar); only survival part varies.\n")

cat(sprintf("\n  B_L range: min=%.3f  max=%.3f\n", min(res$B), max(res$B)))
cat(sprintf("  x_L range: min=%.3f  max=%.3f\n", min(res$x), max(res$x)))
cat(sprintf("  s_L range: min=%.4f  max=%.4f\n", min(res$s), max(res$s)))


# ============================================================
# DEMO 2: Signal recovery (K=1, known loadings)
# ============================================================
sec("DEMO 2: Signal recovery (K=1, true loadings known)")

cat("
Simulate Y = L_true * f' + noise, then call update_L_k.
The posterior means should correlate with the true loading column.
\n")

L_true <- rnorm(n)          # true loading column
F_true <- rnorm(p)          # true factor column
Y <- outer(L_true, F_true) + matrix(rnorm(n * p, sd = 0.5), n, p)

# Set up with correct quantities
Tau_rec <- rep(4.0, p)
EF_rec  <- F_true;  EF2_rec <- F_true^2 + 0.01   # nearly correct factor
EBeta_rec <- 0.0;   EBeta2_rec <- 0.0             # no survival for simplicity
w_rec   <- rep(0.0, n)                             # pure genomics
EL_init <- matrix(0, n, 1); EF_init <- matrix(F_true, p, 1)
R_k_rec <- Y   # with no other factors, R_k = Y
z_no_k_rec <- rep(0, n)

res_rec <- update_L_k(Tau_rec, EF_rec, EF2_rec,
                       w_rec, EBeta_rec, EBeta2_rec,
                       R_k_rec, z_no_k_rec)

rho <- cor(res_rec$mean, L_true)
cat(sprintf("  Correlation with true L: %.4f\n", rho))
cat(sprintf("  Mean posterior SD:       %.4f\n", mean(res_rec$sd)))
cat(sprintf("  Mean |posterior mean|:   %.4f\n", mean(abs(res_rec$mean))))

cat("\n  Top-5 estimated loadings vs true:\n")
top5 <- order(abs(res_rec$mean), decreasing = TRUE)[1:5]
cat("  Sample  True L    Est L     Posterior SD\n")
for (i in top5) {
  cat(sprintf("  %6d  %+7.4f   %+7.4f   %.4f\n",
              i, L_true[i], res_rec$mean[i], res_rec$sd[i]))
}


# ============================================================
# DEMO 3: Genomics-only vs. genomics + survival
# ============================================================
sec("DEMO 3: Genomics-only vs. genomics + survival")

cat("
When survival information is included (w > 0), the precision A_L[i]
increases and the posterior SDs tighten.  The survival term reinforces
the genomics signal when both sources agree.
\n")

# Setup
L_sim  <- rnorm(n)
F_sim  <- rnorm(p) * 0.5
beta_k <- 1.2
Y_sim  <- outer(L_sim, F_sim) + matrix(rnorm(n * p, sd = 0.7), n, p)
tau_sim <- rep(2.0, p)
EF_sim  <- F_sim;  EF2_sim <- F_sim^2 + 0.05
R_k_sim <- Y_sim

# Genomics-only (w=0)
w_zero  <- rep(0, n)
z_sim   <- rep(0, n)
res_gen <- update_L_k(tau_sim, EF_sim, EF2_sim,
                       w_zero, beta_k, beta_k^2 + 0.02,
                       R_k_sim, z_sim)

# With survival signal
w_surv  <- abs(rnorm(n, mean = 1.0, sd = 0.2))
z_surv  <- beta_k * L_sim + rnorm(n, sd = 0.3)
res_both <- update_L_k(tau_sim, EF_sim, EF2_sim,
                        w_surv, beta_k, beta_k^2 + 0.02,
                        R_k_sim, z_surv)

cat("  Comparison (genomics only  vs  genomics + survival):\n\n")
cat(sprintf("  Mean posterior SD:  genomics=%.4f   combined=%.4f\n",
            mean(res_gen$sd), mean(res_both$sd)))
cat(sprintf("  Cor w/ true L:      genomics=%.4f   combined=%.4f\n",
            cor(res_gen$mean, L_sim), cor(res_both$mean, L_sim)))
cat(sprintf("  Mean |A_L|:         genomics=%.2f    combined=%.2f\n",
            mean(res_gen$A), mean(res_both$A)))
cat("\n  Survival tightens posterior (smaller SD) when it aligns with genomics.\n")


# ============================================================
# DEMO 4: Error-in-variables — factor uncertainty inflates A_L
# ============================================================
sec("DEMO 4: Error-in-variables (factor posterior uncertainty)")

cat("
The precision A_L has a genomics part: sum_j(tau_j * E[f^2_jk]).
This uses E[f^2] (the FULL second moment), not f_bar^2.

When we are uncertain about the factor column EF_k, the second moment
E[f^2] = Var_q(f) + f_bar^2 is larger than f_bar^2.
This inflates A_L -> shrinks loading estimates toward zero.
\n")

EF_base <- rnorm(p, sd = 1.0)
R_demo  <- matrix(rnorm(n * p), n, p)
w_demo  <- rep(0, n)   # pure genomics for clarity
tau_demo <- rep(1.0, p)

EF2_low  <- EF_base^2 + 0.01   # nearly certain about factor
EF2_high <- EF_base^2 + 1.00   # highly uncertain about factor

res_low  <- update_L_k(tau_demo, EF_base, EF2_low,
                        w_demo, 0, 0, R_demo, rep(0, n))
res_high <- update_L_k(tau_demo, EF_base, EF2_high,
                        w_demo, 0, 0, R_demo, rep(0, n))

cat(sprintf("  Posterior variance in F:  low=%.2f  high=%.2f\n", 0.01, 1.00))
cat(sprintf("  A_gen (genomics precis):  low=%.1f  high=%.1f\n",
            sum(tau_demo * EF2_low), sum(tau_demo * EF2_high)))
cat(sprintf("  Mean |posterior mean|:    low=%.4f  high=%.4f\n",
            mean(abs(res_low$mean)), mean(abs(res_high$mean))))
cat(sprintf("  Mean posterior SD:        low=%.4f  high=%.4f\n",
            mean(res_low$sd), mean(res_high$sd)))
cat("\n  More uncertain factor -> larger A_L -> more shrinkage in loadings.\n")
cat("  This prevents loadings from chasing poorly estimated factors.\n")


# ============================================================
# DEMO 5: Multi-factor recovery (K=5) with update_L_all
# ============================================================
sec("DEMO 5: Multi-factor recovery (K=5) with update_L_all()")

cat("
update_L_all() loops k=1..K with Gauss-Seidel ordering:
once EL[,k] is updated, the new values are used in compute_R_k
for subsequent factors.

True L matrix: K=5 columns, each of length n=150
\n")

K         <- 5
n5        <- 150
p5        <- 120
L5_true   <- matrix(rnorm(n5 * K), n5, K)
F5_true   <- matrix(rnorm(p5 * K, sd = 0.5), p5, K)
beta5     <- c(1.2, -0.8, 0.5, 0.0, -1.0)
Y5        <- L5_true %*% t(F5_true) + matrix(rnorm(n5 * p5, sd = 0.8), n5, p5)

# Initialize with SVD
sv        <- svd(Y5, nu = K, nv = K)
dk        <- sqrt(pmax(sv$d[1:K], 0))
EL5_init  <- sv$u %*% diag(dk, K, K)
EF5_init  <- sv$v %*% diag(dk, K, K)
EL25_init <- EL5_init^2
EF25_init <- EF5_init^2 + 0.1
Tau5      <- rep(2.0, p5)
w5        <- abs(rnorm(n5, mean = 1.5, sd = 0.3))
z5        <- as.vector(L5_true %*% beta5) + rnorm(n5, sd = 0.4)
EBeta5    <- rep(0, K);  EBeta25 <- rep(0.5, K)

res5 <- update_L_all(Y5, EL5_init, EL25_init, EF5_init, EF25_init,
                      Tau5, w5, z5, EBeta5, EBeta25)

cat("  Column-wise correlations with true L (absolute value):\n\n")
cat("  Factor   cor(EL[,k], L_true[,k])   Notes\n")
cat("  ------   -----------------------   -----\n")

# Compute best permutation by max-abs correlation
cor_mat <- matrix(0, K, K)
for (k in 1:K) for (k2 in 1:K) cor_mat[k, k2] <- abs(cor(res5$EL[, k], L5_true[, k2]))

for (k in 1:K) {
  best_match <- which.max(cor_mat[k, ])
  rho <- cor_mat[k, best_match]
  notes <- if (rho > 0.7) "good recovery" else if (rho > 0.4) "partial" else "weak"
  cat(sprintf("  GEP %-2d   %.4f (matches L_true[,%d])         %s\n",
              k, rho, best_match, notes))
}

cat(sprintf("\n  Second moment check: min(EL2 - EL^2) = %.2e\n",
            min(res5$EL2 - res5$EL^2)))
cat("  (all non-negative, confirming valid posterior variances)\n")


# ============================================================
# Done
# ============================================================
cat("\n")
sep()
cat(" Demo complete.  See code/update_L.R for implementation.\n")
cat(" See derivations/qL/qL_update_derivation.pdf for the math.\n")
sep()
