# =============================================================================
# demos/demo_update_tau.R
#
# Interactive demonstration of code/update_tau.R
#
# Run from project root:
#   Rscript demos/demo_update_tau.R
#
# Five scenarios illustrating the q(tau) update — a CLOSED-FORM MLE
# (no EBNM), with a critical variance correction term.
#
# Related files:
#   code/update_tau.R                              -- the functions demonstrated
#   derivations/qTau/qTau_update_derivation.pdf   -- the math behind them
#   tests/test_update_tau.R                        -- correctness tests
# =============================================================================

source("code/update_L.R")    # for compute_R_k (used in Demo 5)
source("code/update_F.R")    # for update_F_k (used in Demo 5)
source("code/update_tau.R")  # loads compute_var_term, compute_expected_residual_sq, update_tau

sep  <- function() cat(strrep("-", 60), "\n")
sec  <- function(title) { cat("\n"); sep(); cat(title, "\n"); sep() }
val  <- function(label, x, digits = 4) {
  cat(sprintf("  %-38s %s\n", paste0(label, ":"),
              paste(round(x, digits), collapse = "  ")))
}

cat("============================================================\n")
cat(" demo_update_tau.R\n")
cat(" Demonstrations of the q(tau_j) precision update\n")
cat("============================================================\n")

set.seed(42)
n <- 100
p <- 150
K <- 3


# ============================================================
# DEMO 1: Anatomy of Var_Term and R2_bar
# ============================================================
sec("DEMO 1: Anatomy of Var_Term and R2_bar")

cat("
Unlike L, F, and beta (all EBNM problems), tau is updated via a
closed-form MLE.  The expected squared residual is:

  R2_bar[i,j] = (Y[i,j] - sum_k l_bar_ik * f_bar_jk)^2   <- naive
              + Var_Term[i,j]                               <- correction

where the variance correction:
  Var_Term = EL2 %*% t(EF2)  -  EL^2 %*% t(EF^2)

accounts for posterior uncertainty in L and F.  Without it,
tau would be systematically overestimated (noise underestimated).
\n")

EL  <- matrix(rnorm(n * K), n, K)
EL2 <- EL^2 + 0.2       # second moment > squared mean
EF  <- matrix(rnorm(p * K), p, K)
EF2 <- EF^2 + 0.15
Y   <- EL %*% t(EF) + matrix(rnorm(n * p, sd = 0.8), n, p)

VT  <- compute_var_term(EL, EL2, EF, EF2)
naive_resid_sq <- (Y - EL %*% t(EF))^2
R2  <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)

cat(sprintf("  Var_Term dimensions:             %d x %d\n", nrow(VT), ncol(VT)))
cat(sprintf("  Var_Term range:                  min=%.4f  max=%.4f\n",
            min(VT), max(VT)))
cat(sprintf("  All Var_Term >= 0:               %s\n",
            if (all(VT >= -1e-10)) "YES" else "NO"))
cat(sprintf("  Mean (naive resid^2):            %.4f\n", mean(naive_resid_sq)))
cat(sprintf("  Mean Var_Term:                   %.4f\n", mean(VT)))
cat(sprintf("  Mean R2_bar (naive + correction): %.4f\n", mean(R2)))
cat(sprintf("  VarTerm / naive resid^2 ratio:   %.2f%%\n",
            100 * mean(VT) / mean(naive_resid_sq)))
cat("\n  The correction adds approximately that percentage to R2_bar,\n")
cat("  leading to a proportionally lower (more honest) tau estimate.\n")


# ============================================================
# DEMO 2: Known noise recovery
# ============================================================
sec("DEMO 2: Known noise recovery (true tau known)")

cat("
Simulate Y = L * F' + E where E[i,j] ~ N(0, 1/tau_j).
The estimated Tau should be close to the true tau_j.
\n")

tau_true <- runif(p, min = 0.5, max = 5.0)
sigma_j  <- 1 / sqrt(tau_true)

# Use zero posterior variance for clean comparison
EL_rec  <- matrix(rnorm(n * K), n, K)
EF_rec  <- matrix(rnorm(p * K), p, K)
EL2_rec <- EL_rec^2   # zero posterior variance
EF2_rec <- EF_rec^2
noise   <- matrix(rnorm(n * p), n, p) * rep(sigma_j, each = n)  # feature-specific noise
Y_rec   <- EL_rec %*% t(EF_rec) + noise

res_rec <- update_tau(Y_rec, EL_rec, EL2_rec, EF_rec, EF2_rec)

rho_tau  <- cor(res_rec$Tau, tau_true)
rmse_log <- sqrt(mean((log(res_rec$Tau) - log(tau_true))^2))

cat(sprintf("  n=%d, p=%d, K=%d\n", n, p, K))
cat(sprintf("  Cor(estimated tau, true tau):        %.4f\n", rho_tau))
cat(sprintf("  RMSE of log(tau_est) vs log(tau_true): %.4f\n", rmse_log))
cat(sprintf("\n  True tau range:     min=%.3f  max=%.3f\n", min(tau_true), max(tau_true)))
cat(sprintf("  Estimated tau range: min=%.3f  max=%.3f\n",
            min(res_rec$Tau), max(res_rec$Tau)))


# ============================================================
# DEMO 3: Variance correction matters
# ============================================================
sec("DEMO 3: Why the variance correction is critical")

cat("
If we ignore the variance correction (setting EL2=EL^2, EF2=EF^2),
we systematically OVERESTIMATE tau (underestimate noise).

This demo compares two tau estimates for the same data:
  (a) Correct:  uses EL2, EF2 (non-zero posterior variance)
  (b) Biased:   uses EL^2, EF^2 (pretends zero variance)
\n")

EL_c  <- matrix(rnorm(n * K), n, K)
EF_c  <- matrix(rnorm(p * K), p, K)
# Substantial posterior uncertainty
EL2_c <- EL_c^2 + 0.5
EF2_c <- EF_c^2 + 0.3
Y_c   <- EL_c %*% t(EF_c) + matrix(rnorm(n * p, sd = 0.7), n, p)

res_correct <- update_tau(Y_c, EL_c, EL2_c, EF_c, EF2_c)
res_biased  <- update_tau(Y_c, EL_c, EL_c^2, EF_c, EF_c^2)  # no correction

cat(sprintf("  Mean tau (correct, with Var_Term):   %.4f\n",
            mean(res_correct$Tau)))
cat(sprintf("  Mean tau (biased, no Var_Term):      %.4f\n",
            mean(res_biased$Tau)))
cat(sprintf("  Inflation ratio:                     %.3f x\n",
            mean(res_biased$Tau) / mean(res_correct$Tau)))
cat("\n  The biased estimate is systematically HIGHER — it wrongly attributes\n")
cat("  posterior uncertainty in L and F to lower noise.\n")
cat("  Var_Term inflates R2_bar, correctly lowering the tau estimate.\n")


# ============================================================
# DEMO 4: Heteroscedastic noise
# ============================================================
sec("DEMO 4: Heteroscedastic noise recovery")

cat("
Features 1-75 have high precision (tau=5.0, low noise).
Features 76-150 have low precision (tau=0.5, high noise).
The estimated tau should clearly separate the two groups.
\n")

p4       <- 150; n4 <- 120; K4 <- 2
tau_high <- rep(5.0, p4 / 2)
tau_low  <- rep(0.5, p4 / 2)
tau_true4 <- c(tau_high, tau_low)
sigma4    <- 1 / sqrt(tau_true4)

EL4  <- matrix(rnorm(n4 * K4), n4, K4)
EF4  <- matrix(rnorm(p4 * K4), p4, K4)
EL24 <- EL4^2;  EF24 <- EF4^2   # clean case for clarity
noise4 <- matrix(rnorm(n4 * p4), n4, p4) * rep(sigma4, each = n4)
Y4 <- EL4 %*% t(EF4) + noise4

res4 <- update_tau(Y4, EL4, EL24, EF4, EF24)

group_high <- 1:(p4 / 2);  group_low <- (p4 / 2 + 1):p4

cat(sprintf("  True high-precision group (features 1-%d):  tau_true = %.1f\n",
            p4 / 2, tau_high[1]))
cat(sprintf("    Estimated mean tau: %.3f  (SD: %.3f)\n",
            mean(res4$Tau[group_high]), sd(res4$Tau[group_high])))

cat(sprintf("\n  True low-precision group  (features %d-%d): tau_true = %.1f\n",
            p4 / 2 + 1, p4, tau_low[1]))
cat(sprintf("    Estimated mean tau: %.3f  (SD: %.3f)\n",
            mean(res4$Tau[group_low]), sd(res4$Tau[group_low])))

cat(sprintf("\n  Separation ratio (high / low): %.2f x  (true: %.1f x)\n",
            mean(res4$Tau[group_high]) / mean(res4$Tau[group_low]),
            tau_high[1] / tau_low[1]))


# ============================================================
# DEMO 5: ELBO proxy as convergence monitor
# ============================================================
sec("DEMO 5: ELBO proxy increases over manual L-F-tau iterations")

cat("
The tau update also computes the genomics ELBO proxy:
  ELBO_proxy = sum_j [ n/2 * log(tau_j) - tau_j/2 * sum_i R2_bar_ij ]

This should increase (or stay stable) as CAVI converges.
We run 3 manual iterations of L -> F -> tau and track the ELBO.
\n")

n5 <- 100; p5 <- 120; K5 <- 3
L5t <- matrix(rnorm(n5 * K5), n5, K5)
F5t <- matrix(rnorm(p5 * K5, sd = 0.7), p5, K5)
Y5  <- L5t %*% t(F5t) + matrix(rnorm(n5 * p5, sd = 0.6), n5, p5)

# Initialize at SVD
sv5 <- svd(Y5, nu = K5, nv = K5)
dk5 <- sqrt(pmax(sv5$d[1:K5], 0))
EL5  <- sv5$u %*% diag(dk5, K5, K5)
EF5  <- sv5$v %*% diag(dk5, K5, K5)
EL25 <- EL5^2
EF25 <- EF5^2 + 0.1
Tau5 <- 1 / pmax(apply(Y5, 2, var), 1e-8)
w5   <- rep(0, n5);  z5 <- rep(0, n5);  EBeta5 <- rep(0, K5); EBeta25 <- rep(0, K5)

elbo_history <- numeric(4)
res_tau_init <- update_tau(Y5, EL5, EL25, EF5, EF25)
elbo_history[1] <- res_tau_init$elbo_proxy

for (iter in 1:3) {
  # L update
  res_L <- update_L_all(Y5, EL5, EL25, EF5, EF25, Tau5, w5, z5, EBeta5, EBeta25)
  EL5 <- res_L$EL;  EL25 <- res_L$EL2

  # F update
  res_F <- update_F_all(Y5, EL5, EL25, EF5, EF25, Tau5)
  EF5 <- res_F$EF;  EF25 <- res_F$EF2

  # Tau update
  res_tau <- update_tau(Y5, EL5, EL25, EF5, EF25)
  Tau5 <- res_tau$Tau
  elbo_history[iter + 1] <- res_tau$elbo_proxy
}

cat("  Iteration   ELBO proxy    Change\n")
cat("  ---------   ----------    ------\n")
cat(sprintf("  Init        %+11.2f    --\n", elbo_history[1]))
for (iter in 1:3) {
  delta <- elbo_history[iter + 1] - elbo_history[iter]
  improving <- if (delta >= -1e-3) "improving" else "decreasing"
  cat(sprintf("  Iter %d      %+11.2f    %+.2f  (%s)\n",
              iter, elbo_history[iter + 1], delta, improving))
}

cat("\n  ELBO proxy should increase (or plateau) as the algorithm converges.\n")
cat("  A decreasing ELBO would indicate a bug in the update equations.\n")


# ============================================================
# Done
# ============================================================
cat("\n")
sep()
cat(" Demo complete.  See code/update_tau.R for implementation.\n")
cat(" See derivations/qTau/qTau_update_derivation.pdf for the math.\n")
sep()
