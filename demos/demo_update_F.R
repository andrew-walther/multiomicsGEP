# =============================================================================
# demos/demo_update_F.R
#
# Interactive demonstration of code/update_F.R
#
# Run from project root:
#   Rscript demos/demo_update_F.R
#
# Five scenarios illustrating the q(F) CAVI update — a VECTOR EBNM
# problem driven by the genomics likelihood only.
#
# KEY PROPERTY highlighted: tau_j CANCELS in x_j = B_F/A_F, but does
# NOT cancel in s_j = 1/sqrt(A_F).
#
# Related files:
#   code/update_F.R                           -- the functions being demonstrated
#   code/update_L.R                           -- provides compute_R_k
#   derivations/qF/qF_update_derivation.pdf  -- the math behind them
#   tests/test_update_F.R                    -- correctness tests
# =============================================================================

source("code/update_L.R")    # loads compute_R_k
source("code/update_F.R")    # loads update_F_k, update_F_all

sep  <- function() cat(strrep("-", 60), "\n")
sec  <- function(title) { cat("\n"); sep(); cat(title, "\n"); sep() }
val  <- function(label, x, digits = 4) {
  cat(sprintf("  %-36s %s\n", paste0(label, ":"),
              paste(round(x, digits), collapse = "  ")))
}

cat("============================================================\n")
cat(" demo_update_F.R\n")
cat(" Demonstrations of the q(f_{jk}) variational update\n")
cat("============================================================\n")

set.seed(42)
n <- 100   # patients
p <- 200   # genomic features


# ============================================================
# DEMO 1: The tau cancellation in x_j (but not in s_j)
# ============================================================
sec("DEMO 1: The tau cancellation property")

cat("
The F update is pure genomics:
  A_F[j] = tau_j * sum_i(E[l^2_ik])     <- p-vector
  B_F[j] = tau_j * sum_i(R_k[i,j] * l_bar_ik)  <- p-vector

A remarkable property: tau_j CANCELS in the EBNM pseudo-observation:
  x_j = B_F[j] / A_F[j]
      = [tau_j * (...)] / [tau_j * sum(EL2_k)]
      = (...) / sum(EL2_k)           <- tau-FREE!

But tau_j does NOT cancel in the pseudo-noise:
  s_j = 1 / sqrt(A_F[j]) = 1 / sqrt(tau_j * sum_EL2_k)
      -> features with higher tau (more precise) get smaller s_j
      -> less EBNM shrinkage -> posterior mean closer to x_j
\n")

EL_k   <- rnorm(n);   EL2_k  <- EL_k^2 + 0.1
R_k    <- matrix(rnorm(n * p), n, p)

Tau_low  <- rep(1.0, p)
Tau_high <- rep(1000.0, p)

res_low  <- update_F_k(Tau_low,  EL_k, EL2_k, R_k)
res_high <- update_F_k(Tau_high, EL_k, EL2_k, R_k)

cat("  Comparison: Tau=1 vs Tau=1000 (all features)\n\n")

# x_j should be identical
max_x_diff <- max(abs(res_low$x - res_high$x))
cat(sprintf("  Max |x_j(Tau=1) - x_j(Tau=1000)|:  %.2e  <- should be ~0\n",
            max_x_diff))

# s_j should differ
mean_s_low  <- mean(res_low$s)
mean_s_high <- mean(res_high$s)
cat(sprintf("  Mean s_j (Tau=1):                   %.4f\n", mean_s_low))
cat(sprintf("  Mean s_j (Tau=1000):                %.4f  <- 31.6x smaller\n",
            mean_s_high))
cat(sprintf("  Ratio s(Tau=1) / s(Tau=1000):       %.1f\n",
            mean_s_low / mean_s_high))

cat("\n  Larger tau -> smaller s_j -> less shrinkage in posterior:\n")
cat(sprintf("  Mean |posterior mean| (Tau=1):      %.4f\n",
            mean(abs(res_low$mean))))
cat(sprintf("  Mean |posterior mean| (Tau=1000):   %.4f\n",
            mean(abs(res_high$mean))))


# ============================================================
# DEMO 2: Signal recovery (K=1, sparse factor column)
# ============================================================
sec("DEMO 2: Signal recovery (K=1, sparse factor, 30/200 nonzero)")

cat("
Simulate Y = L_true * f_true' + noise, where f_true is sparse.
After one update_F_k call with the correct R_k and L, the posterior
should concentrate on the true nonzero features.
\n")

n_nz     <- 30
f_true   <- rep(0, p)
active   <- sample(p, n_nz)
f_true[active] <- rnorm(n_nz, sd = 1.5)

L_true   <- rnorm(n)
Y_rec    <- outer(L_true, f_true) + matrix(rnorm(n * p, sd = 0.4), n, p)
Tau_rec  <- rep(5.0, p)
EL_rec   <- L_true;   EL2_rec <- L_true^2 + 0.01
R_k_rec  <- Y_rec   # K=1, so R_k = Y

res_rec <- update_F_k(Tau_rec, EL_rec, EL2_rec, R_k_rec)

# Active features (true nonzero)
mean_active   <- mean(abs(res_rec$mean[active]))
mean_inactive <- mean(abs(res_rec$mean[-active]))
rho <- cor(res_rec$mean, f_true)

cat(sprintf("  True nonzero features:       %d / %d\n", n_nz, p))
cat(sprintf("  Correlation with true f:     %.4f\n", rho))
cat(sprintf("  Mean |est| on active feats:  %.4f\n", mean_active))
cat(sprintf("  Mean |est| on null feats:    %.4f\n", mean_inactive))
cat(sprintf("  Ratio (active / null):       %.1f x\n",
            mean_active / max(mean_inactive, 1e-8)))

cat("\n  Top-5 estimated features:\n")
top5 <- order(abs(res_rec$mean), decreasing = TRUE)[1:5]
cat("  Feature  True f    Estimated  In active set?\n")
for (j in top5) {
  in_active <- j %in% active
  cat(sprintf("  %7d  %+7.4f  %+9.4f  %s\n",
              j, f_true[j], res_rec$mean[j],
              if (in_active) "yes" else "NO"))
}


# ============================================================
# DEMO 3: Tau-dependent shrinkage
# ============================================================
sec("DEMO 3: Tau-dependent shrinkage across features")

cat("
Features with different noise precision tau_j experience different
amounts of shrinkage in the posterior, even if they have the same
true signal x_j.

This demo: split 200 features into 100 'clean' (tau=10) and 100 'noisy' (tau=0.5).
Both groups have the same true factor value.
\n")

f_common    <- rep(1.0, p)   # same true factor for all features
L_demo      <- rnorm(n)
Y_demo      <- outer(L_demo, f_common) + matrix(rnorm(n * p, sd = 1.0), n, p)

# Heteroscedastic precision
Tau_het <- c(rep(10.0, p / 2), rep(0.5, p / 2))  # first half clean, second noisy

EL_demo  <- L_demo;   EL2_demo <- L_demo^2 + 0.1
R_demo   <- Y_demo

res_het <- update_F_k(Tau_het, EL_demo, EL2_demo, R_demo)

clean_idx <- 1:(p / 2);   noisy_idx <- (p / 2 + 1):p

cat(sprintf("  Features 1-%d:   tau=10.0 (clean),  mean |est| = %.4f,  mean s = %.4f\n",
            p / 2,
            mean(abs(res_het$mean[clean_idx])),
            mean(res_het$s[clean_idx])))
cat(sprintf("  Features %d-%d: tau=0.5  (noisy),  mean |est| = %.4f,  mean s = %.4f\n",
            p / 2 + 1, p,
            mean(abs(res_het$mean[noisy_idx])),
            mean(res_het$s[noisy_idx])))
cat("\n  Clean features (larger tau) -> smaller s -> less shrinkage -> estimates closer to 1.\n")


# ============================================================
# DEMO 4: Sparsity recovery — confusion matrix
# ============================================================
sec("DEMO 4: Sparsity recovery")

cat("
With a point-normal prior, the posterior can assign exactly zero
mass to null features.  This demo checks how well update_F_k
identifies which features are active vs. null.
\n")

n_active <- 40
f_sparse <- rep(0, p)
active4  <- sample(p, n_active)
f_sparse[active4] <- rnorm(n_active, sd = 2.0)

L4   <- rnorm(n)
Y4   <- outer(L4, f_sparse) + matrix(rnorm(n * p, sd = 0.5), n, p)
Tau4 <- rep(3.0, p)
EL4  <- L4;   EL24 <- L4^2 + 0.05
R_k4 <- Y4

res4 <- update_F_k(Tau4, EL4, EL24, R_k4)

threshold <- 0.05   # treat |estimate| < threshold as 'detected as null'
est_active   <- abs(res4$mean) >= threshold
true_active  <- seq_len(p) %in% active4

TP <- sum(est_active  &  true_active)
FP <- sum(est_active  & !true_active)
FN <- sum(!est_active &  true_active)
TN <- sum(!est_active & !true_active)

cat(sprintf("  True active features: %d / %d (threshold: |est| >= %.2f)\n",
            n_active, p, threshold))
cat(sprintf("  True Positives  (correctly flagged active): %d\n", TP))
cat(sprintf("  False Positives (null flagged as active):   %d\n", FP))
cat(sprintf("  False Negatives (active flagged as null):   %d\n", FN))
cat(sprintf("  True Negatives  (correctly flagged null):   %d\n", TN))
cat(sprintf("  Precision: %.3f   Recall: %.3f\n",
            TP / max(TP + FP, 1), TP / max(TP + FN, 1)))


# ============================================================
# DEMO 5: Multi-factor (K=5) with update_F_all
# ============================================================
sec("DEMO 5: Multi-factor recovery (K=5) with update_F_all()")

cat("
update_F_all() loops k=1..K using Gauss-Seidel ordering.
This demo shows how well the factor columns are recovered
after one pass of update_F_all with near-correct loadings.
\n")

K  <- 5; n5 <- 120; p5 <- 250

# Simulate data
L5_true <- matrix(rnorm(n5 * K), n5, K)
F5_true <- matrix(rnorm(p5 * K, sd = 0.8), p5, K)
Y5      <- L5_true %*% t(F5_true) + matrix(rnorm(n5 * p5, sd = 0.6), n5, p5)

# Near-correct loadings
EL5   <- L5_true + matrix(rnorm(n5 * K, sd = 0.1), n5, K)
EL25  <- EL5^2 + 0.1

# Initialize F at zero
EF5   <- matrix(0, p5, K)
EF25  <- matrix(0, p5, K)
Tau5  <- rep(3.0, p5)

res5 <- update_F_all(Y5, EL5, EL25, EF5, EF25, Tau5)

cat("  Column-wise correlations with true F:\n\n")
cat("  Factor   cor(EF[,k], F_true[,k])   Active features (|est|>0.1)\n")
cat("  ------   -----------------------   ---------------------------\n")

for (k in 1:K) {
  rho      <- cor(res5$EF[, k], F5_true[, k])
  n_active <- sum(abs(res5$EF[, k]) > 0.1)
  cat(sprintf("  GEP %-2d   %+.4f                   %d / %d\n",
              k, rho, n_active, p5))
}

cat(sprintf("\n  Second moment check: min(EF2 - EF^2) = %.2e\n",
            min(res5$EF2 - res5$EF^2)))


# ============================================================
# Done
# ============================================================
cat("\n")
sep()
cat(" Demo complete.  See code/update_F.R for implementation.\n")
cat(" See derivations/qF/qF_update_derivation.pdf for the math.\n")
sep()
