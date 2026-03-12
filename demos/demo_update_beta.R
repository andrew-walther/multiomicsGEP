# =============================================================================
# demos/demo_update_beta.R
#
# Interactive demonstration of code/update_beta.R
#
# Run from project root:
#   Rscript demos/demo_update_beta.R
#
# Shows four scenarios, each illustrating a different aspect of the q(beta)
# CAVI update.  No survival data or full model needed — just vectors.
#
# Related files:
#   code/update_beta.R                       -- the functions being demonstrated
#   derivations/qB/qBeta_update_derivation.pdf -- the math behind them
#   tests/test_update_beta.R                 -- correctness tests
# =============================================================================

source("code/update_beta.R")   # loads update_beta_k, update_beta_all, ebnm

sep  <- function() cat(strrep("-", 60), "\n")
sec  <- function(title) { cat("\n"); sep(); cat(title, "\n"); sep() }
val  <- function(label, x, digits = 4) {
  cat(sprintf("  %-30s %s\n", paste0(label, ":"), paste(round(x, digits), collapse = "  "))
  )
}

cat("============================================================\n")
cat(" demo_update_beta.R\n")
cat(" Demonstrations of the q(beta_k) variational update\n")
cat("============================================================\n")

set.seed(42)
n <- 200   # number of patients used throughout


# ============================================================
# DEMO 1: What do A_k, B_k, x_k, s_k actually represent?
# ============================================================
sec("DEMO 1: Anatomy of a single update_beta_k() call")

cat("
The update reduces to one EBNM problem with two ingredients:
  A_k = sum(w * EL2_k)          <- precision  (how confident we are)
  B_k = sum(w * z_no_k * EL_k)  <- signal     (weighted correlation)
  x_k = B_k / A_k               <- pseudo-observation for EBNM
  s_k = 1 / sqrt(A_k)           <- pseudo-noise for EBNM
\n")

w      <- rep(2.0, n)                   # flat Cox weights
EL_k   <- rnorm(n)                      # loadings for factor k
EL2_k  <- EL_k^2 + 0.1                 # second moment = mean^2 + variance
z_no_k <- 1.5 * EL_k + rnorm(n, 0.3)  # partial working response (true beta=1.5)

res <- update_beta_k(w, z_no_k, EL_k, EL2_k)

val("A_k  (precision)",       res$A)
val("B_k  (signal)",          res$B)
val("x_k  = B/A (naive est)", res$x)
val("s_k  = 1/sqrt(A)",       res$s, digits = 5)
val("Posterior mean",         res$mean)
val("Posterior SD",           res$sd, digits = 5)
val("Second moment (sd^2+mean^2)", res$second)

cat("\n  Note: posterior mean is slightly shrunk toward 0 vs x_k = B/A\n")
cat("  because point-normal prior applies soft thresholding.\n")


# ============================================================
# DEMO 2: Signal recovery — does the estimate find the true beta?
# ============================================================
sec("DEMO 2: Signal recovery (K = 1, true beta = 1.5)")

cat("
We simulate z_no_k = beta_true * l_ik + noise, then call update_beta_k.
The posterior mean should be close to 1.5.
\n")

beta_true <- 1.5
z_no_k    <- beta_true * EL_k + rnorm(n, sd = 0.3)

res_pos <- update_beta_k(w, z_no_k, EL_k, EL2_k)
val("True beta",      beta_true)
val("x_k (naive)",   res_pos$x)
val("Posterior mean", res_pos$mean)
val("Posterior SD",   res_pos$sd, digits = 5)

cat("\n  Repeat with true beta = -1.2:\n\n")
z_no_k    <- -1.2 * EL_k + rnorm(n, sd = 0.3)
res_neg <- update_beta_k(w, z_no_k, EL_k, EL2_k)
val("True beta",      -1.2)
val("Posterior mean", res_neg$mean)


# ============================================================
# DEMO 3: Null factor shrinkage — beta = 0 should stay near zero
# ============================================================
sec("DEMO 3: Null factor shrinkage (true beta = 0)")

cat("
When z_no_k is pure noise (no linear trend with EL_k), the point-normal
prior shrinks the estimate to zero.  Compare:
  (a) z_no_k = signal  (beta=1.5) -> posterior mean far from 0
  (b) z_no_k = noise   (beta=0)   -> posterior mean near 0
\n")

z_signal <- 1.5 * EL_k + rnorm(n, sd = 0.3)
z_noise  <- rnorm(n)                        # pure noise, no trend

res_signal <- update_beta_k(w, z_signal, EL_k, EL2_k)
res_null   <- update_beta_k(w, z_noise,  EL_k, EL2_k)

val("(a) Signal: posterior mean", res_signal$mean)
val("(b) Null:   posterior mean", res_null$mean)
cat("\n  The null factor is shrunk close to 0 by the point-normal prior.\n")


# ============================================================
# DEMO 4: Error-in-variables — posterior uncertainty in L affects beta
# ============================================================
sec("DEMO 4: Error-in-variables effect")

cat("
If we are uncertain about the loadings EL_k, we should be more cautious
about beta_k.  This is the 'error-in-variables' correction:

  A_k = sum(w * EL2_k)   uses the FULL second moment E_q[l^2]
                         = Var_q(l) + mean^2

When posterior variance in L is high, A_k is larger, x_k = B/A is
smaller, and the estimate is more shrunk toward zero.
\n")

z_no_k <- 1.5 * EL_k + rnorm(n, sd = 0.3)

EL2_low  <- EL_k^2 + 0.01   # low posterior variance in L  (certain loadings)
EL2_high <- EL_k^2 + 1.00   # high posterior variance in L (uncertain loadings)

res_low  <- update_beta_k(w, z_no_k, EL_k, EL2_low)
res_high <- update_beta_k(w, z_no_k, EL_k, EL2_high)

cat("  Low  posterior var(L):", 0.01, "\n")
val("    A_k",          res_low$A)
val("    x_k = B/A",    res_low$x)
val("    Posterior mean", res_low$mean)

cat("\n  High posterior var(L):", 1.00, "\n")
val("    A_k",          res_high$A)
val("    x_k = B/A",    res_high$x)
val("    Posterior mean", res_high$mean)

cat("\n  Higher loading uncertainty -> larger A_k -> smaller |x_k| -> more shrinkage.\n")
cat("  This prevents overfitting beta to loadings we are not sure about.\n")


# ============================================================
# DEMO 5: Multi-factor (K = 5) with update_beta_all()
# ============================================================
sec("DEMO 5: Multi-factor recovery (K = 5)")

cat("
update_beta_all() loops k = 1..K using Gauss-Seidel ordering:
each updated beta_k is immediately used when computing z_no_k'
for subsequent factors.

True beta = (1.5, -1.2, 0.8, -0.5, 0.0)
\n")

K         <- 5
beta_true <- c(1.5, -1.2, 0.8, -0.5, 0.0)

EL    <- matrix(rnorm(n * K), n, K)
EL2   <- EL^2 + 0.1
z     <- as.vector(EL %*% beta_true) + rnorm(n, sd = 0.5)
w5    <- abs(rnorm(n, mean = 2, sd = 0.4))
EBeta <- rep(0, K)   # cold start

res_all <- update_beta_all(w5, z, EL, EL2, EBeta)

cat("  Factor  True beta  Estimated  Sign OK?\n")
cat("  ------  ---------  ---------  --------\n")
for (k in 1:K) {
  sign_ok <- if (beta_true[k] == 0) {
    if (abs(res_all$EBeta[k]) < 0.2) "yes (shrunk)" else "NO (not shrunk)"
  } else {
    if (sign(res_all$EBeta[k]) == sign(beta_true[k])) "yes" else "NO"
  }
  cat(sprintf("  GEP %-2d  %+6.2f     %+7.4f    %s\n",
              k, beta_true[k], res_all$EBeta[k], sign_ok))
}

cat("\n  Second moments (sd^2 + mean^2):\n")
for (k in 1:K) {
  cat(sprintf("  GEP %-2d  EBeta2 = %.4f  (sd = %.4f)\n",
              k, res_all$EBeta2[k], res_all$details[[k]]$sd))
}


# ============================================================
# Done
# ============================================================
cat("\n")
sep()
cat(" Demo complete.  See code/update_beta.R for implementation.\n")
cat(" See derivations/qB/qBeta_update_derivation.pdf for the math.\n")
sep()
