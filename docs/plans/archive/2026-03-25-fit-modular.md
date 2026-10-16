# fit_modular.R Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create `code/fit_modular.R` — a reusable `fit_supervised_mf_modular()` function implementing V3 Algorithm 1 (factor-wise Gauss-Seidel CAVI) using the four existing modular update functions — paired with a Quarto companion doc `docs/fit_modular.qmd`.

**Architecture:** The R file contains (1) library/source preamble, (2) `calc_cox_taylor` helper copied verbatim from V2.R, (3) `fit_supervised_mf_modular()` function, and (4) a `DATA_MODE`-gated runner block. The function uses the `_k` variants of each update module inside a k-loop — unlike `run_modular_simulation.R` which uses `_all` block updates. The companion doc matches the style of `docs/update_beta.qmd` (cosmo HTML + xelatex PDF, rendered LaTeX math, consistent section headers).

**Tech Stack:** R, `ebnm`, `survival`, Quarto (xelatex + cosmo theme)

---

## File Map

| File | Action | Responsibility |
|------|--------|---------------|
| `code/fit_modular.R` | **Create** | `fit_supervised_mf_modular()` function + DATA_MODE runner |
| `docs/fit_modular.qmd` | **Create** | Companion doc: math, algorithm walkthrough, function reference, DATA_MODE, distinction from run_modular_simulation.R |
| `docs/fit_modular.pdf` | **Create** (rendered) | PDF output |
| `docs/fit_modular.html` | **Create** (rendered) | HTML output |
| `docs/Makefile` | **Modify** | Add `fit_modular.qmd` to `SRCS` |
| `CLAUDE.md` | **Modify** | Add `fit_modular.R` and `fit_modular.qmd` to quick-reference table |

---

## Task 1: Create `code/fit_modular.R` — preamble, helper, function signature

**Files:**
- Create: `code/fit_modular.R`

- [ ] **Step 1: Write the file preamble and `calc_cox_taylor`**

  Write `code/fit_modular.R` with this opening block. The `calc_cox_taylor` function is copied verbatim from `code/Supervised_Bayesian_MF_V2.R` lines 70–93 per the pattern established in `results/run_modular_simulation.R` lines 47–75. `Supervised_Bayesian_MF_V2.R` cannot be safely `source()`d because doing so executes its top-level simulation and plotting code.

  ```r
  # ==============================================================================
  # TITLE:  fit_supervised_mf_modular — Supervised Bayesian MF CAVI (Modular)
  # AUTHOR: Andrew Walther
  # DATE:   March 2026
  #
  # DESCRIPTION:
  #   Implements V3 Algorithm 1: factor-wise Gauss-Seidel CAVI for the joint model
  #     Genomics:  Y_{n x p} = L_{n x K} F^T_{K x p} + E,  E_{ij} ~ N(0, tau_j^{-1})
  #     Survival:  h(t_i | l_i) = h_0(t_i) exp( sum_k l_{ik} beta_k )   [Cox PH]
  #
  #   Each CAVI iteration runs Steps 1-4 from V3 Algorithm 1 (March 2026 derivation):
  #     Step 1: Cox Taylor expansion (once per outer iter)
  #     Step 2: for k=1..K — update_L_k -> recompute R^{-k} -> update_F_k -> update_beta_k
  #     Step 3: update_tau
  #     Step 4: max-absolute-change convergence check
  #
  #   KEY DISTINCTION FROM results/run_modular_simulation.R:
  #     run_modular_simulation.R uses BLOCK updates (_all variants): all-L -> all-F ->
  #     all-beta -> tau. This script uses FACTOR-WISE updates (_k variants): for each k,
  #     L_k -> F_k -> beta_k. Factor-wise is the true Gauss-Seidel CAVI as derived.
  #     Additionally, convergence here uses max absolute change (V3 Algorithm 1) vs.
  #     mean absolute change in run_modular_simulation.R.
  #     run_modular_simulation.R is a simulation benchmark script with hardcoded DGP and
  #     output figures/tables. This file is a reusable inference function with a
  #     DATA_MODE toggle.
  #
  # USAGE:
  #   source("code/fit_modular.R")          # loads function only (DATA_MODE="real")
  #   Rscript code/fit_modular.R            # runs simulated demo (DATA_MODE="simulated")
  #
  # DEPENDENCIES:
  #   source("code/update_L.R")     # compute_R_k, update_L_k, update_L_all
  #   source("code/update_F.R")     # update_F_k, update_F_all
  #   source("code/update_beta.R")  # compute_z_no_k, update_beta_k, update_beta_all
  #   source("code/update_tau.R")   # update_tau, compute_var_term
  #   Packages: ebnm, survival
  # ==============================================================================

  # ------------------------------------------------------------------------------
  # DATA MODE TOGGLE
  # Set DATA_MODE <- "simulated" to run on synthetic data (default for demo).
  # Set DATA_MODE <- "real"      to supply your own Y, time, status below.
  # ------------------------------------------------------------------------------
  DATA_MODE <- "simulated"

  # [Real data placeholders — fill these in when DATA_MODE = "real"]
  real_Y      <- NULL   # n x p numeric matrix (patients x features)
  real_time   <- NULL   # n-vector: survival/censoring times
  real_status <- NULL   # n-vector: event indicator (1=event, 0=censored)

  # ==============================================================================
  # PART 1: LIBRARIES & SOURCES
  # ==============================================================================

  library(survival)
  library(ebnm)

  source("code/update_L.R")      # compute_R_k, update_L_k, update_L_all
  source("code/update_F.R")      # update_F_k, update_F_all  (needs compute_R_k from L)
  source("code/update_beta.R")   # compute_z_no_k, update_beta_k, update_beta_all
  source("code/update_tau.R")    # compute_var_term, compute_expected_residual_sq, update_tau

  # ==============================================================================
  # Cox Taylor Expansion Helper
  # (Copied verbatim from code/Supervised_Bayesian_MF_V2.R lines 70-93)
  #
  # Transforms the non-conjugate Cox partial likelihood into a locally Gaussian
  # weighted-least-squares form centred at eta_hat = L_bar %*% beta_bar.
  # Working response:  z_i = eta_hat_i + u_i / W_{ii}
  # Weight:            W_{ii} (negative diagonal Hessian, positive)
  # ==============================================================================

  calc_cox_taylor <- function(eta, time, status) {
    n   <- length(time)
    ord <- order(time)
    time_s   <- time[ord]
    status_s <- status[ord]
    eta_s    <- eta[ord]

    theta    <- exp(eta_s)
    risk_sum <- rev(cumsum(rev(theta)))

    h <- status_s / risk_sum
    H <- cumsum(h)

    u_s <- status_s - theta * H
    w_s <- theta * H
    w_s[w_s < 1e-6] <- 1e-6

    u <- numeric(n); w <- numeric(n)
    u[ord] <- u_s;   w[ord] <- w_s
    list(u = u, w = w)
  }
  ```

- [ ] **Step 2: Add the `fit_supervised_mf_modular()` function signature and initialization**

  Append immediately after `calc_cox_taylor`:

  ```r
  # ==============================================================================
  # PART 2: FIT FUNCTION
  # ==============================================================================

  #' Fit Supervised Bayesian Matrix Factorization (Modular CAVI)
  #'
  #' Implements V3 Algorithm 1: factor-wise Gauss-Seidel CAVI using the four
  #' modular update functions. Each iteration: (1) Cox Taylor expansion, (2) for
  #' k=1..K: update_L_k -> update_F_k -> update_beta_k, (3) update_tau,
  #' (4) max-absolute-change convergence check.
  #'
  #' @param Y        n x p numeric matrix: genomics data (patients x features)
  #' @param time     n-vector: survival/censoring times
  #' @param status   n-vector: event indicator (1=event, 0=censored)
  #' @param K        integer: number of latent factors (default 5)
  #' @param max_iter integer: maximum CAVI iterations (default 100)
  #' @param tol      numeric: convergence threshold for max absolute change (default 1e-5)
  #' @param verbose  logical: print iteration summaries every 10 iters (default TRUE)
  #'
  #' @return list with fields:
  #'   EL      n x K posterior means of loadings
  #'   EL2     n x K posterior second moments of loadings
  #'   EF      p x K posterior means of factor weights
  #'   EF2     p x K posterior second moments of factor weights
  #'   EBeta   K-vector posterior means of survival coefficients
  #'   EBeta2  K-vector posterior second moments
  #'   Tau     p-vector noise precision (feature-specific)
  #'   history list(rmse, elbo_proxy, converged, n_iter)
  fit_supervised_mf_modular <- function(Y, time, status,
                                        K        = 5,
                                        max_iter = 100,
                                        tol      = 1e-5,
                                        verbose  = TRUE) {

    n <- nrow(Y); p <- ncol(Y)

    # --------------------------------------------------------------------------
    # Initialization (V3 Algorithm 1, Initialize block)
    # --------------------------------------------------------------------------

    # SVD of Y: deterministic high-variance starting subspace.
    # EL = U sqrt(D), EF = V sqrt(D) so EL %*% t(EF) = Y_rank-K approximation.
    svd_init <- svd(Y, nu = K, nv = K)
    d_k <- sqrt(pmax(svd_init$d[1:K], 0))
    EL  <- svd_init$u %*% diag(d_k, K, K)   # n x K
    EF  <- svd_init$v %*% diag(d_k, K, K)   # p x K

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
    EBeta2 <- EBeta^2

    # Column-specific precision tau_j: initialised from sample variance of Y.
    Tau <- 1.0 / pmax(apply(Y, 2, var), 1e-8)   # p-vector

    # History tracking
    history <- list(
      rmse       = numeric(max_iter),
      elbo_proxy = numeric(max_iter),
      converged  = FALSE,
      n_iter     = max_iter
    )

    if (verbose) {
      cat("=== fit_supervised_mf_modular (V3 Algorithm 1) ===\n")
      cat(sprintf("    n=%d, p=%d, K=%d | max_iter=%d | tol=%.1e\n\n",
                  n, p, K, max_iter, tol))
    }
  ```

- [ ] **Step 3: Verify the file parses without error**

  Run: `Rscript -e "source('code/fit_modular.R')"` from repo root.

  Expected: no error output (the function is incomplete — closing braces are missing — so add a temporary `}` at the end of the file for this check, then remove it).

  Actually, skip the parse-check here — proceed directly to adding the loop in Task 2, then run the full smoke test in Task 3.

---

## Task 2: Add CAVI loop and close `fit_supervised_mf_modular()`

**Files:**
- Modify: `code/fit_modular.R`

- [ ] **Step 1: Append the main CAVI loop body**

  Append the following inside `fit_supervised_mf_modular()`, after the verbose header:

  ```r
    # ==========================================================================
    # Main CAVI Loop  (V3 Algorithm 1, Steps 1-4)
    # ==========================================================================
    for (iter in 1:max_iter) {

      EL_old    <- EL
      EBeta_old <- EBeta

      # Step 1 ----------------------------------------------------------------
      # Compute survival working quantities once per outer iteration.
      # eta_hat_i = sum_k l_bar_{ik} * beta_bar_k
      # z_i = eta_hat_i + u_i / W_{ii}   (Cox working response)
      # w_i = W_{ii}                      (positive Cox weight)
      # -----------------------------------------------------------------------
      eta    <- as.vector(EL %*% EBeta)
      taylor <- calc_cox_taylor(eta, time, status)
      z      <- eta + taylor$u / taylor$w
      w      <- taylor$w

      # Reconstruction RMSE at posterior means (monitoring only, not convergence)
      history$rmse[iter] <- sqrt(mean((Y - EL %*% t(EF))^2))

      # Step 2 ----------------------------------------------------------------
      # Factor-wise coordinate ascent: for k = 1, ..., K
      # Order within each k: L_k -> F_k -> beta_k  (V3 Algorithm 1 Step 2a-c)
      # This is true Gauss-Seidel CAVI: each sub-update uses the most current
      # posterior means, including updates from earlier k' < k in this iteration.
      # -----------------------------------------------------------------------
      for (k in 1:K) {

        # (a) Update q(l_k) — Patient Loadings  (V3 Eqs. 30-33)
        #
        # Partial residual R^{-k}_{ij} = Y_{ij} - sum_{k'!=k} l_bar_{ik'} f_bar_{jk'}
        # compute_R_k uses the CURRENT EL, EF (incorporates updates from k' < k).
        R_k <- compute_R_k(Y, EL, EF, k)

        # Partial survival working response:
        # z^{-k}_i = z_i - sum_{k'!=k} l_bar_{ik'} * beta_bar_{k'}
        #
        # KEY (from update_beta.R docstring lines 50-53):
        # z_no_k does NOT depend on EL[,k] or EBeta[k]. compute_z_no_k computes
        # (EL %*% EBeta - EL[,k]*EBeta[k]), so EL[,k] appears then immediately
        # cancels out. Therefore z_no_k is identical before and after the L_k
        # update. Compute ONCE here and reuse for both step (a) and step (c).
        # (This matches V2.R lines 294-295 and the "z_no_k reuse rationale"
        # in REVISED.tex Sec. 6.)
        z_no_k <- compute_z_no_k(z, EL, EBeta, k)

        res_L   <- update_L_k(Tau, EF[, k], EF2[, k], w, EBeta[k], EBeta2[k],
                               R_k, z_no_k)
        EL[, k]  <- res_L$mean
        EL2[, k] <- res_L$second

        # (b) Update q(f_k) — Biological Factors  (V3 Eqs. 43-46)
        #
        # Recompute R^{-k} using the just-updated EL[,k] before calling update_F_k.
        # R_k DOES depend on EL[,k] (unlike z_no_k), so recomputation is needed.
        # This propagates the fresh l_bar_{ik} into the F update within the same k.
        R_k <- compute_R_k(Y, EL, EF, k)

        res_F   <- update_F_k(Tau, EL[, k], EL2[, k], R_k)
        EF[, k]  <- res_F$mean
        EF2[, k] <- res_F$second

        # (c) Update q(beta_k) — Survival Coefficients  (V3 Eqs. 53-56)
        #
        # Reuse z_no_k from step (a) — it is unchanged by the L_k and F_k updates
        # (see comment above). Uses updated EL[,k] and EL2[,k] for A_k and B_k.
        res_beta <- update_beta_k(w, z_no_k, EL[, k], EL2[, k])
        EBeta[k]  <- res_beta$mean
        EBeta2[k] <- res_beta$second

      }  # end k-loop

      # Step 3 ----------------------------------------------------------------
      # Update tau (V3 Eq. 64): column-specific noise precision.
      # update_tau computes the variance-corrected expected squared residual and
      # returns the ELBO proxy (genomics log-likelihood term).
      # -----------------------------------------------------------------------
      res_tau            <- update_tau(Y, EL, EL2, EF, EF2)
      Tau                <- res_tau$Tau
      history$elbo_proxy[iter] <- res_tau$elbo_proxy

      # Step 4 ----------------------------------------------------------------
      # Convergence check (V3 Algorithm 1 Step 4).
      # Uses MAX absolute change — not mean — consistent with V3 Algorithm 1.
      # Guard iter > 5 prevents early termination on initialization noise.
      # -----------------------------------------------------------------------
      delta_L    <- max(abs(EL - EL_old))
      delta_Beta <- max(abs(EBeta - EBeta_old))

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

    list(
      EL     = EL,
      EL2    = EL2,
      EF     = EF,
      EF2    = EF2,
      EBeta  = EBeta,
      EBeta2 = EBeta2,
      Tau    = Tau,
      history = history
    )

  }  # end fit_supervised_mf_modular
  ```

---

## Task 3: Add DATA_MODE runner block

**Files:**
- Modify: `code/fit_modular.R`

- [ ] **Step 1: Append the simulated-data runner at the bottom of the file**

  ```r
  # ==============================================================================
  # PART 3: DATA_MODE RUNNER
  # ==============================================================================
  # When DATA_MODE = "simulated": generates synthetic genomics + survival data,
  # calls fit_supervised_mf_modular(), and prints a concise summary.
  # When DATA_MODE = "real": uses real_Y, real_time, real_status set at the top.
  #
  # Simulation parameters match run_modular_simulation.R (same seed, n, p, K, B_true)
  # to enable direct comparison of factor-wise vs. blockwise CAVI convergence.
  # ==============================================================================

  if (DATA_MODE == "simulated") {

    # --- Simulation parameters ------------------------------------------------
    set.seed(42)
    n <- 250; p <- 1000; K <- 5

    B_true <- c(1.5, -1.2, 0.8, -0.5, 0.0)   # true survival coefficients

    # --- Generate genomics data -----------------------------------------------
    L_true <- matrix(rnorm(n * K), n, K)
    F_true <- matrix(0, p, K)
    for (k in 1:K) {
      active <- sample(1:p, round(p * 0.05))   # 5% sparse factor structure
      F_true[active, k] <- rnorm(length(active), 0, 5)
    }
    Y <- L_true %*% t(F_true) + matrix(rnorm(n * p), n, p)

    # --- Generate survival data (Weibull baseline) ----------------------------
    eta_true   <- as.vector(L_true %*% B_true)
    raw_times  <- (-log(runif(n)) / (0.01 * exp(eta_true)))^(1 / 1.5)
    cens_times <- rexp(n, rate = 1 / 50)
    time   <- pmin(raw_times, cens_times)
    status <- as.integer(raw_times <= cens_times)

    cat(sprintf("  Simulated data: n=%d, p=%d, K=%d, seed=42\n", n, p, K))
    cat(sprintf("  Censoring rate: %.1f%%\n\n", 100 * mean(status == 0)))

    # --- Fit ------------------------------------------------------------------
    res <- fit_supervised_mf_modular(Y, time, status, K = K, max_iter = 100,
                                     tol = 1e-5, verbose = TRUE)

    # --- Summary --------------------------------------------------------------
    cat("\n=== Results ===\n")
    cat(sprintf("  Converged: %s  (iter %d)\n",
                res$history$converged, res$history$n_iter))
    cat(sprintf("  Final RMSE: %.4f  (expect ~1.0 = true noise SD)\n",
                tail(res$history$rmse, 1)))
    cat("\n  Beta recovery:\n")
    cat(sprintf("    True:      %s\n",
                paste(sprintf("%+.2f", B_true), collapse = ", ")))
    cat(sprintf("    Estimated: %s\n",
                paste(sprintf("%+.2f", res$EBeta), collapse = ", ")))
    cat("    Signs match:", all(sign(res$EBeta[-5]) == sign(B_true[-5])), "\n")

  } else if (DATA_MODE == "real") {

    stopifnot(
      !is.null(real_Y),
      !is.null(real_time),
      !is.null(real_status)
    )
    res <- fit_supervised_mf_modular(real_Y, real_time, real_status,
                                     K = 5, max_iter = 100,
                                     tol = 1e-5, verbose = TRUE)
    cat("Fit complete. Access results via res$EL, res$EF, res$EBeta, etc.\n")

  }
  ```

---

## Task 4: Smoke-test `fit_modular.R`

**Files:**
- Read: `code/fit_modular.R` (verify structure before running)

- [ ] **Step 1: Run the script from repo root**

  ```bash
  cd /Users/ajwalther/GithubProjects/multiomicsGEP/.claude/worktrees/infallible-bassi
  Rscript code/fit_modular.R
  ```

  Expected output (approximate):

  ```
  === fit_supervised_mf_modular (V3 Algorithm 1) ===
      n=250, p=1000, K=5 | max_iter=100 | tol=1.0e-05

    iter  10 | RMSE: ...  | ELBO: ...  | dL: ...  | dB: ...  | beta: [...]
    ...
    Converged at iteration XX  (dL=..., dBeta=...)

  === Results ===
    Converged: TRUE  (iter XX)
    Final RMSE: ~1.0
    Beta recovery:
      True:      +1.50, -1.20, +0.80, -0.50, +0.00
      Estimated: <should match signs: +, -, +, -, ~0>
    Signs match: TRUE
  ```

  If RMSE is far from 1.0 or signs don't match, investigate before proceeding. Common issues:
  - `compute_R_k` called with wrong argument order → check `update_L.R` signature
  - `update_beta_k` called with wrong `z_no_k` → recheck step 2(c) code

- [ ] **Step 2: Verify ELBO proxy is non-decreasing**

  ```bash
  Rscript -e "
  source('code/fit_modular.R')
  if (DATA_MODE == 'simulated') {
    elbo <- res$history$elbo_proxy
    diffs <- diff(elbo)
    cat('ELBO non-decreasing:', all(diffs >= -1), '\n')
    cat('Min ELBO diff:', min(diffs), '\n')
  }
  "
  ```

  Expected: `ELBO non-decreasing: TRUE` (small numerical noise ≤ −1 is acceptable).

- [ ] **Step 3: Commit `fit_modular.R`**

  ```bash
  git add code/fit_modular.R
  git commit -m "Add fit_supervised_mf_modular(): factor-wise CAVI loop (V3 Algorithm 1)

  Implements V3 Algorithm 1 from derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_3_11_26_V3.pdf
  using the four existing modular update functions (update_L_k, update_F_k,
  update_beta_k, update_tau).

  Key differences from results/run_modular_simulation.R:
  - Factor-wise updates (_k variants in k-loop) vs. block updates (_all variants)
  - Max absolute convergence (V3 Algorithm 1) vs. mean absolute convergence
  - Reusable inference function with DATA_MODE toggle vs. hardcoded benchmark script

  DATA_MODE='simulated' verifies sign recovery of beta=(+,-,+,-,0) and RMSE~1.0."
  ```

---

## Task 5: Write `docs/fit_modular.qmd`

**Files:**
- Create: `docs/fit_modular.qmd`

The companion doc follows the exact structure and YAML of `docs/update_beta.qmd`. All math uses Quarto LaTeX (`$...$` inline, `$$...$$` display). Section headers match the convention from the other companion docs.

- [ ] **Step 1: Write the YAML frontmatter and Overview section**

  ```markdown
  ---
  title: "Companion Document: fit_modular.R"
  subtitle: "Full CAVI Loop — Supervised Matrix Factorization (V3 Algorithm 1)"
  author: "Andrew Walther"
  date: today
  format:
    html:
      toc: true
      toc-depth: 3
      toc-location: left
      number-sections: true
      theme: cosmo
      smooth-scroll: true
      self-contained: true
    pdf:
      toc: true
      toc-depth: 3
      number-sections: true
      latex-engine: xelatex
      keep-tex: false
      fig-pos: "htbp"
      include-in-header:
        text: |
          \usepackage{booktabs}
          \usepackage{longtable}
          \usepackage{amsmath}
          \usepackage{amssymb}
          \usepackage{geometry}
          \geometry{margin=1in}
          \renewcommand{\familydefault}{\sfdefault}
  ---

  ## Overview

  `code/fit_modular.R` implements the complete **Coordinate Ascent Variational Inference
  (CAVI)** loop for the Supervised Bayesian Matrix Factorization model, following **V3
  Algorithm 1** from `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_3_11_26_V3.pdf`.

  The public function `fit_supervised_mf_modular()` ties together the four modular update
  scripts — `update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R` — into a complete
  inference procedure. It accepts a genomics matrix $Y$ and paired survival data $(t_i,
  \delta_i)$, and returns variational posterior moments for all model parameters.

  ### Model

  $$Y_{n \times p} = L_{n \times K} F_{K \times p}^\top + E, \qquad
    E_{ij} \overset{\text{iid}}{\sim} \mathcal{N}(0, \tau_j^{-1})$$

  $$h(t_i) = h_0(t_i) \exp\!\left(\sum_k l_{ik} \beta_k\right) \qquad \text{(Cox PH)}$$

  CAVI finds the mean-field variational posterior
  $q(L, F, \beta) = q(L)\,q(F)\,q(\beta)$ that maximises the Evidence Lower Bound (ELBO).
  Each coordinate update is formulated as an **Empirical Bayes Normal Means (EBNM)** problem
  with a point-normal prior, solved via the `ebnm` R package.
  ```

- [ ] **Step 2: Add Mathematical Background section**

  ```markdown
  ## Mathematical Background

  ### Cox Taylor Approximation

  The Cox partial likelihood $\log P(t, \delta \mid L, \beta)$ is non-conjugate. A
  second-order Taylor expansion around the current linear predictor $\hat{\eta}_i =
  \sum_k \bar{l}_{ik} \bar{\beta}_k$ gives

  $$\log P(t, \delta \mid L, \beta) \approx
    \sum_i \left[ -\tfrac{1}{2} W_{ii} (\eta_i - z_i)^2 \right] + C$$

  where $W_{ii} = -\partial^2 \ell_{\text{Cox}} / \partial \eta_i^2 \big|_{\hat{\eta}} \geq 0$
  is the negative diagonal Hessian, and

  $$z_i = \hat{\eta}_i + \frac{u_i}{W_{ii}}, \qquad
    u_i = \frac{\partial \ell_{\text{Cox}}}{\partial \eta_i}\bigg|_{\hat{\eta}}$$

  is the **working response**. This converts the survival term into a Gaussian pseudo-likelihood,
  enabling the same EBNM machinery used for the genomics updates.

  ### Partial Residuals

  For each factor $k$, define the **genomics partial residual**

  $$R_{ij}^{-k} = Y_{ij} - \sum_{k' \neq k} \bar{l}_{ik'} \bar{f}_{jk'}$$

  and the **survival partial working response**

  $$z_i^{-k} = z_i - \sum_{k' \neq k} \bar{l}_{ik'} \bar{\beta}_{k'}$$

  Both remove the contribution of all factors except $k$, isolating the signal for
  factor $k$'s update. Crucially, $z_i^{-k}$ does **not** depend on $l_{ik}$ or
  $\beta_k$, so it can be reused for both the $q(l_k)$ and $q(\beta_k)$ updates within
  the same $k$.

  ### EBNM Coefficient Forms

  Each update reduces to $\text{EBNM}(x, s)$ with point-normal prior. The inputs per
  update are:

  | Parameter | $A$ (precision) | $B$ (signal) | $x = B/A$ | $s = 1/\sqrt{A}$ |
  |-----------|----------------|--------------|-----------|-----------------|
  | $q(l_k)$ | $\sum_j \tau_j \overline{f_{jk}^2} + W_{ii} \overline{\beta_k^2}$ | $\sum_j \tau_j R_{ij}^{-k} \bar{f}_{jk} + W_{ii} z_i^{-k} \bar{\beta}_k$ | $n$-vector | $n$-vector |
  | $q(f_k)$ | $\tau_j \sum_i \overline{l_{ik}^2}$ | $\tau_j \sum_i R_{ij}^{-k} \bar{l}_{ik}$ | $p$-vector | $p$-vector |
  | $q(\beta_k)$ | $\sum_i W_{ii} \overline{l_{ik}^2}$ | $\sum_i W_{ii} z_i^{-k} \bar{l}_{ik}$ | scalar | scalar |

  Note that $\overline{l_{ik}^2} = \text{Var}_q(l_{ik}) + \bar{l}_{ik}^2$ (posterior
  second moment, not squared mean) — this **error-in-variables correction** prevents
  $\beta$ from overfitting to uncertain loadings.

  ### Convergence Criterion

  Following V3 Algorithm 1, convergence uses the **maximum** absolute change:

  $$\Delta_L = \max_{i,k} \left| \bar{l}_{ik}^{\text{new}} - \bar{l}_{ik}^{\text{old}} \right|,
  \qquad
  \Delta_\beta = \max_k \left| \bar{\beta}_k^{\text{new}} - \bar{\beta}_k^{\text{old}} \right|$$

  Stop when $\max(\Delta_L, \Delta_\beta) < \varepsilon$ after at least 5 iterations.
  This is stricter than the mean absolute change used in `run_modular_simulation.R`
  and `Supervised_Bayesian_MF_V2.R`.
  ```

- [ ] **Step 3: Add Algorithm Walkthrough section**

  ```markdown
  ## Algorithm Walkthrough

  The following maps V3 Algorithm 1 to the R implementation step by step.

  ### Initialization

  ```r
  svd_init <- svd(Y, nu = K, nv = K)
  d_k <- sqrt(pmax(svd_init$d[1:K], 0))
  EL  <- svd_init$u %*% diag(d_k, K, K)   # n x K
  EF  <- svd_init$v %*% diag(d_k, K, K)   # p x K
  EL2 <- EL^2;  EF2 <- EF^2
  ```

  This gives $\bar{L}^{(0)} = U \sqrt{D}$, $\bar{F}^{(0)} = V \sqrt{D}$ so that
  $\bar{L}^{(0)} (\bar{F}^{(0)})^\top = Y_{\text{rank-}K}$. Second moments are
  initialised to squared means (zero posterior variance), which is updated at the
  first EBNM call.

  $\bar{\beta}^{(0)}$ is warm-started via `coxph()` on the SVD loadings. $\hat{\tau}_j^{(0)}$
  is initialised from the column-wise sample variance of $Y$.

  ### Step 1: Cox Working Quantities

  Computed once per outer iteration using current $\bar{L}$, $\bar{\beta}$:

  ```r
  eta    <- as.vector(EL %*% EBeta)           # eta_hat_i = sum_k l_bar_ik * beta_bar_k
  taylor <- calc_cox_taylor(eta, time, status)
  z      <- eta + taylor$u / taylor$w         # working response z_i
  w      <- taylor$w                          # Cox weights W_ii
  ```

  ### Step 2: Factor-wise Loop (k = 1, ..., K)

  For each $k$, three sub-updates run in sequence, each using the most current
  posterior means (Gauss-Seidel):

  **(a) Update $q(l_k)$:** Compute `z_no_k` **once** here — it will be reused for step (c):
  ```r
  R_k    <- compute_R_k(Y, EL, EF, k)           # n x p partial residual
  z_no_k <- compute_z_no_k(z, EL, EBeta, k)     # computed once; reused for beta update
  res_L  <- update_L_k(Tau, EF[,k], EF2[,k], w, EBeta[k], EBeta2[k], R_k, z_no_k)
  EL[,k]  <- res_L$mean;  EL2[,k] <- res_L$second
  ```

  **(b) Update $q(f_k)$:** $R^{-k}$ depends on `EL[,k]` (unlike `z_no_k`), so recompute
  it with the fresh loadings before calling `update_F_k`:
  ```r
  R_k   <- compute_R_k(Y, EL, EF, k)            # recomputed: R_k depends on EL[,k]
  res_F <- update_F_k(Tau, EL[,k], EL2[,k], R_k)
  EF[,k]  <- res_F$mean;  EF2[,k] <- res_F$second
  ```

  **(c) Update $q(\beta_k)$:** Reuse `z_no_k` from step (a). It is provably unchanged:
  `compute_z_no_k` computes `EL %*% EBeta - EL[,k]*EBeta[k]`, so `EL[,k]` appears then
  cancels — the result is independent of `EL[,k]`. Uses updated `EL[,k]`, `EL2[,k]`
  from step (a) for $A_k$ and $B_k$:
  ```r
  # z_no_k reused from step (a) — unchanged by L_k and F_k updates
  res_beta <- update_beta_k(w, z_no_k, EL[,k], EL2[,k])
  EBeta[k]  <- res_beta$mean;  EBeta2[k] <- res_beta$second
  ```

  ### Step 3: Update $\hat{\tau}$

  ```r
  res_tau            <- update_tau(Y, EL, EL2, EF, EF2)
  Tau                <- res_tau$Tau
  history$elbo_proxy[iter] <- res_tau$elbo_proxy
  ```

  `update_tau` computes the variance-corrected expected squared residual
  $\overline{R^2}_{ij}$ and the column-specific MLE $\hat{\tau}_j = n / \sum_i
  \overline{R^2}_{ij}$.

  ### Step 4: Convergence

  ```r
  delta_L    <- max(abs(EL - EL_old))
  delta_Beta <- max(abs(EBeta - EBeta_old))
  if (iter > 5 && delta_L < tol && delta_Beta < tol) break
  ```
  ```

- [ ] **Step 4: Add Function Reference, DATA_MODE Toggle, Distinction, and Related Files sections**

  ```markdown
  ## Function Reference

  ### `fit_supervised_mf_modular`

  ```r
  fit_supervised_mf_modular(Y, time, status, K = 5, max_iter = 100,
                             tol = 1e-5, verbose = TRUE)
  ```

  **Arguments:**

  | Argument | Type | Description |
  |----------|------|-------------|
  | `Y` | `n x p` matrix | Genomics data matrix |
  | `time` | `n`-vector | Survival/censoring times |
  | `status` | `n`-vector | Event indicator (1=event, 0=censored) |
  | `K` | integer | Number of latent factors |
  | `max_iter` | integer | Maximum CAVI iterations |
  | `tol` | numeric | Convergence threshold (max absolute change) |
  | `verbose` | logical | Print iteration summaries every 10 iters |

  **Returns:** A list with fields:

  | Field | Dimension | Description |
  |-------|-----------|-------------|
  | `EL` | `n x K` | Posterior means $\bar{L}$ |
  | `EL2` | `n x K` | Posterior second moments $\overline{L^2}$ |
  | `EF` | `p x K` | Posterior means $\bar{F}$ |
  | `EF2` | `p x K` | Posterior second moments $\overline{F^2}$ |
  | `EBeta` | `K`-vector | Posterior means $\bar{\beta}$ |
  | `EBeta2` | `K`-vector | Posterior second moments $\overline{\beta^2}$ |
  | `Tau` | `p`-vector | Noise precisions $\hat{\tau}_j$ |
  | `history` | list | `rmse`, `elbo_proxy` (per iter), `converged`, `n_iter` |

  **Note on naming:** Return fields use the `E`-prefix convention (`EL`, `EF`, `EBeta`)
  matching the internal naming of the modular update functions. This differs intentionally
  from `fit_supervised_mf()` in `Supervised_Bayesian_MF_V2.R`, which returns short names
  (`L`, `F`, `Beta`).

  ## DATA_MODE Toggle

  At the top of `fit_modular.R`, `DATA_MODE` controls what runs when the script is
  executed directly:

  ```r
  DATA_MODE <- "simulated"   # default: run on synthetic data
  # DATA_MODE <- "real"      # uncomment to use real data
  ```

  **Simulated path** (`DATA_MODE = "simulated"`): generates a synthetic dataset with
  $n=250$ patients, $p=1000$ features, $K=5$ factors, and true survival coefficients
  $\beta = (+1.5, -1.2, +0.8, -0.5, 0.0)$. Uses the same seed and DGP as
  `results/run_modular_simulation.R` to enable direct comparison of factor-wise vs.
  blockwise convergence. Prints a summary showing RMSE convergence and $\beta$ sign
  recovery.

  **Real data path** (`DATA_MODE = "real"`): set `real_Y`, `real_time`, and
  `real_status` at the top of the file, then `Rscript code/fit_modular.R`.

  ## Distinction from `run_modular_simulation.R`

  Both files use the four modular update functions, but serve different purposes and
  implement different algorithmic variants:

  | | `results/run_modular_simulation.R` | `code/fit_modular.R` |
  |---|---|---|
  | **Update order** | Block: all-$L$ → all-$F$ → all-$\beta$ → $\tau$ | Factor-wise: for $k$: $L_k$ → $F_k$ → $\beta_k$ |
  | **Update functions** | `_all` variants | `_k` variants |
  | **Convergence** | Mean absolute change | **Max** absolute change (V3 Algorithm 1) |
  | **Purpose** | Simulation benchmark; exports figures and tables | Reusable inference function |
  | **Data** | Hardcoded simulation DGP | `DATA_MODE` toggle (simulated or real) |
  | **Output** | Figures to `results/figures/`, tables to `results/tables/` | Returns a list of posteriors |

  **Why factor-wise is the "correct" CAVI:** In factor-wise updates, after updating
  $q(l_k)$, the new $\bar{l}_{ik}$ is immediately used when computing $R^{-k}$ for
  $q(f_k)$, and the new $\bar{l}_{ik}$ is used in $A_k$ for $q(\beta_k)$. This is
  true Gauss-Seidel coordinate ascent — each sub-update uses the most current
  information — and is precisely what V3 Algorithm 1 specifies. Blockwise updates
  (as in `run_modular_simulation.R`) are a valid but less tight approximation that
  updates all factors of one parameter type before moving to the next.

  ## Related Files

  | File | Role |
  |------|------|
  | `code/update_L.R` | $q(L)$ update: `update_L_k`, `compute_R_k` |
  | `code/update_F.R` | $q(F)$ update: `update_F_k` |
  | `code/update_beta.R` | $q(\beta)$ update: `update_beta_k`, `compute_z_no_k` |
  | `code/update_tau.R` | $\tau$ update: `update_tau` |
  | `code/Supervised_Bayesian_MF_V2.R` | Monolithic V2 implementation (reference) |
  | `results/run_modular_simulation.R` | Block-update simulation benchmark |
  | `docs/update_L.qmd` | Companion doc for $q(L)$ update |
  | `docs/update_F.qmd` | Companion doc for $q(F)$ update |
  | `docs/update_beta.qmd` | Companion doc for $q(\beta)$ update |
  | `docs/update_tau.qmd` | Companion doc for $\tau$ update |
  | `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_3_11_26_V3.pdf` | V3 Algorithm 1 derivation |
  ```

---

## Task 6: Render companion doc and update Makefile

**Files:**
- Modify: `docs/Makefile`
- Modify: `docs/fit_modular.qmd` (if any render errors need fixing)

- [ ] **Step 1: Update `docs/Makefile` to include `fit_modular.qmd`**

  Change line 1 of `docs/Makefile` from:
  ```makefile
  SRCS  := update_beta.qmd update_L.qmd update_F.qmd update_tau.qmd
  ```
  to:
  ```makefile
  SRCS  := update_beta.qmd update_L.qmd update_F.qmd update_tau.qmd fit_modular.qmd
  ```

- [ ] **Step 2: Render PDF and HTML**

  ```bash
  cd /Users/ajwalther/GithubProjects/multiomicsGEP/.claude/worktrees/infallible-bassi/docs
  quarto render fit_modular.qmd --to pdf
  quarto render fit_modular.qmd --to html
  ```

  Expected: `fit_modular.pdf` and `fit_modular.html` created in `docs/`.

  If xelatex errors occur, check:
  - Missing LaTeX packages → run `tlmgr install <package>` via TinyTeX
  - Unescaped `%` or `_` in prose → fix in the `.qmd`
  - Math parse errors → check `$$` delimiters are closed

- [ ] **Step 3: Open and visually inspect HTML output**

  ```bash
  open docs/fit_modular.html
  ```

  Verify:
  - Left sidebar TOC with numbered sections
  - Math renders (not raw LaTeX strings)
  - Comparison table in "Distinction from run_modular_simulation.R" renders correctly
  - Code blocks have `r` syntax highlighting

---

## Task 7: Deprecation note on `results/run_modular_simulation.R`

**Files:**
- Modify: `results/run_modular_simulation.R` (header comment only)

`run_modular_simulation.R` uses block updates (`_all` variants) and is a simulation benchmark
script, not a reusable inference function. Now that `fit_modular.R` exists as the canonical
implementation of V3 Algorithm 1, `run_modular_simulation.R` is no longer aligned with the
project's goals. Add a deprecation notice to its header so future readers understand its status.

- [ ] **Step 1: Add deprecation notice to `results/run_modular_simulation.R` header**

  After the opening comment block (after line ~20, before the `setwd` block), insert:

  ```r
  # ==============================================================================
  # DEPRECATION NOTE (March 2026)
  #
  # This script uses BLOCK coordinate ascent (_all variants): all-L -> all-F ->
  # all-beta -> tau. This ordering was used to showcase the _all helper functions
  # and produce a quick benchmark, but does NOT faithfully implement V3 Algorithm 1,
  # which specifies FACTOR-WISE updates (_k variants): for k=1..K { L_k, F_k, beta_k }.
  #
  # The canonical modular implementation of V3 Algorithm 1 is:
  #   code/fit_modular.R  (fit_supervised_mf_modular)
  #
  # This script is retained as a historical reference for the _all function APIs
  # and for comparison against factor-wise convergence behavior, but is no longer
  # the recommended way to run the supervised MF model.
  # ==============================================================================
  ```

- [ ] **Step 2: Commit the deprecation note**

  ```bash
  git add results/run_modular_simulation.R
  git commit -m "Mark run_modular_simulation.R as deprecated

  Now that fit_modular.R implements V3 Algorithm 1 (factor-wise CAVI) as the canonical
  modular implementation, run_modular_simulation.R's block-update approach is no longer
  aligned with the project goals. Adds a deprecation header noting the distinction and
  pointing to fit_modular.R as the replacement."
  ```

---

## Task 8: Update `CLAUDE.md` and commit everything

**Files:**
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add new entries to the quick-reference table in `CLAUDE.md`**

  In `CLAUDE.md`, add four rows to the Quick Reference table:

  ```
  | Modular CAVI loop (V3 Algorithm 1) | `code/fit_modular.R` |
  | Companion doc for full CAVI loop | `docs/fit_modular.qmd` |
  | Factor-wise simulation runner | `results/run_factor_modular_simulation.R` |
  | Factor-wise simulation report | `results/factor_modular_sim_report.qmd` |
  ```

  Also add a note under the V2 simulation report row that `results/modular_sim_report.qmd`
  is deprecated (blockwise) in favour of `results/factor_modular_sim_report.qmd`.

- [ ] **Step 2: Commit all remaining files**

  ```bash
  git add docs/fit_modular.qmd docs/fit_modular.pdf docs/fit_modular.html \
          docs/Makefile CLAUDE.md
  git commit -m "Add fit_modular.qmd companion doc and update Makefile + CLAUDE.md

  Companion doc for code/fit_modular.R following the same Quarto style as the
  update_*.qmd docs (cosmo HTML, xelatex PDF, rendered LaTeX math, numbered sections).

  Sections cover: Overview, Mathematical Background (Taylor expansion, EBNM
  coefficient forms, max convergence criterion), Algorithm Walkthrough (step-by-step
  mapping of V3 Algorithm 1 to R), Function Reference, DATA_MODE Toggle, Distinction
  from run_modular_simulation.R, Related Files.

  Makefile updated to include fit_modular.qmd in SRCS so 'make all' renders it.
  CLAUDE.md quick-reference table updated with fit_modular.R and fit_modular.qmd."
  ```

---

## Task 9: Write `results/run_factor_modular_simulation.R`

**Files:**
- Create: `results/run_factor_modular_simulation.R`
- Create dirs: `results/figures/factor_modular_sim/`, `results/tables/factor_modular_sim/`

This script mirrors `results/run_modular_simulation.R` in structure and DGP (same seed, n, p, K, B_true), but calls `fit_supervised_mf_modular()` from `code/fit_modular.R` instead of the block-update loop. It exports the same figure and table outputs to `factor_modular_sim/` subdirectories so all three simulation variants are grouped separately:

```
results/figures/full_sim/            ← V2 monolithic (existing)
results/figures/modular_sim/         ← blockwise modular (existing, deprecated)
results/figures/factor_modular_sim/  ← NEW: factor-wise modular
results/tables/full_sim/             ← existing
results/tables/modular_sim/          ← existing, deprecated
results/tables/factor_modular_sim/   ← NEW
```

- [ ] **Step 1: Create output directories**

  ```bash
  mkdir -p results/figures/factor_modular_sim
  mkdir -p results/tables/factor_modular_sim
  ```

- [ ] **Step 2: Write `results/run_factor_modular_simulation.R`**

  Structure mirrors `results/run_modular_simulation.R` exactly:

  ```r
  # ==============================================================================
  # Run Factor-Wise Modular Simulation and Export All Results
  #
  # Uses fit_supervised_mf_modular() from code/fit_modular.R, which implements
  # V3 Algorithm 1: factor-wise Gauss-Seidel CAVI (_k update variants).
  #
  # WHY FACTOR-WISE IS MORE CORRECT THAN BLOCK-WISE:
  #   Block coordinate ascent (as in run_modular_simulation.R) updates all K
  #   columns of L before updating any column of F or beta. Factor-wise CAVI
  #   updates L_k, F_k, and beta_k together for each k before moving to k+1,
  #   immediately propagating each sub-update into subsequent ones. This is
  #   true Gauss-Seidel coordinate ascent — it is the update order specified
  #   in V3 Algorithm 1 and in the EBMF paper — and typically converges in
  #   fewer iterations because each sub-update sees the most current state.
  #   run_modular_simulation.R is retained for reference but is deprecated as
  #   the canonical simulation benchmark.
  #
  # Outputs:
  #   results/tables/factor_modular_sim/  -- CSV tables
  #   results/figures/factor_modular_sim/ -- figures (PDF + PNG)
  #   results/factor_modular_sim_report.qmd/.pdf/.html
  #
  # Run from repo root:
  #   Rscript results/run_factor_modular_simulation.R
  # ==============================================================================

  # [working directory setup identical to run_modular_simulation.R]
  if (Sys.getenv("REPO_ROOT") != "") {
    setwd(Sys.getenv("REPO_ROOT"))
  } else if (file.exists("code/update_L.R")) {
    # already at repo root
  } else if (file.exists("../code/update_L.R")) {
    setwd("..")
  } else {
    stop("Cannot find repo root. Run from project root or set REPO_ROOT env var.")
  }

  library(survival)
  library(ebnm)

  # Source fit_modular.R with DATA_MODE="real" so it only loads the function,
  # does not run the simulated-data block
  DATA_MODE <- "real"
  source("code/fit_modular.R")

  # [same DGP as run_modular_simulation.R: set.seed(42), n=250, p=1000, K=5,
  #  B_true = c(1.5, -1.2, 0.8, -0.5, 0.0), identical L_true/F_true/Y/time/status]

  # [call fit_supervised_mf_modular(Y, time, status, K=K)]

  # [export same 7 CSV tables and 8 figures as run_modular_simulation.R,
  #  written to results/tables/factor_modular_sim/ and results/figures/factor_modular_sim/]
  ```

  Copy the analytics helpers (`get_cindex_comparison`, `get_top_features`,
  `get_factor_summary_table`), DGP block, and figure/table export code verbatim from
  `run_modular_simulation.R`, substituting:
  - The CAVI loop with a single call to `fit_supervised_mf_modular()`
  - Output paths `modular_sim/` → `factor_modular_sim/`
  - Script header comment noting factor-wise vs. block-wise distinction

- [ ] **Step 3: Run the script**

  ```bash
  Rscript results/run_factor_modular_simulation.R
  ```

  Expected:
  - Script completes without error
  - CSV tables written to `results/tables/factor_modular_sim/`
  - Figures written to `results/figures/factor_modular_sim/`
  - Console output shows Converged: TRUE, RMSE ~1.0, β signs match

- [ ] **Step 4: Commit runner script and outputs**

  ```bash
  git add results/run_factor_modular_simulation.R \
          results/tables/factor_modular_sim/ \
          results/figures/factor_modular_sim/
  git commit -m "Add run_factor_modular_simulation.R and factor-wise simulation outputs

  Runner script for factor-wise modular CAVI (V3 Algorithm 1) using
  fit_supervised_mf_modular(). Produces same figures and tables as the deprecated
  run_modular_simulation.R but with factor-wise updates (_k variants).

  Results organized into results/figures/factor_modular_sim/ and
  results/tables/factor_modular_sim/ alongside existing full_sim/ and modular_sim/."
  ```

---

## Task 10: Write and render `results/factor_modular_sim_report.qmd`

**Files:**
- Create: `results/factor_modular_sim_report.qmd`
- Create: `results/factor_modular_sim_report.pdf` (rendered)
- Create: `results/factor_modular_sim_report.html` (rendered)

This report mirrors `results/modular_sim_report.qmd` in structure (same YAML, same sections, same figure/table inclusion pattern) but:
1. Sources `code/fit_modular.R` (factor-wise) instead of the inline block-update loop
2. Adds a brief "Algorithm" section near the top noting why factor-wise CAVI is the correct implementation
3. Does **not** include a side-by-side comparison with blockwise — focuses solely on factor-wise results

- [ ] **Step 1: Write `results/factor_modular_sim_report.qmd`**

  Use the same YAML frontmatter as `results/modular_sim_report.qmd`:

  ```yaml
  ---
  title: "Supervised Bayesian Matrix Factorization"
  subtitle: "Factor-Wise Modular Implementation — Simulation Report"
  author: "Andrew Walther"
  date: today
  format:
    html:
      toc: true
      toc-depth: 3
      toc-location: left
      number-sections: true
      theme: cosmo
      smooth-scroll: true
      self-contained: true
    pdf:
      toc: true
      toc-depth: 3
      number-sections: true
      geometry: "margin=1in"
      mainfont: "Helvetica Neue"
      monofont: "Menlo"
      colorlinks: true
      fig-cap-location: bottom
      tbl-cap-location: top
      keep-tex: false
      fig-pos: "H"
      include-in-header:
        text: |
          \usepackage{booktabs}
          \usepackage{longtable}
          \usepackage{float}
          \usepackage{amsmath}
          \usepackage{amssymb}
          \floatplacement{figure}{H}
          \floatplacement{table}{H}
  execute:
    echo: false
    message: false
    warning: false
    cache: false
  ---
  ```

  Section structure (mirroring `modular_sim_report.qmd`):
  1. **Overview** — model, data parameters, brief note that this uses factor-wise CAVI
  2. **Algorithm Note** — 2–3 sentences on why factor-wise is the correct CAVI order vs. blockwise; reference V3 Algorithm 1
  3. **Convergence** — RMSE trace, ELBO proxy trace (figures)
  4. **Factor Summary** — β estimates, log-rank p-values, sparsity, PVE (table + KM curves)
  5. **Signal Recovery** — estimated vs. true β (figure)
  6. **GEP Heatmap** — top features per factor (figure)
  7. **C-Index** — supervised latent space vs. PCA (table)
  8. **Top Features** — per-factor feature weight tables

  Source pre-computed outputs from `results/figures/factor_modular_sim/` and
  `results/tables/factor_modular_sim/` rather than re-running the simulation inside the .qmd
  (use `knitr::include_graphics()` for figures, `read.csv()` + `kable()` for tables).

- [ ] **Step 2: Render PDF and HTML**

  ```bash
  quarto render results/factor_modular_sim_report.qmd --to pdf
  quarto render results/factor_modular_sim_report.qmd --to html
  ```

  Expected: `factor_modular_sim_report.pdf` and `factor_modular_sim_report.html` created.

- [ ] **Step 3: Commit report and rendered outputs**

  ```bash
  git add results/factor_modular_sim_report.qmd \
          results/factor_modular_sim_report.pdf \
          results/factor_modular_sim_report.html
  git commit -m "Add factor_modular_sim_report: simulation report for factor-wise CAVI

  Quarto report (PDF + HTML) for the factor-wise modular simulation, mirroring the
  structure of modular_sim_report.qmd. Sources pre-computed figures/tables from
  results/figures/factor_modular_sim/ and results/tables/factor_modular_sim/.

  Includes a brief Algorithm Note section explaining why factor-wise updates are the
  correct CAVI implementation vs. the deprecated blockwise approach."
  ```

---

## Updated Expected Final State

```
code/
  fit_modular.R                         ← NEW

docs/
  fit_modular.qmd/.pdf/.html            ← NEW
  Makefile                              ← UPDATED

results/
  run_factor_modular_simulation.R       ← NEW
  factor_modular_sim_report.qmd/.pdf/.html  ← NEW
  run_modular_simulation.R              ← DEPRECATED (header note added)
  figures/
    full_sim/                           ← existing (V2 monolithic)
    modular_sim/                        ← existing (blockwise, deprecated)
    factor_modular_sim/                 ← NEW (factor-wise)
  tables/
    full_sim/                           ← existing
    modular_sim/                        ← existing, deprecated
    factor_modular_sim/                 ← NEW

CLAUDE.md                               ← UPDATED
```
