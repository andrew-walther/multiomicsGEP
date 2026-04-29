# Fixing the L Update: Breaking the β = 0 Cycle

**Created:** 2026-04-29  
**Context:** SSBMF trained on merged TCGA\_PAAD + CPTAC produces β = 0 for all factors.  
**Status:** Diagnosis complete; fix not yet implemented. Start fresh session on this file.

---

## Diagnostic Trail (Summary)

| Experiment | Finding |
|---|---|
| Phase 1 heatmaps | All EBeta = 0 on merged; weak non-zero on single-cohort |
| Lambda sweep (λ ∈ {1,5,10,20} × 3 priors) | β = 0 across all 12 combinations; λ ≥ 5 collapses entire EL matrix to zero |
| EBMF diagnostic | 5/20 unsupervised factors Cox-significant (C-index up to 0.629) → data HAS survival signal |
| Warm-start Exp 1: β-only, EL fixed at EBMF | **β non-zero at iter 1**, 6/20 active → β update is functional |
| Warm-start Exp 2: full CAVI from EBMF init | β collapses to zero in 23 iters → **L update washes out survival signal** |

**Conclusion:** The β CAVI update is correct. The failure is in `update_L_k()`.

---

## The Mechanism

In `update_L_k()`, the EBNM precision for subject *i* on factor *k* is:

```
A_L[i] = (1 - α) * A_gen[i]  +  α * A_surv[i]

where:
  A_gen[i]  = sum_j  τ_j * E[f_jk²]        (sums over p = 2000 genes)
  A_surv[i] = λ * W_ii * E[β_k²]            (scalar; depends on current EBeta[k])
```

**The chicken-and-egg cycle:**

1. At initialization, EBeta ≈ 0 (Cox warm-start on SVD factors of the merged matrix
   is near-zero because SVD is dominated by the RNA-seq vs. proteomics batch factor)
2. A\_surv ≈ 0 → L update is entirely genomics-driven
3. L converges to reconstruction-optimal directions (batch-dominated)
4. β update sees batch-dominated factors → selects nothing → EBeta stays 0
5. Return to step 2

There is also a **scale imbalance**: A\_gen sums over p = 2000 genes, so
`A_gen >> A_surv` even when EBeta is modestly non-zero. The survival gradient
is structurally negligible throughout training, not just at initialization.

**Key evidence from Warm-start Exp 2:** Even initializing EL and EF from the
EBMF posterior means (which carry survival signal), the L update erodes the
structure within 23 iterations. The genomics objective wins at every step.

---

## Candidate Fixes

### Fix 1 — Reorder inner loop: update β before L  *(cheapest)*

**Current order** (per factor k): L_k → F_k → β_k  
**Proposed order**: β_k → L_k → F_k

When β is updated first, the L update for factor k uses the freshest EBeta[k]
from the same inner sweep, rather than the stale value from the previous outer
iteration. On the very first iteration, this means the L update sees whatever
β the Cox warm-start produced (rather than zero).

**Implementation:** Swap the three update blocks inside the k-loop in
`fit_supervised_mf_modular()` (`code/fit_modular.R`, lines ~307–360).

**Risk:** Only helps if the SVD-initialized EL has any survival signal. If the
Cox warm-start on SVD factors also gives EBeta ≈ 0 (likely on the merged
batch-dominated matrix), this alone won't break the cycle.

**Try this first** — it's 3 lines and a valid CAVI reordering.

---

### Fix 2 — β-only burn-in phase before the main CAVI loop  *(simple)*

Run N\_burnin iterations of β-only updates (exactly like Warm-start Exp 1)
using the SVD-initialized EL, *before* entering the joint CAVI loop. This
directly parallels what Exp 1 proved works.

**Implementation:** Add a pre-loop block to `fit_supervised_mf_modular()`:

```r
# β-only burn-in (N_burnin iterations) to break A_surv = 0 cycle
if (N_burnin > 0) {
  EL2_init <- EL^2   # point estimates; zero posterior variance
  for (b in seq_len(N_burnin)) {
    eta    <- as.vector(EL %*% EBeta)
    taylor <- calc_cox_taylor(eta, time, status)
    z      <- eta + taylor$u / taylor$w
    w      <- taylor$w
    for (k in seq_len(K)) {
      z_no_k  <- compute_z_no_k(z, EL, EBeta, k)
      res     <- update_beta_k(w, z_no_k, EL[, k], EL2_init[, k],
                               prior_family = prior_beta, alpha = alpha)
      EBeta[k]  <- res$mean
      EBeta2[k] <- res$second
    }
  }
}
```

A suggested default: N\_burnin = 5–10. Expose as a parameter with default 0
(off) to preserve backward compatibility.

**Risk:** Only works if the SVD loadings carry any survival signal. On the
merged cohort with a dominant batch factor, SVD factors may not, so EBeta may
still be near-zero after the burn-in. **Test empirically** by running with
N\_burnin = 10 on the merged v2 data and checking EBeta before iter 1.

**Combine with Fix 1** for maximum effect.

---

### Fix 3 — Replace multivariate Cox warm-start with ridge Cox  *(moderate)*

The current β warm-start (`fit_modular.R` line ~224) fits a multivariate Cox on
all K SVD loadings simultaneously. With K = 20 and n = 273, this is
borderline overfitted; standard `coxph()` can produce near-zero or degenerate
coefficients on collinear predictors.

**Proposed replacement:** Ridge-penalized Cox via `glmnet`:

```r
# Ridge Cox warm-start — forces non-zero EBeta even on correlated factors
library(glmnet)
cv_fit <- cv.glmnet(EL, Surv(time, status), family = "cox", alpha = 0)
EBeta  <- as.vector(coef(cv_fit, s = "lambda.min"))
```

Ridge does not enforce sparsity, so all K factors get non-zero EBeta. For
initialization purposes this is a feature: any non-zero value breaks the
A\_surv = 0 cycle.

**Risk:** Adds `glmnet` dependency (already likely installed). The ridge EBeta
may not be well-calibrated for the CAVI objective (different loss), but that's
fine — initialization quality, not accuracy, is what matters here.

**When to use:** If Fixes 1+2 fail because SVD loadings genuinely have no
survival signal on the merged cohort.

---

### Fix 4 — Normalize A\_surv and A\_gen to comparable scales in `update_L_k()`  *(most principled)*

Directly addresses the scale imbalance. Instead of raw precision contributions,
normalize each term by its mean across subjects before mixing:

```
A_gen_norm[i]  = A_gen[i]  / mean(A_gen)
A_surv_norm[i] = A_surv[i] / mean(A_surv)   [when mean(A_surv) > 0]

A_L[i] = (1 - α) * A_gen_norm[i]  +  α * A_surv_norm[i]
B_L[i] = (1 - α) * B_gen_norm[i]  +  α * B_surv_norm[i]   [same normalization]
```

This ensures α controls the *fraction of influence*, not just raw scale, so the
survival gradient is never structurally negligible regardless of EBeta's current
value.

**Risk:** Departs from strict ELBO maximization — the normalization introduces
an implicit reweighting of the variational objective that isn't derivable from
the generative model. Must verify ELBO still improves monotonically after the
change. Update `tests/test_update_L.R` to cover the new code path. The
normalization also changes what "α = 0.5" means in practice.

**When to use:** If Fixes 1–3 produce non-zero EBeta at first but β collapses
after a few CAVI iterations (i.e., the scale imbalance is the dominant issue,
not just initialization).

---

## Recommended Investigation Sequence

### Step 0 — Instrument first, fix second

Before implementing any fix, add logging to `fit_supervised_mf_modular()` to
capture at iteration 1:

1. **What is EBeta after the Cox warm-start (line ~224)?** If it's already
   non-zero but β still collapses, the problem is the L update scale. If it's
   zero, the initialization is the entry point.

2. **What is A\_surv / A\_gen at iter 1, factor k=6 (EBMF6 equivalent)?**
   This quantifies the scale imbalance directly.

```r
# Add inside the CAVI loop, iter == 1 only:
if (iter == 1 && verbose) {
  A_gen_k  <- sum(Tau * EF2[, k])
  A_surv_k <- mean(w) * EBeta[k]^2 * lambda
  cat(sprintf("    [iter1, k=%d] A_gen=%.2e  A_surv=%.2e  ratio=%.2e\n",
              k, A_gen_k, A_surv_k, A_surv_k / A_gen_k))
}
```

### Step 1 — Try Fix 1 + Fix 2 together

Reorder inner loop (β → L → F) AND add a β-only burn-in (N\_burnin = 10).
Check if EBeta is non-zero after burn-in and whether it stays non-zero through
the full CAVI.

### Step 2 — Add Fix 3 if Step 1 fails

Replace the Cox warm-start with ridge Cox. Re-run with N\_burnin = 0 first
(just the better init) to isolate the effect.

### Step 3 — Implement Fix 4 if collapse still occurs

If non-zero EBeta is achieved at initialization but consistently collapses
during CAVI, the scale imbalance is confirmed as the root cause. Implement
A\_surv / A\_gen normalization in `update_L_k()` and re-run the full test suite.

---

## Key Files

| File | Role |
|---|---|
| `code/fit_modular.R` lines ~197–235 | Initialization block; Cox warm-start; CAVI loop header |
| `code/fit_modular.R` lines ~307–360 | Inner k-loop with L → F → β order |
| `code/update_L.R` | `update_L_k()` — where A\_gen and A\_surv are computed and combined |
| `tests/test_update_L.R` | Must pass after any change to `update_L_k()` |
| `results/benchmark_sim/run_ebmf_warmstart.R` | Exp 1 (β-only) and Exp 2 (full warm-start) for re-running diagnostics |
| `results/benchmark_sim/outputs/ebmf_warmstart/tables/` | Warm-start results to compare against after fixes |

---

## Success Criteria

After the fix, on the merged TCGA\_PAAD + CPTAC v2-preprocessed training set:

- At least **1 factor with |EBeta| > 0.05** (any prior)
- Full CAVI **does not collapse EL** (max|EL| stays > 0.01)
- ELBO improves monotonically after the fix
- All **171/171 tests pass**
- External C-index on ≥1 cohort improves vs. SVD-initialized baseline
