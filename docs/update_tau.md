# update_tau.R — Companion Document

**Module:** `code/update_tau.R`
**Project:** Supervised Bayesian Matrix Factorization (multiomicsGEP)
**Role:** Closed-form MLE update for per-feature noise precisions q(τ)

---

## Table of Contents

1. [Overview](#overview)
2. [Mathematical Background](#mathematical-background)
3. [Function Reference](#function-reference)
   - [compute_var_term](#compute_var_term)
   - [compute_expected_residual_sq](#compute_expected_residual_sq)
   - [update_tau](#update_tau)
4. [Key Properties](#key-properties)
5. [Test Suite Summary](#test-suite-summary)
6. [Demo Scenarios](#demo-scenarios)
7. [Pipeline Context](#pipeline-context)
8. [Related Files](#related-files)

---

## Overview

`update_tau.R` implements the variational update for τ = (τ₁, …, τ_p), the per-feature noise precisions in the model

```
Y = L F' + E,    E_ij ~ N(0, 1/τ_j)
```

This module is **unique among the four update modules** in two important ways:

1. **No EBNM dependency** — the update is a closed-form maximum likelihood estimate, not an empirical Bayes normal means problem. The module uses pure base R with no external packages.
2. **No per-k loop** — all K factors are handled simultaneously through matrix operations. The entire update is computed in a handful of matrix multiplications.

The module exports three functions:

| Function | Returns | Purpose |
|---|---|---|
| `compute_var_term` | n × p matrix | Posterior variance correction term |
| `compute_expected_residual_sq` | n × p matrix | Expected squared residual under q |
| `update_tau` | list (Tau, R2_bar, Var_Term, elbo_proxy) | Full τ update + ELBO proxy |

---

## Mathematical Background

### The CAVI Objective for τ_j

Taking the variational expectation of the complete-data log-likelihood with respect to q(L) and q(F), the optimal q(τ_j) is Gamma-distributed. Under a flat (improper) prior on τ_j, the CAVI fixed-point reduces to a column-wise MLE:

```
τ_j = n / Σᵢ R̄²_ij
```

where R̄²_ij is the **expected squared residual** under the current posterior q(L) ⊗ q(F).

### Expected Squared Residual

The expected squared residual cannot be computed naively as (Y - EL EF')². Posterior uncertainty in L and F inflates the true expectation:

```
R̄²_ij = E_q[(Y_ij - Σ_k L_ik F_jk)²]
        = (Y_ij - Σ_k EL_ik · EF_jk)²  +  Var_Term_ij
```

The second term, `Var_Term`, accounts for the fact that E[(LF')²] ≠ (EL · EF')² when L and F have nonzero posterior variance.

### Variance Correction (Var_Term)

Under the mean-field factorization q(L) ⊗ q(F) with independent factors:

```
Var_Term_ij = Σ_k E[L²_ik] E[F²_jk]  −  Σ_k (EL_ik)² (EF_jk)²
            = (EL2 · EF2')_ij  −  (EL² · (EF²)')_ij
```

where EL2_ik = E[L²_ik] and EF2_jk = E[F²_jk] are posterior second moments (not squared means).

This simplification uses independence across factors and across the q(L) ⊗ q(F) factorization. By Jensen's inequality, E[X²] ≥ (E[X])² for any random variable X, so:

```
Var_Term_ij ≥ 0    for all i, j
```

Var_Term is exactly zero only when **both** EL2 = EL² **and** EF2 = EF² — i.e., all posterior variances are identically zero (a degenerate point-mass posterior).

### Why Ignoring the Correction Is Dangerous

If one naively substitutes EL2 → EL² and EF2 → EF² (ignoring posterior uncertainty), the effective R̄²_ij is **underestimated**, causing τ_j to be **overestimated** (noise underestimated). In practice this bias is 3–5× in typical simulations (see Demo 3).

### ELBO Proxy

The module also computes an ELBO proxy (the τ-dependent terms of the evidence lower bound):

```
ELBO_proxy = Σ_j [ n/2 · log(τ_j)  −  τ_j/2 · Σᵢ R̄²_ij ]
```

At the MLE solution τ_j = n / Σᵢ R̄²_ij, this is maximized over τ. Monitoring ELBO_proxy across iterations provides a convergence diagnostic: it should be non-decreasing. A decrease signals a bug in the implementation.

---

## Function Reference

### compute_var_term

```r
compute_var_term(EL, EL2, EF, EF2)
```

**Arguments**

| Argument | Shape | Description |
|---|---|---|
| `EL` | n × K | Posterior means E[L] |
| `EL2` | n × K | Posterior second moments E[L²] |
| `EF` | p × K | Posterior means E[F] |
| `EF2` | p × K | Posterior second moments E[F²] |

**Returns** an n × p matrix.

**Formula**

```r
(EL2 %*% t(EF2)) - (EL^2 %*% t(EF^2))
```

**Notes**
- All entries are ≥ 0 by Jensen's inequality.
- Entries are zero when the corresponding column of EL (or EF) has zero posterior variance, i.e., EL2[,k] = EL[,k]^2.
- This function is a pure computation — no side effects, no floor applied.

---

### compute_expected_residual_sq

```r
compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
```

**Arguments**

| Argument | Shape | Description |
|---|---|---|
| `Y` | n × p | Observed genomics matrix |
| `EL` | n × K | Posterior means E[L] |
| `EL2` | n × K | Posterior second moments E[L²] |
| `EF` | p × K | Posterior means E[F] |
| `EF2` | p × K | Posterior second moments E[F²] |

**Returns** an n × p matrix R̄².

**Formula**

```r
(Y - EL %*% t(EF))^2 + compute_var_term(EL, EL2, EF, EF2)
```

**Notes**
- The first term is the squared residual at posterior means (what naive MLE would use).
- The second term inflates it to account for posterior uncertainty.
- All entries are ≥ 0 because both summands are non-negative.

---

### update_tau

```r
update_tau(Y, EL, EL2, EF, EF2, tau_floor = 1e-8)
```

**Arguments**

| Argument | Shape | Default | Description |
|---|---|---|---|
| `Y` | n × p | — | Observed genomics matrix |
| `EL` | n × K | — | Posterior means E[L] |
| `EL2` | n × K | — | Posterior second moments E[L²] |
| `EF` | p × K | — | Posterior means E[F] |
| `EF2` | p × K | — | Posterior second moments E[F²] |
| `tau_floor` | scalar | 1e-8 | Minimum τ value (prevents Tau = Inf) |

**Returns** a named list:

| Element | Shape | Description |
|---|---|---|
| `Tau` | p-vector | Estimated noise precisions τ_1, …, τ_p |
| `R2_bar` | n × p | Expected squared residuals |
| `Var_Term` | n × p | Variance correction term |
| `elbo_proxy` | scalar | ELBO proxy (τ-dependent terms) |

**Algorithm**

```r
Var_Term <- compute_var_term(EL, EL2, EF, EF2)
R2_bar   <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)
RSS_j    <- colSums(R2_bar)                          # length-p vector
Tau      <- n / pmax(RSS_j, n * tau_floor)           # MLE with floor
elbo_proxy <- sum(n/2 * log(Tau) - Tau/2 * RSS_j)
```

**Numerical stability**

`pmax(RSS_j, n * tau_floor)` prevents division by near-zero RSS when the reconstruction is nearly perfect. Without this floor, τ_j → ∞ (infinite precision), which can destabilize downstream L and F updates. The floor corresponds to a maximum τ of 1/tau_floor = 1e8.

**No per-k loop**

Unlike `update_beta_all`, `update_L_all`, and `update_F_all`, this function processes all K factors in a single pair of matrix multiplications. The sum over k in the Var_Term formula collapses into matrix products.

---

## Key Properties

| Property | Statement |
|---|---|
| Update type | Closed-form MLE (NOT EBNM) |
| Loop structure | No per-k loop; fully vectorized |
| External dependencies | None — pure base R |
| Output positivity | Tau > 0 guaranteed by tau_floor |
| Var_Term sign | >= 0 everywhere (Jensen's inequality) |
| ELBO monotonicity | elbo_proxy non-decreasing at convergence |
| Perfect reconstruction | Tau hits ceiling 1/tau_floor when RSS ≈ 0 |
| Heteroscedastic noise | Each τ_j estimated independently per column |

---

## Test Suite Summary

**File:** `tests/test_update_tau.R`
**Total:** 27 tests across 9 groups (all passing)

### T1: Mathematical Identities (5 tests)

Verify that each formula is implemented exactly as derived.

| Test | Assertion |
|---|---|
| T1.1 | `Var_Term == EL2 %*% t(EF2) - EL^2 %*% t(EF^2)` (to machine precision) |
| T1.2 | `all(Var_Term >= 0)` |
| T1.3 | `R2_bar == (Y - EL %*% t(EF))^2 + Var_Term` |
| T1.4 | `Tau == n / colSums(R2_bar)` (when above floor) |
| T1.5 | ELBO proxy formula matches manual computation |

### T2: Perfect Reconstruction (3 tests)

Edge case where Y is exactly reproduced by posterior means and posteriors collapse to point masses.

| Test | Assertion |
|---|---|
| T2.1 | Var_Term = 0 when EL2 = EL^2 and EF2 = EF^2 |
| T2.2 | R2_bar ≈ 0 when Y = EL %*% t(EF) and no posterior variance |
| T2.3 | Tau hits ceiling (1/tau_floor) for zero-residual case |

### T3: Known Noise Recovery (3 tests)

Simulates Y = L F' + N(0, 1/τ_j) with zero posterior variance; checks that estimated τ recovers truth.

| Test | Assertion |
|---|---|
| T3.1 | cor(Tau_est, Tau_true) ≥ 0.9 (n = 200) |
| T3.2 | Heteroscedastic columns separated correctly |
| T3.3 | Larger n gives better recovery (consistency check) |

### T4: Variance Correction (3 tests)

Demonstrates that the Var_Term matters quantitatively.

| Test | Assertion |
|---|---|
| T4.1 | Var_Term inflates R2_bar relative to naive (Y - pred)^2 |
| T4.2 | Ignoring correction overestimates mean(Tau) |
| T4.3 | Larger posterior variance → larger Var_Term |

### T5: ELBO Proxy (3 tests)

| Test | Assertion |
|---|---|
| T5.1 | elbo_proxy matches manual formula |
| T5.2 | elbo_proxy is finite (no NaN/Inf) |
| T5.3 | MLE Tau maximizes elbo_proxy (ELBO_MLE ≥ ELBO_perturbed) |

### T6: Dimension & Shape (2 tests)

| Test | Assertion |
|---|---|
| T6.1 | dim(Var_Term) == c(n, p) |
| T6.2 | length(Tau) == p |

### T7: Numerical Stability (4 tests)

| Test | Assertion |
|---|---|
| T7.1 | Y = 0 produces finite, positive Tau |
| T7.2 | 1×1×1 case (n=1, p=1, K=1) does not error |
| T7.3 | Large entries (1e6 scale) do not overflow |
| T7.4 | tau_floor prevents Tau = Inf |

### T8: Pipeline Interaction (2 tests)

| Test | Assertion |
|---|---|
| T8.1 | Tau is positive, finite, length p (valid downstream input) |
| T8.2 | Tau can be used in A_L computation without error |

### T9: V2.R Consistency (2 tests)

Compares output against the corresponding lines in `code/Supervised_Bayesian_MF_V2.R`.

| Test | Assertion |
|---|---|
| T9.1 | Var_Term and Tau match V2.R lines 374–376 to within 1e-12 |
| T9.2 | elbo_proxy matches V2.R line 379 to within 1e-12 |

---

## Demo Scenarios

**File:** `demos/demo_update_tau.R`

Each demo is a narrative scenario with printed output. There are no PASS/FAIL assertions — these are exploratory illustrations, not unit tests.

### Demo 1: Anatomy of Var_Term and R2_bar

**Purpose:** Decompose R2_bar into its two additive components and show that Var_Term is non-negligible.

**Setup:** Small n=20, p=10, K=3 example with realistic posterior variances.

**Key output:** Mean ratio Var_Term / R2_bar; illustrates that ignoring the correction would systematically bias τ upward.

### Demo 2: Known Noise Recovery

**Purpose:** Validate that update_tau recovers true per-feature precisions under controlled conditions.

**Setup:** Y = L F' + noise with τ_j ~ Uniform(0.5, 5); point-mass posteriors (EL2 = EL^2, EF2 = EF^2); n=200.

**Key output:** cor(estimated Tau, true Tau) ≈ 0.9+; scatter plot of estimated vs. true.

### Demo 3: Variance Correction Matters

**Purpose:** Quantify the bias introduced by ignoring EL2 and EF2.

**Setup:** Same Y as Demo 2, but run update_tau twice — once with true second moments, once substituting EL2 ← EL^2 and EF2 ← EF^2.

**Key output:** mean(Tau_naive) / mean(Tau_correct) ≈ 3–5×; shows systematic overestimation of precision (underestimation of noise) when the correction is omitted.

### Demo 4: Heteroscedastic Noise

**Purpose:** Confirm that the column-wise MLE correctly separates features with different noise levels.

**Setup:** p=150 features split into two groups — high-precision (τ=5, low noise) and low-precision (τ=0.5, high noise).

**Key output:** Estimated Tau clearly bimodal; group means close to 5 and 0.5 respectively.

### Demo 5: ELBO as Convergence Monitor

**Purpose:** Show that elbo_proxy is a valid convergence diagnostic that increases (or plateaus) across iterations.

**Setup:** 3 manual CAVI iterations — update L → update F → update τ — on a small example.

**Key output:** elbo_proxy at each iteration; should be non-decreasing. A decrease would indicate a bug in L or F update.

---

## Pipeline Context

### Data Flow

```
update_L_all  →  EL, EL2  ─┐
                             ├─→  update_tau  →  Tau (p-vector)
update_F_all  →  EF, EF2  ─┘
```

**Upstream:** `update_tau` consumes posterior moments from `update_L.R` and `update_F.R`. Specifically, it requires both the posterior means (EL, EF) and the posterior second moments (EL2, EF2).

**Downstream:** The output `Tau` feeds back into:
- `update_L_k` via the `A_L` precision accumulation: `A_L <- A_L + Tau[j] * EF2[j, k]`
- `update_F_k` via the `A_F` precision accumulation: `A_F <- A_F + Tau[j] * EL2[i, k]`

### Position in CAVI Loop

```
Repeat until convergence:
  1. update_beta_all(...)       # q(β) update — scalar EBNM per factor
  2. update_L_all(...)          # q(L) update — vector EBNM, dual-source
  3. update_F_all(...)          # q(F) update — vector EBNM, genomics only
  4. update_tau(...)            # q(τ) update — closed-form MLE  <-- THIS MODULE
  5. Check convergence (ΔL + Δβ < tol, ELBO non-decreasing)
```

### V2.R Reference Lines

The logic in this module corresponds to the following lines in `code/Supervised_Bayesian_MF_V2.R`:

| V2.R Lines | Content |
|---|---|
| 374–376 | Var_Term, R2_bar, Tau computation |
| 379 | ELBO proxy accumulation |

Tests T9.1 and T9.2 verify exact agreement (to 1e-12) between this module and those V2.R lines.

---

## Related Files

| File | Role |
|---|---|
| `code/update_tau.R` | Source implementation |
| `tests/test_update_tau.R` | 27-test suite (9 groups) |
| `demos/demo_update_tau.R` | 5 narrative demo scenarios |
| `derivations/qTau/qTau_update_derivation.tex/.pdf` | Full mathematical derivation (variance correction proof) |
| `code/Supervised_Bayesian_MF_V2.R` | Main algorithm (lines 374–376, 379) |
| `code/update_L.R` | Upstream: produces EL, EL2 |
| `code/update_F.R` | Upstream: produces EF, EF2 |
| `code/SupervisedMF_Context.md` | AI/code quick-reference for the full R implementation |
| `PROJECT_STATUS.md` | Project-level documentation and session log |

---

## Quick Reference Card

```
# Minimal usage
source("code/update_L.R")   # provides EL, EL2
source("code/update_F.R")   # provides EF, EF2
source("code/update_tau.R")

result     <- update_tau(Y, EL, EL2, EF, EF2)
Tau        <- result$Tau        # p-vector of noise precisions
elbo_proxy <- result$elbo_proxy # convergence monitor

# Internals
Var_Term <- compute_var_term(EL, EL2, EF, EF2)         # n × p, >= 0
R2_bar   <- compute_expected_residual_sq(Y, EL, EL2, EF, EF2)  # n × p

# MLE formula (closed form, no EBNM, no per-k loop)
Tau[j]   <- n / sum(R2_bar[, j])

# ELBO proxy
elbo_proxy <- sum(n/2 * log(Tau) - Tau/2 * colSums(R2_bar))
```

**Critical reminders**
- Var_Term uses **second moments** EL2, EF2 — not (EL)², (EF)². Confusing these is the most common implementation error and causes 3–5× overestimation of τ.
- tau_floor is applied to RSS (not to Tau directly): `pmax(RSS_j, n * tau_floor)`.
- ELBO_proxy should be **non-decreasing** across CAVI iterations. If it drops, check the L or F update first.
- This is the **only module** with no EBNM call and no `ebnm` package dependency.
