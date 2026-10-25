# Design: `fit_modular.R` — Modular CAVI Loop + Companion Doc

**Date:** 2026-03-25
**Author:** Andrew Walther
**Status:** Approved

---

## Overview

Implement V3 Algorithm 1 ("Supervised MF CAVI Update Procedure") as a reusable R function
`fit_supervised_mf_modular()` in `code/fit_modular.R`, using the four existing modular update
functions (`update_L_k`, `update_F_k`, `update_beta_k`, `update_tau`). Pair it with a
companion document `docs/fit_modular.qmd` in the same style as the existing `update_*.qmd` docs.

---

## Deliverable 1: `code/fit_modular.R`

### Public API

```r
fit_supervised_mf_modular(
  Y,             # n×p genomics matrix
  time,          # n-vector of survival/censoring times
  status,        # n-vector of event indicators (1=event, 0=censored)
  K       = 5,   # number of latent factors
  max_iter = 100, # maximum CAVI iterations
  tol     = 1e-5, # convergence threshold (max absolute change)
  verbose  = TRUE  # print iteration summaries
)
# Returns: list(EL, EL2, EF, EF2, EBeta, EBeta2, Tau, history, converged, n_iter)
```

### DATA_MODE Toggle

At the top of the script, a `DATA_MODE` flag controls whether to run on synthetic or real data:

```r
DATA_MODE <- "simulated"   # "simulated" or "real"
```

- **`"simulated"`**: generates synthetic Y (n×p, from true L, F with Gaussian noise) and
  survival data (time/status from Cox model with true β), then calls
  `fit_supervised_mf_modular()` and prints a summary. Parameters (n, p, K, seed, true β) set
  as named constants at the top.
- **`"real"`**: user sets `real_Y`, `real_time`, `real_status` and the function runs on them.
  (Mirrors the DATA_MODE convention from `Supervised_Bayesian_MF_V2.R`.)

Running `Rscript code/fit_modular.R` executes the simulated path by default.

### Algorithm (V3 Algorithm 1, factor-wise Gauss-Seidel)

**Initialization (b=0):**
1. Rank-K SVD of Y → L̄^(0), F̄^(0); set L̄²^(0) = L̄^(0)², F̄²^(0) = F̄^(0)²
2. Warm-start β̄^(0) via Cox PH on initial loadings; set β̄²^(0) = β̄^(0)²
3. τ̂^(0) from column-wise sample variance of Y

**Loop (b = 1, 2, … until convergence):**

- **Step 1** — Compute survival working quantities once:
  - η̂_i = Σ_k l̄_ik β̄_k
  - (u_i, W_ii) = Cox score and negative Hessian at η̂
  - z_i = η̂_i + u_i / W_ii

- **Step 2** — For k = 1, …, K (factor-wise Gauss-Seidel):
  - **(a) Update q(l_k):** call `update_L_k`; store updated l̄_ik, l̄²_ik
  - **(b) Recompute R̄^{-k}** using the just-updated l̄_ik, then call `update_F_k`; store f̄_jk, f̄²_jk
  - **(c) Update q(β_k):** compute z̄_i^{-k} via `compute_z_no_k(z, EL, EBeta, k)` using the
    live `EL` matrix — which already contains the updated l̄_ik from step (a). `compute_z_no_k`
    excludes factor k from the sum so the updated l̄_ik does not create a circular dependency.
    Call `update_beta_k`; store β̄_k, β̄²_k.

- **Step 3** — Update τ̂: call `update_tau`; store τ̂_j and elbo_proxy

- **Step 4** — Check convergence (max, not mean); only checked after iter > 5 to avoid
  early termination on initialization noise:
  - Δ_L = max_{i,k} |l̄_ik^new − l̄_ik^old|
  - Δ_β = max_k |β̄_k^new − β̄_k^old|
  - Stop if iter > 5 AND max(Δ_L, Δ_β) < ε

**Returns:** posterior means L̄, F̄, β̄; second moments L̄², F̄², β̄²; noise precision τ̂;
history (rmse, elbo_proxy per iteration); converged flag; n_iter.

Return fields use the `E`-prefix naming convention (`EL`, `EF`, `EBeta`, etc.) matching the
internal storage names used throughout the modular update functions. This intentionally differs
from `fit_supervised_mf()` in V2.R which uses short names (`L`, `F`, `Beta`). The E-prefix
convention is more consistent with the module API and avoids name collision with base R.

### Dependencies

Sources (in order): `update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R`
Packages: `ebnm`, `survival`

`compute_R_k` is defined in `update_L.R` and used here for both L_k and F_k steps. The Cox
Taylor helper (`calc_cox_taylor`) is copied verbatim from `code/Supervised_Bayesian_MF_V2.R`
lines 70–93 with an attribution comment — following the pattern in `run_modular_simulation.R`
lines 47–75. `Supervised_Bayesian_MF_V2.R` cannot be safely `source()`d as a dependency because
doing so executes its top-level simulation and plotting code.

### Key Distinction from `run_modular_simulation.R`

| | `run_modular_simulation.R` | `fit_modular.R` |
|---|---|---|
| **Update order** | Block (_all variants): all-L → all-F → all-β → τ | Factor-wise (_k variants): for k: L_k → F_k → β_k |
| **Convergence** | Mean absolute change | Max absolute change (V3 Algorithm 1) |
| **Purpose** | Simulation benchmark with plots | Reusable inference function |
| **Data** | Hardcoded simulation DGP | DATA_MODE toggle (simulated or real) |
| **Output** | Figures + tables written to `results/` | Return value (list of posteriors) |

---

## Deliverable 2: `docs/fit_modular.qmd`

### Style

Matches `docs/update_beta.qmd` / `update_L.qmd` / `update_F.qmd` / `update_tau.qmd`:

```yaml
---
title: "Companion Document: fit_modular.R"
subtitle: "Full CAVI Loop — Supervised Matrix Factorization"
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
    ...
---
```

### Sections

1. **Overview** — what the script does; relationship to V3 Algorithm 1 and the four update modules
2. **Mathematical Background** — model equations, objective function, Taylor approximation (rendered LaTeX); initialization; convergence criterion (max vs. mean distinction)
3. **Algorithm Walkthrough** — step-by-step mapping of V3 Algorithm 1 to R code; factor-wise Gauss-Seidel explained
4. **Function Reference** — `fit_supervised_mf_modular()` signature, arguments, return value
5. **DATA_MODE Toggle** — how to switch between simulated and real data; parameter constants
6. **Distinction from `run_modular_simulation.R`** — explicit comparison table + prose explanation of why factor-wise > blockwise for faithful CAVI
7. **Related Files** — links to the four update module docs, V2.R, derivation PDFs

### Math conventions

All mathematical expressions use Quarto LaTeX: `$...$` inline, `$$...$$` display,
`$\begin{aligned}...\end{aligned}$` for multi-line. Consistent with the other companion docs.

---

## Out of Scope

- No new tests (existing 105 tests cover the individual modules; integration correctness
  verified by the simulated-data run recovering true β signs and RMSE ≈ 1)
- No new demos directory entry (DATA_MODE toggle serves the demo purpose)
- No changes to `run_modular_simulation.R`

---

## File Summary

| File | Action |
|------|--------|
| `code/fit_modular.R` | **New** — function + DATA_MODE runner |
| `docs/fit_modular.qmd` | **New** — companion doc |
| `docs/fit_modular.pdf` | **New** — rendered output |
| `docs/fit_modular.html` | **New** — rendered output |
| `docs/Makefile` | **Update** — add fit_modular.qmd to SRCS |
| `CLAUDE.md` | **Update** — add fit_modular.R to quick-reference table |
| Memory / PROJECT_STATUS | **Update** — after implementation |
