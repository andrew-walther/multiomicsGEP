# SupervisedMF — Code Context File

**For AI sessions and future collaborators picking up this codebase.**

This file lives in `code/` and is the fast-access reference for the R
implementation. For broader project context (repo map, derivation review,
session log), see `PROJECT_STATUS.md` at the repo root.

---

## The Active Script

**`Supervised_Bayesian_MF_V2.R`** is the current implementation.
`code/legacy/Supervised_Bayesian_MF.R` is the original V1 — archived for reference only.

---

## Model in One Line

```
Y (n×p) = L (n×K) × F' (K×p) + E        [genomics]
h(t_i) = h₀(t_i) exp(L_i · β)            [survival / Cox PH]
```

Inference: **CAVI** (Coordinate Ascent Variational Inference).
Each update reduces to an **EBNM** (Empirical Bayes Normal Means) sub-problem.

---

## Variable Dictionary

| R Variable | Math Symbol | Dimension | Meaning |
|------------|-------------|-----------|---------|
| `EL` | $\bar{L}$, $\bar{l}_{ik}$ | n × K | Posterior **means** of patient loadings |
| `EL2` | $\overline{l^2_{ik}}$ | n × K | Posterior **second moments** = Var + mean² |
| `EF` | $\bar{F}$, $\bar{f}_{jk}$ | p × K | Posterior means of factor weights |
| `EF2` | $\overline{f^2_{jk}}$ | p × K | Posterior second moments |
| `EBeta` | $\bar{\beta}_k$ | K | Posterior means of survival coefficients |
| `EBeta2` | $\overline{\beta^2_k}$ | K | Posterior second moments |
| `Tau` | $\tau_j$ | p | Feature-specific noise **precision** (1/variance) |
| `w` | $W_{ii}$ | n | Cox neg-diagonal Hessian (working weights) |
| `z` | $z_i$ | n | Cox working response = η̂ + u/W |
| `eta` | $\hat{\eta}_i$ | n | Linear predictor = `EL %*% EBeta` |
| `R_k` | $R^{-k}_{ij}$ | n × p | Genomics partial residual (factor k removed) |
| `z_no_k` | $z^{-k}_i$ | n | Survival partial working response (factor k removed) |
| `A_L` | $A_{ik}$ | n | EBNM precision for loadings update |
| `B_L` | $B_{ik}$ | n | EBNM signal numerator for loadings update |
| `A_F` | $A_{jk}$ | p | EBNM precision for factors update |
| `B_F` | $B_{jk}$ | p | EBNM signal numerator for factors update |
| `A_Beta` | $A_k$ | scalar | EBNM precision for β update |
| `B_Beta` | $B_k$ | scalar | EBNM signal numerator for β update |
| `Var_Term` | $\sum_k[\overline{l^2}\,\overline{f^2} - \bar{l}^2\bar{f}^2]$ | n × p | Variance inflation for τ update |
| `R2_bar` | $\bar{R}^2_{ij}$ | n × p | Expected squared residual |

> **Critical:** `EL2 ≠ EL^2` whenever there is posterior uncertainty.
> `EL2[i,k] = res_L$posterior$sd^2 + res_L$posterior$mean^2`

---

## Code Walkthrough (V2.R)

### Part 1 — Libraries & Helpers (lines 1–145)

| Function | Purpose |
|----------|---------|
| `calc_cox_taylor(eta, time, status)` | Cox 2nd-order Taylor → returns `list(u, w)` (score & neg-Hessian). Sorts by time internally. Floor on `w` at `1e-6`. |
| `get_cindex_comparison(EL, data)` | Harrell C-index: supervised L vs top-5 PCA |
| `get_top_features(EF, n_top)` | Ranks features by absolute weight per factor |
| `get_factor_summary_table(res, data)` | Per-factor: β, log-rank p, NonZero_Pct (% features with nonzero weight), PVE % |

### Part 2 — `fit_supervised_mf()` (lines 171–424)

**Signature:**
```r
fit_supervised_mf(Y, time, status,
                  K = 5, max_iter = 100, tol = 1e-5,
                  orthogonalize = FALSE, refresh_taylor = FALSE,
                  verbose = TRUE)
# Returns: list(L, F, Beta, Beta2, Tau, history)
```

**Initialization (lines 183–226):**
- SVD of Y → `EL = U√D`, `EF = V√D` (rank-K approximation)
- `EL2 = EL^2`, `EF2 = EF^2` (zero posterior variance at start)
- Warm-start `EBeta` via `coxph()` on initial loadings
- `Tau` = 1 / column variance of Y

**Outer CAVI loop (lines 237–421):**

```
for iter in 1:max_iter:
  1. Cox Taylor expansion → (z, w)           [lines 254–257]
  2. for k in 1:K:
       compute R_k, z_no_k                   [lines 290–295]
       (a) update L[:,k]  via EBNM(A_L, B_L) [lines 310–321]
       (b) update F[:,k]  via EBNM(A_F, B_F) [lines 334–340]
       (c) update Beta[k] via EBNM(A_β, B_β) [lines 356–361]
  3. update Tau (column-specific MLE)         [lines 374–376]
  4. track ELBO proxy                         [line 379]
  5. [optional] orthogonalise EL             [lines 389–395]
  6. check dual convergence (ΔL & Δβ < tol) [lines 400–420]
```

**Returns:**
```r
list(
  L      = EL,      # n×K posterior mean loadings
  F      = EF,      # p×K posterior mean factors
  Beta   = EBeta,   # K   posterior mean survival coefficients
  Beta2  = EBeta2,  # K   posterior second moments (for uncertainty)
  Tau    = Tau,     # p   noise precision per feature
  history = list(rmse, elbo_proxy, converged, n_iter)
)
```

### Part 3 — Visualisation (lines 430–509)

| Function | Output |
|----------|--------|
| `plot_gep_heatmap(res, n_features=50)` | Red-blue heatmap of top features × factors |
| `visualize_dashboard(res, data)` | Sequential plots: RMSE trace, ELBO proxy, heatmap, KM curves, signal recovery |

### Part 4 — Simulation & Execution (lines 515–583)

`sim_data_fn(n=250, p=1000, k=5)` generates:
- `L ~ N(0,1)`, `F` 5%-sparse with `N(0,25)` active weights
- `Y = LF' + N(0,1)`, Weibull survival with `β_true = (1.5, -1.2, 0.8, -0.5, 0)`

---

## Math ↔ Code Quick-Reference

### L update (per factor k)

| Math | R Code |
|------|--------|
| $A_{ik} = \sum_j \tau_j \overline{f^2_{jk}} + W_{ii}\overline{\beta^2_k}$ | `sum(Tau * EF2[,k]) + w * EBeta2[k]` |
| $B^{\text{gen}}_{ik} = \sum_j \tau_j R^{-k}_{ij} \bar{f}_{jk}$ | `R_k %*% (Tau * EF[,k])` |
| $B^{\text{surv}}_{ik} = W_{ii} z^{-k}_i \bar{\beta}_k$ | `w * z_no_k * EBeta[k]` |
| EBNM: $x_i = B_{ik}/A_{ik}$, $s_i = 1/\sqrt{A_{ik}}$ | `ebnm(x = B_L/A_L, s = 1/sqrt(A_L), ...)` |

### F update (per factor k)

| Math | R Code |
|------|--------|
| $A_{jk} = \tau_j \sum_i \overline{l^2_{ik}}$ | `Tau * sum(EL2[,k])` |
| $B_{jk} = \tau_j \sum_i R^{-k}_{ij} \bar{l}_{ik}$ | `Tau * t(R_k) %*% EL[,k]` |
| EBNM: $x_j = B_{jk}/A_{jk}$, $s_j = 1/\sqrt{A_{jk}}$ | `ebnm(x = B_F/A_F, s = 1/sqrt(A_F), ...)` |

### β update (per factor k)

| Math | R Code |
|------|--------|
| $A_k = \sum_i W_{ii} \overline{l^2_{ik}}$ | `sum(w * EL2[,k])` |
| $B_k = \sum_i W_{ii} z^{-k}_i \bar{l}_{ik}$ | `sum(w * z_no_k * EL[,k])` |
| EBNM: $x_k = B_k/A_k$, $s_k = 1/\sqrt{A_k}$ | `ebnm(x = B_Beta/A_Beta, s = 1/sqrt(A_Beta), ...)` |
| `z_no_k` **reused** from L update — no recomputation needed | same `z_no_k` variable |

### τ update

| Math | R Code |
|------|--------|
| $\bar{R}^2_{ij} = (Y - \bar{L}\bar{F}')^2 + \sum_k[\overline{l^2}\,\overline{f^2} - \bar{l}^2\bar{f}^2]$ | `(Y - EL%*%t(EF))^2 + (EL2%*%t(EF2) - EL^2%*%t(EF^2))` |
| $\hat{\tau}_j = n / \sum_i \bar{R}^2_{ij}$ | `n / colSums(R2_bar)` |

---

## Modular Update Functions (Sessions 3–5)

All four CAVI parameters now have standalone modular update files.
Each module follows the same pattern: R code + test file + demo + LaTeX derivation.

```r
source("code/update_beta.R")    # β update (Session 3)
source("code/update_L.R")       # q(L) update (Session 5)
source("code/update_F.R")       # q(F) update (Session 5) — requires update_L.R for compute_R_k
source("code/update_tau.R")     # q(τ) update (Session 5) — base R only, no EBNM
```

### update_L.R — Vector EBNM, dual-source

```r
compute_R_k(Y, EL, EF, k)
# Returns: n×p partial residual R^{-k} = Y - EL%*%t(EF) + outer(EL[,k], EF[,k])
# Shared between update_L.R and update_F.R

update_L_k(Tau, EF_k, EF2_k, w, EBeta_k, EBeta2_k, R_k, z_no_k,
           prior_family = "point_normal", A_floor = 1e-10)
# Key: A_L[i] and B_L[i] are n-VECTORS (patient-specific via Cox weights W_ii)
#   A_gen  = sum(Tau * EF2_k)              [SCALAR — same for all i]
#   A_surv = w * EBeta2_k                  [n-vector — varies via W_ii]
#   A_L    = pmax(A_gen + A_surv, floor)   [n-vector]
#   B_gen  = R_k %*% (Tau * EF_k)         [n-vector]
#   B_surv = w * z_no_k * EBeta_k         [n-vector]
#   ebnm(x = B_L/A_L, s = 1/sqrt(A_L))  ← n-observation EBNM call
# Returns: list(mean, second, sd, A, B, B_gen, B_surv, x, s, ebnm_result)

update_L_all(Y, EL, EL2, EF, EF2, Tau, w, z, EBeta, EBeta2, ...)
# Returns: list(EL, EL2, details)
```

### update_F.R — Pure genomics, τ cancellation

```r
update_F_k(Tau, EL_k, EL2_k, R_k,
           prior_family = "point_normal", A_floor = 1e-10)
# KEY PROPERTY: τ_j CANCELS in x_j = B_F[j]/A_F[j], does NOT cancel in s_j
#   sum_EL2_k = sum(EL2_k)                [SCALAR]
#   A_F = pmax(Tau * sum_EL2_k, floor)    [p-vector]
#   B_F = Tau * t(R_k) %*% EL_k           [p-vector]
#   x_j = B_F[j]/A_F[j] = ... / sum_EL2_k [tau-FREE]
#   s_j = 1/sqrt(Tau * sum_EL2_k)         [tau-dependent]
#   ebnm(x = x_F, s = s_F)              ← p-observation EBNM call
# Returns: list(mean, second, sd, A, B, x, s, sum_EL2_k, ebnm_result)

update_F_all(Y, EL, EL2, EF, EF2, Tau, ...)
# Returns: list(EF, EF2, details)
```

### update_tau.R — Closed-form MLE, no EBNM

```r
compute_var_term(EL, EL2, EF, EF2)
# Returns: n×p matrix V = EL2%*%t(EF2) - EL^2%*%t(EF^2)  [all entries >= 0]
# This is the posterior variance correction preventing tau overestimation

update_tau(Y, EL, EL2, EF, EF2, tau_floor = 1e-8)
# NOT an EBNM problem — closed-form MLE with NO per-k loop
#   Var_Term = compute_var_term(...)       [n×p]
#   R2_bar   = (Y - EL%*%t(EF))^2 + Var_Term  [n×p]
#   Tau      = n / pmax(colSums(R2_bar), n*tau_floor)  [p-vector]
#   elbo_proxy = sum(n/2*log(Tau) - Tau/2*colSums(R2_bar))  [scalar]
# Returns: list(Tau, R2_bar, Var_Term, elbo_proxy)
```

### β update (Session 3)

`code/update_beta.R` extracts the β update from V2.R into testable, decoupled functions.
Source it independently: `source("code/update_beta.R")`.

### Function signatures

```r
# Partial working response (excludes factor k's contribution)
compute_z_no_k(z, EL, EBeta, k)
# Returns: length-n vector z^{-k}_i = z_i - Σ_{k'≠k} l̄_{ik'} β̄_{k'}

# Single-factor β update — core EBNM call
update_beta_k(w, z_no_k, EL_k, EL2_k,
              prior_family = "point_normal", A_floor = 1e-10)
# Returns: list(mean, second, sd, A, B, x, s, ebnm_result)
#   A      = sum(w * EL2_k)       — precision (error-in-variables)
#   B      = sum(w * z_no_k * EL_k)  — signal
#   x, s   = B/A, 1/sqrt(A)       — EBNM inputs
#   mean   = posterior$mean        — replaces EBeta[k]
#   second = sd^2 + mean^2         — replaces EBeta2[k]

# Full K-factor update with Gauss-Seidel ordering
update_beta_all(w, z, EL, EL2, EBeta,
                prior_family = "point_normal", A_floor = 1e-10)
# Returns: list(EBeta, EBeta2, details)
#   EBeta  = updated posterior means (K-vector)
#   EBeta2 = updated posterior second moments (K-vector)
#   details = list of K update_beta_k results (for diagnostics)
```

### Key properties
- `update_beta_k` is a **drop-in replacement** for V2.R lines 356–361 (verified by T9.1)
- `update_beta_all` propagates `EBeta_curr[k]` immediately (Gauss-Seidel), matching V2.R A1
- `A_floor = 1e-10` matches V2.R A3 precision floor
- `z_no_k` computed inside `update_beta_all` but also returned by `compute_z_no_k`
  for reuse in the L update (when extracting `update_L.R`)

### Test suite

```r
# Run all tests (24 tests, 9 groups):
Rscript tests/run_tests.R

# Test groups:
# T1: Mathematical identities (A_k, B_k, x_k, s_k, second moment)
# T2: WLS limit (EL2 = EL^2 → no posterior uncertainty → standard WLS)
# T3: K=1 signal recovery
# T4: Multi-factor K=5 (signs correct, null factor zeroed)
# T5: Null factor shrinkage (z_no_k pure noise → β≈0)
# T6: Error-in-variables (higher EL2 → larger A_k → smaller |x_k| → more shrinkage)
# T7: Numerical stability (w=0, extreme weights, n=1, floor test)
# T8: Gauss-Seidel ordering (updated β_1 propagates to z_no_k for β_2)
# T9: V2.R consistency (matches V2.R lines 356-361 exactly to 1e-12)
```

---

## Known Issues & Design Decisions

### Why `orthogonalize = FALSE` (default)
Orthogonalisation applies SVD rotation to EL every 10 iterations for factor
identifiability. But it resets `EL2 = EL^2` and `EF2 = EF^2`, **dropping all
posterior variance** until the next EBNM call. This causes:
- Transient underestimate of posterior uncertainty
- Temporarily reduced shrinkage on β (A_Beta biased low)
- Temporarily inflated τ (variance inflation term omitted)

Point-Normal EBNM already promotes sparsity and factor distinctness without
rotation, so orthogonalisation is off by default.

### Why `refresh_taylor = FALSE` (default)
`(z, w)` are computed once per outer iteration (standard IRLS-within-VI).
Within the inner k-loop, as EL and EBeta update, the true linear predictor
shifts but z and w stay stale. The outer loop restores exactness. Setting
`refresh_taylor = TRUE` recomputes at every k (more accurate, slower).

### Gauss-Seidel vs Jacobi [A1]
V1 had a subtle bug: `z_no_k` was computed from the linear predictor
**at the start of the outer iteration**, so it didn't reflect updates to
earlier factors within the same iteration (Jacobi-style). V2 recomputes
`eta_no_k` from the current `EL` and `EBeta` inside the k-loop
(Gauss-Seidel), which is the correct CAVI update.

### Single-observation EBNM for β
Each β_k update passes exactly one pseudo-observation (x_k, s_k) to EBNM.
This works because:
- Point-normal family has only 2 parameters (π₀, σ²)
- Shrinkage toward zero occurs naturally when |x_k|/s_k is small
- The prior is re-estimated each iteration as loadings stabilise

### `z_no_k` reuse is mathematically valid
$z^{-k}_i = z_i - \sum_{k'\neq k} \bar{l}_{ik'}\bar{\beta}_{k'}$ depends
only on k'≠k quantities. Since the L and β updates for factor k don't change
any k'≠k values, the same `z_no_k` computed before the L update is valid for
the β update of the same factor k.

### Numerical floors [A3]
All EBNM precision inputs are floored at `1e-10` via `pmax(...)`. The τ
denominator is floored at `n * 1e-8`. This prevents:
- Division by zero in `s = 1/sqrt(A)`
- Infinite pseudo-observations when A ≈ 0
- Unbounded prior variance in EBNM

---

## Backlog

### High Priority
- [ ] **Real data run:** Set `DATA_MODE <- "real"` in V2.R and supply:
  - `real_genomics_mat` — n×p numeric matrix (rows=patients, cols=features)
  - `real_clinical_df` — data frame with `time` and `status` columns
  - Adjust `real_time_col` / `real_status_col` if column names differ
  - Consider pre-processing: log-normalise RNA-seq; scale features
- [ ] **Select K:** Implement cross-validated C-index or held-out log-likelihood
  to choose K objectively. Currently hardcoded to 5.

### Medium Priority
- [ ] **Full ELBO:** Track the complete ELBO including the Cox survival
  Taylor approximation term and all KL divergences. Currently only
  the genomics log-likelihood is tracked as `elbo_proxy` (line 379).
- [ ] **Alternative priors:** Try `prior_family = "point_laplace"` or
  `"normal_scale_mixture"` in the EBNM calls for different sparsity.
- [ ] **Scalability:** For large p (>10,000 features), the main bottleneck
  is the n×p `Var_Term` matrix and `R_k` construction — consider
  batching over features or using sparse representations.

### Lower Priority
- [ ] **Manuscript:** `paper/multiomicsGEP_manuscript.qmd` — update to
  reference V2 results (RMSE=0.998, C-index comparison) and corrected
  derivations.
- [ ] **Bibliography:** Add `\bibliography{refs}` to REVISED.tex to resolve
  the `\citep{wang2022}` undefined-reference warning in the compiled PDF.
- [x] **CI/testing (β):** 24 tests for β update — all passing.
- [x] **Modular L/F/τ updates:** `update_L.R`, `update_F.R`, `update_tau.R` with
  28+26+27 tests (105 total passing) and derivations in `derivations/qL/qF/qTau/`.
- [ ] **Factor density / sparsity tuning:** Simulation shows 0% sparsity (all
  1,000 features have nonzero weight — `NonZero_Pct = 100`). The point-normal
  prior did not induce exact zeros. For cleaner GEP interpretation on real data,
  experiment with tighter EBNM priors or post-hoc thresholding of small-weight
  features.

---

## How to Run

```r
# Simulated benchmark (default)
source("code/Supervised_Bayesian_MF_V2.R")

# Real data
DATA_MODE         <- "real"
real_genomics_mat <- your_matrix        # n × p numeric matrix
real_clinical_df  <- your_dataframe     # must have columns: time, status
real_time_col     <- "time"
real_status_col   <- "status"
source("code/Supervised_Bayesian_MF_V2.R")
```

**Expected runtime** (simulated, n=250, p=1000, K=5): ~30–90 seconds
depending on convergence.

**Key output to check:**
1. ELBO proxy non-decreasing ✓
2. RMSE converging near 1.0 ✓
3. β signs match expectation ✓
4. Supervised C-index > PCA C-index ✓

---

## Derivation Documents

| File | Purpose |
|------|---------|
| `derivations/qB/qBeta_update_derivation.tex/.pdf` | **(Session 3):** Self-contained q(β) derivation (11 pages). |
| `derivations/qL/qL_update_derivation.tex/.pdf` | **(Session 5):** q(L) derivation — dual-source (genomics+survival), vector EBNM. |
| `derivations/qF/qF_update_derivation.tex/.pdf` | **(Session 5):** q(F) derivation — τ cancellation proof, pure genomics. |
| `derivations/qTau/qTau_update_derivation.tex/.pdf` | **(Session 5):** q(τ) derivation — variance correction, closed-form MLE vs EBNM. |
| `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_REVISED.tex/.pdf` | Corrected derivations (R1–R8 fixed) **+ full step-by-step algebra** (21 pages). Three tcolorbox types: gray `derivbox` = algebra steps, red `correctionbox` = R1–R8 fixes, blue `ebnmbox` = EBNM problem statements. |
| `derivations/MF_UpdateDerivations/MF_V2_Companion.tex/.pdf` | Section-by-section math → V2.R code mapping (17 pages) |
| `derivations/MF_UpdateDerivations/MF_Derivations_UpdateAlgo_2_12_26.pdf` | **Original Feb 12 document — contains R1–R8 errors. Reference only.** |

## Simulation Results (results/)

Benchmark: n=250 patients, p=1000 features, K=5 factors.
True β = (1.5, −1.2, 0.8, −0.5, 0), 5% feature sparsity, τ=1 noise.

| Metric | Value |
|--------|-------|
| Final RMSE | 0.9978 (converged near true noise σ=1) |
| Convergence | 45 iterations |
| C-index (PCA top-5) | 0.827 |
| C-index (Supervised MF) | 0.828 |
| PH test p-value | 0.258 (no proportional-hazards violation) |
| β recovery | Signs correct on GEPs 1–4; GEP5 (true β=0) correctly zeroed |
| Factor density | GEP1–5: NonZero_Pct = 100 (all 1000 features have nonzero weight; 0% sparse) |

Full report: `results/simulation_report.pdf` (23 pages) or `.md` source.
Figures: `results/figures/fig1_rmse_trace.png` through `fig8_tau_distribution.png`.
To reproduce: `source("results/run_simulation.R")`.
