# Plan: Factor Classification & K-Selection Validation

## Context

The recommended D4 model (YFB, K=7 CV-selected) shows K_eff=2 survival-active factors.
Two methodological gaps need to be addressed before the manuscript:

1. **K was selected by cross-validated C-index** — which uses survival outcomes to choose
   model structure, conflates structure selection with fitting, and was advised against.
   The alternative (Method 2) is to start at a large K and let the ARD prior prune.

2. **Factor classification is binary** — current code only distinguishes "active" (|β|>threshold)
   vs. "shrunk". There is no accounting of factors that are *genomics-active but survival-silent*
   (real gene expression programs that don't associate with survival) vs. *dead* (fully pruned).
   The manuscript needs: "X total programs, Y are survival-associated (Z adverse, W protective),
   X−Y are genomics-active without prognostic value."

3. **Simulations always set K_fit = K_true** — ARD's ability to recover the correct K when
   over-specified has never been validated. The one archived data point (K_init=10, K_true=5)
   gave K_eff=4, not 5.

4. **Signal-ratio sweep used K_fit = K_true = 6** — the comparison against the two-step
   EBMF→Cox baseline should be re-run with K_init ≫ K_true, using ARD pruning, to confirm
   the conclusion holds under the corrected K-selection approach.

---

## What Exists (Reuse These)

| Utility | Location | Purpose |
|---|---|---|
| `compute_pve(res, Y)` | `code/select_K.R:28` | Per-factor PVE (genomics activity), K-vector |
| `auto_prune_K(res, Y, ...)` | `code/select_K.R:69` | Active flag via OR(|β|>thresh, PVE>thresh) |
| `generate_multicohort_data(...)` | `results/multi_cohort_sim/generate_multicohort_data.R` | DGP with shared (survival-active) and specific (genomics-only, β=0) factors; returns `$factor_labels` |
| `history$factor_pve` | model object | Per-iteration PVE matrix (n_iter × K) |
| `desurv_comparison_fits.rds` | `results/benchmark_sim/outputs/desurv_comparison/` | Current D4 fit at K=7 — starting reference |
| `globals.yml` thresholds | `config/globals.yml` | `beta_threshold: 0.001`, `pve_threshold: 0.01` |

**No existing utility classifies factors into the three-way survival-active / genomics-only / dead
table. This is the central missing piece.**

---

## Implementation Plan

### Step 1 — Add `classify_factors()` to `code/select_K.R`

New function added after `auto_prune_K()`. Reuses `compute_pve()` internally.

```r
classify_factors <- function(res, Y,
                              beta_thresh = 0.001,   # from globals.yml
                              pve_thresh  = 0.01) {
  pve     <- compute_pve(res, Y)          # reuse existing function
  ab_beta <- abs(res$EBeta)
  surv_active <- ab_beta > beta_thresh
  geno_active <- pve     > pve_thresh
  data.frame(
    factor      = seq_len(ncol(res$EL)),
    EBeta       = res$EBeta,
    abs_EBeta   = ab_beta,
    PVE         = pve,
    surv_active = surv_active,
    geno_active = geno_active,
    category    = ifelse(surv_active, "survival_active",
                  ifelse(geno_active, "genomics_only",
                                      "dead"))
  )
}
```

Also fix the hardcoded `beta_thresh = 0.05` in `run_phase1_diagnostics.R:52`
to read from `globals.yml` (consistency with the rest of the codebase).

**Test:** Add 3–4 unit tests to `tests/test_select_K_cv.R` covering each category
(survival_active, genomics_only, dead) and the edge case where both thresholds
are borderline.

---

### Step 2 — Analysis A: K_init Stability Sweep on Real PDAC

**New script:** `results/benchmark_sim/run_k_init_sweep.R`

**What it does:**
- Loads TCGA+CPTAC with D4 preprocessing (DeSurv-aligned gene selection, per-platform z-std)
  by sourcing existing preprocessing from `run_desurv_comparison.R`
- Fits YFB (`fit_cox_on_yf`) at K_init ∈ {7, 10, 15, 20} — includes current CV K=7 as baseline
- For each fit, calls `classify_factors()` and records:
  - K_total, K_survival_active, K_genomics_only, K_dead
  - Per-factor EBeta and PVE
  - External C-index across 5 held-out cohorts (reuses existing `predict_cox_on_yf` + `concordance`)
- Saves: `results/benchmark_sim/outputs/k_init_sweep/k_init_sweep_results.csv`

**Key question answered:** Is K_eff_survival always 2, and K_eff_genomics always N,
regardless of starting K? If yes → ARD pruning is stable and Method 2 is valid.

---

### Step 3 — Analysis B: ARD K-Recovery Simulation

**New script:** `results/multi_cohort_sim/run_k_recovery_sim.R`

**DGP** (via `generate_multicohort_data()` with `specific_prognostic=FALSE`):

| Condition | K_shared (survival-active) | K_specific per cohort (genomics-only) | K_true_total |
|---|---|---|---|
| A | 1 | 2 | 5 |
| B | 2 | 2 | 6 |
| C | 2 | 3 | 8 |

For each condition, fit at K_init ∈ {K_true+5, K_true+10, K_true+15}, 5 seeds.

**Metrics per fit:**
- K_eff_survival (from `classify_factors()`) vs. K_true_shared (ground truth)
- K_eff_genomics vs. K_true_specific (ground truth)
- Factor recovery: max-cor of estimated F vs. true F for shared programs separately from specific
- β RMSE for survival-active factors

**Key question answered:** Does the model correctly identify *both* the number of survival-active
factors AND the number of genomics-only factors when K_init ≫ K_true?

**Output:** `results/multi_cohort_sim/outputs/k_recovery_sim_results.csv`

---

### Step 4 — Analysis C: Re-run Signal-Ratio Sweep with K_init ≫ K_true

**Modify:** `results/multi_cohort_sim/run_signal_ratio_sweep.R`

**Change:** Add a `K_INIT` parameter (default 20, separate from K_FIT=6=K_true) passed through to
each model arm's fit call. The existing `K_FIT` variable is retained for data generation only;
each arm's model fit uses `K_INIT`.

This re-runs the sweep over `a_specific` ∈ {12, 18, 24, 36, 48} with the corrected
K-selection approach, confirming the original finding (YFB outperforms EBMF→Cox when
survival signal is strong) holds under ARD pruning from over-specified K.

**Output:** `results/multi_cohort_sim/outputs/signal_ratio_sweep_results_kinit20.csv`
(new file, preserving the original)

---

## Files Created / Modified

| Action | File |
|---|---|
| **Modify** | `code/select_K.R` — add `classify_factors()` after `auto_prune_K()` |
| **Modify** | `tests/test_select_K_cv.R` — add unit tests for `classify_factors()` |
| **Fix** | `results/benchmark_sim/run_phase1_diagnostics.R` — replace hardcoded `beta_thresh=0.05` with `cfg$k_selection$beta_threshold` |
| **New** | `results/benchmark_sim/run_k_init_sweep.R` — Analysis A |
| **New** | `results/multi_cohort_sim/run_k_recovery_sim.R` — Analysis B |
| **Modify** | `results/multi_cohort_sim/run_signal_ratio_sweep.R` — Analysis C (add K_INIT param) |
| **Update** | `DECISIONS.md` — record K-selection method change and findings |
| **Update** | `ROADMAP.md` — mark Analysis A/B/C as in-progress |

---

## Execution Order

1. `code/select_K.R` + tests (prerequisite for all scripts)
2. Analysis A (`run_k_init_sweep.R`) — quick, runs on existing real data; informs B and C
3. Analysis B (`run_k_recovery_sim.R`) — validates ARD in simulation
4. Analysis C (`run_signal_ratio_sweep.R` update) — re-run only after A+B confirm the approach

---

## Verification

- `Rscript tests/run_tests.R` — all 374 tests still passing + new classify_factors tests
- Analysis A output: K_eff_survival == 2 at all K_init (expected); C-index stable
- Analysis B output: K_eff_survival == K_true_shared, K_eff_genomics ≈ K_true_specific
  across seeds and K_init values
- Analysis C output: signal-ratio curve shape unchanged from original sweep
