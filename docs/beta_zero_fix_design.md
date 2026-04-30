# Design Document: Fixing β=0 in SSBMF on Multi-Platform Data

**Created:** 2026-04-29
**Status:** Design phase — both clusters documented, Cluster A approved for implementation
**Companion docs:** `docs/update_L_fix.md` (Cluster A debugging guide), `derivations/qL/qL_update_derivation.tex` (reference for Cluster B math)
**Related plan file:** `~/.claude/plans/here-s-the-combined-comprehensive-gleaming-lightning.md`

---

## §1. Context & Problem

### 1.1 The observed failure

Training the SSBMF model (`code/fit_modular.R`) on the merged TCGA_PAAD (RNA-seq) + CPTAC (proteomics) cohort produces all `EBeta = 0`. No factor is selected as prognostic, so the model provides no survival prediction beyond a constant hazard.

### 1.2 Diagnostic chain (4 steps)

| # | Experiment | Finding |
|---|---|---|
| 1 | Lambda sweep, λ ∈ {1,5,10,20} × 3 priors | β = 0 in all 12 combinations; λ ≥ 5 collapses entire EL matrix → not a tuning issue |
| 2 | EBMF unsupervised diagnostic | 5/20 factors Cox-significant (C-index up to 0.629) → **survival signal exists in the data** |
| 3 | Warm-start Exp 1: β-only CAVI with EL fixed at EBMF loadings | β non-zero at iteration 1, 6/20 factors active → **β update is functional** |
| 4 | Warm-start Exp 2: full CAVI from EBMF init | β collapses to zero in 23 iterations → **L update washes out survival signal** |

**Conclusion:** Neither the data nor the β update is the failure mode. The L update structurally erodes survival signal whenever the joint CAVI is run.

### 1.3 Root cause: chicken-and-egg + scale imbalance

In `update_L_k()` (`code/update_L.R` lines 140–144), the EBNM precision for subject *i* on factor *k* is:

```
A_L[i] = (1 - α) * A_gen[i]  +  α * A_surv[i]

where:
  A_gen     = sum_j  τ_j * E[f²_jk]        (scalar; sums over p = 2000 features)
  A_surv[i] = λ * W_ii * E[β_k²]            (n-vector; depends on current EBeta[k])
```

**Cycle:**

1. At initialization, the Cox warm-start fits a multivariate Cox on K SVD loadings of the merged matrix. The dominant SVD direction is the RNA-seq vs. proteomics batch axis → SVD loadings carry minimal survival signal → EBeta ≈ 0.
2. EBeta ≈ 0 → A_surv ≈ 0 → L update is entirely genomics-driven.
3. L converges to reconstruction-optimal directions (batch-dominated).
4. β update sees batch-dominated factors → selects nothing → EBeta stays 0.
5. Return to step 2.

**Scale imbalance, separate from the cold-start:** Even if EBeta were modestly non-zero, A_gen sums over p = 2000 genes while A_surv contributes a scalar per subject. The ratio A_surv / A_gen is structurally small at all times, not just at initialization. Warm-start Exp 2's collapse-in-23-iterations is the empirical signature of this imbalance.

### 1.4 Multi-platform-specific considerations

The merged Y matrix concatenates RNA-seq and proteomics features column-wise after v2 preprocessing (intersect → quantile normalize → log2 transform). Three properties of the merged matrix exacerbate the failure:

- **Dominant batch factor in SVD:** The first SVD component captures platform identity, not biology. SVD-warm-started L therefore reflects platform, not survival.
- **Heterogeneous noise scales:** RNA-seq counts and protein abundances have different post-preprocessing variance structures. The shared τ_j vector treats them identically; some platform-specific scaling may be informative.
- **Feature-space asymmetry:** Some genes are measured on both platforms, producing correlated columns of Y that can inflate A_gen (the sum-over-features genomics precision) without adding biological information.

Single-cohort runs do not exhibit this failure (e.g., TCGA_PAAD alone converges with non-zero β), confirming the multi-platform setting is the trigger.

---

## §2. Goals

In priority order:

1. **Recover non-zero, stable EBeta** when training on merged multi-platform data
2. **Generalize to external cohorts** (single-platform, e.g., RNA-seq-only test sets)
3. **Preserve the Bayesian framework** (EBNM priors, posterior uncertainty) where possible
4. **Maintain the existing test infrastructure** (171/171 tests passing)
5. **Document train/test alignment** so external validation is principled, not ad hoc

Non-goals:

- Achieving state-of-the-art C-index on every cohort (we want correct mechanism, not benchmark performance)
- Replacing CAVI with gradient-based optimization unless both Cluster A and Cluster B fail
- Building platform-specific submodels unless evidence emerges that shared τ_j is the binding constraint

---

## §3. Approach Landscape

Three clusters of approaches, ordered by effort and degree of structural change:

| Cluster | Strategy | Effort | Risk | Recommended Order |
|---|---|---|---|---|
| **A** | In-model fixes: change initialization, scheduling, and optimization order without altering the generative model | 1–3 days | Medium — does not address the structural scale imbalance | **First** |
| **B** | Cox-on-YF reformulation: change the survival linear predictor from `η = Lβ` to `η = (YF)β` | 1–2 weeks (incl. derivations) | Low conceptual risk; medium implementation risk | **Second (if A insufficient)** |
| **C** | Architectural rewrites: platform-specific noise, gradient-based ELBO, supervised NMF replacing EBNM | 1–3 months | High implementation risk; loses some Bayesian properties | **Deferred** |

The cluster choice ultimately depends on **what fails after Cluster A**:

- Full success → Cluster B becomes optional research; ship Cluster A
- EBeta non-zero but unstable → Cluster B's structural argument becomes load-bearing
- EBeta still zero → Either Cluster A's Fix 3 (ridge Cox) or Cluster B's reformulation

---

## §4. Cluster A Design — In-Model Fixes

This cluster is documented in detail in [`docs/update_L_fix.md`](update_L_fix.md). The implementation plan below is the actionable summary.

### 4.1 Step 1 — Instrument first

**Why first:** We have a hypothesis (chicken-and-egg + scale imbalance) but no quantitative measurement of the imbalance during a real fit. Logging answers a critical interpretive question without changing implementation order:

- If EBeta is **already non-zero** after the Cox warm-start (line ~256 of `fit_modular.R`), the cold-start is not the entry point — the scale imbalance is the more likely dominant issue. Fix 4 (§4.8) is more likely to be needed; Fixes 1+2 may be insufficient on their own.
- If EBeta ≈ 0 after warm-start, the cycle starts at initialization, and Fixes 1+2 (§4.2–4.3) targeting the cold-start are the natural first attempt.

In either case the implementation order is the same: instrument → Fix 1 → Fix 2 → tests → smoke fit → decision gate (§4.6). The instrumentation calibrates expectations and lets us interpret the smoke fit results correctly.

**Where:** `code/fit_modular.R`

**What:**

```r
# After the Cox warm-start block (after line 256):
if (verbose) cat(sprintf("    [init] EBeta range: [%.3e, %.3e]\n",
                         min(EBeta), max(EBeta)))

# Inside the inner k-loop, at iter == 1 only (after Taylor expansion):
if (iter == 1 && verbose) {
  A_gen_k  <- sum(Tau * EF2[, k])
  A_surv_k <- mean(w) * EBeta2[k] * lambda
  cat(sprintf("    [iter1, k=%d] A_gen=%.2e  A_surv=%.2e  ratio=%.4f\n",
              k, A_gen_k, A_surv_k, A_surv_k / (A_gen_k + 1e-30)))
}
```

### 4.2 Step 2 — Fix 1: Reorder inner k-loop (β → L → F)

**Current order:** `compute R_k → compute z_no_k → L_k update → recompute R_k → F_k update → β_k update`

**Proposed order:** `compute R_k → compute z_no_k → β_k update → L_k update → recompute R_k → F_k update`

**Why this is safe (no extra recomputation):** Expanding `z_no_k`:

```
z_no_k = z − (EL %*% EBeta − EL[,k] * EBeta[k])
       = z − Σ_{k' ≠ k} EL[,k'] * EBeta[k']
```

`z_no_k` depends on EL[,k'] and EBeta[k'] for k' ≠ k, but **not on EL[,k] or EBeta[k]** — both cancel by the same algebra. The current code's reuse rationale (Companion.tex Sec. 6) is "EL[,k] cancels"; the symmetric statement "EBeta[k] cancels" also holds. After the β update changes EBeta[k], `z_no_k` for the same k is unchanged, so the L update can use it directly.

Similarly, `R_k = Y − Σ_{k' ≠ k} EL[,k'] * EF[,k']'` depends on EL[,k'] and EF[,k'] for k' ≠ k — neither EL[,k] nor EF[,k] for the current k. The β update doesn't touch EL or EF at all, so R_k from before the β update is still valid for the L update.

**Implementation cost:** A 3-block rearrangement. Same number of `compute_R_k()` calls (2: before L, before F) and `compute_z_no_k()` calls (1) as the current code. The benefit comes purely from the L update seeing the freshly-updated EBeta[k] in its A_surv and B_surv terms.

### 4.3 Step 3 — Fix 2: β-only burn-in phase

Add `N_burnin` parameter to `fit_supervised_mf_modular()` (default 0 = backward-compatible). Pre-CAVI loop:

```r
if (N_burnin > 0) {
  EL2_init <- EL^2   # point estimates; zero posterior variance during burn-in
  for (b in seq_len(N_burnin)) {
    eta_b    <- as.vector(EL %*% EBeta)
    taylor_b <- calc_cox_taylor(eta_b, time, status)
    z_b      <- eta_b + taylor_b$u / taylor_b$w
    w_b      <- taylor_b$w
    for (k in seq_len(K)) {
      z_no_k_b <- compute_z_no_k(z_b, EL, EBeta, k)
      res_b    <- update_beta_k(w_b, z_no_k_b, EL[, k], EL2_init[, k],
                                prior_family = prior_beta, alpha = alpha)
      EBeta[k]  <- res_b$mean
      EBeta2[k] <- res_b$second
    }
  }
}
```

This directly replicates Warm-start Exp 1 (which proved β is functional when EL is fixed). First trial: `N_burnin = 10`.

### 4.4 Step 4 — A2: Progressive α schedule

Add `alpha_schedule` parameter to support curriculum learning:

```r
# alpha_schedule = list(warmup_iters = 5, ramp_iters = 10)
# At outer iteration `iter`:
if (iter <= alpha_schedule$warmup_iters) {
  alpha_iter <- 0   # pure genomics
} else if (iter <= alpha_schedule$warmup_iters + alpha_schedule$ramp_iters) {
  prog <- (iter - alpha_schedule$warmup_iters) / alpha_schedule$ramp_iters
  alpha_iter <- alpha * prog
} else {
  alpha_iter <- alpha   # target reached
}
```

This lets L settle into biologically meaningful directions before survival pressure is applied. Cheap to implement; pairs well with Fix 2 (burn-in primes β; ramped α prevents L from being yanked toward survival before genomics structure exists).

### 4.5 Step 5 — Run tests

```bash
Rscript tests/run_tests.R
```

Expected: 171/171. Any failure in `test_update_L.R` or `test_update_beta.R` blocks merge. The reorder (Fix 1) does not change the EBNM math, so existing tests should still pass; the burn-in adds new behavior behind a default-off parameter.

### 4.6 Step 5 decision gate

Run a smoke fit on merged TCGA_PAAD + CPTAC v2-preprocessed data with `N_burnin = 10`, verbose:

| Outcome | Action |
|---|---|
| EBeta non-zero AND stable through CAVI | **Done.** Validation phase (Phase 3). |
| EBeta non-zero but collapses during CAVI | Implement **Fix 4** (A_surv/A_gen normalization in `update_L_k()`). Scale imbalance is confirmed as the dominant issue. |
| EBeta still ~ 0 after burn-in | Implement **Fix 3** (ridge Cox warm-start via `glmnet`). SVD loadings carry no survival signal. |

### 4.7 Fix 3 (conditional) — Ridge Cox warm-start

Replace the `coxph()` block in `fit_modular.R` lines 244–257 with:

```r
library(glmnet)
tryCatch({
  cv_fit <- cv.glmnet(EL, survival::Surv(time, status),
                      family = "cox", alpha = 0)  # alpha=0 → pure ridge
  EBeta <- as.vector(coef(cv_fit, s = "lambda.min"))
}, error = function(e) {
  EBeta <<- rep(0, K)
})
EBeta2 <- EBeta^2
```

Ridge enforces non-sparsity → all K factors get non-zero EBeta. Any non-zero value breaks the A_surv ≈ 0 cycle. `glmnet` is likely already installed — verify before implementing.

### 4.8 Fix 4 (conditional) — Normalize A_surv/A_gen

Modify `code/update_L.R` lines 140–158:

```r
A_gen  <- sum(Tau * EF2_k)                                # scalar
A_surv <- lambda * w * EBeta2_k                           # n-vector

A_gen_norm  <- A_gen  / (mean(A_gen)  + 1e-30)            # = 1
A_surv_norm <- A_surv / (mean(A_surv) + 1e-30)            # n-vector; guard mean=0
A_L <- pmax((1 - alpha) * A_gen_norm + alpha * A_surv_norm, A_floor)

B_gen  <- as.vector(R_k %*% (Tau * EF_k)) / (mean(A_gen)  + 1e-30)
B_surv <- lambda * w * z_no_k * EBeta_k    / (mean(A_surv) + 1e-30)
B_L    <- (1 - alpha) * B_gen + alpha * B_surv
```

**Caveat:** This departs from strict ELBO maximization (the normalization is not derivable from the generative model). Must verify ELBO still improves monotonically and update `tests/test_update_L.R` to cover the new code path. Edge case: when EBeta_k = 0 exactly, `mean(A_surv) = 0` — the `+ 1e-30` guard prevents division by zero but the survival contribution evaluates to 0/ε ≈ 0, which is the desired behavior.

### 4.9 Cluster A success criteria

On the merged TCGA_PAAD + CPTAC v2-preprocessed training set:

- ≥ 1 factor with `|EBeta| > 0.05` (any prior)
- Full CAVI does not collapse EL (`max|EL| > 0.01` throughout)
- ELBO improves monotonically
- 171/171 tests pass
- External C-index on ≥ 1 cohort improves vs. SVD-initialized baseline

---

## §5. Cluster B Design — Cox-on-YF Reformulation

### 5.1 The proposal

Replace the survival linear predictor:

```
Current:    h_i(t) = h_0(t) * exp(l_i · β)              where l_i is the i-th row of L
Proposed:   h_i(t) = h_0(t) * exp(y_i · F · β)          where y_i is the i-th row of Y
```

Equivalently, `η = (YF)β` — the Cox model is supervised by the projection of observed Y onto the factor matrix F, rather than by the latent loadings L.

### 5.2 Why this addresses the chicken-and-egg

`YF` depends on observed data Y (fixed) and learned F (initialized via SVD). Even when SVD loadings carry no survival signal, `YF` carries the signal that's actually present in Y projected onto the SVD subspace. The Cox warm-start on `YF` therefore produces non-zero EBeta from iteration 1, and the F update's survival precision is non-zero from the start.

### 5.3 Train/test alignment argument

The current model's train/test mismatch:

| Phase | Linear predictor | Latent or observed? |
|---|---|---|
| Training | `η_i = l_i · β` | l_i is latent EBNM posterior |
| Test (external cohort) | `η_test = (Y_test · F · (F'F)⁻¹) · β` | OLS projection of observed Y_test onto F |

The training-time `l_i` is not the OLS projection — it's pulled toward zero by the EBNM prior and modulated by survival via `A_surv`. So the model is trained on one quantity but evaluated on a different one. Cluster B closes this gap: training and test both use the projection `(YF(F'F)⁻¹)β`.

A useful reparameterization absorbs `(F'F)⁻¹` into β:

```
β̃ = (F'F)⁻¹ β    →    η = (YF) β̃
```

At prediction time, `β̃` is the coefficient on the projection scores `Y · F`. Numerical caveat: when `F'F` is near-singular, the inverse hides instability rather than solving it; a small ridge `(F'F + ε I)⁻¹` is recommended for stability.

### 5.4 Update structure under Cluster B

| Update | Current role | Cluster B role |
|---|---|---|
| q(L) | Dual-source (genomics + survival) | **Genomics only** — simpler than current q(F) |
| q(F) | Genomics only | **Dual-source (genomics + survival)** — inherits structure of current q(L) |
| q(β) | Working response z_no_k = z − (Lβ − L[,k]β_k) | Working response defined via `(YF)β`: z_no_k = z − ((YF)β − (YF)[,k]β_k) |
| q(τ) | Closed-form MLE | **Unchanged** |

**Symmetry observation:** L and F play structurally symmetric roles in the bilinear term `Y = LF'`. Cluster B is a "swap" of which latent variable carries the survival burden. The math derivations for the new q(F) supervised update should be analogous to the existing q(L) derivation, with the n-dimensional and p-dimensional roles swapped.

### 5.5 Required derivations (Phase 4 work)

Deliverables in `derivations/qF_supervised/`:

| File | Purpose |
|---|---|
| `qF_supervised_derivation.tex` | New q(F) dual-source EBNM update under η = (YF)β |
| `qBeta_YF_derivation.tex` | Updated q(β) with z_no_k redefined via YF |
| `ELBO_YF_derivation.tex` | Full ELBO under the reformulation; KL terms re-expressed |
| `qL_unsupervised_derivation.tex` | Reduced q(L) update (genomics only) — simpler, smaller change |

**Critical analytical checks during derivation:**

1. **Scale imbalance audit for the new F update.** A_gen for q(F) sums over n samples; the new A_surv contribution is roughly `Σᵢ y_{ij}² · W_ii · E[β_k²]` (sketch — needs careful derivation). Verify whether this ratio exhibits the same imbalance that crippled q(L), shifted to a different domain.
2. **Handling of (F'F)⁻¹.** Decide between explicit pseudoinverse (carry through every CAVI step) vs. β reparameterization (absorb into β̃) vs. ridge regularization. Each has different numerical and interpretability consequences.
3. **Identifiability.** Under η = (YF)β with both F and β learned, is there a rotation/scaling ambiguity that needs anchoring? The current model has L and F sharing a scale ambiguity that's broken by the EBNM prior; the new model may have additional ambiguities through the YF product.
4. **Train/test consistency check.** Confirm analytically that `η_test = Y_test · F · β̃` matches what's used during training, i.e., that no scaling factor (e.g., the missing F'F term in the reparameterization) leaks in only at test time.

### 5.6 Risks specific to Cluster B

- **Scale imbalance may not vanish, just shift.** The fundamental issue is that genomics has p×n data points and survival has n events. Whichever latent variable bridges to survival inherits this asymmetry. Cluster B reframes the problem; it doesn't prove the imbalance disappears.
- **Loss of L's biological interpretability.** Currently L is the variable with prognostic loadings (per-patient survival risk decomposition). After Cluster B, L is purely a reconstruction quantity — interpretation shifts to "Y · F · β̃ = patient survival score" which uses gene-level programs F directly.
- **Implementation surface area.** All four update modules (`update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R` indirectly via Tau coupling) and their tests need revision. Test count will change from 171 to whatever the new structure requires.
- **Backward incompatibility.** Existing fits saved with the current model cannot be directly compared to Cluster B fits — the β̃ coefficients are on a different scale. External validation pipelines that consume `EBeta` directly need updating.

### 5.7 Cluster B success criteria

Same as Cluster A (§4.9), evaluated on the same merged training set and external cohorts. Additional requirement specific to Cluster B:

- Train/test predictions are computed by the same formula (`Y · F · β̃`), with no special-case projection logic at inference time.

---

## §6. Cluster C — Deferred Architectural Changes

Documented for completeness; not pursued unless both Cluster A and Cluster B fail.

### 6.1 Platform-specific noise + shared/private factor structure

Generative model: separate τ_p per platform; some factors are platform-private (technical/batch), others are shared (biological); only shared factors couple to survival. Requires a full ELBO rewrite and updates for ~6 parameter types. Most principled for the multi-platform case but largest implementation cost.

### 6.2 Gradient-based ELBO maximization

Replace CAVI with auto-diff gradient ascent (e.g., PyTorch). Avoids local optima — the current β=0 is a local stationary point of the ELBO under coordinate updates; full gradients can escape it. Loses much of the EBNM machinery and posterior uncertainty quantification.

### 6.3 Supervised NMF / penalized factorization

Drop EBNM in favor of explicit L1/L2 penalties; reformulate as a convex (or biconvex) optimization with a Cox loss term. Connects to existing literature (e.g., supervised PCA, survival-supervised NMF). Loses Bayesian inference but gains optimization tractability.

---

## §7. Sequencing & Verification

### 7.1 The 5 phases

| Phase | Session | Deliverable | Decision Gate |
|---|---|---|---|
| 1 | Current | This design doc | User approval |
| 2 | New | Cluster A code on `fix-L-update-beta-cycle` branch | 171/171 tests pass + smoke fit |
| 3 | New | Validation report on merged + external cohorts | Cluster A success criteria met? |
| 4 | New | Cluster B derivations in `derivations/qF_supervised/` | Algebraic review |
| 5 | Conditional | Cluster B implementation on `cox-on-YF-reformulation` branch | Same as Phase 3 |

### 7.2 Session boundaries (explicit handoffs)

🔄 **SESSION BOUNDARY 1 — between Phases 1 and 2.**
*Pause when:* design doc is approved.
*Reason:* the design phase and implementation phase are different cognitive modes. A fresh session for implementation gets a clean context window dedicated to R code, test running, and git branch management — without the design discussion crowding the context.
*Resume prompt for new session:* "Implement Cluster A from `docs/beta_zero_fix_design.md` §4 on a new branch `fix-L-update-beta-cycle`. Start with Step 1 (instrumentation) and proceed through the decision gate at §4.6."

🔄 **SESSION BOUNDARY 2 — between Phases 2 and 3.**
*Pause when:* Cluster A is implemented and 171/171 tests pass.
*Reason:* validation involves long-running R fits across multiple cohorts, plot inspection, and write-ups. The implementation session's context is dominated by code-level details that aren't useful for validation analysis.
*Resume prompt for new session:* "Validate Cluster A per `docs/beta_zero_fix_design.md` §4.9. Run on merged TCGA+CPTAC + external cohorts. Update `docs/update_L_fix.md` and `DECISIONS.md` with results."

🔄 **SESSION BOUNDARY 3 — between Phases 3 and 4.**
*Pause when:* Cluster A validation is complete and the decision gate (§4.6, §4.9) is resolved.
*Reason:* derivation work is math-heavy LaTeX and benefits from a session dedicated to algebra, not entangled with R debugging or validation results.
*Resume prompt for new session:* "Sketch Cluster B derivations in `derivations/qF_supervised/` per `docs/beta_zero_fix_design.md` §5.5. Use existing `derivations/qL/qL_update_derivation.tex` as a structural template — F now plays the role L currently plays."

🔄 **SESSION BOUNDARY 4 — between Phases 4 and 5.**
*Pause when:* derivations are complete and reviewed.
*Reason:* implementation of Cluster B depends on whether Cluster A succeeded (Phase 3) and on derivation findings (Phase 4); a fresh session is the right place to make the go/no-go call with both results in hand.
*Resume prompt for new session (if go):* "Implement Cluster B per derivations in `derivations/qF_supervised/`. New branch `cox-on-YF-reformulation`. Update `code/update_F.R`, `code/update_L.R`, `code/update_beta.R`, and create new tests."

### 7.3 Verification per phase

| Phase | Success Criteria |
|---|---|
| 1 | Design doc covers both clusters; user approves; self-review finds no placeholders or contradictions |
| 2 | Code compiles; 171/171 tests pass; instrumentation prints expected values on a smoke run |
| 3 | At least 1 \|EBeta\| > 0.05 on merged training set; ELBO monotone; max\|EL\| > 0.01 throughout; external C-index ≥ SVD baseline on ≥ 1 cohort |
| 4 | Derivations algebraically verified (peer review or independent rederivation); train/test alignment confirmed |
| 5 | Same as Phase 3 success criteria, on the new model |

### 7.4 Rollback points

- After Phase 2: if 171/171 tests cannot be made green, revert the branch to `main` and reconsider.
- After Phase 3: if Cluster A produces non-zero β only on single cohorts (not merged), proceed to Phase 4. If Cluster A fully succeeds on merged, Phases 4–5 become optional research.
- After Phase 4: if derivations reveal Cluster B has the same scale imbalance shifted to a different domain, escalate to Cluster C (treated as out-of-scope research planning, not implementation).

---

## §8. Open Questions

To resolve during implementation:

1. **Cluster A — Fix 4 hyperparameter:** The normalization `mean(A_surv) + 1e-30` uses the empirical mean across subjects. Should we use a more robust scale (median, trimmed mean) or expose the normalization scale as a tunable parameter?
2. **Cluster A — Progressive α schedule:** Is the linear ramp the right shape, or should it be sigmoidal / cosine-annealed? This is small but worth a sensitivity check.
3. **Cluster B — Identifiability:** Under η = (YF)β, what anchors the rotation/scaling of F? Does the EBNM prior on F suffice, or do we need an additional constraint (e.g., diag(F'F) = 1)?
4. **Cluster B — F'F regularization:** If we keep the explicit (F'F)⁻¹ form (vs. β-reparameterization), what value of ε in `(F'F + ε I)⁻¹` is principled? Should it be a hyperparameter or set adaptively from the spectrum of F'F?
5. **Cluster B — Backward compatibility:** Should we keep the current model as `cox-on-L` and the new model as `cox-on-YF` as a `model_type` parameter in `fit_supervised_mf_modular()`, or is a clean break preferable?
6. **Multi-platform: per-platform τ:** Cluster A leaves τ_j shared across platforms. If Cluster A fully succeeds, does the residual variance pattern suggest τ_j is the next bottleneck, or is single-τ adequate? (Inspect `τ` distribution post-fit to decide.)
7. **External validation pipeline:** Does the existing `predict_supervised_mf()` in `code/predict.R` need updates for Cluster A (no — it already uses the projection at test time), or for Cluster B (yes — coefficient interpretation changes from β to β̃)?

These do not block Phase 1. They are checkpoints for Phases 2–5.
