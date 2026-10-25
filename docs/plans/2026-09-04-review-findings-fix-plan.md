# Fix Plan — Codex Review Findings on the 2026-09-04 Chapter

> **Branch:** `fix/2026-09-04-review-findings` (created off `main` at the commit `main` was at when
> this plan was written — `main` itself is left untouched by this work). **Why a separate branch:**
> the 2026-09-04 progress-book chapter and its underlying benchmark outputs are needed, as-is, for
> an upcoming lab meeting. Doing this fix work on a branch means the meeting can show the current
> `main` version of the Quarto book regardless of whether these fixes are finished in time — nothing
> here needs to land before the meeting. Merge to `main` only once Steps 1-6 below are complete and
> reviewed, not partially.
>
> **Provenance:** this plan was produced from a two-round Codex code review
> (`docs/reviews/2026-09-04_progress_notebook_review.md`, reviewed at commit `a650491` then updated
> through `6196b59`) plus direct independent verification of the highest-priority findings against
> the current code, done in a Claude Code session on 2026-09-04. See that session's transcript or
> `DECISIONS.md` (2026-09-04 entries) for the full narrative. This file is the actionable version.

## Ordering constraint (do not reorder or parallelize)

Per Codex's explicit follow-up: implement and validate Breslow tied-event handling FIRST. Then
define the orientation convention from training data only and use it consistently everywhere.
Then fix `EBeta_pooled`. Only after all three are done and tested, re-run downstream analyses.

---

## Step 1: Breslow tied-event-time Cox likelihood (first, before anything else)

`calc_cox_taylor_yf()` in `code/fit_cox_on_yf.R` computes one risk-set denominator per sorted row.
For tied event times, Breslow's method requires the SAME denominator for every event at that time.

**Verified bug:** a 3-observation synthetic check gives logPL = -1.864706 from this implementation
vs. -2.102889 from `coxph(..., ties="breslow")`; permuting the order of two tied rows changes the
implementation's answer to -1.244592 (should be invariant — it isn't). Real training data has 9
TCGA + 12 CPTAC event rows with times tied to another row, so this is not a hypothetical edge case.
This function is used in fitting itself (CAVI u/w weights), `code/compute_bic.R`, and
`code/compute_cv_loglik.R` — everything downstream is stale until this is fixed.

**Tasks:**
- Implement correct Breslow-tied score, diagonal Hessian, and log partial likelihood.
- Cross-check against `survival::coxph(..., ties="breslow")` on multiple synthetic tied-time cases.
- Add a regression test: permuting rows with identical event times must not change the result.
- Existing `strata` support in `calc_cox_taylor_yf` (added 2026-09-04) must keep working — ties are
  computed within each stratum independently.
- Do NOT re-run any downstream benchmark yet — that happens in Step 4.

---

## Step 2: Orientation convention, defined from training data only, used everywhere

Two different orientation hacks currently exist and disagree with each other:

**(a) Phase C** in `fit_cox_on_yf.R` (~line 842-863) flips `EBeta`'s sign based on training
concordance, but the `concordance()` call omits `reverse=TRUE` — `survival::concordance()`'s
formula method assumes larger predictor = LONGER survival by default, the opposite of a Cox risk
score. This makes the flip decision backwards. (Already documented as a "KNOWN ISSUE" in the code
comments.)

**(b) `oriented_cindex()`** in `results/benchmark_sim/run_cohort_beta_comparison.R` (line 149-152)
instead computes `max(c_raw, 1-c_raw)` using the EXTERNAL test cohort's own outcomes — i.e., it
picks the orientation by looking at the answer it's trying to score. This is circular: it's not
evaluating a prospectively-oriented predictor, it's reporting a sign-invariant association
statistic.

**Context on why (b) exists** (Andrew's note, 2026-09-04): likely introduced because early
cross-validated C-index results came back very small (well below 0.5), and `max(c, 1-c)` was a
quick way to avoid reporting nonsensical numbers. It's possible the root cause (the Phase C
orientation bug, or some other bug since fixed) is what produced those small numbers, and the hack
may no longer be necessary — worth checking `DECISIONS.md` history for exactly when/why it was
added before assuming it's still needed.

**Tasks:**
- Decide ONE orientation convention: establish sign/direction from TRAINING data only (with correct
  `reverse=TRUE` semantics for a Cox risk score), freeze it, and apply the same frozen orientation
  to every external cohort and every bootstrap replicate. No evaluator may look at the data it's
  scoring to decide its own sign.
- Remove or replace `oriented_cindex()`'s `max(c_raw, 1-c_raw)` with scoring under the frozen
  training-derived orientation. If the resulting external C-index comes back near or below 0.5 for
  some configuration, that is itself a valid (if disappointing) finding — report it, don't mask it.
- Apply this consistently across `code/concordance_ci.R`, `run_cohort_beta_comparison.R`, and any
  other benchmark runner using an orientation hack (grep for `"reverse"`, `"max(c"`, `"1 - c"`,
  `"1-c"`).
- Add a test verifying the orientation is fixed at fit time and does not vary by which data it's
  later applied to.

---

## Step 3: `EBeta_pooled` — coherent pre/post-Phase-C state

In `fit_cox_on_yf.R` (~line 903), `compute_pooled_beta(w, z, ZF, rowMeans(EBeta), ...)` is called
AFTER the Phase C sign-flip block, so `rowMeans(EBeta)` is post-flip, but `w`/`z`/`ZF` are whatever
the main CAVI loop last left them (pre-Phase-C — Phase C never recomputes these). Concrete impact
on the cached D4 fit: Program 7's pooled beta is 0.0404 with Phase C disabled vs. 0.0204 with it
enabled — a 2x difference from this inconsistency alone, not from Phase C's actual intended effect.
`compute_pooled_beta()` (`code/update_beta_cohort.R:263`) is a single `update_beta_all()` call —
one Gauss-Seidel sweep, not iterated to convergence — so a bad initializer doesn't wash out.

**Tasks:**
- Per Codex's guidance: implement the PRE-Phase-C coherent-state fallback (i.e., compute
  `EBeta_pooled` from the same pre-flip `w`/`z`/`ZF`/`EBeta` snapshot, all mutually consistent)
  UNLESS a separately-converged shared-beta refit (iterating `update_beta_all` to convergence, not
  one sweep) can be justified as the actually-intended target — if so, implement that instead and
  document why.
- Add a regression test asserting that a global sign flip of the fitted cohort-beta matrix changes
  `EBeta_pooled` by a global sign ONLY (not magnitude) — the test that would have caught this bug.
- Also fix the smaller, related Phase C bug: line ~861 sets `EBeta2 <- EBeta^2` on a sign flip,
  which destroys posterior variance (a sign flip should leave `E[beta^2]` unchanged). Fix: negate
  only `EBeta`, leave `EBeta2` untouched. Add a test that posterior second moments are invariant to
  a sign correction.

---

## Step 4: Fix `compute_cv_loglik.R`'s two held-out scoring gaps, THEN re-run downstream

Both in `code/compute_cv_loglik.R`'s `cv_survival_loglik()`, around line 257-279:

**(a)** The held-out scoring call
`calc_cox_taylor_yf(pred$risk_scores, time[test_idx], status[test_idx])` omits
`strata = strata_id[test_idx]`, even though `calc_cox_taylor_yf` supports a `strata` arg and
training-fold fitting correctly uses it. A model fit under stratified risk sets is being scored
under one pooled risk set. **Fix:** pass held-out strata through. Add a test that the stratified
test-fold likelihood equals the sum of cohort-specific likelihoods.

**(b)** When `beta_cohort_id` is set, EVERY fold scores with `EBeta_pooled`, even though ordinary
within-cohort CV folds have KNOWN held-out cohort labels. This conflates two different targets:
within-cohort CV (should score with `fit$EBeta` + held-out `beta_cohort_id`) vs. unseen-cohort
generalization (should use the repaired `EBeta_pooled` fallback from Step 3). Decide which is
wanted for which report, implement both paths distinctly, and label them clearly in output.

**Only after Steps 1-4 are complete and tested, re-run:**
- `results/benchmark_sim/run_cohort_beta_comparison.R` (5-arm external C comparison)
- `results/benchmark_sim/run_cohort_beta_bootstrap_ci.R` (bootstrap CIs — also fix the estimand
  mismatch: the chapter reports the mean of 5 within-cohort C-indices, but the pooled bootstrap
  concatenates all 616 patients across cohorts, giving larger cohorts more weight and including
  between-cohort pairs, so it answers a different question — 0.0135 mean-of-cohorts vs. 0.0221
  pooled-patient difference on the current, pre-fix numbers. Bootstrap PAIRED WITHIN-COHORT
  differences and average with the same weights as the headline metric; state this is conditional
  on the 5 fixed cohorts, not a claim about generalization to a new cohort.)
- `results/benchmark_sim/run_cohort_beta_supplementary.R` (strata_only, held-out LL per arm)
- `results/benchmark_sim/run_k_init_sweep.R`'s CV/BiCV likelihood columns (tied-time fix changes
  these; use `--reuse-cache` for the ELBO/BIC parts that don't depend on ties, don't reuse the CV
  cache)
- `results/multi_cohort_sim/run_training_set_subanalysis.R` (pooling-vs-single-cohort) is NOT
  affected by ties/orientation in the same way, but its "pooling clearly helps" claim should be
  softened regardless — it changes both sample size AND gene selection simultaneously, so it
  doesn't isolate pooling's effect. Reword rather than re-run.

---

## Step 5: Regenerate chapter + `DECISIONS.md` from the new outputs

`docs/progress_book/chapters/2026-09-04.qmd` currently has THREE CONFIRMED factual errors
independent of the above (verify these still hold against regenerated CSVs, they may or may not
change with the fixes above):

1. §1 table claims in-sample joint log-likelihood prefers K_init=3 — the CSV actually shows K=15
   highest (-750479.9 vs -768147.9 at K=3). ELBO and BIC do prefer K=3; log-likelihood doesn't.
2. §2 claims retained-factor count "stays flat around both dips [K=11,13] under every threshold" —
   `threshold_vs_kinit.csv` shows K_survival_active=3 (not flat at 2) at K=11 for `beta_thresh`
   0.001 AND 0.01, and at K=13 for `beta_thresh=0.001`.
3. §4 delta-figure caption claims "flat or negative below K=6" — the hybrid scenario is strongly
   POSITIVE below K=6, peaking at K=4 (+0.157), larger than at K=6.

Also fix §3's overstated equivalence claim: the pooled CI for baseline minus `beta_cohort_id` was
[-0.0093, 0.0505] (pre-fix) — a CI including zero does not demonstrate equivalence, especially when
the upper bound represents a meaningful loss. Reword to something like: "the point estimate favored
the baseline; this analysis did not detect a difference, while the interval remained compatible
with a meaningful loss" — then regenerate with the corrected estimand from Step 4.

Also narrow (no new computation needed, just reword): partial-pooling description (zero-centered
shared shrinkage, not shrinkage toward an estimated nonzero common mean), the `alpha_F=0`
"mathematically identical" claim (survival can still affect which iterate is returned via the
stopping rule), the top-2 factor-merge diagnostic (label as single-seed exploratory unless
aggregated across all 15 seeds), gene-set overlap p-values (descriptive/nominal, not inferential —
genes selected from the same expression matrix, correlated selection events).

---

## Step 6: Full test suite

Run `Rscript tests/run_tests.R` (currently 443/443) — must stay green. Add all the targeted tests
named above. Update `CLAUDE.md`'s Quick Reference test count if the total changes.

**Do not skip or reorder Steps 1-4.** Ask before starting Step 5 (chapter rewrite) if anything from
Steps 1-4 changes the substantive conclusions (not just the numbers) — e.g., if the corrected
orientation convention reveals the cohort-beta model is NOT roughly equivalent to the pooled model
after all, that's a finding to discuss, not just silently fold into a rewritten table.

---

## Merge criteria

Merge `fix/2026-09-04-review-findings` back to `main` only when Steps 1-6 are all complete, tests
are green, and the chapter/DECISIONS.md reflect the corrected numbers — not partially. Until then,
`main` keeps showing the current (pre-fix) 2026-09-04 chapter, which is what the upcoming lab
meeting will present.
