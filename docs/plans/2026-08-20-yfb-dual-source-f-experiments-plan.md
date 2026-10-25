# YFB Dual-Source F Update — Experiment Plan

**Date:** 2026-08-20
**Source:** Review of `docs/notes/YFB_derivation_05_08_26.pdf` (Section 13 caveat) against
`DECISIONS.md` 2026-05-04, 2026-05-22, 2026-07-15, and the DeSurv comparison
([[reference-desurv-paper]]). Not yet implemented — this is a plan only.

**Branching:** one branch per step, cut from `main`, per standing repo convention (TDD for new
functions, `superpowers:code-reviewer` dispatch before each commit, full test suite run,
`DECISIONS.md`/`ROADMAP.md` updated with real numbers regardless of outcome, no push without
confirmation).

---

## Background

The recommended YFB model (η = (YF)β̃) fixes `alpha_F=0` in `update_F_surv_YFB_k`, meaning F is
learned from genomics only and survival enters solely through β acting on the fixed `YF` scores.
This was a deliberate fix for a positive-feedback instability (`DECISIONS.md` 2026-04-30), and a
2026-05-22 retest of `alpha_F ∈ {0.1, 0.3, 0.5}` reproduced the same failure: RMSE jumped from
~290 to 750-800 and β still collapsed to 0.

**DeSurv** (Young et al., PNAS 2026) uses a structurally analogous dual-source loading matrix `W`
(survival gradient acts on `W` directly, mixing weight `α=0.334`, tuned jointly with `k` via
Bayesian optimization) and reports it working, with survival concentrated in one factor. That rules
out "dual-source F is inherently unstable" as an explanation and points to two specific mechanical
differences instead:

- **H1 — shrinkage-spike floor.** Our β update is an EBNM point-normal/point-exponential posterior
  mean, which can land exactly at its shrinkage spike (β=0). DeSurv's β is an elastic-net Cox
  regression estimate, which stays small-but-nonzero under weak signal. F's survival term scales
  with `EBeta_k²`/`EBeta_k` (`code/update_F_surv_YFB.R`), so an exact-zero β provides no gradient at
  all, while a small-but-nonzero one still nudges F every iteration.
- **H2 — thin search over the escape mechanism.** DeSurv jointly tunes `(k, α, λ)` via Bayesian
  optimization across many candidate configurations; we tried 3 fixed `alpha_F` values once each,
  from a cold start, with no warm-up/restart strategy.

These two hypotheses are independent and are tested as separate experiments below, then combined
in Step 4.

**Guardrails (apply to every step):**
- No change to `config/globals.yml` defaults, `update_beta.R`'s default behavior, or
  `update_F_surv_YFB.R`'s `alpha_F=0` default until Step 4 concludes with a clear win — all new
  behavior is opt-in via explicit new parameters, so the existing 374/374 + 88/88 test suites stay
  green unchanged.
- LB (Cluster A) is out of scope. YFB/Cluster B only.
- Reuse existing infrastructure where it fits (`results/benchmark_sim/run_yfb_beta_fix_diagnostic.R`
  for diagnostics; `code/select_k_alpha_bo.R`'s `rBayesianOptimization` wrapper and
  `pick_trustworthy_bo_winner` validity-gate pattern for the BO piece) rather than duplicating it.
  Note: `select_k_alpha_bo.R`'s `alpha` is the Phase 1a objective-normalization weight passed to the
  β update — a different parameter from `alpha_F`. Step 3 extends the same BO *pattern* to a new
  `alpha_F` search dimension; it does not reuse that script's existing search space as-is.

---

## Step 0 — Baseline reconfirmation

**Branch:** `yfb-dualF-baseline-check`

**Goal:** confirm the May 2026 failure (RMSE blowup, β collapse at `alpha_F ∈ {0.1,0.3,0.5}`) still
reproduces under the current codebase before investing in fixes — three months of Phase 1-3 changes
(objective normalization, beta_threshold, cohort_id, K=7 selection) sit on top of the original
diagnostic run.

**Implementation:** rerun `results/benchmark_sim/run_yfb_beta_fix_diagnostic.R` unchanged against
the current recommended D4 config (YFB, per-platform z-std, DeSurv gene selection, K=7,
`alpha_F ∈ {0, 0.1, 0.3, 0.5}`). Record RMSE trajectory and final β magnitudes for each.

**Decision rule:** if the failure no longer reproduces (e.g. Phase 1-3's fixes incidentally
resolved it), stop here, update `DECISIONS.md`, and re-scope Steps 1-4 — they assume the failure is
still live. Otherwise proceed to Steps 1 and 2.

---

## Step 1 — Experiment A: non-degenerate β estimator (attacks H1)

**Branch:** `yfb-dualF-ridge-beta`

**Goal:** test whether the hard shrinkage-spike floor, not dual-source F itself, is what blocks
`alpha_F>0` from converging.

**Implementation:**
1. Add a new, explicitly experimental β-update variant — a ridge/L2-regularized posterior-mode
   estimate (small fixed penalty, no spike-and-slab component) that is structurally guaranteed to
   return a strictly nonzero point estimate under weak signal, matching what an elastic-net Cox
   step would produce. Keep this isolated (e.g. a new function in
   `results/benchmark_sim/` or a clearly-labeled experimental module) — do **not** change
   `update_beta.R`'s default `EBNM`-based update or its call sites in the production YFB path.
2. Wire this alternate estimator into a copy of the CAVI loop used only for this experiment (do not
   add a permanent branch inside `fit_cox_on_yf.R` at this stage), and rerun the `alpha_F ∈
   {0.1, 0.3, 0.5}` sweep from Step 0 with it in place of the EBNM β update.
3. TDD the new estimator function in isolation (fixed inputs → known ridge-regression output);
   dispatch code-reviewer; run the full test suite to confirm nothing in the production path moved.

**Decision rule (mechanical):** for each `alpha_F` tested, check (a) RMSE stays within 2× of the
`alpha_F=0` baseline (no blowup) and (b) β does not decay back toward ~0 by the final iteration
(final |β_max| stays within an order of magnitude of its post-warm-up value). If **both** hold for
at least one `alpha_F>0` value → H1 supported, proceed to full external validation for that
configuration in Step 4. If **neither** holds for any `alpha_F>0` → H1 rejected; record in
`DECISIONS.md` and rely on Step 3's result alone.

---

## Step 2 — Experiment B: β warm-up schedule + broader α_F search (attacks H2)

**Branch:** `yfb-dualF-alpha-schedule`

**Goal:** test whether giving β room to move away from zero *before* the survival term enters F's
update — combined with a properly searched (not hand-picked) `alpha_F` — lets dual-source F
converge, operationalizing the frozen-F preconditioning idea from `DECISIONS.md` 2026-05-22.

**Implementation:**
1. Add an `alpha_F_schedule = NULL | list(warmup_iters, ramp_iters)` parameter to
   `fit_cox_on_yf.R`, mirroring the existing `alpha_schedule` mechanism (`code/fit_cox_on_yf.R:205`,
   line ~298): hold `alpha_F=0` for `warmup_iters` (β updates freely off the fixed genomics-optimal
   `YF`), then ramp `alpha_F` linearly up to its target over `ramp_iters`. Default `NULL` preserves
   current behavior exactly (`alpha_F` fixed at whatever is passed, no schedule) — existing tests
   unaffected.
2. Extend the BO pattern in `code/select_k_alpha_bo.R` (new function, e.g.
   `select_k_alphaF_bayesopt()`, not a modification of the existing `alpha`-search function) to
   search jointly over `K` and `alpha_F_target` (with `warmup_iters`/`ramp_iters` either fixed at a
   reasonable default or included as a third search dimension if the two-dimensional search proves
   cheap enough), using mean CV external C-index as the objective — same fold-fitting machinery,
   same `pick_trustworthy_bo_winner`-style validity gate adapted to flag a degenerate
   `alpha_F_target≈0` winner (which would just reconfirm the current default, not a genuine result).
3. TDD the schedule ramp function and the new BO wrapper (mirroring `tests/test_select_k_alpha_bo.R`
   structure); dispatch code-reviewer; run full test suite.

**Decision rule:** does the BO-selected `(K, alpha_F_target)` achieve mean external C-index within
1 SE of the current baseline (0.627, margin ≈0.598, matching the project's existing 1-SE
convention) **and** land at an `alpha_F_target` meaningfully greater than 0 (not a degenerate
near-zero winner)? If yes → H2 supported, proceed to Step 4. If no → H2 rejected; record in
`DECISIONS.md`.

---

## Step 3 — Synthesis and final comparison

**Branch:** `yfb-dualF-synthesis` (only if Step 1 and/or Step 2 supported their hypothesis;
otherwise skip and just write up the negative result)

**Goal:** produce one clean, apples-to-apples verdict on whether dual-source F is viable for YFB on
this data, run through the exact same external validation protocol as the current recommended
config (same 5 held-out cohorts, same preprocessing).

**Implementation:**
1. If both H1 and H2 were supported, test them combined (ridge-style β estimator + scheduled
   `alpha_F`) as the strongest candidate. If only one was supported, carry that one forward alone.
2. Run the candidate through `results/benchmark_sim/run_YFB_benchmark.R`'s external validation path
   (or an equivalent script matching its exact preprocessing/gene-selection/cohort protocol) against
   all 5 external cohorts.
3. Update `DECISIONS.md` with the full comparison table (candidate vs. current `alpha_F=0`
   recommended config) regardless of outcome — win, tie, or loss.

**Decision rule for adopting a change to the default:** only change `config/globals.yml` /
`update_F_surv_YFB.R`'s default if the candidate **clearly beats** (not just ties) the current
recommended config's mean external C-index by more than 1 SE. A tie is recorded as "an equally
valid alternative, current default kept for reproducibility/simplicity" — consistent with how the
K=4/K=5 warm-start ties were handled relative to K=7 (`DECISIONS.md` 2026-07-13).

---

## What this plan will settle

- Whether "F only sees genomics, survival enters solely through β" is a genuine model limitation we
  should fix, or a stable, already-validated design choice that happens to have a plausible but
  ultimately non-beneficial alternative.
- If a fix works, whether it's worth adopting given the added complexity (new schedule parameter,
  new estimator, no CAVI-derivation-purity guarantee for the ridge β estimator) versus the
  reproducibility-simplicity of the current single-fresh-fit `alpha_F=0` default.
- If no fix works, a documented, DeSurv-informed explanation for *why* — not just "we tried alpha
  and it broke" — suitable for stating honestly in the manuscript/talk if the dual-source question
  comes up.
