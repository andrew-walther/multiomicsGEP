# SSBMF — K-Parsimony Follow-Up & Validation Plan

**Date:** 2026-07-13
**Source:** Phase 3 of `docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md`
(`DECISIONS.md` 2026-07-13) found that the recommended config's K=7 is not free to shrink in a
single-seed comparison — but flagged that K=2/K=4's poor performance is suspicious (fast convergence
to near-zero β), consistent with the CAVI factor-collapse vulnerability documented in Phase 2
(`DECISIONS.md` 2026-07-12/13). This plan resolves that open question and produces a validated,
final answer plus a research-team-facing executive summary.

**Branching:** one branch per step (Step 1 = `phase3-followup-warmstart`), cut from `main`. Merge
each step's branch into `main` once its own tests pass and the independent code-reviewer agent
finds no unresolved issues — do not wait for a separate go-ahead to merge each step (the user has
pre-approved this plan); **never push to remote without asking**, per standing repo convention.

**Automation intent:** this plan is written so a session can execute Steps 1–6 end to end with
minimal back-and-forth. Steps 2 and 3 are mutually exclusive and their selection is a mechanical
decision rule evaluated from Step 1's own results (see Step 1's "Decision rule," not a judgment
call) — do not stop to ask the user which branch to take; follow the rule and state which branch
you took and why. Standing repo conventions still apply within each step (TDD for new functions,
`superpowers:code-reviewer` dispatch before each commit, full test suite run, `DECISIONS.md`/
`ROADMAP.md`/`PROJECT_STATUS.qmd` updated with real numbers, no push without confirmation). Only
pause for the user on a genuine blocker (missing `PDAC_DATA_ROOT`, a test failure you can't resolve,
an ambiguity this plan doesn't resolve) — not to re-confirm something this plan already decided.

---

## Step 1 — Multistart + warm-start re-check (always run)

**Branch:** `phase3-followup-warmstart`, cut from `main`.

**Goal:** determine whether Phase 3's K=2/K=4 underperformance
(`results/benchmark_sim/run_k_parsimony_curve.R`, `DECISIONS.md` 2026-07-13) reflects a genuine
capacity ceiling or the CAVI factor-collapse artifact from Phase 2.

**Implementation:**
1. Warm-start K=2,3,4,5 fits using `fit_cox_on_yf(..., init_method="custom", EL_init=, EF_init=)`
   (parameter already supported), seeded from the already-converged K=7 fit's top-K columns ranked
   by proportion-of-variance-explained (not just |β| — K serves both reconstruction and prediction).
2. Best-ELBO multistart at the same K values: extend `code/fit_modular_multistart.R` (currently
   LB-only) to support YFB. An earlier attempt to test this hit an unrelated setup error (a
   `!is.null(real_Y)` assertion from `fit_cox_on_yf.R`'s standalone-runner block, likely a missing
   `tryCatch(source(...))` wrap) — root-cause and fix that, don't just work around it.
3. Re-run external validation (same 5 held-out cohorts) for both the warm-started and multistart
   fits at each K, alongside the original fresh-SVD numbers for direct comparison.
4. TDD any new reusable code (e.g., a PVE-ranked column-extraction helper); dispatch code-reviewer;
   run full test suite; update `DECISIONS.md` with the full comparison table regardless of outcome.

**Decision rule (mechanical, not a judgment call):** for each K in {2,3,4,5}, compute the best
external mean C-index achieved by either warm-start or multistart at that K.
- If **any** K < 7 achieves external mean C-index within 1 SE of K=7's 0.627 (margin 0.607, per
  the existing 1-SE convention) → outcome is **OPTIMIZATION-LIMITED** → proceed to **Step 2**,
  skip Step 3.
- If **no** K < 7 reaches that margin even with both improved-optimization strategies → outcome is
  **CAPACITY-LIMITED** → proceed to **Step 3**, skip Step 2.
- Record which outcome occurred and the exact numbers in `DECISIONS.md` before proceeding — the
  next step's prompt depends on reading this record, not on re-deriving it.

---

## Step 2 — Deflation-init fix for CAVI factor collapse (only if OPTIMIZATION-LIMITED)

**Branch:** `cavi-deflation-init`, cut from `main` (after Step 1 merges).

**Goal:** make the optimization improvement from Step 1 permanent and general, not a one-off
warm-start hack — fixes the root cause diagnosed in Phase 2 (near-symmetric, tied-amplitude
factors causing SVD-init CAVI to hit a degenerate fixed point).

**Implementation:**
1. Add a deflation-style initialization option to `fit_cox_on_yf`/`fit_supervised_mf_modular`:
   seed factor 1 from a rank-1 SVD of Y, subtract its contribution, take a rank-1 SVD of the
   residual for factor 2, and so on through K factors — then run the existing joint CAVI updates
   unchanged.
2. Verify this fixes the `sparse_synthetic` collapse scenario in
   `results/multi_cohort_sim/run_survival_strength_sweep.R` (Phase 2's `diagnose_factor_collapse.R`
   is the fastest way to check: dead-factor count should drop from the current ~1-2/4 for YFB,
   ~2.4/4 average for LB, toward 0).
3. Confirm no regression on the full test suite (276+ tests) and on the existing real-data D4
   benchmark (`run_desurv_comparison.R`).
4. TDD, code-reviewer, `DECISIONS.md`/`ROADMAP.md` update with before/after numbers.

---

## Step 3 — Joint (K, α, penalty) Bayesian optimization (only if CAPACITY-LIMITED)

**Branch:** `joint-k-alpha-bayesopt`, cut from `main` (after Step 1 merges).

**Goal:** if better optimization alone doesn't unlock a smaller K, test DeSurv's own approach —
jointly tuning K, α, and an elastic-net-style penalty on β via Bayesian optimization, rather than
our current fixed-α, no-penalty, K-only CV.

**Implementation:** follow `docs/plans/joint_k_alpha_bayesopt_plan_07_12_2026.md` (already drafted,
not yet implemented). Present the concrete implementation plan in-session before writing code (this
one genuine judgment call — an untried, larger piece of new methodology — is worth a quick sanity
check even under the "automate the sequence" intent above, since scope/complexity here could
reasonably shift once the code is in front of you). TDD, code-reviewer, `DECISIONS.md`/`ROADMAP.md`
update.

---

## Step 4 — Final validation (always run)

**Branch:** whichever branch Step 2 or Step 3 produced; no new branch needed.

**Goal:** produce the conclusive, validated K-vs-external-performance answer.

**Implementation:** re-run `results/benchmark_sim/run_k_parsimony_curve.R` (or an updated version)
using whichever fitting procedure emerged from Step 2 or Step 3, across the same K grid and 5
held-out cohorts. Update `DECISIONS.md`/`ROADMAP.md`/`CLAUDE.md`'s "Current model status" line with
the conclusive result — either a smaller K is adopted as the new recommended configuration (update
`config/globals.yml`'s `k_merged_yfb_desurv` accordingly), or K=7 is reconfirmed as necessary under
two independent optimization strategies (multistart/warm-start and, if run, joint BO) — a
stronger, doubly-verified version of Phase 3's original finding either way.

---

## Step 5 — Pathway enrichment (always run, after Step 4)

**Branch:** `pathway-enrichment-plan` (already exists, commit `d7aee46`, plan at
`docs/plans/pathway_enrichment_plan.md`).

**Goal:** biological interpretation of the final recommended config's active factors.

**Implementation:** before starting, re-verify the plan's target factors (originally Factor 7
adverse / Factor 3 protective) are still correct for whatever config Step 4 concluded with — the K
and active-factor identities may have changed. If they changed, update the plan's targets rather
than proceeding on stale factor numbers. Then execute the plan's Step 1 onward as written.

---

## Step 6 — Executive summary / progress report (always run, last)

**Goal:** a research-team-facing summary of everything since the last lab meeting, ready to present.

**Implementation:** synthesize Phase 1 (objective normalization correction), Phase 2 (joint-vs-
two-step validation, including the CAVI collapse finding), Phase 3 (K-parsimony curve), Steps 1–4
above (the follow-up validation and whichever fix path was taken), and Step 5 (pathway enrichment)
if complete. Follow the project's documentation-audience convention (`CLAUDE.md`'s "Documentation
audience" note; biostatistician collaborators reading cold, model/prior/metric/result framing, not
internal phase-label/session narration — avoid "Phase 1/2/3" and step numbers as the primary
structure, use them only as lookup keys if needed). Check `ROADMAP.md`'s flagged to-do item and the
`project_executive_summary_todo` memory entry for full context. Present a draft for review before
finalizing — this is the one deliverable meant for a human audience outside this repo, worth a
real look before it goes out.
