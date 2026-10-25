# SSBMF — Item 2: Progress Report & Executive Summary (implementation plan)

**Date:** 2026-07-15
**Purpose:** Execute Item 2 of the consolidated remaining-work plan
(`docs/plans/ssbmf_progress_consolidation_and_remaining_work_07_14_2026.md`) — a two-tier written
synthesis of everything accomplished since the 6/18/2026 lab meeting. Written to be executed end to
end by a fresh session. **Items 1, 3, 4, 5 of the consolidation plan are already complete
(2026-07-15); only Items 2 and 6 remain. Item 6 is gated behind Item 2.**

**Read first (do not skip):**
- `docs/plans/ssbmf_progress_consolidation_and_remaining_work_07_14_2026.md` — the authoritative plan;
  its Part A is the orientation summary this report synthesizes, its Item 2 is the spec this expands.
- `DECISIONS.md` entries dated **2026-07-12, 2026-07-13, 2026-07-15** — the substance of the arc.
- `docs/reports/pathway_enrichment_report_07_15_26.{qmd,pdf}` — Item 1's result (an input here).
- Memory: `project_executive_summary_todo`, `project_html_visual_summary_todo`.

---

## Current state (all merged to LOCAL `main`; nothing pushed)

- **Recommended configuration D4:** YFB (η = (YF)β), DeSurv-aligned gene selection (combined_rank,
  top-3000 per cohort before per-platform z-standardization, 2064 genes after intersection), K=7, no
  cohort indicator. **Mean external C = 0.627 across 5 held-out PDAC cohorts, K_eff = 2.**
- **Tests:** `Rscript tests/run_tests.R` → **374/374**; `Rscript tests/test_real_data_loading.R` →
  **88/88** (local-only, needs `PDAC_DATA_ROOT`).
- **Item 1 (pathway enrichment, done):** Program 7 = adverse = basal-like/squamous/MET-EGFR;
  Program 3 = protective = classical/differentiated-epithelial — confirmed by 5 independent methods
  (gene-set enrichment, PurIST subtype concordance, external-cohort survival HR, SBMF-vs-DeSurv gene
  overlap, ORA cross-check). Report: `docs/reports/pathway_enrichment_report_07_15_26`.
- **Item 3 (stratified Cox baseline hazard, done):** optional `strata_id=` (stratified partial
  likelihood, study-specific baseline, no parametric h0) in both fits, **off by default**,
  **performance-neutral** on D4 (mean external C 0.6267 vs 0.6263). DECISIONS.md 2026-07-15.
- **Item 4 (manuscript scaffold, done & pushed):** separate repo `~/GithubProjects/SSBMF-paper`
  (Quarto scaffold + drafted intro + reference library). Do not touch unless asked.
- **~47 commits on local `main` are NOT pushed** (origin stuck at `33fee99`, pre-6/18-followup).

---

## Deliverables (two tiers — related, not duplicate)

### A. Long-form comprehensive progress report → `docs/reports/`
- **Format:** `.qmd` rendering to **PDF + HTML**. Dated filename following the existing convention,
  e.g. `docs/reports/ssbmf_progress_report_MM_DD_26.qmd` (use the implementing session's date).
- **Audience:** biostatistician collaborators reading cold. Per `CLAUDE.md`'s documentation-audience
  convention: frame everything by **model / prior / training set / metric / result** — **NOT** by
  internal "Phase N" / "Step N" narration (those labels may appear only as parenthetical lookup keys,
  never as the structure of the prose).
- **This is the primary deliverable.** Present a draft for review before finalizing/rendering final.

### B. Short meeting-ready executive summary → `docs/progress_report/`
- **Format:** `.qmd` → PDF/HTML, ~1–2 pages, distilled from Deliverable A. Match the existing folder
  naming convention (`SSBMF_Status_Update_MM_DD_26.{qmd,pdf}`), or
  `SSBMF_Executive_Summary_MM_DD_26` if clearer — confirm the name with the user.
- **Purpose:** prep for a follow-up lab meeting the user will schedule (Item 5 "help me prep" call).
  Tight, presentable, headline-first.

### C. (Later, gated) Item 6 — visual HTML executive summary / demo
- **Do NOT start until Deliverables A and B are done and reviewed.** Noted here only so it is not
  forgotten; it draws on A/B as input. See the consolidation plan's Item 6 and the
  `project_html_visual_summary_todo` memory.

---

## Pre-work — final whole-branch review (do BEFORE writing Deliverable A)

The consolidation plan (Item 5) and the original plan's Global Verification call for **one final
whole-branch review** of the entire arc, beyond the per-item reviews already done. Dispatch an
independent `superpowers:code-reviewer` over the whole arc (post-lab-meeting Phases 1–3 +
K-parsimony follow-up Steps 1–4 + Item 1 pathway enrichment + Item 3 strata) against the
**net-benefit gate**: (1) all tests green (374/374 core, 88/88 real-data), and (2) external mean
C-index not worse than the 0.636 pre-work baseline (it is 0.627 — a documented, understood correction
from a preprocessing-bug fix, not a regression; confirm the reviewer agrees this is accounted for).
The report should rest on a verified arc. Record the review outcome (and any findings addressed) in
`DECISIONS.md`.

---

## Prior-report & figure staleness audit (REQUIRED — do while gathering material)

The report must reflect the **current** state of the method with **all relevant results**, not
reproduce figures/numbers from prior reports that recent corrections have superseded. Several prior
reports predate the 2026-07-12 corrections (external mean C 0.636 → 0.627; the erroneous
`boost_beta=TRUE` K_eff 4 → 2 figure; the K-parsimony conclusions) and the 2026-06-16 program-
direction correction. **Audit every prior report and every figure it would draw on; carry forward
only what is still valid, and refresh/replace what is stale — do not silently inherit outdated
figures.**

Concrete inventory (verified 2026-07-15 by grepping the `.qmd` sources):

**Carry the superseded `0.636` external C figure (current is 0.627) — treat as STALE sources:**
- `docs/reports/desurv_alignment_report_05_27_26.qmd` — the "DeSurv record" report (named in
  `CLAUDE.md`); its external-C figures are pre-correction.
- `docs/reports/multicohort_sim_proposal_06_14_26.qmd`
- `docs/reports/desurv_factor_diagnostics_05_27_26.qmd`
- `docs/reports/pathway_enrichment_overview.qmd`
- `docs/progress_report/SSBMF_Status_Update_4_29_26.qmd`
- `docs/progress_report/SSBMF_Status_Update_5_27_26.qmd`
- `docs/progress_report/SSBMF_Status_Update_5_28_26.qmd`

**Mention `K_eff=4` / "4 active" (current is K_eff=2) — VERIFY each (some may be false positives):**
`docs/reports/desurv_factor_diagnostics_05_27_26.qmd`, `docs/reports/ssbmf_summary_report_04_29_26.qmd`,
`docs/progress_report/SSBMF_Status_Update_5_28_26.qmd`.

**Already current (0.627) — safe to reuse figures/framing from:**
`docs/reports/joint_vs_twostep_sweep_07_12_2026.qmd` (Phase 2 value-add sweep),
`docs/reports/pathway_enrichment_report_07_15_26.qmd` (Item 1 biology).

**Also check (beyond text numbers):** any figure image files or figure-generating scripts that
produced plots (C-index bars, KM curves, factor/loading heatmaps, K-curves) from **pre-correction
fits** — regenerate from current fits or the current stored results
(`results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_results.csv`, current;
`docs/reports/pathway_enrichment_report_07_15_26` assets, current) rather than embedding an old image.
When in doubt about a figure's provenance, regenerate it. The authoritative current numbers live in
`DECISIONS.md` (2026-07-12/13/15), `CLAUDE.md` (Current model status), and the current results CSV.

The report is meant to be a **comprehensive review of the current state, results, and potential value
of the method** — so the audit is not just error-correction; it is the mechanism for pulling every
still-valid result forward into one current, coherent picture.

## Content spec for Deliverable A (what the report must cover)

**Honest headline (the through-line):** the recommended configuration itself did not change since the
6/18 deck, but its numbers were corrected (mean external C 0.636 → 0.627, attributable to a
train/test preprocessing fix) and its K=7 choice was rigorously stress-tested from multiple
independent angles — all converging on the same answer: **the model uses 2 of its 7 factors for
survival prediction (K_eff = 2), and this is robust and real, not a tuning artifact.** K=7 is kept
for one-step reproducibility; K=4/K=5 are validated statistically-equivalent (but two-step
warm-start) alternatives.

Sections to include (adapt structure to the audience, not to this list literally):
1. **Model overview** — the YFB parameterization (η = (YF)β), how it differs from the LB variant
   (η = Lβ), and why YFB is recommended (closes the train/test formula mismatch; L pure-genomics, β
   pure-survival).
2. **Objective on a principled footing** — the genomics/survival normalization work and retirement of
   the redundant λ (α is the sole mixing weight, DeSurv-aligned); the honest finding that this gave no
   YFB performance benefit and the corrected `boost_beta=FALSE` default; the train/test preprocessing
   fix that moved 0.636 → 0.627.
3. **Value-add over a two-step baseline** — joint YFB vs. PCA+Cox and EBMF+Cox across a
   survival-signal sweep: advantage emerges and grows with signal in 3/4 DGP scenarios; the 4th
   reversed and was root-caused to a genuine, understood CAVI factor-collapse vulnerability (not a
   survival-coupling bug).
4. **Factor parsimony (K=7 vs DeSurv's 3; K_eff=2)** — total-K vs active-K_eff distinction; the
   three independent stress tests (warm-start, deflation-init, joint (K,α) Bayesian optimization) that
   all confirm K_eff=2 and a statistical tie across K∈{4,5,7,8}; why K=7 is kept as default.
5. **External validation** — mean C=0.627 across the 5 held-out PDAC cohorts (RNA-seq, microarray,
   proteomics); per-cohort numbers; the essential role of per-platform z-standardization (10/12
   non-per-platform configs collapse to β=0).
6. **Biology of the two active programs (from Item 1)** — Program 7 adverse/basal-like/squamous/
   MET-EGFR; Program 3 protective/classical; confirmed by 5 independent methods; the marginal-vs-joint
   β sign subtlety (`feedback-yfb-km-sign-suppression` — label factors empirically).
7. **Honest open items / shortfalls** — the YFB `point_normal` K-CV C=0.5 bug (open, ROADMAP.md); the
   stratified-baseline-hazard neutral result (Item 3, available but off by default); the K=7-vs-K=4
   reproducibility-vs-parsimony judgment call.
8. **Potential value / positioning (make this explicit)** — where the method stands and why it
   matters: SBMF as the Bayesian counterpart to DeSurv (same design/data/scoring), differentiated by
   empirical-Bayes prior shrinkage, posterior uncertainty, error-in-variables E[L²], and
   heteroscedastic τ; the joint model recovering prognostic structure a two-step misses (only when
   survival carries it); interpretable, externally-validated gene programs with a clear biological
   story. Frame the head-to-head with DeSurv as a *pending* matched comparison, not a superiority
   claim. This is the "so what" of the report — state it plainly.

Pull exact numbers from `DECISIONS.md` (2026-07-12/13/15), `CLAUDE.md` (Current model status), and
`ROADMAP.md` — do not re-derive or re-run fits for the report; cite the stored results.

---

## Standing conventions

- Independent review before any merge; `Rscript tests/run_tests.R` must stay **374/374** if any code
  is touched (the report work should be prose/`.qmd` only — no model-code changes expected).
- Do report/plan work on a branch off `main`; merge once the draft is approved.
- **Never push to remote without explicit user confirmation.**
- Update `DECISIONS.md`/`ROADMAP.md`/relevant memory when Item 2 completes.

## Push gate (the ~47 unpushed commits)

The user wants to push the backlog "soon." Agreed sequence: **after the whole-branch review confirms
the arc is clean**, propose the push; execute only on the user's explicit "push now" (a plain
`git push origin main`). Flag the current unpushed count when raising it.

## Done criteria

- Deliverable A (`.qmd` + rendered PDF + HTML) in `docs/reports/`, reviewed and approved by the user.
- Deliverable B (`.qmd` + rendered) in `docs/progress_report/`, distilled from A.
- Whole-branch review completed and its outcome recorded in `DECISIONS.md`.
- `DECISIONS.md`/`ROADMAP.md`/memory updated; tests still 374/374.
- Item 6 explicitly left for a subsequent session (gated).
