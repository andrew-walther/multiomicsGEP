# SSBMF — Post-Lab-Meeting Action Plan

**Date:** 2026-07-08
**Source feedback:** `docs/plans/rashid_lab_meeting_notes_06_18_2026.md` (Rashid lab, 6/18/2026)
**Branching:** one branch per phase (Phase 1 = `objective-normalization`); see integration strategy
below — **do not touch `main` until each branch's net-benefit gate is confirmed.**

---

## Context

On 6/18/2026 the Rashid lab reviewed the Supervised Bayesian Matrix Factorization (SSBMF) results.
The feedback raised concerns about the objective specification, preprocessing, factor count, and —
most importantly — asked for evidence that the *joint* supervised model adds value over the standard
*2-step* (unsupervised factorization → Cox) procedure. This plan triages that feedback and sequences
the work so that (a) the objective is put on a principled, interpretable footing first, (b) the
"value-add over 2-step" story is then demonstrated on top of a correct objective, and (c) the
manuscript is scaffolded once the model work stabilizes.

The current recommended configuration is **D4** (YFB, DeSurv-aligned preprocessing, K=7, mean external
C-index = 0.636 across 5 held-out PDAC cohorts). This is the baseline every change is measured against.

### Triage summary

| Feedback item | Code reality (verified 2026-07-08) | Disposition |
|---|---|---|
| α default = 1/2 | Already default (`config/globals.yml:140`) | ✓ confirm & document |
| h₀(t) specification | Breslow Cox **partial** likelihood — h₀ cancels (`fit_modular.R:98-101`) | ✓ document; the lab's own question answers "yes" |
| Normalize by n and p | **Absent** — genomics ELBO is O(n·p), survival is O(n) (`update_tau.R:169`, `compute_elbo.R:150`) | **Phase 1** (lead) |
| λ scale / (1−λ)·L_gen | Our λ is an ad-hoc survival multiplier over (0,∞), **not** DeSurv's β-penalty; overlaps α | **Phase 1** — retire, use α |
| Z-std then platform-correct | Per-platform z-std happens before merge (`preprocess_desurv.R:314-317`) | ✓ mostly in place |
| Study-specific baseline hazard `strata(study)` | **Absent** — single h₀(t); cohort_id affects genomics only | **Phase 4** (extension) |
| External datasets normalized? / rank transform | **Train/test mismatch**: external is rank-transformed & not per-platform z-std'd; training is the reverse (`run_desurv_comparison.R:290-296`) | **Phase 1** (fold-in) |
| Why 7 factors (DeSurv=3)? K shrink if survival off? | K=7 by CV+1-SE+floor≥3; K_eff=2 active; no survival-off diagnostic | **Phase 3** |
| Joint vs 2-step equivalence when survival has no signal | EBMF 2-step exists in `multi_cohort_sim`; no PCA+Cox baseline; no survival-strength sweep in main narrative | **Phase 2** (the value-add story) |
| Factor interpretation (basal-like adverse / protective) + GO | Pathway-enrichment plan already written, not implemented (branch `pathway-enrichment-plan`) | **Phase 5** |
| Paper repo / Quarto draft | Not started (`github.com/andrew-walther/SSBMF-paper`) | **Phase 6** |
| Post-commit Codex review hook; RAG lit DB | Not started | **Backlog** |

### Why α, not a second λ (DeSurv-grounded)

DeSurv's objective is `(1−α)·L_NMF(W,H) − α·L_Cox(W,β)`: **α is the single mixing weight**
(α=0 ⇒ plain NMF; PDAC selected α=0.334), and DeSurv's **λ is an elastic-net penalty on β** — a
regularizer, a completely different role. SBMF's α already plays DeSurv's α role, but SBMF's current
λ is an ad-hoc survival-scale multiplier that duplicates α. SBMF gets its shrinkage from
empirical-Bayes priors, so it needs no β-penalty λ. **Decision:** normalize the two terms (÷np and ÷n)
so α is interpretable, use α∈(0,1) as the sole mixing knob, and retire the current λ. This also makes
SBMF's α directly comparable to DeSurv's α in the head-to-head.

---

## Branch & integration strategy

**One branch per phase/workstream** (not a single mega-branch) — so each merge gate is small,
independently reviewable, and `main` only ever gains work that has been confirmed a net benefit.

| Phase | Branch | Cut from | Depends on |
|---|---|---|---|
| 1 — Objective normalization | `objective-normalization` | current `main` | — |
| 2 — Joint vs 2-step value-add | `validation-two-step` | `main` after Phase 1 merges (else off `objective-normalization`) | Phase 1 (scaled α) |
| 3 — Factor parsimony / K | `factor-parsimony-k` | `main` after Phase 1 merges | Phase 1 |
| 4 — Study-specific baseline hazard | `cohort-strata-hazard` | `main` after Phase 1 merges | Phase 1 |
| 5 — Factor interpretation / enrichment | reuse existing `pathway-enrichment-plan` (rebase on updated `main`) | `main` after Phases 1–3 | Phases 1–3 |
| 6 — Manuscript | separate repo `SSBMF-paper` (not a branch here) | — | model work stable |

- **Sequencing note:** Phase 1 is the shared dependency for 2–4. If Phase 1 is merged before starting
  the next phase, cut the new branch off updated `main`; if Phase 1 is still under review, cut off the
  `objective-normalization` branch and rebase once it lands. Do not run Phases 2–4 on the *old* objective.
- **Checkpoint-commit within each branch** (between sub-steps and phases) so any regression can be
  `git restore`d without losing completed work.
- **`main` is not touched** until a branch clears its net-benefit gate.
- **Net-benefit gate (per-branch merge criterion):** merge only when, on the corrected pipeline,
  (1) all tests pass (`246/246` core, `77/77` real-data local), and (2) external mean C-index is **not
  worse** than the 0.636 D4 baseline — ideally improved — with any trade-offs documented in
  `DECISIONS.md`. If the gate is not met, keep the branch open and report findings; do not merge.
- Record the pre-work `main` commit hash as a backtrack checkpoint in `MEMORY.md`/`PROJECT_STATUS.md`.

---

## Independent review protocol (every checkpoint)

After each sub-step / phase — **before** committing — spin up an **independent reviewer agent**
(`superpowers:code-reviewer`, or the `requesting-code-review` skill) with fresh context. The reviewer
does not inherit the implementer's assumptions; it is given the phase's stated goal from this plan, the
diff since the last checkpoint, and the repo coding standards (`CLAUDE.md`), and is asked to check for:

1. **Correctness / bugs** — especially that normalization propagates into the CAVI *update
   coefficients*, not just the ELBO monitor (the most likely place for a silent error).
2. **Consistency** — with the plan's stated intent, existing patterns, and the math in
   `derivations/`; no drift between the ELBO objective and what the updates actually optimize.
3. **Test integrity** — assertions verify *correct behavior*, not just that they pass; any test whose
   expected value changed must be justified, not merely re-baselined to green.
4. **Simplification opportunities** — dead parameters (e.g. retired λ), redundant code paths, or
   places where the change could be smaller/clearer.

The reviewer reports findings; the implementer addresses or explicitly defers each (with rationale)
**before** the checkpoint commit. Surface disagreements rather than averaging them. Record any deferred
items in the phase's commit message or `ROADMAP.md`.

---

## Phase 1 (LEAD) — Principled, interpretable objective

**Goal:** Put the two likelihood terms on a comparable per-element scale (DeSurv's `(1−α)/α`
structure), remove the redundant λ, and fix the train/test preprocessing mismatch.

**What we hope to achieve — success criteria (define before implementing):**
- **Interpretable α.** Report the ratio of the two ELBO-term magnitudes at α=0.5, before vs after
  normalization. *Success:* the ratio moves from ≈O(p) apart to ≈O(1) — the terms are now genuinely
  comparable, so α=0.5 means "balanced."
- **α behaves sensibly (monotone control).** At α=0 the fit is survival-blind and numerically matches an
  unsupervised factorization (loadings correlation ≈1, β≈0); at α=1 the factorization is driven by
  survival. *Success:* these two limits reproduce cleanly — the primary evidence the normalization
  reached the update coefficients, not just the ELBO monitor.
- **Fewer, clearer parameters.** λ retired with no test regression; the parameter set is α (+ K, priors).
- **No performance loss from the correction.** *Success gate:* `246/246` tests pass, and the corrected
  train/test pipeline gives external mean C-index **≥ 0.636** (report the delta honestly either way; a
  drop is still informative but blocks merge).

### 1a. Normalize the objective by n and p
- **Genomics term** ÷ (n·p); **survival term** ÷ n.
- Critical: normalization must propagate into the **CAVI update coefficients**, not just the ELBO
  monitor. The `A`/`B` precision-and-signal terms in `update_L.R:146-160,191-192`,
  `update_beta.R:140-148`, `update_F.R:114-126`, and `update_tau.R` carry the same α weighting and must
  be rescaled consistently, or the posterior will not optimize the normalized objective.
- Update ELBO assembly at `fit_modular.R:613-615` and `fit_cox_on_yf.R:519-521`.
- Files: `code/update_tau.R`, `code/update_L.R`, `code/update_beta.R`, `code/update_F.R`,
  `code/compute_elbo.R`, `code/fit_modular.R`, `code/fit_cox_on_yf.R`.
- **Verify:** α=0 must reduce to survival-blind factorization; α=1 drives β purely from survival.
  `Rscript tests/run_tests.R` (update tests that assert on unnormalized ELBO magnitudes). **Independent review** (spin up a reviewer agent per the protocol above), then **commit.**

### 1b. Reconcile α / λ → single α (DeSurv-aligned)
- Retire the `lambda` survival-scale multiplier (`update_L.R:146`, `config/globals.yml:19-22`); α is
  the sole mixing knob over (0,1), selected by CV + 1-SE (`code/select_alpha_cv.R`).
- Document that DeSurv's λ (elastic-net β penalty) is a different role filled by SBMF's EB priors.
- **Decision point:** fully remove λ from signatures (recommended for clarity) vs freeze at 1.0 with a
  deprecation note. Remove once tests confirm no regression. **Independent review** (spin up a reviewer agent per the protocol above), then **commit.**

### 1c. Fix train/test preprocessing mismatch
- Make external-validation preprocessing identical to training: external cohorts are currently
  rank-transformed and *not* per-platform z-std'd, while D4 training is the reverse
  (`code/preprocess_desurv.R:90-115`, `results/benchmark_sim/run_desurv_comparison.R:290-296`).
- **Verify:** re-run `results/benchmark_sim/run_desurv_comparison.R`; report the 5-cohort external
  C-index delta vs 0.636 honestly (up, down, or flat) — a correctness fix regardless of direction.
  **Independent review** (spin up a reviewer agent per the protocol above), then **commit.**

---

## Phase 2 — Joint vs 2-step value-add (the "sell the method" story)

**Goal:** Show in simulation that when survival carries signal the joint model beats the 2-step, and
when survival has no signal they are equivalent. Depends on Phase 1 (α must be correctly scaled for
α→0 to mean "survival off").

**What we hope to achieve — success criteria (define before implementing):**
- **Equivalence at zero signal (the honesty check).** At survival effect = 0, joint (tuned α) vs both
  2-step baselines (PCA+Cox, EBMF+Cox) agree in C-index (overlapping CIs / |ΔC| < ~0.01) and in factor
  recovery. *Success:* no spurious advantage when there is nothing to gain — this is what makes the
  positive result credible.
- **Separation that scales with signal (the value-add).** As the prognostic effect grows 0 → large,
  joint C-index (and β / factor-recovery accuracy) exceeds the 2-step, and the gap widens monotonically.
  *Success:* a clearly separated curve with a one-sentence takeaway defensible to the lab —
  "the joint model recovers prognostic structure the 2-step misses, and only when survival carries it."
- **α=0 reproduces the unsupervised 2-step exactly** (internal sanity control).
- *Deliverable:* a short Quarto report in `docs/reports/` with the sweep figures.

- **Survival-strength sweep:** extend the DGP (`results/multi_cohort_sim/generate_multicohort_data.R`,
  `config/globals.yml` synthetic block) to vary the prognostic effect 0 → large on 1–2 factors,
  holding genomics fixed.
- **Add a PCA+Cox 2-step baseline** alongside the existing EBMF 2-step
  (`results/multi_cohort_sim/run_multicohort_sim.R:163-191`): `fit_pca_cox()` = `prcomp` →
  `coxph(Surv ~ scores)`.
- **Equivalence at signal = 0:** joint (α tuned) ≈ 2-step in C-index and factor recovery; divergence =
  objective bug. **Separation at signal > 0:** joint exceeds 2-step, gap grows with signal.
- **α=0 internal control:** joint at α=0 reproduces the unsupervised 2-step exactly.
- **Deliverable:** short Quarto diagnostic report in `docs/reports/` with sweep curves (use the
  `dataviz` skill for figures). **Independent review** (spin up a reviewer agent per the protocol above), then **commit.**

---

## Phase 3 — Factor parsimony & why K=7 vs DeSurv's 3

**Goal:** Explain/reduce the factor count and answer "does K shrink when survival is off?"

**Framing — total K vs active K_eff (critical to keep separate):** DeSurv retains **3 total** factors
with survival concentrated in **one**; SBMF D4 retains **7 total** with **K_eff=2** active (|β̂|>1e-3).
The 7-vs-3 gap is in *total* factors, not survival-carrying factors — on the latter axis both methods
agree survival lives in ~1–2 factors. **Leading hypothesis for the gap:** the unnormalized objective
under-weights survival, so SBMF behaves like weakly-supervised NMF — and DeSurv's paper reports standard
NMF needs k≈7 to match DeSurv's k=3. If so, Phase 1 normalization alone should pull K down toward 3.
This is a hypothesis to test, *not* established — which is why Phase 3 runs only on the corrected objective.

**What we hope to achieve — success criteria (define before implementing):**
- **A principled answer to "why 7 vs DeSurv's 3."** Produce a K_eff-vs-α curve. *Success:* either (a) on
  the normalized objective the CV-selected total K drops from 7 toward DeSurv's 3 with external C-index
  maintained (≥0.636) — the desired outcome — or (b) a documented, defensible reason K stays higher
  (e.g. genomic variance structure requires it), so the number is no longer an open question.
- **Report total K and K_eff separately at every K.** *Success:* we can see whether the extra total
  factors earn their keep (improve held-out C) or are reconstruction-only passengers.
- **Survival-off K behavior characterized.** Report selected K at α=0 vs α=tuned. *Success:* a clear
  statement of whether/how much supervision concentrates signal into fewer active factors (K_eff should
  rise above the current 2 if supervision is working as intended).
- *Success gate:* K, K_eff, and external C-index reported side-by-side, before vs after normalization.

- Diagnostic refitting across α ∈ {0, 0.5, tuned}, reporting K_eff (active by |β̂|>1e-3 or PVE>1%)
  as a function of α (`code/select_K.R:69-98`).
- Hypothesis: once the objective is normalized, survival supervision concentrates prognostic signal
  into fewer factors (DeSurv: k=3, signal in one factor D1).
- Re-run merged K-CV (`results/benchmark_sim/run_merged_kcv.R`) on the normalized objective; keep the
  1-SE rule but **relax the hard K≥3 floor — test K=2 explicitly** (a protective + adverse pair may
  suffice). The floor was imposed to match DeSurv's k=3–4 back when 1-SE selected K=2; if CV on the
  corrected objective prefers K=2, do not force a third factor. Caveat: K=2 *total* forces both factors
  to carry reconstruction *and* survival (leaner than DeSurv's 3-total/1-active) — decide empirically.
- **Verify:** report K, K_eff, external C-index before/after normalization side by side. **Independent review** (spin up a reviewer agent per the protocol above), then **commit.**

---

## Phase 4 — Study-specific baseline hazard (modeling extension)

**Goal:** Let baseline survival differ by cohort (`+ strata(study)`), since cohort_id currently only
absorbs *genomic* platform offsets, not differing baseline risk.

**What we hope to achieve — success criteria (define before implementing):**
- **Correct handling of cross-cohort survival heterogeneity.** The stratified Breslow partial
  likelihood forms risk sets *within* each study, so a cohort with systematically shorter survival no
  longer distorts β. *Success:* new unit tests for the stratified likelihood pass, and the fit is
  unchanged when all samples share one stratum (reduces to the current model — a correctness anchor).
- **Measurable effect on generalization.** *Success:* external mean C-index with cohort strata is
  **≥** without (or documented as neutral); β estimates are stable/more defensible across cohorts.
  If strata do not help, that is a reportable finding, not a failure.

- Extend the Breslow risk-set computation to be **stratified by cohort** (separate risk sets per study)
  in `fit_cox_on_yf.R:74-102` / `fit_modular.R:71-102`. Baseline still cancels *within* each stratum —
  standard stratified Cox partial likelihood, no parametric h₀ introduced.
- **Higher effort / decision point:** touches the CAVI survival derivation and needs new unit tests.
  Sequenced after Phases 1–3 because benefit is only measurable on a corrected objective.
- **Verify:** new tests for the stratified partial likelihood; external C-index with vs without cohort
  strata. **Independent review** (spin up a reviewer agent per the protocol above), then **commit.**

---

## Phase 5 — Factor interpretation & pathway enrichment

**Goal:** Characterize the prognostic factors biologically (adverse = basal-like; protective per
DeSurv), answering the GO/gene-ID feedback.

**What we hope to achieve — success criteria (define before implementing):**
- **Biological identity for each active factor.** *Success:* the adverse factor is enriched for
  basal-like / known-aggressive PDAC programs and the protective factor is concordant with DeSurv's,
  with GO/pathway terms reported and factor signs labeled *empirically* (not assumed — the joint β sign
  can invert vs marginal KM direction; memory `feedback-yfb-km-sign-suppression`).
- **Concordance with a published comparator.** *Success:* factor–gene loadings correlate with DeSurv's
  reported programs, giving external corroboration of what the model found.

- Execute the already-written plan on branch `pathway-enrichment-plan`
  (`docs/plans/pathway_enrichment_plan.md`); gene lists in
  `presentation/.../assets/active_factor_genes.csv`.
- Mind the sign/suppression gotcha (memory `feedback-yfb-km-sign-suppression`): joint β sign can invert
  vs marginal KM direction — label factors empirically.
- Cross-reference DeSurv's gene-program correlations for concordance.

---

## Phase 6 — Manuscript track (deferred; after model work stabilizes)

**Goal:** Stand up paper infrastructure and begin the introduction.

**What we hope to achieve — success criteria (define before implementing):**
- **A working manuscript scaffold.** *Success:* `SSBMF-paper` repo renders a Quarto skeleton;
  benchmark scripts and figure generators are organized into a `paper/` layout so results accrue in
  place rather than being reconstructed later.
- **A cited, defensible introduction draft.** *Success:* an intro draft with a reference list
  (supervised/guided MF, prognostic gene programs, PDAC subtyping, Bayesian EBMF) and an explicit
  positioning of SBMF as the Bayesian counterpart to DeSurv — framed as a *pending* head-to-head, not a
  superiority claim.

- Scaffold `github.com/andrew-walther/SSBMF-paper`: Quarto manuscript skeleton; organize benchmark
  scripts and figure generators into a `paper/` layout.
- **Background research for the intro:** prompt the user for their existing reference list, then run
  deep literature research (`deep-research` skill) to fill gaps — supervised/guided matrix
  factorization, prognostic gene programs, PDAC subtyping (PurIST, DeCAF, basal/classical), Bayesian
  EBMF. DeSurv is the primary comparator (memory `reference-desurv-paper`).
- Positioning: SBMF as the Bayesian counterpart to DeSurv (same design/data/scoring), differentiated by
  EB-prior shrinkage, posterior uncertainty, error-in-variables E[L²], heteroscedastic τ. A matched
  same-protocol head-to-head is a *pending* result, not yet a superiority claim.

---

## Backlog (low priority)

- **Post-commit Codex review hook:** global git `post-commit` hook running an external Codex review of
  the current-vs-prior diff. Dev-experience only.
- **RAG literature database:** download PDFs as encountered, build a retrieval index for future lit
  search. Useful once Phase 6 research volume grows.
- **Follow-up lab meeting:** schedule once Phase 1–2 results are in hand.

---

## Key decisions / open questions

1. **λ removal vs freeze (Phase 1b):** recommend removing `lambda` from signatures; fallback is
   freezing at 1.0 with a deprecation note if downstream scripts depend on it.
2. **Normalization constant convention:** divide by exactly n·p and n (per-element mean) — the
   interpretation that makes α match DeSurv's mixing. Confirm before touching the update math.
3. **Phase 4 scope:** stratified partial likelihood is the minimal, correct route (no parametric
   baseline); confirm we do not want a fully parametric (e.g. Weibull) baseline instead.

## Global verification

- `Rscript tests/run_tests.R` → 246/246 after each phase.
- `Rscript tests/test_real_data_loading.R` → 77/77 (local, requires `PDAC_DATA_ROOT`).
- Re-run `results/benchmark_sim/run_desurv_comparison.R` after Phases 1 and 3; track external C-index
  vs the 0.636 baseline at every step and report deltas honestly.
- Update `DECISIONS.md` (α/λ reconciliation, normalization convention, strata extension) and
  `ROADMAP.md` (milestones) as each phase completes.
- **Independent review at every checkpoint** (see protocol) — plus a **final whole-branch review**
  against the net-benefit gate before any merge is proposed.
- **Merge to `main` only when the net-benefit gate is met.**
