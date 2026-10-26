# Parallel work on the current method and the rebuild: plan

> **For agentic workers:** One coordinating session executes this plan. It starts session B (the rebuild) first, then does the merge, starts session A, and produces the meeting chapter. Each session reads this file, the plan it executes, and `DECISIONS.md` before starting. Andrew approves every **Gate**. Within a task, use superpowers:executing-plans. If the coordinating session cannot start other sessions, it tells Andrew immediately and runs session B's work through subagents instead, ahead of everything except the chapter.

**Priority if the work cannot all run in parallel:** session B's answer to "does the rebuilt model fix supervision and convergence?" comes first, then the meeting chapter, then session A.

**Goal:** Make progress on the current single- and multi-modality methods (supervision of F, multimodal convergence, preprocessing, feature experiments) while a separate session rebuilds the method as the `ssbmf` package, so that whichever route produces a working model first is ready, and each informs the other.

**North star (from the advisors):** a working model comes first. A model is *working* when it meets the operational criteria O1–O8 of `2026-10-22-unified-method-rebuild-plan.md` (monotone objective, convergence, recovery of a planted low-variance prognostic program, external C no worse than the benchmark, coefficient uncertainty). Experiments that depend on convergence wait for it.

## Decisions taken in this plan (approve with the plan)

1. **Merge `round2/supervised` into `main` today, after tagging it.** All round2 model changes are opt-in with defaults unchanged (684/684 tests on round2), the current-code work needs them (feature-sensitivity helpers, log-space Cox arithmetic, CpG tools), and two diverged lines make parallel sessions harder. Cost: round2's diagnostic CSV traces come with it; trimming them is later documentation work.
2. **The round2 close-out plan is not a blocker except for its Task 1 (tag) plus the merge.** Its Task 2 (findings report) becomes one section of the 10/26 chapter; Task 3 (port vetted code) is replaced by the merge; Tasks 4–5 (document restructure, pinned chapter) happen after the meeting.
3. **Three sessions, each with its own branch and its own folder** (git worktrees beside the main checkout, so sessions never switch each other's branches):

| Session | Branch | Folder | Writes only |
|---|---|---|---|
| Coordinator | `main` | `multiomicsGEP/` | `docs/`, `DECISIONS.md`, `ROADMAP.md`, `PROJECT_STATUS.qmd`, `CLAUDE.md`, `results/feature_selection/`, `results/benchmark_sim/run_single_log_check.R` |
| A: current method | `fix/current-method` | `../multiomicsGEP-current/` | `code/`, `tests/`, `results/benchmark_sim/`, `results/multimodal_real/`, its own `docs/findings/` files |
| B: rebuild | `rebuild/unified-ssbmf` | `../multiomicsGEP-rebuild/` | `ssbmf/`, `derivations/unified_ssbmf/`, `results/ssbmf/`, its own `docs/findings/` files |

   Only the coordinator edits the shared living documents; A and B report results in their own findings files and messages, and the coordinator records decisions. A merges into `main` at its gates; B stays separate until rebuild Gate 6b. Session B starts from the current `main` right away (it writes only new folders and does not need the round2 merge); session A starts after the merge.
4. **Division of work between A and B.** A runs *short, decisive tests* on the current code that answer questions the rebuild also needs (Does a survival weight of the right size make supervision of F work? Do the bookkeeping fixes make the multimodal fit converge?). A does not attempt open-ended repair: two failed attempts on the same problem end that line and the finding passes to B. B is the long-term route to a unified, packaged method. Both report against the same criteria.
5. **Expression-only experiments may run now on the current code.** The single-modality model converges (all 24 feature-count fits did), so the survival-screen comparison is valid on it. CpG and joint feature-count curves wait for multimodal convergence (session A Phase A2 or the rebuild).

---

## Coordinating session: before the 4 pm advisor meeting (10/26)

### P0. Start session B (first action)
- [ ] Start session B in its own worktree on branch `rebuild/unified-ssbmf` from the current `main`, with the brief in the Session B section below. Then continue with P1.

Estimated fit cost: about 3 s per single-modality fit at 2,064 genes (feature-expansion runtimes), so P3–P4 take minutes of compute.

### P1. Tag and merge round2 (about 1 h)
- [ ] `git tag archive/round2-supervised-2026-10 round2/supervised`; verify `git show archive/round2-supervised-2026-10:docs/plans/Round2_attempts_log.md` resolves.
- [ ] `git merge --no-ff round2/supervised` on `main`. Resolve: `CLAUDE.md` and `ROADMAP.md` (keep main's reorganized structure, carry round2's content); move round2 documents into the 10/25 naming (`docs/reports/*` → `docs/findings/2026-10-16-*`, round2 plans/handoffs/attempt log → `docs/plans/2026-10-15-*`/`2026-10-16-*`, `docs/plans/archive/*` → plans folder marked superseded), adding rows to each folder README; do not commit rendered `_book/` files if main ignores them.
- [ ] Replace `CLAUDE.md`'s "Start here (2026-10-16)" with a pointer to this plan and the rebuild plan.
- [ ] Verify: `Rscript tests/run_tests.R` passes on merged `main` (expect the round2 count, 684); `shasum -a 256 -c results/benchmark_snapshot_2026-10-14/MANIFEST.sha256` passes. Commit; ask Andrew before pushing.

### P2. Preprocessing audit (rebuild Task 13 / 10-16 Task A1; about 1.5 h)
- [ ] Write `docs/findings/2026-10-26-preprocessing-audit.md`: for each cohort and data type, the source file, transformations in order, missing values, screening (rule, count, which patients), scaling, alignment of validation data, and the function that implements each step; side-by-side single-modality vs multimodal table.
- [ ] Record the confirmed double log₂(x + 1) on TCGA, Dijk and PACA-AU RNA-seq (`load_pdac_raw()` then `preprocess_merged_cohorts()` with `platform_log_transform`), single log on CPTAC and microarrays.
- [ ] Quantify on the frozen D4 training matrix (training only): share of negative entries after z-scoring, share of reconstruction error on negative entries (nonnegative L and F cannot reconstruct below-average values).
- [ ] End with proposed alignments (one choice per step), marked for Andrew's decision.

### P3. Effect of the double log (rebuild Task 14; about 30 min)
- [ ] `results/benchmark_sim/run_single_log_check.R`: D4 configuration (seed 42, K = 7, α_F = 0) on single-logged data. Report training 5-fold CV C (same folds), overlap of the screened genes with the frozen 2,064, and per-cohort external C with a donor-level paired bootstrap against the frozen risk scores. External scoring of this pipeline correction is part of approving this plan; nothing is selected from it.
- [ ] Outputs to `results/benchmark_sim/outputs/single_log_check/` (summaries committed; fits ignored). The benchmark snapshot is not modified.

### P4. Survival-screen experiment, expression only, training CV (feature-selection plan Tasks E1–E4; about 1 h)
- [ ] Run `2026-10-16-feature-selection-experiment-plan.md` E1–E4 on `main` after the merge (its helpers come with round2). Gate E2 defaults, approved with this plan: N ∈ {250, 500, 1000} matched to U(N)'s shared counts; screen statistic = Wald p-value of `coxph(Surv ~ y + strata(cohort))` on fold-training patients; no survival-only third arm.
- [ ] Pipeline: the benchmark's current (double-log) preprocessing, so arms compare with D4; if P3 shows a material change, rerun on the corrected pipeline later.
- [ ] Report CV C ± SE by arm and size, convergence, programs retained, and program agreement with U2064. **External scoring (E5) waits for Andrew's review after the meeting.**

### P5. The 10/26 progress-book chapter (about 1.5 h)
Fill `docs/progress_book/chapters/2026-10-26.qmd` (currently a placeholder).

**Audience and style.** Write for a reader outside the project (advisors, lab members): what was done, what was found, what failed, what comes next. Organize it as discrete sections, one per issue, experiment or piece of progress, each able to start a discussion. Do not use planning labels or internal names: no "round 2", task or phase codes (P2, A1, Gate…), "two-track", session names, assistant or reviewer names ("Gemini", "Codex", "Claude"), or branch names. Describe methods and results as a statistician would: model, data, metric, result, limitation. Follow the repository's no-mannered-prose rule.

**Sections (rename to describe their content):**
- [ ] Summary: the state of the method in a few sentences, the main findings since 10/16, and what is being decided.
- [ ] Problems found in the current fitting algorithm: survival has almost no influence on the loadings (≈ 10⁻⁷ of the loading-update precision); the multimodal fit does not converge and removes no programs on real data; inconsistencies in how the fitting objective is tracked; survival coefficients that become stuck at zero. One short paragraph each, with the evidence.
- [ ] Fixes tried on the current code and their outcomes: a table (log-space risk-set arithmetic; fixed coefficient prior; profiled Cox curvature; fixed loading prior; protected and posterior-preserving pruning) with what each fixed and what remained; the candidates not yet tried.
- [ ] Data preprocessing: the audit across both analyses, the double log transform and its measured effect on the benchmark (P2, P3).
- [ ] Survival-associated gene screening: design, training cross-validation results, what it means, what is next (P4).
- [ ] Rebuilding the method as an R package: why, the model and fitting principles (one model for one or more data types; an objective that every update must improve), the planned interface (fit, predict, coefficients with intervals), and the criteria that define a working model. Include session B's progress and any early results available by 2 pm.
- [ ] Work continuing on the current method: the tests planned for supervision of the loadings and multimodal convergence, and which experiments wait for a working model.
- [ ] Questions for the advisors: the criteria and thresholds for a working model; whether a survival weight is acceptable if the model cannot set the balance itself; external scoring of the screening experiment.
- [ ] Render the book, check the chapter, and publish only with Andrew's OK (the book deploys to `gh-pages`). Target: rendered by 2 pm for Andrew's review.

If session A's Phase A1 or session B has results by 2 pm, add them as preliminary and say so.

---

## Session A: the current method (`fix/current-method`)

Starts after P1 (branch from merged `main`). Each phase ends at a **Gate** with keep / revise / drop; two failed attempts on the same problem end that line. A writes `docs/findings/2026-10-2x-current-method-<phase>.md` and messages the coordinator.

### A1. Supervision of F, single modality (start today)
Question: is the survival contribution simply too small? The round2 audit measured the survival share of the F-update precision at about 10⁻⁷ (multimodal) and recorded small shares in the single-modality fits.
- [ ] On the existing simulation with a planted low-variance prognostic program, fit `fit_cox_on_yf()` with survival in the F update (the `w_S` survival contribution as implemented after the merge; first record whether that update includes the within-factor coupling term, since the older `alpha_F` update did not) over w_S ∈ {1, 10², 10⁴, 10⁶} and with supervised initialization on/off. Measure, with the loading-balance audit helpers: survival share of the F precision; recovery of the planted program (|r| ≥ 0.7, survival-active, correct sign); convergence; held-out C; reconstruction loss.
- [ ] Repeat the best two settings on real training data (TCGA + CPTAC, 5-fold CV only).
- [ ] **Gate A1:** does any w_S make supervision of F work without breaking convergence or reconstruction? The answer goes to session B's Part 4 either way.

### A2. Multimodal convergence and β collapse
- [ ] Add opt-in controls (defaults unchanged) to `fit_multimodal_yfb()`: KL recomputed against the refitted loading prior; the β prior variance carried through canonicalization; a floor on the Cox working weights (as in the single-modality code); a lower bound π_β ≥ 1/K on the shared β prior (the untried R2). Tests: defaults reproduce the saved objects exactly; each control changes only what it names.
- [ ] Run the existing multimodal simulation panel and the TCGA K = 7 fit with all four controls on, then one at a time if they help. Measure convergence within 500 sweeps, objective decreases per sweep, β collapse.
- [ ] **Gate A2:** if the fit converges on real data, A3 and A5 open; if not after two attempts, stop and hand the diagnosis to session B.

### A3. Multimodal programs without deletion
- [ ] On converged A2 fits, classify programs by role on one fit (as `classify_factors()` does for single-modality) instead of deleting and refitting; report roles across K_init ∈ {3, 5, 7, 10}. **Gate A3.**

### A4. Preprocessing alignment in the current code
- [ ] After Andrew decides the P2 alignments: single log for every RNA-seq cohort; the same expression scaling in both analyses; one TCGA expression source. Implement as configuration in `config/globals.yml` with the old setting reproducible. If P3 showed a material change, rerun the D4 benchmark on the corrected pipeline. **Gate A4.**

### A5. Feature-count curves (after Gate A2 passes)
- [ ] CpG-count curve (10k → 305k) and joint expression × CpG curves with the converged multimodal configuration; expression curve rerun if A4 changed preprocessing. Goal: the smallest feature set that keeps performance and program stability. **Gate A5.**

### A6. Survival-screen follow-up
- [ ] External scoring of the P4 experiment once (E5) after Andrew's review; if the multimodal fit converges, the same comparison for CpGs. **Gate A6.**

---

## Session B: the rebuild (`rebuild/unified-ssbmf`)

Executes `2026-10-22-unified-method-rebuild-plan.md`.

### Gate 0 (needed before B starts)
Andrew confirms D1–D11 of the rebuild plan. Recommended answers, accepted if he approves this plan without changes: fresh package `ssbmf/`; criteria O1–O8 as tabled; projection predictor; Cox strata for cohorts; no missing-block prediction in version 1; two-step comparator everywhere, DeSurv later; Longleaf for panels; install `testthat` and `devtools`; package `ssbmf` with Andrew as author until release; features × samples input by default; package inside this repository until release.

### Unattended run (starts now, continues overnight 10/26 → 10/27)
**Goal of the run:** answer whether the rebuilt model fixes the problems of the current implementation, in simulation and on real data, for one data type (expression) and two (expression + methylation). Because B works only on its own branch and new folders, the rebuild plan's review gates (1a–1c, 3, 4a–4b, 6a–6b) become *review after the fact*: B makes the provisional choice the rebuild plan recommends, records it with its reason, and continues. B still stops where a step fails its own check.

- [ ] Task 0: scaffold; `R CMD check` passes on the empty package.
- [ ] Tasks 1–3: specification and derivations; choose the Cox treatment provisionally; run `verify_derivation.R`. **Stop and report** if the bound or stationarity checks fail.
- [ ] Tasks 4–8: package implementation with tests.
- [ ] Tasks 9–10: simulators and a reduced local panel (5 datasets per cell instead of 20; both the one-block and two-block realistic designs and the low-variance prognostic scenario), evaluated against O1–O4 and O7. If O1–O3 fail, **stop and report**; do not tune parameters to pass.
- [ ] Part 4, only if O4 fails at w = 1: try the model-based balance options in the order of the rebuild plan (effective-feature block weight; then a non-absorbing β prior if not already the default). A survival weight chosen by training cross-validation is the last option and is reported as such. Every choice is made on training data or simulation only.
- [ ] Part 6 on real data (local data: `PDAC_DATA_ROOT` set to the OneDrive `PDAC_data` folder; matched TCGA/ICGC cache `data/processed_multiomics_tcga_icgc.rds`), using the current benchmark preprocessing unless the coordinator has recorded a change:
  - Expression only, TCGA + CPTAC training: `ssbmf_cv()` over K_init ∈ {5, 7, 10, 15}; freeze; score the five external cohorts once; compare with the frozen benchmark (O5) and report convergence (O2), program roles and coefficient intervals.
  - Expression + methylation, TCGA training (n = 144): joint and expression-only fits through the same `ssbmf_fit()`, K_init by training CV; freeze; score ICGC primary (n = 50) once; paired donor-level comparison (O6).
- [ ] Do not submit Longleaf jobs, push, or change `main`. Fits stay git-ignored under `results/ssbmf/outputs/`.
- [ ] Report: `docs/findings/2026-10-27-rebuild-first-results.md`, written for a reader outside the project: a table of O1–O8 (met / not met / not tested, with numbers), what the rebuilt model does differently from the current one, every provisional choice to review, and test counts. Message the coordinator with a short summary as soon as the simulation verdict (O1–O4) is known, so it can reach the 4 pm chapter if ready.

### After the review
Andrew reviews the provisional choices in one sitting. Accepted choices are recorded in `DECISIONS.md`; the full simulation panel then runs on Longleaf, and the rebuild plan continues from there. If the rebuilt model meets the criteria, the experiments in session A's A5–A6 move to the package.

---

## Coordination rules

- Information between sessions goes through messages to the coordinator and the findings files; the coordinator records decisions in `DECISIONS.md` and keeps `ROADMAP.md` current.
- Session A's Gate A1 and A2 results go to session B before B's Part 4 and Task 2 choices respectively.
- Shared criteria: both tracks report against O1–O8. A result from either track that meets them goes to Andrew as a candidate to replace the benchmark (snapshot rule: at least as good on every goal and better on one).
- External cohorts are scored once per frozen configuration and never choose anything.
- Neither A nor B pushes; the coordinator asks Andrew before every push.
- When the rebuild passes Gate 6b, a separate cleanup plan decides which current code and documents to archive.
