# 10/30 realignment plan: cleanup, then survival supervision of F (YFB)

Written 2026-10-30. Approved by Andrew 2026-10-30. Execution prompt: [`../prompts/2026-10-30-codex-execution-prompt.md`](../prompts/2026-10-30-codex-execution-prompt.md).

## Context

The project has several branches of attempted fixes (round2, rebuild, fix/current-method). None has made survival supervise the factorization, which is the method's contribution (Naim, Yusha). Andrew wants to get realigned, clean up, and work the biggest issues in priority order instead of running more experiments on a model that does not yet do what it claims. Claude builds the plan and the Codex handoff. Codex executes the plan and writes a new 10/30 progress-book chapter.

**Decisions made in this session (2026-10-30):**
- Keep η = (YF)β and implement the Gemini specification in `docs/plans/2026-10-30-Unified-Architecture-Alignment.md`. Yusha's suggestion to switch to η = Lβ with fold-in prediction is recorded as the fallback, not adopted.
- Codex does Naim's "before your next step" cleanup items first, on their own branch.
- `rebuild/unified-ssbmf` stays as it is, unmerged. It is used as a reference: it has a correct KL, 216 tests, and an M-block engine. A new clean M-modality implementation comes later (priority 4).

## Where things stand (summary for the chapter)

| Thread | State |
|---|---|
| round2/supervised | Merged and tagged `archive/round2-supervised-2026-10`. The pre-registered supervision test was null (ΔC = −0.003). The writeup does not report this yet. |
| fix/current-method (A1) | Merged to main (45f52d1). No single w_S supervises F at p ≈ 2,000 genes. Survival supplies 10⁻⁶–10⁻⁴ of the F precision, and the w_S ≈ 10⁴ that gives it a 0.5 share breaks convergence. At p = 50, w_S = 1 recovers the planted program in 18–20 of 20 datasets. A2–A6 were not started. |
| rebuild/unified-ssbmf | Not merged or pushed. It meets O1–O3, O5 and O6 as defined, but not O4 (low-variance recovery: 13–20%, against 27% for the two-step method). The joint fit is worse than expression alone (C 0.554 vs 0.707), and it still uses the double log. |
| Naim item 1, KL normalizer | Still on main: `code/compute_elbo.R:82` with the scalar `A_L` from `code/update_L_surv_YFB.R:79`. |
| Naim item 2, double log | Still active: `config/globals.yml:211-218`. Measured effect is C 0.625 vs 0.627 (`docs/findings/2026-10-26-single-log-check.md`). |
| Naim item 3, BO-T12 on Linux | No record. It needs a Longleaf run. |
| Preprocessing, 11 alignments | All awaiting decision (`docs/findings/2026-10-26-preprocessing-audit.md:203-213`). Recommendations are below. |

**The full Gemini thread** is at https://share.gemini.google/aHdNyfCXQhle (it redirects to gemini.google.com/share/2c28b7fef714). I opened it on 10/30. Its final turns reproduce `docs/plans/2026-10-30-Unified-Architecture-Alignment.md` word for word, so that file is the working specification. The thread also recommends **point-normal as the default β prior** ("Use Point-Normal as your default"), with point-Laplace only for very large K. Record the link in `docs/reviews/README.md` next to the 10/17 transcript.

**A caution on the Gemini γ, which Codex must test before building on it.** At D4 dimensions (n = 273, p = 2,064, 139 events), γ = α/(1−α)·p/N_event is about 15 at α = 0.5. That is inside the range A1 already tested. A per-gene precision count under the spec's own scaling (s_k = 1/(√p d_k), (1/n)Σl² = 1) gives a survival share of about β²α / ((1−α) n τ_j d_k²), which is roughly 1/n, not O(1). Balancing the two likelihoods as totals (n·p entries against N_event terms) needs γ ≈ 10⁴, and at that value A1's fit failed. So the real question is this: with the other Gemini pieces in place (damping ρ, a non-collapsing β prior, the hybrid start, point-exponential sparsity, the √p projection scaling), does the fit stay stable at a γ that gives survival a meaningful share? Step S0 checks this before anything else is built. DeSurv solved this same balance problem for NMF + Cox, so the problem can be solved. S0 therefore also reads how DeSurv does it. The DeSurv package is in the isolated library the tests use. Codex should find its loss normalization, the scale of its α, and how it keeps the factorization stable at large p, then compare these with γ. If the problem persists, Codex either finds a fix by diagnosis or writes a standalone problem statement that Andrew can take to Gemini or Codex: what was tried, the numbers, and the exact point where it fails.

## Phase C: cleanup (branch `cleanup/before-next-step` from main; merge to main when tests pass)

1. **KL normalizer (Naim 1).** In `compute_ebnm_kl` (`code/compute_elbo.R:82`), replicate a scalar A to `length(x)`. Add a test that compares the KL to its exact value with scalar and vector A. Check every caller (`code/fit_cox_on_yf.R:1088, 1257-1259, 1446`). → verify: the new test fails before the fix and passes after; `Rscript tests/run_tests.R` passes.
2. **Single log (preprocessing #1, Naim 2).** Set the RNA-seq double-log flags in `config/globals.yml` so every cohort is logged once. → verify: the loader-plus-wrapper output for TCGA, Dijk and PACA-AU equals log2(x+1) of the raw values once (a small test).
3. **Recompute the frozen D4 YFB real-data benchmark once with items 1 and 2** (`results/benchmark_sim/run_YFB_benchmark.R`, run in the background). Write new outputs alongside the old ones. Do not overwrite the tag or `results/benchmark_snapshot_2026-10-14/`. Report old and new external C per cohort, and say whether pruning still removes the planted program (the KL feeds the stopping and pruning steps).
4. **Orientation in comparator scripts (Naim 5).** Replace `max(c_raw, 1 − c_raw)` in `run_ebmf_cox_regularized.R:51-53` and `run_ebmf_cox_external.R:94-97` with orientation fixed on the training data (`frozen_reverse_cindex()`, `code/concordance_ci.R:34`). Rerun them, and report the fair comparison (EBMF at K = 7 + Cox) with the corrected numbers.
5. **Test portability (Naim 8c).** Check whether any test still runs `git show feature/supervised-F:` (grep on main finds none, so confirm). If any does, store its baselines as fixtures under `tests/fixtures/`.
6. **Yusha's points 1–2.** Check in code whether the single-block risk score has the documented sign (adverse program β > 0 under hazard ∝ e^η), and whether stored β variances are zero when A.5 implies they should be positive. Fix any real bug and add a test for it. If the behaviour is a convention, document it.
7. **Branch inventory (read only).** Make a table of every local branch: merged into main or not, last commit date, purpose, recommended keep/archive. **Delete nothing.** Deletion is Andrew's decision.
8. **BO-T12 on Linux.** This is not pursued, because all work runs locally. It goes on the tracker as a known platform-specific risk, and it is revisited only if work ever moves to Linux.

Phase C merges to main once the tests pass, because every item fixes a problem that exists on main. Its items are the "progress made" list in the 10/30 chapter.

## Phase S: survival supervision of F under YFB (branch `fix/yfb-supervision` from the merged cleanup)

**S0. Precision-share audit (no fitting code).** Using D4 dimensions and a current fit's L, F, β and τ, compute the per-gene survival share of A_F under the Gemini scaling for α ∈ {0.5, 0.9, 0.99, 0.999}, analytically and numerically. Report the γ needed for a share of 0.1 and of 0.5. If the spec's γ formula gives a share far below 0.1, use γ that targets the share directly (γ chosen so the median survival share equals a target ρ_S ∈ {0.05, 0.2, 0.5}) as the grid in S2. Record this as a deviation from the spec in DECISIONS.md.

**S1. Implement the specification at M = 1** as one new opt-in fitter, `code/fit_yfb_supervised.R`. It holds the modality blocks in a list so that M = 2 is a later extension, not a rewrite. It follows spec §4 exactly:
- nonnegative Y (within-cohort division by the SD, no centering), with centered Ỹ in the projection;
- s_k = 1/(√p d_k) with d_k fixed at the start of each sweep;
- γ and the curvature floor w_i = max(a_i, 10⁻⁶), using the existing log-space Breslow risk sets (`control$cox_arithmetic = "logspace"`);
- point-exponential priors on L and F;
- Gauss–Seidel F update over genes with the coupling term and damping ρ (reuse or adapt `multimodal_yfb_update_F_mk()`);
- a non-collapsing β prior: `ebnm_point_normal` with the slab proportion bounded to π_β ∈ [1/K, 0.5] and slab SD ≥ 0.05, which is the Gemini default. The inverse-gamma Normal is a fallback, used only if the bounded point-normal cannot be fitted through `ebnm`;
- the hybrid start (K−2 NMF/SVD factors plus ± Cox-score columns, followed by a β warm-up);
- RMS canonicalization of L;
- KL computed directly from the `ebnm` return value;
- no `prune_nullcheck`; `classify_factors()` is applied afterwards.

Do not change `fit_cox_on_yf()` defaults. Read the ssbmf package on `rebuild/unified-ssbmf` (`ssbmf/R/updates.R`, `cox.R`, `priors.R`) for code to reuse, but do not merge it.
- Tests (`tests/test_fit_yfb_supervised.R`):
  - The KL equals its exact value.
  - Z̃, β and g_β do not change when L and F are rescaled by c_k.
  - The β prior variance stays > 0 when the data carry no survival signal.
  - With α = 0, the fit equals the same engine with the survival term removed (an unsupervised EBMF).
  - The objective does not decrease beyond tolerance, or a decrease is reported, on a small fixture.
  - Fail loudly: non-convergence, NaN η and an all-zero L/F are returned as flags, never as silent convergence (Naim item 3's failure mode).

**S2. Decision test in simulation.**
- Extend the A1 test bed (`results/benchmark_sim/run_supervision_diagnostics.R`) to p ∈ {50, 2,000}. Use the calibrated bed and a planted program with variance share in {0.25%, 1%, 5%}.
- Use 10 seeds per cell to start, and use the A1 metrics (recovered, false active, survival share, convergence, validation C).
- Arms:
  - α = 0 (unsupervised two-step);
  - the S0 γ grid;
  - the current `fit_cox_on_yf` canonical path with w_S = 1, as the baseline;
  - a survival-only lasso, as the ceiling.
- Register the arms and the gate in `docs/plans/2026-10-16-round2-attempts-log.md` before running.
- **Gate S (pre-registered):** supervision "works" if, at p = 2,000 and some γ:
  - recovery exceeds α = 0 by ≥ 3/10 datasets in at least two variance cells;
  - mean false-active ≤ 1;
  - ≥ 90% of fits converge.
- If the gate passes, run a one-at-a-time ablation at the best γ, removing the hybrid start, damping and the β prior in turn, to show which pieces matter.
- If it fails, do not tune further. Write diagnostics instead:
  - per-sweep traces of survival share, β, ‖L‖ and the objective;
  - the sweep where instability starts;
  - which piece fails first.

  Record the failure, and stop Phase S. Andrew then decides whether the η = Lβ fallback is next.

**S3. Real-data training CV, only if Gate S passes.** Adapt `results/benchmark_sim/run_a1_real_ws_cv.R` to the new fitter on the single-log D4 data. Use training data only, and do not score any external cohort (Naim item 6: keep at least one cohort untouched for the final claim).

## Phase P: preprocessing decisions (recommendations; Andrew approves this plan or overrides an item)

| # | Recommendation |
|---|---|
| 1 | Single log for every RNA-seq cohort (done in C2). |
| 2 | One TCGA expression source; keep the one D4 uses and record which it is. |
| 3 | The Gemini option: X/sd within each cohort, nonnegative, with centering only inside the projection. This is a new option (d) in the audit, and it applies to the new fitter only. |
| 4 | Scale each validation cohort with its own within-cohort SD, consistent with #3. |
| 5 | Gene names from training data only; keep the DeSurv per-cohort top-3,000 screen with intersection. Survival pre-screening only as an ablation (Gemini §7). |
| 6, 7 | Choose one missing-feature rule and one duplicate-symbol rule, namely the ones D4 uses now, and record them. |
| 8 | Keep β-values for now; M-values are deferred to the multimodal phase. |
| 9 | Keep KNN imputation and state the donor counts (50 vs 67) explicitly. |
| 10 | Apply the cross-reactive and SNP probe filters when multimodal work resumes. |
| 11 | One record per donor in pooled analyses (568 unique donors), with a donor-level or cohort-stratified bootstrap. |

Codex records these in DECISIONS.md and in the audit's §5. Only #1 changes code in this session. #3 and #4 apply inside the new fitter. The rest are recorded for when the pipeline is next touched.

## Later, not in this session (listed on the tracker)

These items proceed as soon as Gate S passes, in this order:
- Priority 2: multimodal convergence (Gemini issues 5 and 8: modality weight v_m, √p_m scaling), as the M = 2 case of the same fitter.
- Priority 3:
  - Survival pre-screening, as an ablation only.
  - A feature-count sensitivity curve: refit as the screened feature set grows, and find the smallest set at which C-index and the loading correlations with the largest set are about the same.
- Priority 4: a unified ssbmf package, M-modality with M = 1 as the single-block case, on a new branch.
- Writeup:
  - Naim 4: report the 10/15 null result.
  - Naim 6: list every use of the external cohorts; donor-level resampling.
  - Naim 7: PurIST, DeSurv, elastic-net and clinical comparators.
  - Naim 8: writeup computed from saved results with R chunks; fits on /proj/rashidlab.
- CpG promoter annotation; compressing the CSV traces (Gemini 13).

Later code consolidation follows the Gemini layout: a few core scripts that become the R package, with other scripts used only to produce results and biology analyses.

## Phase D: documentation (Codex, throughout)

- **New chapter `docs/progress_book/chapters/2026-10-30.qmd`**, added to `_quarto.yml`. This chapter is the main deliverable for the upcoming advisor meeting. It opens with a short list of clear takeaways (what was fixed, what was learned, what is still blocked and why) that can be shared as is. **Chapter style rules:**
  - Name no reviewer or advisor; describe issues, solutions and progress only.
  - No branch names, file paths or commit hashes. Write "we tried a unified implementation in a separate line of work", not a branch name.
  - Plain statistical prose.

  Do not edit the 10/26 chapter; carry content over from it only where needed. Sections:
  1. Where the project stands: the thread table above and the branch inventory.
  2. **Issue tracker:** one row per issue (the 13 specification issues, the 8 items from the 10/26 code review, and the 6 points from the March-vs-October derivation comparison, merging duplicates), with issue, diagnosis, proposed solution, status and evidence. The internal CSV may keep source tags; the rendered table shows no names. Store it as a CSV, `docs/progress_book/data/2026-10-30-issue-tracker.csv`, and render it with an R chunk so later chapters update the same file.
  3. Cleanup results: old vs new benchmark numbers, read from the saved CSVs by R chunks, not typed in (Naim 8).
  4. The S0 precision audit, and the S1/S2 results or the diagnosed failure.
  5. Preprocessing decisions.
  6. How the model changed since the March derivation: changes are tracked in DECISIONS.md, the March-to-October comparison is summarized, and the choice to keep the projection predictor is recorded with the latent-L predictor as its fallback.
  7. Questions for advisors.
- Update DECISIONS.md, ROADMAP.md, PROJECT_STATUS.qmd, the `docs/` folder READMEs and the CLAUDE.md "Start here" section (it points to this plan). Name documents `YYYY-MM-DD-topic-type`.

## Handoff mechanics

Claude does Phase C first (2026-10-30), then hands the rest to Codex. Progress on Phase C is recorded at the end of this file.

1. Copy this plan to `docs/plans/2026-10-30-realignment-plan.md` and add its row to the README.
2. Write `docs/prompts/2026-10-30-codex-execution-prompt.md`: the order of phases (C → P records → S0 → S1 → S2 → gate → S3 or diagnostics → D), the stop conditions, and the project rules (attribution, commit style, plain prose, tests after every change to an update script, report failures).
3. Update the CLAUDE.md "Start here" section, commit on main and push. Claude does no other execution.

## Verification (end to end)

- `Rscript tests/run_tests.R` passes after each phase, with the new KL, single-log and fitter tests included.
- The benchmark rerun completes and its outputs are saved. The chapter renders with `quarto render docs/progress_book`, and its tables come from CSV/RDS files.
- Gate S is reported against the criteria registered before the run, including when it fails.

## Stop conditions for Codex

Stop and write a status note in the chapter, without working around the problem, if:
- tests fail and the cause is unclear;
- the benchmark rerun changes external C by more than 0.02 without explanation;
- Gate S fails, after the diagnosis described in S2 and a written problem statement (S0 is a check, not a stop: it changes the γ grid);
- any step would touch external cohorts in a way that informs model choices.

## Phase C progress (updated 2026-10-30)

Work is on `cleanup/before-next-step`, branched from `main` at 69ece3b.

| Item | Status | Evidence |
|---|---|---|
| C1 KL normalizer | **done** (31db892) | An exact-KL test fails before the fix (−28.24 vs −22.82) and passes after. Supervision-T6 now compares the two paths after the same 400 sweeps; they agree to 4×10⁻⁷. |
| C2 single log | **done** (6a0b96f) | All `platform_log_transform` flags are false. Real-data test T8.1 passes 89/89. The loader maxima (TCGA 22.6, PACA-AU RNA-seq 17.9) show the test fails under the old flags. DECISIONS.md has a 2026-10-30 entry. |
| C3 benchmark rerun | **running / see below** | `run_single_log_check.R --out results/benchmark_sim/outputs/kl_fix_rerun_2026-10-30` refits D4 under both log settings with the KL fix: 5-fold training CV plus one external scoring. The KL fix matters for default fits too, because every path stops on `elbo_full` (`code/fit_cox_on_yf.R:1325`). |
| C4 orientation | **done** (b85d390) | Both comparator scripts use `frozen_reverse_cindex()`. Rescoring the saved risks changes none of the 20 cohort-arm C values (`outputs/ebmf_cox_external/2026-10-30-orientation-rescore.csv`), so the reported comparator numbers stand. |
| C5 test portability | **done** (1703bc6) | The baselines are in `tests/fixtures/feature-supervised-F/`, read with `read_baseline_source()`. No test calls git. 711/711 pass. |
| C6 sign and β variance | **open, for Codex** | Starting points: Phase C's global flip in `fit_cox_on_yf()` (around lines 1362–1420) can override the partial-likelihood sign; `update_beta.R:185` stores variance plus mean². Check the saved D4 fit's `EBeta`/`EBeta2` and the sign of the adverse program. |
| C7 branch inventory | **done** | `docs/progress_book/data/2026-10-30-branch-inventory.csv`. 27 of 32 local branches are merged into main. Nothing was deleted. |
| C8 BO-T12 on Linux | **dropped** | All work runs locally. It goes on the tracker as a known platform risk. |
