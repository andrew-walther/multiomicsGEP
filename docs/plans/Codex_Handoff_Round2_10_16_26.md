# Codex handoff, round 2: integrate, fix, diagnose, compare, expand features (2026-10-16)

This replaces the open items of `Codex_Handoff_M6_M7_10_15_26.md`. Read, in this order:
- `CLAUDE.md`;
- the 2026-10-15 entries at the top of `DECISIONS.md` (on `feature/supervised-F`);
- `docs/plans/Codex_to_Claude_M6_M7_10_15_26.md` and
  `git show m6/robustness:docs/reports/Codex_M6_M7_Summary_10_15_26.md` (round-1 results);
- `docs/plans/Meeting_Notes_Followup_10_2_26.md`: "What the meeting concluded", Workstream B and
  Workstream C.

## What Andrew decided (2026-10-16)

1. **K selection (round-1 Prompt E) is resumed with two corrections:**
   - each model is judged by its own convergence rule: the ELBO rule for single modality, and for
     the multimodal model the output-change rule it has used since September;
   - protection must not shield dead programs. A program is protected only if it is credibly
     prognostic: |E[β]| ≥ 1.96·SD(β) **and** SD(β) > 0. Programs that explain neither survival nor
     reconstruction must be removable.
2. **The collapse of the supervised multimodal model on real data is the top priority.** The 10/1
   framework still works on real data (ICGC primary C 0.70 in round 1; 0.685 on 10/2). The M4
   additions (supervised start and the shared point-Laplace β prior) produced β ≡ 0 on TCGA →
   ICGC. Diagnose on TCGA training data only. Remedies are chosen in simulation and by training CV,
   never against ICGC.
3. **DeSurv is installed from the authors' repositories** (approved by Andrew):
   https://github.com/rashidlab/DeSurv (package) and https://github.com/rashidlab/DeSurv-paper
   (paper pipeline).
4. **Feature filtering (advisors, 10/2):** the Bayesian model should be able to use the full gene and
   CpG sets instead of the stringent screening inherited from NMF/DeSurv. Approach: start from the
   current screened set and grow to the full set, warm-starting each step from the previous fit.
5. **Simulation comparators** are rerun with per-arm error isolation and pre-specified handling of
   very small EBMF fits.

Standard for every decision (Andrew): valid, robust, defensible, and performing well. In practice:
- write the design and any adoption rule into the script header before running;
- never change a rule after seeing results;
- report failures as failures, and never select or tune on external validation data.

## Rules for every task

- **Branch.** Work on one cumulative branch, `round2/supervised`, created in Task 0. Make one or
  more commits per task, each message starting with the task number. Never push, never merge into
  `feature/supervised-F` or `main`, and never edit those branches. Run `git branch --show-current`
  before every commit and stage only your files.
- **Do not edit a script while it is running.** (Round 1's D exited 1 for this reason.) Copy it, or
  wait.
- **Defaults stay unchanged.** Default fits of `fit_cox_on_yf()` and `fit_multimodal_yfb()` must
  stay bit-for-bit unchanged unless a task explicitly changes a default. New behavior is opt-in,
  with a test that the default path is identical to the base branch
  (`git show feature/supervised-F:<file>`).
- **Tests.** Run `Rscript tests/run_tests.R` before every commit and report the count. New
  functions get tests that check correct behavior, not only output shape.
- **Errors.** Isolate errors per arm and per fit, so that one failure never discards other arms'
  results. Record every failure with its message.
- **Reporting.** Record results in `docs/progress_book/chapters/2026-10-15.qmd` (new sections) and
  every adoption or rejection in `DECISIONS.md`, in the repo's plain style. Report convergence
  (and iteration caps) for every fit.
- **Longleaf.** Note memory and runtime in the header of any script whose runs exceed about 1 hour
  locally (Task 5), and write a SLURM array wrapper for it. Do not submit jobs; Andrew does.
- **Stops.** If a task cannot be done as specified, stop that task, write down why in
  `docs/plans/Codex_to_Claude_Round2.md`, and continue with the next task whose prerequisites
  are met.

---

## Task 0: integrate round 1

From `feature/supervised-F`, create `round2/supervised` and merge, with `--no-ff` and in this order:
`m5/k-selection`, `m7/sim-comparators`, `m7/multimodal-comparators`, `m6/robustness`,
`m7/desurv-h2h`. Resolve conflicts in `DECISIONS.md`, the progress-book chapter, `PROJECT_STATUS.qmd`
and `tests/run_tests.R` by keeping every entry, in date order. Run the full suite; every test
from every branch must pass. Commit the round-1 handoff note
`docs/plans/Codex_to_Claude_M6_M7_10_15_26.md` if it is not already tracked.

## Task 1: protection fix and the K_init panel (round-1 Prompt E, resumed)

1. **Fix the protection rule** in both pruning nullchecks: protected ⇔ |E[β]| ≥ 1.96·SD(β) and
   SD(β) > 0. Tests:
   - a dead program (β point mass at 0, L and F zero) is removed under protection;
   - a credibly prognostic program is kept even when its removal raises the ELBO;
   - with protection off, behavior is identical to before.
2. **Run `run_k_selection_panel.R`** with the convergence rule clarified: single modality, ELBO
   rule (`converge_on = "elbo"`, the default); multimodal, its own output-change rule. For single
   modality, also record whether the output criterion was met at the last iteration, but do not
   adopt `converge_on = "outputs"` here. Panel as specified in round 1: calibrated beds, noise 0.25
   and 1, 20 seeds, K_init ∈ {5, 7, 10, 15}.
   - Single-modality arms: the selected configuration with `prune = TRUE` and protection; the
     α_F = 0 reference.
   - Multimodal arms: the selected configuration with and without protection, and the 10/1
     framework as reference.
3. **Adoption rule** (unchanged from round 1 except the convergence rule per model): adopt
   protected pruning for a model if, at every K_init and both noise levels:
   - recovery is ≥ 16/20;
   - mean null survival-active programs are ≤ reference + 0.5;
   - every fit converges under that model's rule;
   - median K_final is 3.

   If adopted for single modality, run `run_k_selection_real_1se.R` (Amber's 1-SE rule on 5-fold
   training CV, K_init ∈ {3, 5, 7, 9, 11, 13}). If single modality fails only on convergence, also
   report the panel under the output rule. That second report is for Andrew's decision only; do not
   adopt it.

## Task 2: why the supervised multimodal model collapses on real data (TCGA only)

Do not load ICGC in this task until step 4.

1. **Trace.** Add an opt-in `control$trace = TRUE` to `fit_multimodal_yfb()` that records, per
   outer iteration:
   - the shared β prior (π, Laplace scale);
   - E[β_k], SD(β_k), and the pseudo-observation and its precision (x_k, A_k);
   - the number of live programs;
   - the ELBO components (reconstruction, survival, KL);
   - the relative output changes used for convergence.

   Test that the trace does not change the fit.
2. **Diagnose on TCGA** (the matched expression + methylation training data of
   `results/multimodal_real/run_supervised_multimodal_comparison.R`, K_init = 7). Arms:
   - (a) selected configuration;
   - (b) the same with the β prior held fixed at its starting value (no refit);
   - (c) shared Normal β prior with warm-up;
   - (d) the 10/1 framework (per-factor empirical-Bayes β);
   - (e) the selected configuration without the supervised start;
   - (f) the same arms for the expression block alone (single modality, `fit_cox_on_yf()` with the
     selected configuration, on TCGA only).

   Report when and how β goes to zero, and whether that coincides with the shared prior's slab
   weight π reaching 0. Test the leading hypothesis (the refit shared prior becomes a point mass,
   after which β ≡ 0 is absorbing). Also report why every joint fit reaches the 500-iteration cap:
   which quantity keeps changing.
3. **Remedies** (decide the candidates from step 2's evidence, write them and the adoption rule into
   the header before running, and keep them to the smallest change that addresses the cause).
   Plausible candidates, if the hypothesis holds:
   - R1: a fixed shared β prior (point-Laplace, π = 0.5, scale 1, never refit);
   - R2: refit the shared prior but bound π below at 1/K;
   - R3: per-factor priors as in the 10/1 framework, combined with the supervised start only.

   Evaluate each candidate on (i) the calibrated simulation beds of both models with the round-1
   criteria (recovery ≥ 16/20 at both noise levels, false positives ≤ reference + 0.5, convergence
   under each model's rule) and (ii) 5-fold cross-validated C on TCGA training data (folds
   stratified by event).

   Adoption rule: among candidates passing (i), choose the simplest whose mean training-CV C is
   not below the 10/1 framework's by more than one paired SE. If none passes, report that and stop.
4. **One external check.** Only for the adopted remedy, fit on all of TCGA and score ICGC (primary
   n = 50 and all n = 67) against the 10/1 framework and the expression-only model, exactly as in
   `run_supervised_multimodal_comparison.R`, with the paired bootstrap. This is a report, not a
   selection step.

## Task 3: DeSurv head-to-head (single modality, real data)

1. Install DeSurv with `remotes::install_github("rashidlab/DeSurv")`. Pin the commit SHA in the
   script header and in `DECISIONS.md`. Clone or read https://github.com/rashidlab/DeSurv-paper for
   the documented pipeline: preprocessing, rank and penalty selection (Bayesian optimization), and
   how a fitted model scores new samples. Do not modify the package.
2. Train DeSurv on the same patients as D4 (TCGA_PAAD + CPTAC), using **DeSurv's own preprocessing
   and tuning** on training data only. DeSurv uses NMF on nonnegative input, so it cannot share D4's
   centered matrix. A fair comparison uses each method's own pipeline on the same patients. If the
   paper's tuning procedure is too expensive to run locally, use the hyperparameters the paper
   reports for this training set, say so, and add a Longleaf wrapper for the full tuning.
3. Score the 5 external cohorts with DeSurv's own projection, with orientation fixed on the training
   data. Report:
   - per-cohort C and mean-of-cohorts C;
   - paired differences with the cohort-stratified bootstrap (B = 2000) against D4 and against the
     selected supervised single-modality model (`run_supervised_real_check.R` arms; refit with
     `set.seed(42)`).

Script: `results/benchmark_sim/run_desurv_h2h.R`; outputs in `results/benchmark_sim/outputs/desurv_h2h/`.

## Task 4: simulation comparators, rerun cleanly (round-1 Prompt B)

Rerun `run_supervision_comparators.R` on both noise levels, 20 seeds each, after these fixed changes:
- errors are isolated per arm;
- when EBMF returns 1 factor, the two-step methods use ordinary Cox only (lasso is recorded as
  not applicable);
- when EBMF returns 0 factors, the risk score is constant: the dataset counts as not recovered,
  with C undefined (NA), and is kept in the denominator;
- report how many datasets gave 0, 1 and ≥ 2 EBMF factors, since that is itself a result.

Add DeSurv as an arm only if Task 3 produced a working fitter and DeSurv's preprocessing applies to
the simulated data (it needs nonnegative input). Otherwise state why it is omitted.

## Task 5: feature sets, from screened to full (Workstream B, advisor suggestion)

**Prerequisites:** Task 1 and Task 2 settled the configurations. If Task 2 found no working
supervised multimodal remedy, use the 10/1 framework for the multimodal part and say so.

1. **Continuation fit.** A helper that fits a model on a feature set, then fits on a larger set by
   warm start. Existing features keep their loadings. New features get loadings initialized from
   the current scores: F_new = Y_newᵀ L (LᵀL)⁻¹, the least-squares loadings. Their posterior second
   moments and loading priors are re-estimated in the first sweep. CAVI then continues on the full
   set. It must fit the same model, only reached by a different path. Test on a small example that:
   - the warm-started fit's final ELBO is at least the cold-started fit's, or within tolerance;
   - default cold fitting is unchanged.
2. **Single modality** (TCGA + CPTAC, D4 preprocessing except screening): gene sets of the top
   2,064 (current) → 5,000 → 10,000 → all genes expressed in both cohorts, each step warm-started
   from the previous. Also cold-fit the largest set for comparison.
3. **Multimodal** (TCGA): vary one modality at a time, as the meeting notes ask.
   - Expression: current → all genes, with CpGs fixed at 10,000.
   - CpGs: 10,000 → 30,000 → 100,000 → all available (about 305,000), with genes fixed. If all
     CpGs exceed local memory or time, run up to the largest feasible set locally and write the
     Longleaf wrapper for the rest.
4. **Per step, report:**
   - 5-fold training-CV C (single modality) or training C with a TCGA CV (multimodal);
   - retained and survival-active programs;
   - convergence;
   - each program's loading correlation with the previous step's programs (factor stability).

   The external cohorts are scored once, at the end, for the current and the largest feature set.
5. **No adoption here.** This is a sensitivity analysis. Report whether the larger sets change C,
   the programs or K, and what fitting them costs.

Scripts: `results/benchmark_sim/run_feature_expansion.R`, `results/multimodal_real/run_feature_expansion_multimodal.R`.

## Task 6 (optional, last): clean rerun of round-1 D

Rerun `run_supervised_robustness.R` into a fresh output folder, so the documented process completes
with exit 0. Compare it with the saved tables, which should match exactly.

---

## Order and handback

Order: 0 → 1 → 2 → 3 → 4 → 5 → 6.
- Tasks 3 and 4 do not depend on Task 2. If Task 2 stops, continue with them.
- Task 5 depends on Tasks 1 and 2.

At the end, write `docs/plans/Codex_to_Claude_Round2.md`. For each task give:
- status;
- commits and test count;
- key tables;
- anything that failed or looked wrong.

Claude Code reviews `round2/supervised` before anything is merged.
