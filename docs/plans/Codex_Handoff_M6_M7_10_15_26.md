# Codex handoff: M6 robustness and M7 comparators (2026-10-15)

> **Historical Round1 specification — do not execute these prompts as a restart.**
> Its branching instructions and selected configurations describe the original
> experiment, not the current approved state. Round2 superseded the open work.
> Start with [current status](Codex_to_Claude_Round2.md) and the
> [active plan](Advisor_Priorities_and_Execution_10_15_26.md).
> Current work stays on `round2/supervised`; no new Round1 branches, pushes,
> protected-branch edits or merges. Original prompts below are preserved.

Context for every prompt below: survival supervision of the loadings (Workstream C) is implemented
on branch `feature/supervised-F`. Read first:
- `CLAUDE.md`;
- the 2026-10-15 entries at the top of `DECISIONS.md`;
- `docs/plans/Meeting_Notes_Followup_10_2_26.md`, Workstream C and "Sequencing and handoffs";
- `docs/plans/Comparator_Inventory_10_15_26.md`.

**Selected configurations** (do not change them; flag a problem instead):
- **Single modality:** `fit_cox_on_yf(canonicalize = TRUE, w_S = 1, init_supervised = TRUE,
  beta_prior_mode = "shared", prior_beta = "point_laplace", beta_warmup = 10L,
  prior_LF = "point_laplace")`, **without** `prune` (ELBO pruning removes the planted program in
  single modality; the K-selection rule is an open M5 decision). About 1 in 20 fits oscillates
  stably at the 300-iteration cap (DECISIONS.md, KL-sign-fix entry); record `history$converged` and
  report non-converged fits, but do not change the stopping rule.
- **Multimodal:** `fit_multimodal_yfb(..., control = list(intercept = TRUE, prior_update = "ebnm",
  prior_F = list(expression = "point_laplace", methylation = "point_laplace"), prune = TRUE,
  survival_elbo = "plugin", init_supervised = TRUE, beta_prior = "shared_point_laplace",
  beta_warmup = 10L))`. Confirmed by M4 (DECISIONS.md 2026-10-15).
- **Comparator / reference:** the reported single-modality model (α_F = 0, point-exponential), with
  cached risk scores in `results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_riskscores.rds`
  (`[["D4"]]`). Use frozen orientation risk = −η (DECISIONS.md 2026-10-02).

**Rules for every task:**
- One writer in the checkout at a time; Claude Code is idle while Codex runs.
- Work on a new branch off `feature/supervised-F` named in the prompt. Run `git branch --show-current`
  before every commit, stage only your files, and do not push.
- Do not edit `code/fit_cox_on_yf.R`, `code/fit_multimodal_yfb.R` or `code/multimodal_yfb_*.R`. If a
  change is needed there, stop and report it.
- Fix arms, seeds and metrics in the script header before running anything. Do not change them
  after seeing results.
- Run `Rscript tests/run_tests.R` before committing and report the count.
- Each script gets the repo header (Script, Purpose, Author: Andrew Walther, Created,
  Dependencies, Run) and roxygen comments on its functions. Write outputs to its own folder under
  `results/benchmark_sim/outputs/`.
- Record results in a short section of `docs/progress_book/chapters/2026-10-15.qmd` and any
  decision in `DECISIONS.md`, in the repo's plain style (no mannered prose, no session labels).

---

## Prompt E — M5: protected pruning, the K_init panel, and the 1-SE rule (run first)

Agreed with Andrew on 2026-10-15. Background: ELBO pruning (`fit_cox_on_yf(prune = TRUE)`,
`prune_yf_nullcheck()`; multimodal `multimodal_yfb_nullcheck()`) keeps the true 3 programs in the
multimodal simulation, but removes the planted survival program in single modality. Two fits of the
selected single-modality configuration also oscillate stably at the iteration cap
(DECISIONS.md 2026-10-15). The advisors asked for: the model choosing K, a replicated K_init panel
with K_eff ≈ true K and survival programs kept, and Amber's rule (best K over a grid, then the
most parsimonious K within 1 SE of it on held-out C).

```text
On a new branch m5/k-selection off feature/supervised-F, read
docs/plans/Codex_Handoff_M6_M7_10_15_26.md (context and rules; for this
prompt only, you MAY edit code/fit_cox_on_yf.R, code/fit_multimodal_yfb.R
and code/multimodal_yfb_updates.R as described below, keeping every default
fit bit-for-bit unchanged).

1. Protected pruning. Add an option (fit_cox_on_yf(prune_protect_survival =
   FALSE); multimodal control$prune_protect_survival, default FALSE) under
   which a program whose coefficient is credibly nonzero
   (|E[beta_k]| >= 1.96 * sqrt(E[beta_k^2] - E[beta_k]^2)) is never removed
   by the ELBO nullcheck: among the remaining programs, remove the one with
   the largest ELBO gain if that gain is >= 0, refit, repeat. Pass the option
   through the pruning refits. Tests (tests/test_fit_yf_supervision.R and
   tests/test_multimodal_yfb_supervision.R): a survival-active program is
   kept when its removal would raise the ELBO; with the option FALSE the
   pruning sequence is identical to before; default fits identical to the
   current branch (compare against git show feature/supervised-F:<file>).

2. Stopping rule, reported alongside, not adopted. Add
   fit_cox_on_yf(converge_on = c("elbo", "outputs")), default "elbo".
   "outputs" stops when the relative change of L F' and of eta = Y F beta
   are both < 1e-4, using multimodal_yfb_relative_change() as
   fit_multimodal_yfb() does. In the runs below record convergence under
   both rules for every fit (fit with "elbo"; also report whether the
   "outputs" criterion was met at the last iteration).

3. K_init panel in simulation (pre-specified; write it in the script
   header before running). Calibrated beds as in
   results/benchmark_sim/run_supervision_diagnostics.R (n_train = 400,
   beta3 = 3.3) and results/benchmark_sim/run_multimodal_supervision.R,
   noise 0.25 and 1, 20 seeds, K_init in {5, 7, 10, 15}. Arms:
   single modality: selected configuration (DECISIONS.md) + prune = TRUE +
   prune_protect_survival = TRUE; reference alpha_F = 0 (for false
   positives). Multimodal: selected configuration with and without
   prune_protect_survival. Metrics as in those runners plus K_final.
   Adoption rule (fixed now): adopt protected pruning for a model if, at
   every K_init and both noise levels, recovery >= 16/20, mean null
   survival-active <= reference + 0.5, every fit converges (ELBO rule), and
   the median K_final is 3. Record the outcome in DECISIONS.md either way;
   if any part fails, say which and stop there (do not tune).

4. Real data, single modality (only if step 3 adopts protected pruning for
   single modality): Amber's 1-SE rule. With the selected configuration
   plus protected pruning, run 5-fold cross-validated C on the training data
   for K_init in {3, 5, 7, 9, 11, 13}, using the folds and code of
   results/benchmark_sim/run_supervised_train_cv.R (training data only; do
   not load external cohorts). Take the K_init with the highest mean fold
   C, then the smallest K_init whose mean is within one SE (SD of fold C /
   sqrt(5)) of it. Report the table, the chosen K_init, and K_final and
   survival-active programs of the full-data fit at that K_init.

Scripts: results/benchmark_sim/run_k_selection_panel.R and
results/benchmark_sim/run_k_selection_real_1se.R; outputs under
results/benchmark_sim/outputs/k_selection/. Update the progress book and
DECISIONS.md. Run the full test suite before each commit.
```

After Prompt E, Prompts A–D use the configuration that E records in DECISIONS.md (with protected
pruning and the chosen K_init if adopted; otherwise as listed above).

---

## Prompt A — M7: DeSurv head-to-head (single modality, real data)

```text
On a new branch m7/desurv-h2h off feature/supervised-F, read
docs/plans/Codex_Handoff_M6_M7_10_15_26.md (context and rules) and
docs/plans/Comparator_Inventory_10_15_26.md.

No head-to-head DeSurv run exists. Build one:
1. Find how DeSurv (Young et al., PNAS 2026; see memory/reference notes and
   paper/ssbmf-refs.bib) is fitted. If an R implementation is installed or
   vendored in the repo or data folders, use it; if not, stop and report what
   is needed rather than reimplementing DeSurv.
2. Train DeSurv on the same training data and preprocessing as D4 (copy the
   loading block of results/benchmark_sim/run_supervised_real_check.R), with
   its rank and penalty chosen by its own documented procedure using only the
   training data.
3. Score the 5 external cohorts exactly as run_supervised_real_check.R does
   (same gene intersection and per-platform standardization), with the
   orientation fixed on the training data.
4. Report per-cohort C, mean-of-cohorts C, and paired differences against
   (a) the reported D4 model and (b) the supervised single-modality fit
   (refit it with the selected configuration and set.seed(42)), each with the
   cohort-stratified paired bootstrap (B = 2000) used in
   run_supervised_real_check.R.
Script: results/benchmark_sim/run_desurv_h2h.R; outputs in
results/benchmark_sim/outputs/desurv_h2h/. Ask before committing.
```

## Prompt B — M7: comparators on the planted simulation bed

```text
On a new branch m7/sim-comparators off feature/supervised-F, read
docs/plans/Codex_Handoff_M6_M7_10_15_26.md.

The calibrated planted test bed is make_testbed() in
results/benchmark_sim/run_supervision_diagnostics.R with n_train = 400 and
beta3 = 3.3 (noise 0.25 and 1, 20 seeds each, seeds 20261015 + 1:20). Add
these comparators, scored with the same metrics where they apply (planted
recovery |cor| >= 0.7 with correct sign, validation C; for two-step methods
the "program" is the factor with the largest |Cox coefficient z|):
1. unsupervised EBMF (flashier, point-Laplace, greedy + backfit, K up to 7)
   alone: recovery only;
2. two-step EBMF -> Cox (Cox on the EBMF scores; also lasso-Cox);
3. PCA -> Cox with K = 7;
4. the selected supervised model (from DECISIONS.md) for reference;
5. DeSurv, if Prompt A found a usable implementation.
Script: results/benchmark_sim/run_supervision_comparators.R; outputs in
results/benchmark_sim/outputs/supervision_comparators/. Use a socket
cluster as run_supervision_diagnostics.R does. Ask before committing.
```

## Prompt C — M7: multimodal comparators and the constant two-step C

```text
On a new branch m7/multimodal-comparators off feature/supervised-F, read
docs/plans/Codex_Handoff_M6_M7_10_15_26.md.

1. In results/multimodal_real/outputs/real_sweep_summary.csv the two-step
   EBMF -> Cox C is identical (0.418) at starting ranks 6, 7 and 8. Find out
   why (is the stage-1 rank actually varying?). Report the cause with file
   and line; fix it only in results/multimodal_real/ scripts.
2. Rerun the multimodal comparison on TCGA -> ICGC with the selected
   multimodal configuration (DECISIONS.md, M4 entry) against: the 10/1
   pruned framework without supervision, two-step EBMF -> Cox (stacked),
   expression-only selected single-modality model, at starting rank 7.
   Report validation C, number of retained programs, survival-active
   programs, convergence.
Outputs in results/multimodal_real/outputs/supervised/. Ask before committing.
```

## Prompt D — M6: robustness of the selected single-modality model

```text
On a new branch m6/robustness off feature/supervised-F, read
docs/plans/Codex_Handoff_M6_M7_10_15_26.md.

Using the calibrated planted bed (make_testbed(), n_train = 400, beta3 = 3.3,
noise 0.25 and 1, 20 seeds) and the selected single-modality configuration:
1. Jittered starts: 20 starts per dataset, each perturbing the SVD start with
   Gaussian noise of SD 0.1 times each column's SD (EL_init/EF_init via
   init_method = "custom", keeping init_supervised's program by
   constructing it yourself as in fit_cox_on_yf()). Keep the start with the
   highest final history$elbo_full. Report recovery, false survival-active
   programs, and how often the best-ELBO start differs in recovery from the
   default start.
2. Loading prior families: prior_LF in {point_laplace, normal} (the
   point-exponential prior collapses on centered data; see DECISIONS.md).
3. Rerun results/benchmark_sim/run_supervised_train_cv.R unchanged (its
   saved outputs predate a fix to the monitored objective; every fold hit
   the iteration cap). Save the old outputs alongside as *_prefix.csv and
   report whether the C comparison changes.
4. Feature count (real data, Workstream B): top 1000, 2000, 3000 (D4's
   setting; 2064 genes after the cross-cohort intersection) and 5000 genes
   per cohort with D4 preprocessing otherwise; report
   training cross-validated C (5-fold, folds as in
   results/benchmark_sim/run_supervised_train_cv.R) and the number of
   survival-active programs. Do not score the external cohorts.
Script: results/benchmark_sim/run_supervised_robustness.R; outputs in
results/benchmark_sim/outputs/supervised_robustness/. Ask before committing.
```
