# Codex follow-up prompt for tasks D and E

Date: 2026-11-12. Author: Andrew Walther.

This follows Prompt 2 of the [C, D and E prompts](2026-11-12-codex-c-d-e-prompts.md).
That prompt produced `2bfd61b5` (D steps 1 and 2) and `153ac5ea` (E code). An
independent review of both is recorded in the 2026-11-12 notes under Task D and Task E
of the [week plan](../plans/2026-11-09-week-plan.md). Part 1 runs E's fits. Part 2 is the D
fitter work for the option Andrew chooses. Run Part 1 first, or run Parts 1 and 2 in
two separate Codex sessions, because they touch different files.

**Before pasting:** delete the Part 2 options Andrew did not choose. Part 1's priors are
settled (A `bdbd9989`, A2 `b37783c7`, B `b0c73b1c`).

---

## Part 1: Task E fixes and fits

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or change
branches. Read AGENTS.md. Then read these sections of the week plan: "Reporting rules for
all tasks", "Task A" and "Task E" (including its 2026-11-12 notes). Also read the
[E code findings](../findings/2026-11-12-lb-yfb-study-code-findings.md) and the 2026-11-12
entries at the top of DECISIONS.md.

**Fixed choices (do not change without asking Andrew).**

- Input: the z-scored D4 genes, pooled TCGA + CPTAC, no `strata_id` and no
  `center_projection`. Scaling is learned on each training fold.
- L and F prior: `point_laplace` (task A).
- β prior: per-program Normal (task B, `b0c73b1c`), which is `normal` in this runner.
- K = 7 (`cfg$benchmark$k_merged_yfb_desurv`), the config's `max_iter` and `tol`, no
  pruning.
- Genes: the current D4 set in the prepared input.
- Training w_S = 1 (S-15). Simulation w_S ∈ {1, 10}.
- Folds: the five saved S-11 folds stored in `data/lf_prior_real_2026-11-12.rds`.

**Conditional on task F.** Task F (gene screen, scaling, rank) is running now and may
change the gene set or K. These fits are valid for D4 and K = 7 only. Before starting,
check the top of DECISIONS.md for a task F entry. If F has changed the gene set or K,
stop and ask Andrew instead of running.

**Step 1. Fixes before any scientific fit** (code in `code/lb_yfb_study.R` and
`results/benchmark_sim/run_lb_yfb_study.R`, tests in `tests/test_lb_yfb_study*.R`):

1. Write each cohort's held-out concordance (TCGA, CPTAC), with its comparable-pair count,
   to the CV and simulation aggregate tables, next to the within-cohort and pooled C.
   Test it against a hand-computed two-cohort case.
2. Make the training w_S a named runner option (`--w_S`, default 1) and write it to
   `settings.txt`. Test that the setting is recorded.
3. Use `timefix = TRUE` in `lb_yfb_study_concordance()` so it matches task A's
   `survival::concordance()`. If an existing test depends on `FALSE`, explain why in
   the findings and ask before keeping it.
4. Assert that YFB's recorded starting state equals the supplied common start, as is
   already done for LB. Make the mutation check fail when the YFB start is altered.
5. Run the two focused test files, then `Rscript tests/run_tests.R` (816 before this
   change). Commit with a detailed message and push before fitting.

**Step 2. Fits.** Run the two modes at the same time as two separate background
processes, each single-threaded. That is two CPU workers in total. Check first that task
C's methylation runs are not using more than two cores. If they are, run the modes one
after the other.

```sh
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
E_SCALING=zscore
E_PRIOR_LF=point_laplace
E_PRIOR_BETA=normal     # per-program Normal (task B)
RUN_DATE=$(date +%Y-%m-%d)

nohup Rscript results/benchmark_sim/run_lb_yfb_study.R \
  --mode=simulation --scaling="$E_SCALING" \
  --prior_LF="$E_PRIOR_LF" --prior_beta="$E_PRIOR_BETA" \
  --out=results/benchmark_sim/outputs/lb_yfb_simulation_$RUN_DATE \
  > data/lb_yfb_simulation_$RUN_DATE.log 2>&1 &

nohup Rscript results/benchmark_sim/run_lb_yfb_study.R \
  --mode=training --scaling="$E_SCALING" \
  --prior_LF="$E_PRIOR_LF" --prior_beta="$E_PRIOR_BETA" --w_S=1 \
  --input=data/lf_prior_real_2026-11-12.rds \
  --out=results/benchmark_sim/outputs/lb_yfb_training_$RUN_DATE \
  > data/lb_yfb_training_$RUN_DATE.log 2>&1 &
```

What runs:

| Mode | Cells | Fits | Expected time (serial) |
|---|---|---|---|
| Simulation | S-9 shares 0.25%, 1%, 5% and S-14 0.05%, seeds 1–10, w = 1 and 10 (80 cells). 2026-11-06 matched generator, amplitudes 0.5 and 1, seeds 1–6, w = 1 and 10 (24 cells) | 208, plus LB fold-in of each validation set | 45–75 min at 10–15 s per fit |
| Training | 5 folds × (LB, YFB), then one full-data pair on all 273 patients | 12, plus fold-in of each held-out fold | 20–35 min at 1–2 min per fit |

The wall time with both running together is about 1 to 1.5 hours. If a mode has run
for twice its upper estimate, report the progress from the checkpoints and ask before
continuing. A failed cell or fold stays in the tables as an explicit failure. Do not
rerun it with different settings.

**Step 3. Findings note** `docs/findings/<run date>-lb-yfb-study-findings.md` (add a row
to `docs/findings/README.md`). Report:

1. Simulation, per cell and model: recovery successes out of the seeds, mean
   false-active programs, number converged, and mean validation C. Apply the plan's only
   adoption rule. If one model fails A's recovery (≥ 8/10), false-activity (≤ 1) or
   convergence (≥ 9/10) criterion in a cell where the other passes, the failing model is
   not carried forward. The 2026-11-06 cells have no binary recovery gate, so report
   their matched score correlations descriptively.
2. Training CV: for each model and fold, the within-cohort C, each cohort's C and the
   pooled C. Give the five-fold means, the paired YFB − LB difference and its SE, and the
   numbers of converged folds, programs with nonzero loadings and survival-active
   programs.
3. Full-data pair: the correlation of the two linear predictors, the per-program
   correlations of L_k with (YF)_k, the top 20 genes per program, their overlap between
   the two models, and their overlap with DeSurv's programs. Use the existing helpers in
   `code/pathway_enrichment.R`.
4. A short interpretation section suitable for the progress book. Say which score
   (L or YF) each number refers to and what a hazard ratio per unit of each means (week
   plan, "Interpretation to aim for").
5. Limits: the w = 10 cells use survival block weighting, which is not the historical
   S-9/S-14 weighting at w ≠ 1. Fold-internal z-scoring means CV C does not exactly
   reproduce S-11. External scoring is not part of this prompt.

Commit only the aggregate CSVs and settings under the two output folders, the code and
test changes, and the findings note. Patient-level fits and scores stay under ignored
`data/`.

**Done when.** The fixes are tested and pushed, both modes have finished (or stopped with
explicit failures), and the findings note with all five items is pushed.

---

## Part 2: Task D fitter (option-keyed; keep only Andrew's choice)

Read the "Task D" section of the week plan, including its 2026-11-12 notes, and the
[M = 1 findings](../findings/2026-11-12-multimodal-m1-equivalence-findings.md). Also read
the "Extension to M Modalities" section of
`derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`. The
working configuration is the S-11 call in
`results/benchmark_sim/run_yfb_supervised_real.R` (`s11_fit`): point-Laplace L and F,
z-scored input, `canonicalize = TRUE`, w_S = 1, per-program Normal β, pooled.

### Option C (recommended): generalize `fit_ssbmf_yfb()` with feature groups

The derivation's model (shared L, F_m per modality, feature-specific τ_m,
Z = Σ_m Y_m F_m) is the single-modality model fitted to the column-concatenated
[Y_1 … Y_M]. The only difference is that the F prior is fitted for each (modality,
program) pair. The A_l, B_l, Z, V^Z and τ formulas are already sums over features.

1. Tests first (`tests/test_feature_groups_yfb.R`, added to `tests/run_tests.R`):
   - With `feature_group = NULL` and with a single group label, the fit is identical
     to the current fit (`identical()` on EL, EF, EBeta and the ELBO history). Use two
     fixtures: the pooled working configuration and one with `strata_id` and
     `center_projection = TRUE`.
   - For two groups with very different sparsity, each group's fitted F prior equals
     `ebnm()` run on that group's pseudo-observations for the same program. The two
     fitted priors differ.
   - The last monitored objective, recomputed from the returned moments and the
     per-group priors, matches the recorded value.
   - `predict_ssbmf_yfb()` on concatenated new data equals Σ_m Y_m,new F_m β.
2. Code: add `feature_group = NULL` to `fit_ssbmf_yfb()`. In the F update, fit the prior
   separately within each group (`update_F_surv_YFB_k()` and its non-survival
   counterpart), sum the F KL over groups, store the priors by group, and rescale each
   group's prior in canonicalization. Error loudly if `feature_group` is combined with an
   option it does not support (for example `cohort_id`). Do not change
   `code/fit_multimodal_yfb.R`. It becomes historical.
3. Step 3 simulation: write `results/benchmark_sim/run_feature_group_simulation.R` with
   the `low_variance_prognostic_factor` scenario of `code/simulate_multimodal_yfb.R`,
   seeds 1–10. z-score each modality on training rows and fit with the working
   configuration and `feature_group`. Report recovery of the shared prognostic program
   (≥ 8/10), mean false-active programs (≤ 1), convergence (≥ 9/10), and, as a check of
   methylation dominance, each modality's share of the L precision Σ_j τ_mj E[f²_mjk].
   Repeat at a 1:5 expression-to-methylation feature ratio. Modality weights need a
   derivation first, so do not add them. Report dominance if it appears.
4. Step 4 (only if task C's matched-patient preparation exists and step 3 passed): five-fold
   CV with task C's folds for expression alone, methylation alone (10k CpGs) and both.
   Report the CV C and the programs that appear only in the joint fit. This is descriptive.

**Done when.** The tests pass, the full suite passes, step 3 is reported against its
criteria, and step 4 is either reported or marked not run with the reason. Findings go in
`docs/findings/<date>-feature-group-fitter-findings.md`.

### Option A: extend `code/fit_multimodal_yfb.R` to mirror the single sweep

Add point-Laplace L, `strata_id` risk sets, `center_projection`, w_S with block
weighting, the β → L → F order per program with canonicalization after the sweep, a
`posterior_init` that is shared with `fit_ssbmf_yfb()`, the full-ELBO stopping rule, the
dead-program rule, and one KL sign convention. Then rewrite the M = 1 test: same start,
fixed sweep count (`tol = 0`), all multimodal defaults pinned, the working configuration
plus one stratified fixture, and the same objective on both sides at relative tolerance
10⁻⁶. **Done when** both fixtures pass and the full suite passes. Steps 3 and 4 then
follow as in Option C.

### Option B: shared sweep extracted from `fit_ssbmf_yfb()`

Refactor the single-modality iteration into sweep functions that take a list of
modality blocks. Have both `fit_ssbmf_yfb()` and a new multimodal entry point call them.
Prove before anything else that `fit_ssbmf_yfb()` is unchanged: the full suite plus
reproduction of the saved S-9, S-14 and S-11 refits. **Done when** those hold, the M = 1
test passes as in Option A, and steps 3 and 4 follow.

---

## Rules for both parts

- Script headers and roxygen per AGENTS.md, author Andrew Walther, no tool names.
- Do not edit DECISIONS.md, PROJECT_STATUS.qmd, ROADMAP.md, the week plan or
  `docs/plans/2026-10-16-round2-attempts-log.md`. Another session owns them. Put proposed
  text for them in a final "Text for the shared status files" section of each findings
  note, and add each findings note to `docs/findings/README.md`.
- Do not touch the task C files (`results/multimodal_real/run_methylation_single_modality.R`
  and its outputs) or `results/benchmark_sim/run_lf_prior_comparison.R`.
- Stage only your own paths (`git add <paths>`, never `git add -A`). Commit each step
  separately with a detailed message and no tool attribution, then push.
- At most 2 CPU workers. Run a new modular-update change against
  `Rscript tests/run_tests.R` before declaring it done.
- Get an independent review of the code changes and the findings before the final
  commit. Report its result in the findings note.
- Stop and ask if a fixed choice above turns out to be impossible, or if a result would
  require changing the design.
