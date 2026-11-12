# Codex prompts for tasks C, D and E

Date: 2026-11-12. Author: Andrew Walther.

Two prompts, run **concurrently in two separate Codex sessions**. Prompt 1 (C) is
long-running compute. Prompt 2 runs D and then E **sequentially in one session**,
since both are code and tests with no long fits. Task A runs at the same time in a
Claude session, which writes `results/benchmark_sim/run_lf_prior_comparison.R` and
owns the shared status files.

---

## Prompt 1: Task C, methylation alone (Codex session 1)

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or
change branches. Read AGENTS.md, then in `docs/plans/2026-11-09-week-plan.md` the
"Reporting rules for all tasks", "Task A" and "Task C" sections, then the
2026-11-12 entry at the top of DECISIONS.md, which sets the design below and
replaces the plan's fixed 10,000-CpG screen.

**Goal.** Fitted to methylation alone, does the single-modality model find
prognostic programs that expression does not, and how do performance and run time
change with the number of CpGs?

**Design (fixed, do not change without asking Andrew).**

- Patients: the matched TCGA patients with expression, methylation and survival
  (about 144). ICGC is used only to score fixed models, never to choose anything.
- CpGs: the prepared universe in `data/round2_tcga_feature_universe.rds`
  (built by `results/multimodal_real/prepare_feature_universe.R`; check its contents
  and the mask helpers in `code/cpg_feature_sensitivity.R` first). Keep CpGs whose
  mean β over TCGA training patients is in [0.1, 0.9] (expected 182,934). Rank them
  by variance. Counts: top 10k, 25k, 50k, 100k and all. Both steps are outcome-free
  and computed once before the folds. Verify the 182,934 count and report it.
- Values: β values, not M-values.
- Model: the single-modality `fit_ssbmf_yfb()` with the working-model call in
  `results/benchmark_sim/run_yfb_supervised_real.R` (`s11_fit`): K = 7, no pruning,
  `canonicalize = TRUE`, `w_S = 1`, `prior_beta = "normal"`, the config's `max_iter`
  and `tol`. Two setups at every count:
  1. PL: z-score each CpG, `prior_LF = "point_laplace"`.
  2. PE: divide each CpG by its SD (no centering, stays nonnegative),
     `prior_LF = "point_exponential"`, `center_projection = TRUE`.
  One cohort only, so no `strata_id`.
- Comparator: expression with the D4 genes on the same matched patients and the same
  folds, under both setups.
- Folds: five-fold, balanced by event, one fixed seed, saved and shared by every arm
  and count.

**Order.**

1. Write `results/multimodal_real/run_methylation_single_modality.R` (script header
   per AGENTS.md conventions, roxygen on every function, author Andrew Walther, no
   tool names). Smoke-test it on 2,000 CpGs and one fold.
2. Fit one full-data fit per count and setup, smallest first, recording wall time,
   iterations and convergence. Report a projected time for the CV at each count.
   If a count's projected five-fold CV for both setups exceeds 3 hours, do not run
   its CV now: record it as deferred with the timing and continue with the others.
3. Five-fold CV for every count that converged and fits the time limit, plus the
   expression comparator.
4. Describe each converged full fit: programs with nonzero loadings, survival-active
   programs (|E[β_k]| ≥ 1.96 posterior SD), and for each methylation program the
   largest absolute correlation of its patient scores with any expression-program
   score. A program is "methylation-only" if that correlation is below 0.5 and it is
   survival-active in the full fit and in at least three of five folds.
5. Score each fixed full fit on the matched ICGC set, as description only. Never use
   these results to choose a count or setup.

**Reporting.** Primary metric: CV concordance per fold and mean, with the paired
difference against expression and its SE. Also convergence, run time and memory by
count, nonzero and survival-active program counts, methylation-only programs, and
ICGC concordance labeled as external description. At n ≈ 144 everything is
descriptive. If a fit fails to converge, report and save it; do not tune settings.

**Compute.** At most 3 parallel R workers, each with `OMP_NUM_THREADS=1`. Run the
100k and all-CpG fits one at a time. Task A is using 5 cores on this machine.

**Files and rules.**

- Patient-level matrices, fits and scores go under the ignored `data/` directory.
  Commit only aggregate CSVs and figures under
  `results/multimodal_real/outputs/methylation_single_modality_2026-11-12/`.
- Write `docs/findings/2026-11-12-methylation-single-modality-findings.md` and add
  a row to `docs/findings/README.md`.
- Do not edit DECISIONS.md, PROJECT_STATUS.qmd, ROADMAP.md, the week plan or
  `docs/plans/2026-10-16-round2-attempts-log.md`. Another session owns them. Put the
  proposed attempts-log entry and DECISIONS text in a final section of your findings
  note titled "Text for the shared status files".
- Stage only your own paths (`git add <paths>`, never `git add -A`). Commit with a
  detailed message and no tool attribution, then push.
- Run `Rscript tests/run_tests.R` once at the end if you changed anything under
  `code/`, and report the count.

**Done when.** Every count has a full-fit row (converged, capped, failed or
deferred with timing), CV is complete for every count that qualified, the findings
note answers the goal with the caveat that β is Normal pending task B, and the work
is committed and pushed. Stop and report; do not start other tasks.

---

## Prompt 2: Task D steps 1 and 2, then task E code (Codex session 2, sequential)

Work in `~/GithubProjects/multiomicsGEP` on `main`. Do not create a worktree or
change branches. Read AGENTS.md, then in `docs/plans/2026-11-09-week-plan.md` the
"Reporting rules for all tasks", "Task D" and "Task E" sections, and the 2026-11-12
entry at the top of DECISIONS.md (task A fits pooled training data, no strata, with
within-cohort concordance as the primary CV metric). Do D first, then E. No long
fits in this session.

### Part 1: Task D, steps 1 and 2

1. Read `code/fit_multimodal_yfb.R`, `code/multimodal_yfb_updates.R` and
   `code/fit_ssbmf_yfb.R`, and the "Extension to M Modalities" section of
   `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`.
2. Write the test before any fitter change: `tests/test_multimodal_m1_equivalence.R`.
   With M = 1, the multimodal fitter must match `fit_ssbmf_yfb(canonicalize = TRUE)`
   under identical fixed options (including `center_projection`) in E[L], E[F], E[β]
   and the final ELBO, to relative tolerance 10⁻⁶, on two small fixtures: one pooled
   and one with `strata_id`. Run it against the current code and record the result.
3. If it fails, find where the two fitters differ (update order, canonicalization,
   w_S, centering, KL terms, stopping rule) and list each difference. Small
   alignments to the multimodal fitter are in scope. **If matching would need a
   rewrite around the single-modality sweep, stop and write up the options for
   Andrew** instead of rewriting.
4. Step 2: check the multimodal fitter's KL helpers against the two November 6
   single-modality fixes (the scalar-precision normalizer and the exactly flat β
   likelihood; see DECISIONS.md 2026-11-06 and `tests/test_empty_program_kl.R`).
   Either show by test that they agree or replace them with the shared helpers.
5. Run `Rscript tests/run_tests.R` and report the count.

Write `docs/findings/2026-11-12-multimodal-m1-equivalence-findings.md`. Do not run
D steps 3 and 4; they wait for tasks A to C.

### Part 2: Task E, code only

Build what task E needs so its fits can run once A and B choose the priors.

1. A five-fold CV loop for LB that scores each held-out fold by fold-in
   (`predict_ssbmf_lb_foldin()` in `code/predict_ssbmf_lb.R`), and for YFB by
   projection. Fits are pooled (no `strata_id`). Scoring reports within-cohort
   concordance (pairs only within TCGA or within CPTAC, pooled by pair count, as
   `concordancefit()` with `strata` computes) as primary and pooled concordance
   beside it. Orientation is fixed on training data as in `frozen_reverse_cindex()`.
2. A runner, extending `results/benchmark_sim/run_lb_yfb_matched_comparison.R` or a
   new `results/benchmark_sim/run_lb_yfb_study.R`, that takes the input scaling,
   `prior_LF` and `prior_beta` as arguments. It builds the input once, passes the same
   object to both fitters, and stops if the two inputs, folds or starting values
   differ. It also records the correlation of the two linear predictors and the
   per-program correlations between L_k and (YF)_k.
3. Tests for the CV loop and the identical-input check (a test that the check stops
   on a mismatch), then a smoke run on a tiny synthetic input only.
4. Run `Rscript tests/run_tests.R` and report the count.

Write `docs/findings/2026-11-12-lb-yfb-study-code-findings.md` with the exact
command to run E's simulation and training-CV fits once A and B are decided.

### Rules for both parts

- Script headers and roxygen per AGENTS.md, author Andrew Walther, no tool names.
- Do not edit DECISIONS.md, PROJECT_STATUS.qmd, ROADMAP.md, the week plan or the
  attempts log. Put proposed text for them in a final "Text for the shared status
  files" section of each findings note. Add each findings note to
  `docs/findings/README.md`.
- Another Codex session is writing `results/multimodal_real/run_methylation_single_modality.R`
  and a Claude session `results/benchmark_sim/run_lf_prior_comparison.R`. Do not
  touch them. Stage only your own paths, commit each part separately with a
  detailed message and no tool attribution, and push.
- At most 2 CPU workers for tests and smoke runs.

**Done when.** D's equivalence test exists and its result is reported (passing,
fixed, or stopped with options), the KL check is resolved, E's CV code, runner and
tests are committed with a passing full suite, and both findings notes are pushed.
Stop and report.
