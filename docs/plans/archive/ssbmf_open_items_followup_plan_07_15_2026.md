# SSBMF — Follow-up Plan for Open Items (post-Item-2)

**Date:** 2026-07-15
**Purpose:** A saved "next steps" backlog for every open item surfaced by Item 2's progress report
(`docs/reports/ssbmf_progress_report_07_15_26.qmd`, §7 "Honest open items and shortfalls") and its
preceding review passes (`DECISIONS.md`, 2026-07-15 entries). Not yet approved for execution — this
is a menu to pick from in a future session, not an in-progress plan. Nothing here should be started
without picking it up explicitly and confirming scope first (per the project's standing convention:
TDD for new functions, independent review before merge, tests green after any code change, one
branch per item, never push without asking).

**How to use this document:** each item has a goal, a concrete approach, an effort estimate, and a
success criterion, so a future session (or a fresh Claude Code instance) can pick one up without
re-deriving the plan. Items are grouped by priority tier; within a tier, order is not significant
unless a dependency is noted.

---

## Priority 1 — validity gaps most likely to be asked about next

These two are the most likely first questions from a statistically literate advisor. Neither
requires re-fitting the model — both are analyses layered on top of already-fitted results.

### 1a. Per-cohort concordance confidence intervals — **COMPLETE 2026-07-16**

**Done, branch `external-cindex-cis`.** Implemented as a percentile bootstrap (option 2 below), not
the analytic option — the user also asked for a paired significance test against the two-step
baseline, which needs the bootstrap machinery anyway (the analytic concordance variance doesn't
directly give a paired covariance between two different risk scores on the same patients). Also
surfaced and fixed a second stale-baseline finding along the way (the EBMF→Cox two-step comparator's
cached output predated the Phase 1c preprocessing fix; refreshed 0.564→0.581). New tested utility
`code/concordance_ci.R` (12 tests, TDD); results in
`results/benchmark_sim/outputs/desurv_comparison/{external_cindex_ci.csv,external_paired_diff_ci.csv}`.
Full results, methodology, and the paired-comparison finding (recommended model significantly more
concordant than the two-step baseline, pooled across cohorts: +0.042, 95% CI 0.013–0.071):
DECISIONS.md 2026-07-16. Folded into both progress-report documents.

**Original goal (for reference):** replace the bare point-estimate C-index per external cohort with a proper interval, and a
test against chance (0.5) — especially for the weakest cohort (Moffitt microarray, C=0.549), whose
separation from chance is not currently established.

**Approach:** the recommended fit's per-patient risk scores for each external cohort are already
computed inside `results/benchmark_sim/run_desurv_comparison.R` (Section 6, external validation) —
this is a scoring step on existing data, not a re-fit. Two options, in order of preference:
1. **Analytic:** `survival::concordance()` returns a variance estimate directly (`$var`); use it to
   build a Wald CI and a one-sided test against 0.5 for each cohort. Fastest, no new dependency.
2. **Bootstrap (if the analytic variance looks unstable at small cohort n, e.g. Moffitt):**
   resample patients within each external cohort with replacement (fixed seed, e.g. 1000 reps),
   recompute concordance each time, report the empirical 2.5/97.5 percentiles. Needs `boot` or a
   manual loop; more defensible for small/skewed cohorts, more compute (still trivial — no CAVI
   refitting involved, purely re-scoring cached risk predictions).

**Effort:** Small (a half-day). No new model code; extends the existing benchmark runner's Section 6
and the results CSV schema (add `c_index_lo`, `c_index_hi` columns).

**Success criterion:** every cohort in `desurv_comparison_results.csv` has a CI; new unit test
confirms the CI construction against a known `survival::concordance()` reference case; the progress
report's §5 table and figure are updated to show error bars instead of bare bars.

### 1b. Factor-stability under resampling

**Goal:** show that the two survival-active gene programs (Program 7 adverse, Program 3 protective)
are not an artifact of one particular training draw — that refitting on a bootstrap or
cross-validation split of the training cohort recovers essentially the same two programs.

**Approach:** refit the recommended configuration (YFB, DeSurv-aligned gene selection, K=7, no
cohort indicator) on B bootstrap resamples of the pooled training cohort (TCGA-PAAD + CPTAC,
n=273), or on the existing 5-fold CV splits already used for K-selection (cheaper, reuses existing
fold-assignment code in `code/select_K.R`). For each resample: (1) identify the survival-active
factor(s) by the same $|\hat\beta|$-thresholding rule already used (`beta_threshold` in
`config/globals.yml`), (2) compute the gene-weight correlation between each resample's active
factor(s) and the reference fit's Program 7/Program 3 columns (Pearson or cosine similarity on
$\mathbf{F}$), (3) report the distribution of these correlations and the fraction of resamples where
exactly 2 factors are recovered as active.

**Effort:** Medium (1-2 days) — this is genuinely a set of new model fits (B$\approx$50-100 refits
at the existing K=7, not a cheap re-scoring step like 1a), though each individual fit is fast
(the existing pipeline converges in well under a minute per fit at this n/p).

**Success criterion:** report what fraction of resamples recover 2 active factors matching Program
7/3 (by gene-weight correlation above some threshold, e.g. $|r|>0.7$); if instability is found,
document it honestly rather than treating it as a negative result to hide — either outcome is
reportable, per the project's existing "neutral/negative is an acceptable outcome" convention.

---

## Priority 2 — documented robustness/compliance gaps, not currently blocking

### 2a. YFB `point_normal` (spike-and-slab) K-CV collapse

**Goal:** either fix or formally close out the collapse to $\hat\beta=0$ under the point_normal
prior in cross-validation (C=0.5 for every K).

**Approach (per the existing ROADMAP.md note):** try a $\hat\beta$-only warm-up (a few iterations
holding $\mathbf{L}$/$\mathbf{F}$ fixed while $\hat\beta$ moves under the point_normal prior before
the spike component can lock in), or an adaptive/scheduled spike weight (softening the mixture
prior's spike probability early in CAVI, the same `alpha_schedule`-style ramp idiom already used
elsewhere in this codebase). If neither works within a reasonable effort budget, document as a
formally closed, understood limitation of point_normal specifically for YFB's reduced-fold training
size, rather than leaving it as a perpetually open question.

**Effort:** Medium. **Priority:** Medium-Low — the recommended configuration already uses the normal
prior and is unaffected; this only matters if point_normal's sparsity is wanted for a future
manuscript emphasis on the point_normal prior.

### 2b. General (non-warm-start) fix for the CAVI factor-collapse vulnerability

**Goal:** a from-scratch initialization or damping strategy that avoids the degenerate near-zero
fixed point under near-symmetric, disjoint-support factor loadings, without depending on a
pre-existing higher-$K$ fit (the current escape, warm-start, requires one).

**Approach:** deflation-init was tried and analytically shown not to help (it matches plain SVD-init
whenever singular values are non-degenerate — DECISIONS.md 2026-07-13, Step 2). Candidate directions
not yet tried: (1) a trust-region / damped update on the EBNM shrinkage step specifically in the
early iterations (limiting how far a factor's precision can move per iteration, giving borderline
factors more chances to recover before being permanently shrunk); (2) an explicit tie-breaking
perturbation when two factors' updated precisions are within some tolerance of each other (directly
targeting the "near-tied-amplitude" degenerate condition, rather than changing the initial subspace).

**Effort:** Large — this touches the CAVI update derivation directly and needs new theory + TDD, not
a parameter tweak. **Priority:** Low-Medium — no evidence this affects the real PDAC fits reported
here; it is currently an occasional issue only in adversarial/lower-SNR synthetic settings for LB.
Revisit if it is ever observed on real data, or before a manuscript claims general robustness.

### 2c. Joint $(K,\alpha)$ Bayesian-optimization re-run with a fixed acquisition/search issue

**Goal:** re-run the Step 3 joint-tuning search (`code/select_k_alpha_bo.R`) with the identified fix
(avoid the integer-bounds convention that pre-rounds $K$ before the Gaussian process sees it — e.g.
a continuous relaxation of $K$ rounded only inside the objective function, not in `BayesianOptimization()`'s
own bounds) and a larger `init_points` budget, to more thoroughly rule out whether a smaller $K$ is
reachable via joint tuning (the current run's search collapsed to ~9 unique evaluated points, not
23, per DECISIONS.md 2026-07-13 Step 3).

**Effort:** Small-Medium (mostly a re-run with a different bounds convention plus more init points;
existing runner and objective function need only a small modification, not new derivation).
**Priority:** Low — this is a completeness/rigor question about the search, not a known problem with
the current recommended $K=7$; relevant mainly if a manuscript needs a stronger claim that smaller
$K$ was thoroughly ruled out (vs. the current honest "not found by this specific, limited run").

### 2d. Review of the DeSurv paper's own simulation construction ("Amber's simulation")

**Goal:** close the 6/18 feedback item asking for a review of how the DeSurv authors' own simulation
was constructed, to compare against this project's independently-built multi-cohort simulation
(`results/multi_cohort_sim/`) — a request never formally actioned (DECISIONS.md 2026-07-15,
compliance-review entry).

**Approach:** read the DeSurv paper's methods/supplement simulation section (`reference-desurv-paper`
memory has the paper's location) and, if available, the DeSurv GitHub repo's simulation code; write
a short comparison note (not new code) covering what data-generating assumptions the two simulations
share vs. differ on, and whether any differences bear on the interpretation of this project's own
`sparse_synthetic` collapse finding or the general joint-vs-two-step conclusions.

**Effort:** Small (a focused reading + writeup session, no code). **Priority:** Medium — this is a
direct, specific ask from the meeting that has gone unaddressed the longest of anything on this list.

---

## Priority 3 — infrastructure / bigger future work (not urgent)

### 3a. HPC / Longleaf cluster workflow

**Goal:** close the 6/18 "figure out HPC cluster workflow" ask — establish how the existing
benchmark runners (fits, K-CV, the multi-cohort simulation, pathway enrichment) would actually be
submitted and run on Longleaf (SLURM), including memory/parallelization notes per this project's
global HPC conventions.

**Approach:** audit which existing runner scripts are the ones that would actually need to run at
scale on Longleaf (likely the K-CV grid searches and the multi-cohort simulation sweeps, which are
the most compute-heavy, embarrassingly-parallel-by-seed/fold jobs); write SLURM array-job wrappers
for those; document memory requirements based on observed local run behavior.

**Effort:** Medium. **Priority:** Low unless/until an actual compute need (e.g. the resampling study
in 1b at scale, or a much larger simulation grid) makes local runtime a bottleneck.

### 3b. Matched, same-protocol DeSurv head-to-head

**Goal:** the comparison this project has consistently and correctly declined to claim yet — running
DeSurv's own joint $(k,\alpha,\lambda)$ Bayesian-optimization tuning procedure and this project's
model on identical data with identical hyperparameter-search rigor, to make an apples-to-apples
performance claim rather than "the Bayesian counterpart, comparison pending."

**Effort:** Large — likely the single largest remaining piece of work in this project's arc, closer
to a manuscript-scale undertaking than an "item" on this list. **Priority:** revisit when a
manuscript timeline requires it; not needed for the current progress-reporting purpose.

### 3c. Dev-tooling backlog (unchanged from the original 6/18 ask)

- Post-commit Codex review hook (automated review triggered on every commit).
- RAG literature database (download+index papers for retrieval instead of repeated web search).

**Effort:** Small-Medium each. **Priority:** Low — genuinely backlog, not blocking any scientific
conclusion; action only if/when the user wants to invest in the workflow tooling itself.

---

## Summary table

| Item | Priority | Effort | Blocks anything now? |
|---|---|---|---|
| 1a. Per-cohort concordance CIs | **DONE 2026-07-16** | Small | — |
| 1b. Factor-stability under resampling | **High** | Medium | No, but a natural reviewer ask |
| 2d. Review DeSurv's own simulation | Medium | Small | No — longest-outstanding unaddressed ask |
| 2a. point_normal K-CV collapse | Medium-Low | Medium | No — normal prior is the default |
| 2c. Re-run joint (K,α) BO with fixed search | Low | Small-Medium | No |
| 2b. General CAVI-collapse fix | Low-Medium | Large | No evidence it affects real fits |
| 3a. HPC/Longleaf workflow | Low | Medium | Only if compute need arises |
| 3c. Dev-tooling backlog | Low | Small-Medium | No |
| 3b. Matched DeSurv head-to-head | Revisit at manuscript stage | Large | No |
