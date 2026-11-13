# Codex prompt: task D weight path (find a methylation weight between dominance and collapse)

Date: 2026-11-12. Author: Andrew Walther.

This continues [the modality-weights prompt](2026-11-12-codex-d-modality-weights-prompt.md).
Its results (`e51e03c3` to `1161999b`) are in
[the modality-weights findings](../findings/2026-11-12-modality-weights-findings.md).
Andrew approved this plan on 2026-11-12. The "Rules for both parts" section at the end of
[the D/E follow-up prompt](2026-11-12-codex-d-e-followup-prompt.md) applies unchanged.

**This is the last task D run for now.** Andrew is pausing the multimodal work after it,
to prepare for the advisor meeting and the prelim submission. Report the results and stop.
Do not start follow-up runs or propose new studies inside this task. Put suggested next
steps in the findings note only.

## What the audit of the last run found

The implementation matches the specification. The weight multiplies τ only in the
reconstruction coordinates, and the survival score still uses Σ_m Y_m F_m. No run required
by the prompt is missing. The scaling test scope is adequate.

**The 1/q weight overcorrects.**

- *What happened.* At 50:250 with methylation weight 0.2, both rules set almost every
  methylation block to zero: the methylation precision share was 0 in every program. The
  weight lowers methylation's total influence on L, as intended. But it also lowers each
  CpG's evidence for its own loading by the same factor, so empirical Bayes refits the
  methylation prior to a point mass and the block empties.
- *Why recovery read 0/10.* The true prognostic direction is mostly methylation features.
  With an empty methylation block, the correlation with the full truth cannot reach 0.7.
  The 0/10 follows from the collapse. It does not mean the expression part was missed.
- *What an unweighted fit does instead.* At weight 1, methylation dominates L.

**Block B had no full-pattern success in any arm,** but the per-program results favor
separate priors. At 50:50, where both arms are unweighted, grouped beat pooled on loading
plus modality pattern:

| Program | Grouped | Pooled |
|---|---|---|
| Expression-only | 5/10 | 1/10 |
| Methylation-only | 2/10 | 0/10 |
| Joint | 7/10 | 5/10 |

The methylation-only program (β = 0.35) is hard to recover at n = 120 in every arm.

## Question

Is there a methylation weight between 1 (dominance) and 1/q (collapse) at which methylation
keeps non-empty loadings on the true programs while no longer drowning out expression?

## Part A: simulation weight path

Use `results/benchmark_sim/run_modality_weight_simulation.R` and its fixed settings:
grouped priors, working configuration, K = 5, seeds 1–10, 120 training patients, the same
survival regeneration, signed-SVD starts and recovery rules.

- **Expression weight:** 1.
- **Methylation weight w:** 1, 0.7, 1/√5 (0.4472), 0.3 and 0.2. These are fixed values, not
  rules computed from the data.
- **Cells:**
  - Block A (`both_informative`): 50:250 with 15% methylation nonzero, and 50:250 with 3%;
  - Block B (the expression-only, methylation-only and joint scenario): 50:250 with 15%.
- **Total:** 3 cells × 5 weights × 10 seeds = 150 fits.
- **Reproducibility check:** the w = 1 and w = 0.2 arms repeat the previous unweighted and
  count-weighted grouped arms. Their aggregates must match the earlier run exactly. If they
  do not, stop and report before interpreting anything.

Report for every cell and weight:

- primary recovery (Block A), mean false activity and convergence, as before;
- **block-wise correlations** for the program matched to each true program: the fitted
  expression block against the expression part of the truth, and the fitted methylation
  block against the methylation part, each renormalized. This shows whether a fit found the
  expression part even when its methylation part collapsed;
- the methylation precision share of each matched true program, and the number of seeds in
  which its methylation block is empty or under 5% of its precision;
- final empty-block counts by modality;
- Block B: loading recovery, modality pattern, and loading plus pattern for each of the
  three programs.

**Registered rule, set before any fit.** A weight w < 1 is **usable** if all of the
following hold:

1. *No loss in Block A.* In each Block A cell, its primary recovery is at least the w = 1
   count, and its mean false activity is at most 1.
2. *Methylation stays in.* In each Block A cell, the matched primary program's methylation
   block carries at least 5% of its precision in at least 8 of 10 seeds.
3. *No loss in Block B.* Its loading-plus-pattern count for the expression-only program
   and for the methylation-only program is each at least the w = 1 count.
4. *A gain somewhere.* It is strictly better than w = 1 in at least one of: summed Block A
   primary recovery, or summed Block B loading plus pattern over the three programs.

If several weights are usable, select the one with the largest summed Block B loading plus
pattern. Break ties with the larger w, the one closer to unweighted.

If none is usable, say so. Then describe where along the path the methylation blocks
collapse, using the precision shares and empty-block counts. The next step would be a design
decision for Andrew, such as weighting methylation only in the L update. Do not implement any
such alternative in this task.

Use new output directories, for example `feature_group_simulation_2026-11-12_weight_path`.
Commit only aggregates.

## Part B: TCGA full-data weight path, descriptive

Andrew allows these descriptive full-data fits regardless of Part A's outcome. They choose
nothing. Label every result **not validated by simulation** unless Part A selects a usable
weight.

- **Data and settings:** task C's 144 matched patients, the 10k-CpG point-Laplace setup and
  D4 genes, read-only, scaled on all 144 patients as in the earlier full fits. K = 5, the
  working configuration, grouped priors.
- **Weights:** methylation weight w ∈ {0.7, 1/√5, 0.3, 0.22}. The value 0.22 is 2,218 / 10,000,
  the count ratio here. Reuse the unweighted full-data joint fit from the next-steps run as
  w = 1.
- **Scope:** four new full-data fits. No cross-validation and no external scoring.

Report for each w:

- each modality's precision share per live program, and empty blocks by modality;
- programs that are expression-exclusive or methylation-exclusive (the other block is
  empty or under 5% of the precision), and joint programs (both blocks at least 5%);
- each program's best absolute patient-score correlation with the expression-only and
  methylation-only full fits, flagging programs below 0.5 with both;
- for each exclusive or new program, the top-20 genes or CpGs, with CpG island and
  promoter annotation;
- convergence, iterations and run time.

Patient-level data stay under ignored `data/`. Commit only aggregates.

## Order and commits

1. Runner changes for the weight list, block-wise correlations and the reproducibility
   check, with tests first. Full suite. Commit.
2. Part A smoke run (one seed per cell and weight), then the full run. Commit aggregates.
3. Part B fits. Commit aggregates.
4. Findings in a new note, `docs/findings/2026-11-12-weight-path-findings.md`, with a
   "Text for the shared status files" section and a row in `docs/findings/README.md`.
   Independent review of the code and findings, then the final commit and push.

**Done when.** The reproducibility check passes. Every cell and weight is reported, and the
usable-weight rule is applied. Part B is reported for all four weights. The full suite
passes, and the review is recorded.
