# Codex prompt: task D rerun with the empty-block fix and two diagnostics

Date: 2026-11-12. Author: Andrew Walther.

This continues [the D finishing prompt](2026-11-12-codex-d-finish-prompt.md), whose
registered simulation failed every cell (`b90d73ae`,
[findings](../findings/2026-11-12-feature-group-fitter-findings.md)). Andrew approved
this plan on 2026-11-12. The "Rules for both parts" section at the end of
[the D/E follow-up prompt](2026-11-12-codex-d-e-followup-prompt.md) applies unchanged.

The failed run had two separate problems. This rerun fixes the first and adds the
measurements needed to diagnose the second. It keeps the scenario, the cells, the seeds,
the working configuration, the survival regeneration and the gate from the finishing
prompt, so the results are directly comparable.

## Problem 1: the grouped F update crashes when a program is empty in one modality

In `yf_update_F_groups()` (`code/feature_groups_yfb.R`), the sequential sweep computes a
block posterior q under the incoming prior and then refits the prior. When the refit
returns a pure point mass at zero while q still has slab mass, KL(q || g_new) is
infinite and the fit stops. This caused 38 of 40 grouped failures. A program that is
active in expression and absent in methylation is exactly what separate priors are for,
so this case must be handled, not avoided.

**Decision (Andrew, 2026-11-12).** When a group's refitted prior for program k is a pure
point mass at zero, treat that block as empty for program k, as the existing dead-block
branch already does:

1. Remove the block's contribution from E[Z_k] and V[Z_k] using the moments the sweep
   just returned, before visiting the next group.
2. Set E[F] and E[F²] (and VF, if carried) to zero on that block.
3. Store the point-mass prior for that group and program. The block's loading KL is 0,
   since q and g are then the same point mass.
4. Count these events per fit and program and return the count in the fit object
   (for example `fit$empty_blocks`), so they are reported rather than hidden.

Keep the explicit error for any other nonfinite KL. Do not change the NULL or one-group
path: the frozen bitwise fixtures must still pass. If the fix cannot be made without
changing them, stop and ask.

**Tests first:**

- A two-group fixture in which one group has no signal for one program. The fit
  completes, that block's loadings are exactly zero, its stored prior is a point mass,
  the other group's loadings for the same program stay nonzero, and the empty-block count
  is positive.
- The monitored objective recomputed from the returned moments and priors matches the
  recorded value on that fixture (extend the existing objective check).
- E[Z_k] after the update equals Y F_k from the returned moments (no stale contribution
  from the zeroed block).
- The existing 14 fitter checks and the bitwise fixtures still pass.

## Problem 2: the pooled comparator rarely recovers the program

The pooled arm, which is the single-modality model on the concatenated data, recovered
2/10 and 4/10 at 50:50 and 0/10 at 50:250. Two explanations fit this: the scenario is too
weak to recover from (low-variance program, 120 training patients), or the extra
methylation features take over L as the feature count grows. Add two measurements:

1. **True-risk C.** For every seed, the validation C of the true linear predictor η used
   to draw the outcomes. This bounds what any fit can reach.
2. **Expression-only arm.** Fit the expression block alone (`feature_group = NULL`) on the
   same seed, the same scaled expression columns, the same survival draw and the
   matching start columns. Score recovery against the expression part of g, renormalized
   to unit norm, with the same rule (|cor| ≥ 0.7, survival-active, correct sign). For the
   grouped and pooled arms, also report the matched program's correlation with this
   expression part, so all three arms are compared on the same target.

## Runs

4 cells (feature ratio 1:1 and 1:5 × methylation nonzero 15% and 3%) × 3 arms (grouped,
pooled, expression-only) × seeds 1–10 = 120 fits, at most 2 workers. Write to new output
and private directories (for example `feature_group_simulation_2026-11-12_rerun`). Do not
touch the first run's directories. Report per cell and arm: recovery out of 10, mean false
activity, converged out of 10, errors, empty-block events, mean validation C and mean
true-risk C (with defined counts), and the precision shares as before.

## Registered rules, set before any fit

- **Gate:** unchanged (recovery ≥ 8/10, mean false activity ≤ 1, convergence ≥ 9/10).
- **Which prior goes to step 4:** unchanged. Grouped, unless pooled passes the gate in a
  cell where grouped fails. Step 4 (the matched-TCGA comparison in the finishing prompt)
  runs only if the selected arm passes the gate in at least one cell.
- **Reading problem 2, descriptive and stated in the findings:**
  - If expression-only passes recovery in a cell where both joint arms fail, adding
    methylation hurts recovery. The next step would be a derivation of modality weights
    on the reconstruction terms. Do not add weights in this task.
  - If expression-only also fails, the scenario is too weak at this sample size for the
    gate to separate the arms. The next step would be choosing a scenario where it can
    pass. Do not change the scenario in this task.

## Order and commits

1. Tests for the fix and the new arm (fail first). Commit.
2. Fitter fix, runner changes, one-seed smoke run per cell, full suite. Commit.
3. Full 120-fit run. Commit aggregates.
4. Step 4 only if the rule allows. Commit aggregates.
5. Add a dated rerun section to the existing findings note, with a revised "Text for the
   shared status files". Independent review of the code and findings, then the final
   commit and push.

**Done when.** The fix is tested and the bitwise fixtures pass; every cell and arm is
reported against the gate; the problem-2 reading is stated; step 4 is reported or marked
not run with the reason; the full suite passes; the review is recorded.
