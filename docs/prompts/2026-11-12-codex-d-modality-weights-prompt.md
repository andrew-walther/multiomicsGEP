# Codex prompt: task D, balancing expression and methylation with modality weights

Date: 2026-11-12. Author: Andrew Walther.

This continues [the D next-steps prompt](2026-11-12-codex-d-next-prompt.md). Its results
are in [the TCGA findings](../findings/2026-11-12-multimodal-tcga-comparison-findings.md)
and the [feature-group findings](../findings/2026-11-12-feature-group-fitter-findings.md).
Andrew approved this plan on 2026-11-12. The "Rules for both parts" section at the end of
[the D/E follow-up prompt](2026-11-12-codex-d-e-followup-prompt.md) applies unchanged.

## Why

**What the joint model is for.** Andrew's aims for the multimodal model are interpretive,
not a gain in prediction:

- find which genes and CpGs load on the same programs;
- find programs that load on one modality only (methylation-exclusive or
  expression-exclusive) as well as joint programs;
- use a large set of CpGs (the 182,934 with mean β in [0.1, 0.9]), not a set cut down to
  the size of the gene list.

**The problem is that the larger modality takes over.** With equal feature counts (50:50
in simulation), the joint fits recovered the prognostic program as well as or better than
expression alone (9 vs 5 of 10 in one cell). With five times as many methylation features
(50:250), they did worse (4 vs 7 of 10, false activity 1.6–1.7). On TCGA with 10,000 CpGs
against 2,218 genes, methylation supplied 83–97% of the L precision in three of five
programs, and the joint CV C was 0.06 below expression alone.

**Task C showed why a large CpG set is needed.** Methylation-alone CV C does not rise
beyond 10k CpGs (0.545 at 10k, 0.50–0.52 from 25k to all). But the clearest
methylation-specific factor, an immune-related program that expression does not carry,
forms as a single factor only from about 50k CpGs. It has R² 0.06 on the expression
factors at 100k. Andrew chose 10k as the first joint fit and 100k for this factor. At
100k, methylation has about 45 times the features of expression.

**Why the F prior cannot fix this.**

- *The prior is estimated, not set.* Each modality's F prior is estimated from its own
  data by empirical Bayes. On TCGA both fitted priors were dense (point-mass weight about
  10⁻⁷), so they shrink nothing. Fixing the methylation prior by hand would abandon
  empirical Bayes and would set sparsity, not balance.
- *Rescaling the data cancels out.* Because τ is feature-specific and estimated, scaling
  a modality's data by c only rescales τ by 1/c² and F by c. The fit is unchanged.

**The fix is a power weight on each modality's reconstruction likelihood.** This is the
same device as w_S on the survival term: log p(Y_m | L, F_m, τ_m) is multiplied by w_m.
It is also the multimodal analogue of block weighting in multiple factor analysis, where
each block is scaled so that no block dominates the leading dimensions.

## Part 1: derivation, before any code

Add a "Modality weights" subsection to `\section{Extension to $M$ Modalities}`
(line ~1062 of `derivations/SSBMF-YFB-model-derivations/SSBMF-YFB-model-derivation-full.tex`).
Model it on the `w_S` keybox near line 232, and cover the following.

1. **The objective.** Each reconstruction term is weighted by w_m, while the priors and KL
   terms keep weight 1. It is the bound for a power posterior, not the original posterior.
2. **The updated CAVI equations.**
   - *L:* the precision becomes Σ_m w_m Σ_j τ_mj E[F²_mjk], and the mean changes to match.
   - *F_m:* the likelihood precision w_m τ_mj Σ_i E[L²_ik] and the matching signal
     enter the ebnm pseudo-observations.
   - *τ_mj:* state the update, and whether w_m cancels from its point estimate.
   - *Survival terms:* unchanged. Say explicitly whether Z = Σ_m Y_m F_m needs any change.
     It should not, since the weights act on reconstruction only.
   - *Canonicalization:* show that the rescaling c_k leaves the weighted objective
     invariant.
3. **Two weight rules, both relative to expression (w_expression = 1):**
   - *feature count:* w_m = p_expression / p_m. If methylation has q times as many
     features as expression, its weight is 1/q. This is Andrew's proposal, and it aims to
     reproduce the balance of the 50:50 simulation;
   - *leading eigenvalue (multiple factor analysis):* w_m = λ₁(expression) / λ₁(m), where
     λ₁ is the largest eigenvalue of the training-fold z-scored block's X′X / n.
     Methylation features are highly correlated, so this gives a different weight from
     the count rule.
4. **The known failure mode.** Dividing the genomics term by p collapsed L and F in LB
   (DECISIONS.md, 2026-07-12; derivation line ~1002). A small w_m makes the prior
   relatively stronger, so a methylation block may refit to a point mass and empty out.
   Say how that shows up (empty-block counts) and that it must be reported, not
   prevented.
5. **Checks.** Add checks for the weighted updates to `verify_derivation.R` in that folder
   and run it.

Get an independent review of the derivation (statistical logic) before Part 2. Commit.

## Part 2: implementation

- Add `modality_weight = NULL` to `fit_ssbmf_yfb()`: a named positive numeric vector keyed
  by the `feature_group` labels. NULL or all ones must reproduce the current grouped fit
  bit for bit. It requires `feature_group`; error otherwise.
- Add a helper that computes both weight rules from training-fold data only.

Tests first:

- all-ones weights equal the unweighted grouped fit bit for bit;
- the L precision equals the derived formula on a fixture;
- the recomputed weighted objective matches the recorded value;
- multiplying a block's data by a constant leaves the weighted fit unchanged (τ absorbs
  it);
- prediction is unchanged in form.

Then run the full suite. Commit.

## Part 3: simulation

Use the runner from the next-steps prompt, its survival regeneration, recovery rule
(S-9 convention, option 1 on false activity), working configuration, K = 5, seeds 1–10
and 120 training patients.

**Block A, balance (`both_informative`).** Three cells: 50:250 with 15% and with 3%
methylation nonzero (where the joint fits failed), and 50:50 with 3% as a control where
they passed. Four arms:

1. grouped, unweighted;
2. grouped, feature-count weights;
3. grouped, eigenvalue weights;
4. expression only.

That is 120 fits.

**Block B, joint and modality-exclusive programs (the main test of Andrew's aims).** Add a
generator option (default off, so existing scenarios are unchanged) that plants:

- program 1 on expression only, by zeroing its methylation loadings;
- program 2 on methylation only, by zeroing its expression loadings;
- program 3 on both.

Use β = (0.70, 0.35, 0) as in `both_informative`, so the expression-only and
methylation-only programs are prognostic. Two cells, 50:50 and 50:250 with 15% nonzero,
and four arms, for 80 fits:

1. grouped, feature-count weights;
2. grouped, eigenvalue weights;
3. pooled, feature-count weights;
4. grouped, unweighted.

For each true program and arm, report:

- *Loading recovery:* the matched fitted program, with |cor| ≥ 0.7 against its standardized
  loading over all features. Program 3 is not prognostic, so this is loading recovery only.
- *Modality pattern:* whether the matched program has the right modality pattern. For
  programs 1 and 2, the absent modality's block is empty or carries under 5% of the
  program's L precision. For program 3, both blocks carry at least 5%.
- *Survival activity and sign* for programs 1 and 2.

A seed counts as a **full pattern success** when all three programs are recovered with the
right modality pattern. Report full pattern successes out of 10 per arm and cell. This is
the case where grouped and pooled priors should differ, because only a grouped prior can
shrink one modality's block for one program without shrinking the other.

**Registered rules, set before any fit.**

- *Gate:* unchanged (recovery ≥ 8/10, mean false activity ≤ 1, convergence ≥ 9/10).
- *When a weight rule counts as effective:* it passes the gate in at least one 50:250
  cell where the unweighted grouped arm fails, and it still passes the 50:50 control.
- *Selecting a rule:* if only one rule is effective in Block A, select it. If both are,
  select the one with more full pattern successes in Block B, summed over its two cells.
  Break ties first with program-1 recoveries over the two 50:250 cells in Block A, then
  with the feature-count rule. If neither rule is effective, report that and do not run
  Part 4.
- *Separate priors:* report whether grouped beats pooled on full pattern successes under
  feature-count weights. That is the evidence on whether separate priors serve their
  purpose.

Report empty-block counts and precision shares in every arm. Use new output directories,
for example `feature_group_simulation_2026-11-12_weights`. Commit the aggregates.

## Part 4: TCGA at 10k and 100k CpGs, descriptive

Run only if a weight rule was selected. Use task C's 144 matched patients, folds, nested
CpG ranking (10k and 100k), point-Laplace setup and D4 genes, read-only, with training-fold
scaling. The weights are computed on each training fold only. K = 5, working
configuration. Three configurations, each with five folds plus a full-data fit, for 18
fits:

1. grouped, selected weights, 10k CpGs;
2. grouped, selected weights, 100k CpGs;
3. grouped, unweighted, 100k CpGs. The unweighted 10k fit already exists from Part A of
   the next-steps prompt.

The readouts are listed in order of importance. Andrew's interest is in what the
methylation data add to the programs, not in prediction. Label the results **not validated
by simulation** where the gate was not passed for this setting.

- *Balance:* the weights used, each modality's L-precision share per live program, and
  empty-block counts by modality and program.
- *Modality-exclusive programs:* programs whose expression or methylation block is empty
  or under 5% of the precision.
- *Programs new to the joint fit:* programs that match no program in the expression-only
  or methylation-only full fits (|correlation| < 0.5). Also whether the 100k joint fit
  recovers task C's methylation-specific factor (point-Laplace factor 2 at 100k; patient
  score correlation).
- *Genes and CpGs together:* for each full-data joint program, the top-20 genes and the
  top-20 CpGs with their island and promoter annotation. Count the top CpGs annotated to
  a promoter of one of the program's top-50 genes, and the sign of each such pair's
  loadings. Expression and methylation are expected to be anti-correlated at promoters.
- *Prediction:* CV C per fold and the mean, with paired differences against expression
  alone (from the next-steps run) and SEs. This is descriptive only. Prediction is not
  expected to improve.
- *Cost:* convergence, iteration caps and run time. Task C's 100k methylation fits took
  7–15 minutes each.

Patient-level data stay under ignored `data/`. Commit only aggregates.

## Order and commits

1. Derivation and its checks, then independent review. Commit.
2. Tests (failing first), then implementation and the full suite. Commit each.
3. Simulation smoke run, full run and aggregates. Commit.
4. TCGA (if allowed). Commit aggregates.
5. Findings in a new note, `docs/findings/2026-11-12-modality-weights-findings.md`, with
   a "Text for the shared status files" section and a row in `docs/findings/README.md`.
   Independent review of the code and findings, then the final commit and push.

**Done when.** The derivation is reviewed and its checks pass. The weighted fitter is
tested, and the full suite passes. Both simulation blocks are reported against their
rules. Part 4 is reported or marked not run with the reason. The review is recorded.
