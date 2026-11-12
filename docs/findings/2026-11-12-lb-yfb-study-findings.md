# Matched LB and YFB comparison at K = 5

Date: 2026-11-12. Author: Andrew Walther.

## Scope and registered choices

This study implements Part 1 of the [follow-up prompt](../prompts/2026-11-12-codex-d-e-followup-prompt.md).
The fixed inputs are D4 genes, per-cohort training-fold z-scores, point-Laplace
L and F, per-program Normal β, K = 5 and pooled fits without Cox strata.
Training uses w_S = 1. Simulation uses w_S = 1 and 10. No pruning is used.
The maximum iteration count and tolerance come from `config/globals.yml`.
The saved S-11 folds in `data/lf_prior_real_2026-11-12.rds` are retained.

## Pre-fit software corrections

`code/lb_yfb_study.R` returns within-cohort, pooled and individual-cohort
concordance, with each comparable-pair count. `survival::concordancefit()` uses
`reverse = TRUE` and `timefix = TRUE`. Larger scores always mean higher hazard.
Counts sum concordant, discordant and tied-score pairs, including every stratum.
Constant scores remain undefined under this study's reporting convention, with
pair counts retained. Near-equal event times are treated as ties, matching task A's
current concordance calculation. No existing test required keeping `timefix = FALSE`.
The hand-computed unequal-cohort case has six concordant pairs in one cohort
and one discordant pair in the other, giving within-cohort C = 6/7.

The runner accepts a positive integer `--K` and a positive finite `--w_S`.
Rank defaults to the historical configuration value in every mode. Training
weight defaults to 1. The scientific commands explicitly set K = 5. Simulation
retains its registered weights 1 and 10. Actual settings are written once to
`settings.txt`. The smoke test checks the saved models as well as this file.

When a full posterior start is supplied, `code/fit_ssbmf_yfb.R` now saves its
actual moments, precisions and persistent priors immediately before the first
sweep. Ordinary fit objects are unchanged. The paired-study assertion compares
this snapshot with the common start using exact identity. A mutation test changes
the recorded β start and requires the assertion to reject it. Andrew explicitly
authorized this small addition to the fitter shared with Part 2 before E's fits.

Failed folds and cells retain aggregate rows with explicit error messages and
undefined metrics. An incomplete five-fold comparison has an undefined mean and
paired SE. Patient-level fits and scores remain under ignored `data/`.

## Verification and independent review

The two focused files pass 10/10 and 5/5 checks. Independent source review found
one mismatch with the specification: smoke mode still defaulted to rank 2.
That exception was removed so all modes use the configuration default unless
`--K` is supplied. Review found no other defect in concordance, pair counts,
initialization checks, weight propagation or failure schemas. The first full-suite run picked up Part 2's temporarily registered feature-group tests
and reported 821 passed and seven failed. All seven failures were the intentionally
unimplemented feature-group checks. Andrew authorized a temporary copy of the committed
manifest at `19f910bd`, with current Part 1 code and tests. That full suite passed
819 tests with zero failures. It emitted the existing bounded β optimizer
profile-likelihood fallback warning. Scientific fits have not started.

## Text for the shared status files

The matched LB/YFB runner has cohort-specific concordance and pair counts,
explicit K and survival-weight options, and an exact check of both actual starts.
Scientific fitting and external scoring remain pending this software checkpoint.
The registered recovery, false-active and convergence gate remains in force.
