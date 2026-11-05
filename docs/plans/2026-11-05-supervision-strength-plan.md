# Attempt S-13: setting the survival supervision strength without losing stability

Written 2026-11-05. **Status: registered 2026-11-05 after Andrew's go-ahead, before any grid fit. Andrew asked that A2 run as well, so both arms run. The registration details are under "Registration" at the end.**
Follows S-9 (corrected), S-11 and S-12. Andrew asked for survival to have real pull on the
loadings, set from the data with as little manual tuning as possible, in a Bayesian form rather
than a copy of DeSurv's optimizer.

## What is known

1. **At w_S = 1 supervision is too weak.** At 2,000 genes the validation calibration slope is
   1.45 to 3.23 (exact event times), and the risk score has a third to a half of the true
   spread. On real data survival supplies 10⁻⁸ to 10⁻⁶ of the loading precision, and the loadings
   equal the unsupervised ones (|r| = 1.000).
2. **At w_S = 10 calibration is close to 1** in the two stronger simulated cells (0.71, 1.74), and
   every fit converges.
3. **At w_S ≥ 40 the fits drift (diagnosed 2026-11-05 from the saved S-9 traces).** The objective
   rises every sweep (two decreases in 300) but by about 4 × 10⁻⁵ per sweep, above the stopping
   tolerance. In the matched program:
   - survival's share of the loading precision reaches 1.0, so expression no longer constrains it;
   - the program's activity norm falls from 580 to 15, so it stops explaining expression;
   - E β grows steadily from −1.7 to −25 with no plateau.

   This is the pattern of a Cox fit on a separating predictor: with 2,000 genes the program can
   be turned into a direction that orders the training patients almost perfectly, and the Cox
   likelihood keeps rising as β grows. Validation slopes of 0.07 to 0.17 confirm the overfit.
4. **Under the current objective w_S also tempers β.** The fitter maximizes
   E log p(Y) + w_S · E log PL(η) − KL(L) − KL(F) − KL(β). The β update therefore sees the Cox
   likelihood w_S times against one prior. The β posterior SD shrinks as 1/√w_S, which by itself
   raises false-active counts with w_S (independent review, item 2) and removes the prior's
   resistance to the drift in point 3.
5. **DeSurv** (source at the pinned commit, read 2026-11-05) normalizes both losses to order one
   (reconstruction / ‖X‖², Cox / events), so its α = 0.33 is nominally a very large weight. Its
   updates then cap the Cox push on W at α times the NMF gradient norm, clamp the multiplicative
   step, backtrack on the loss, rescale W columns to unit norm, and put an elastic net on β on the
   standardized scale. The realized balance is set by these step rules, not by α. Translated to
   our scale, DeSurv's nominal α would be w_S ≈ 4,000, which is far into the region where our fit
   overfits. **So we do not adopt DeSurv's normalization as a weight.** Its safeguards map to
   Bayesian counterparts below.

## The proposed objective change (arm A1)

Weight the whole survival block, the Cox likelihood and the β prior together:

  L_w(q) = E log p(Y | L, F, τ) − KL(L) − KL(F) + w_S · [E log PL(η) − KL(q(β) ‖ g_β)].

- **L, F, τ updates are unchanged.** L and τ do not involve the Cox term. In the F update the
  survival precision and signal are multiplied by w_S as now (A_F = A_gen + w_S A_surv).
- **The β update becomes untempered.** For q(β) the objective is w_S · [E log PL − KL(β)], whose
  maximizer does not depend on w_S. β is fitted to the Cox likelihood at weight 1 against its
  empirical-Bayes prior, as in an ordinary Bayesian Cox model on the current program scores.
- **Reading.** This is a generalized-Bayes (power) posterior in which the survival block
  (likelihood plus prior on β) is raised to the power w_S as a whole. Equivalently, dividing by
  w_S, the genomics block (likelihood and the L and F priors) is tempered by κ = 1/w_S with Cox at
  weight 1. That is the downweighted-genomics objective of S-10 with the KL option "(b) temper the
  priors with the likelihood" applied to the genomics block, and the Gemini guide's Issue 1
  proposal. At w_S = 1 it is the current posterior exactly.
- **Why it should help with the drift:** the β prior keeps its full weight against the Cox
  likelihood, so β cannot grow without paying the prior cost, and the F update sees E[β²] from an
  untempered posterior. Whether this alone stops the separation drift is the first thing tested.
- **Code change.** Opt-in argument (for example `survival_block_weight = TRUE`) in
  `fit_ssbmf_yfb()` under `canonicalize = TRUE`: the β update receives Cox weight 1, and the
  objective multiplies KL(β) by w_S. Defaults unchanged. Derivation written first
  (`derivations/yfb_survival_block_weight/`), then tests, then code.

## Arms

| Arm | Objective | β posterior |
|---|---|---|
| A0 | current (S-9) | mean-field, tempered by w_S |
| A1 | survival block weighted | mean-field, untempered |
| A2 | A1 | joint K × K Gaussian q(β) (Gemini Issue 4), so correlated program scores share uncertainty |

A2 runs alongside A1 (Andrew, 2026-11-05: run it if it is easy, for a fuller picture).

## Setting w_S from the data

The weight is chosen by five-fold **training CV partial log-likelihood** over a grid
w_S ∈ {1, 2, 5, 10, 20, 50}, not by C. C ignores the size of risk differences, which is
exactly what goes wrong at both ends of the range. CV partial log-likelihood rewards
calibrated risk, so it is a direct, frequentist analogue of choosing the learning rate by
predictive performance as SafeBayes does. If one weight is chosen in nearly every simulated
dataset and on real data, it can become the default and the grid search optional.

## Stability safeguards (Bayesian counterparts of DeSurv's)

| DeSurv | Here |
|---|---|
| unit-norm W columns | scale canonicalization already in the fitter |
| elastic net on standardized β | untempered empirical-Bayes prior on β (A1); point-normal prior as an option |
| capped Cox step, backtracking | exact coordinate ascent of one fixed objective; ELBO decreases are counted, and backtracking is added only if they appear |
| multiple random starts | `fit_ssbmf_yfb_multistart()` for the final real-data fit |

Each fit also reports, per program, survival's share of the loading precision and the program's
share of expression variance. A survival-active program explaining under 1% of expression
variance is labelled a survival-only direction, not a gene program (Gemini Issue 5).

## Test beds

1. **The S-9 cells** at p = 2,000 (planted variance 0.25%, 1%, 5%; seeds 1–10).
2. **A weak-signal cell matched to the real data:** p = 2,000, planted variance 1%, true
   coefficient scaled so that the true predictor's validation C is about 0.68. The S-9 bed is
   much stronger than the real data (C 0.94–0.99), so this is the cell that says whether the
   result transfers.
3. **The null:** the S-11 permuted-outcome datasets.

## Registered criteria (simulation, at the CV-chosen weight, per arm)

1. Gate S as registered (recovery, false activity, convergence), unchanged.
2. Median calibration slope (exact times) between 0.8 and 1.25 in at least two of the three S-9
   cells.
3. At least 9 of 10 fits converge in every cell.
4. Null: mean false-active at most 1 per cell.
5. Weak-signal cell: recovery and C against the w_S = 0 two-step fit, reported without a threshold.

An arm meeting 1–4 goes to real data: refit on the D4 genes at the CV-chosen weight and report
how far the loadings move from the unsupervised ones, survival share, programs, CV C and CV
partial log-likelihood. External scoring then gets its own registration, as S-12 did.

## Registration (2026-11-05)

- **Code:** `fit_ssbmf_yfb(canonicalize = TRUE, w_S, survival_block_weight = TRUE,
  beta_posterior = "meanfield" | "joint")`, derivation
  `derivations/yfb_survival_block_weight/2026-11-05-survival-block-weight-derivation.qmd`, tests
  `tests/test_survival_block_weight.R` (full suite 761/761). All other settings are the S-9 call.
- **Runner:** `results/benchmark_sim/run_yfb_supervision_strength.R`, output
  `results/benchmark_sim/outputs/yfb_supervision_strength_2026-11-05/`.
- **Data:** p = 2,000, seeds 1–10. Cells: the S-9 cells (planted variance 0.25%, 1%, 5%) and the
  weak cell (variance 1%, true predictor scaled by 0.045, giving oracle validation C ≈ 0.68; the
  scale was set on seeds 101–112 from the generator alone, before any fit). Null: the three S-9
  cells with training outcomes permuted (seed + 9000), full fits only.
- **Weights:** w_S ∈ {1, 2, 5, 10, 20, 50}. Folds: five, stratified by event (seed + 7000).
- **Selection rule:** for each dataset and arm, the weight with the largest summed five-fold CV
  partial log-likelihood (Verweij–van Houwelingen) among weights whose five folds all fit; ties
  go to the smaller weight. Criteria 1–5 above are read at the selected weight. Criterion 1 uses
  the S-9 Gate S definition and the S-7 reference. Criterion 2 uses the exact-times slope.
  Criterion 4 is read at the weight selected most often across the signal cells for that arm.
- **Also reported, not gated:** C and slope at every weight, CV C, per-program roles (prognostic
  program, survival-only direction, genomics-only program, inactive), how often the matched
  program is a survival-only direction, and A1 versus A2.

