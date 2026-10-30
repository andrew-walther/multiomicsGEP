# Problem statement: loading recovery without controlled survival activity

## Question

Can a nonnegative Bayesian factorization with normalized projection predictor
η = centered(Y) F β / ‖F_k‖ learn a low-variance prognostic program more often
than a reconstruction-only factorization, while controlling unmatched
survival-active programs and converging reliably?

The registered October 30 grid does not establish that claim. Numerical
posterior failures were diagnosed and repaired independently, but the stable
updates still fail the registered scientific gate. This statement is intended
for a methods discussion, not as evidence that every possible projection
model or weight must fail.

## What was fixed and what was held fixed

- Log-space Breslow risk sets, curvature floor and coupled gene-wise updates.
- Nonnegative cohort-SD scaling with centering only inside the Cox projection.
- Loading normalization and activity RMS canonicalization with prior pushforwards.
- Hybrid supervised start; the unsupervised reference disables that start.
- Point-normal coefficient prior with bounded nonzero slab proportion and SD.
- Numerical normal-means scaling, stable exponential moments/KL checked by
  independent integration, and a warned same-model prior-optimizer fallback.
- No factor deletion/refitting, no external cohort model choice, no changed gate.

The seventh registered attempt corrects the shared-loading units and uses a paired design extended from the calibrated A1 bed:
10 datasets per cell; p=50 or 2000; raw centered planted variance share 0.25%,
1% or 5%; K=7; training n=400 and independent synthetic validation n=120.
The grid targets precision shares 0.05, 0.2 and 0.5 at the paired reference,
using its original damping denominator, then applies reduced F-message
damping 0.02. β retains the full γ weight. Each γ is fixed through its fit.

## Saved numerical evidence

At the saved weak pseudo-message, independent quadrature yields posterior
variance about 8.54; the direct ebnm second moment instead loses essentially
all of it. Canceling the common Gaussian density before computing the
spike/slab normalizer restores the same mathematical posterior and its KL.
The independently checked repair is isolated from the historical kernel.
The earlier failed example then converges at all four diagnosed damping
values, while retaining several unmatched survival-active factors.

The earlier grids used raw-expression risk directions with SD-divided
reconstruction loadings and therefore cannot establish shared-loading
recovery. The corrected grid generates risk from normalized F/SD_train,
using the same training direction for validation. The coefficient now has
prepared-projection units. Raw and prepared variance shares and oracle C are
saved separately; predictive strength is not comparable across attempts.

The registered failures and repairs are retained separately. A valid result
cannot be obtained by dropping failed fits or pooling repeated attempts.

## Final registered results

The corrected shared-loading grid completed all **360/360** fits, with
**0 raw numerical-error checkpoints**,
**0 lost-variance records**, and
**153 nonconverged fits**. **Gate S failed:** none of the
three share-target arms qualifies in the required two large-feature cells.
The original criteria are unchanged. Earlier grids are retained as numerical
and design diagnostics because their risk and reconstruction directions used
different units.

| attempt | saved | missing | task_level_errors | raw_checkpoint_errors | variance_loss_records | sweep_cap_records | nonconverged | undefined_C |
|---|---|---|---|---|---|---|---|---|
| attempt1 | 80 | 280 | 56 | 4 | 0 | 1 | 285 | 284 |
| attempt2 | 360 | 0 | 0 | 149 | 0 | 12 | 161 | 149 |
| attempt3 | 360 | 0 | 0 | 155 | 0 | 6 | 161 | 0 |
| attempt4 | 360 | 0 | 0 | 0 | 53 | 115 | 140 | 0 |
| attempt5 | 90 | 270 | NA | 0 | 4 | 22 | 293 | 270 |
| attempt6 | 360 | 0 | 0 | 0 | 0 | 129 | 129 | 0 |
| attempt7 | 360 | 0 | 0 | 0 | 0 | 153 | 153 | 0 |

Missing files, fitting exceptions, variance loss and caps are distinct and can
overlap in the conservative nonconvergence total. Task-level status is
unavailable for the deliberately interrupted fifth driver; its stop reason
and complete missing-record denominator are saved. Undefined C is counted
explicitly. Partial-state C does not establish a successful fit.

| attempt | all_zero_returned_states | nonfinite_returned_eta | inactive_factor_total | orphan_factor_total | unavailable_factor_states | warning_records | objective_decreases |
|---|---|---|---|---|---|---|---|
| attempt1 | 0 | 0 | 261 | 0 | 228 | 4 | 479 |
| attempt2 | 0 | 0 | 431 | 45 | 149 | 60 | 1325 |
| attempt3 | 0 | 0 | 882 | 302 | 0 | 60 | 2951 |
| attempt4 | 0 | 0 | 781 | 267 | 0 | 124 | 13708 |
| attempt5 | 0 | 0 | 238 | 82 | 225 | 31 | 3318 |
| attempt6 | 0 | 0 | 844 | 305 | 0 | 112 | 16600 |
| attempt7 | 0 | 0 | 823 | 293 | 0 | 145 | 16862 |

Warning-record counts are cumulative within datasets in the original driver;
they are not independent warning-event counts. Zero/nonfinite counts concern
returned states; unavailable states are explicit. Objective decreases concern
the monitored approximation and cannot be described as monotone exact ascent.

| variance | arm | recovered | baseline_recovered | recovery_gain | mean_false | converged | qualifies |
|---|---|---|---|---|---|---|---|
| 0.0025 | share0.05 | 0 | 0 | 0 | 6.8 | 2 | FALSE |
| 0.01 | share0.05 | 0 | 0 | 0 | 6.8 | 1 | FALSE |
| 0.05 | share0.05 | 0 | 1 | -1 | 7 | 1 | FALSE |
| 0.0025 | share0.2 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.01 | share0.2 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.05 | share0.2 | 0 | 1 | -1 | 6.8 | 1 | FALSE |
| 0.0025 | share0.5 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.01 | share0.5 | 0 | 0 | 0 | 7 | 0 | FALSE |
| 0.05 | share0.5 | 0 | 1 | -1 | 7 | 0 | FALSE |

The exact cell and arm verdicts are in `attempt7/gate_cells.csv` and
`attempt7/gate.csv`. All fit-level failures and derived penalties are in the
attempt-specific `metrics.csv`; raw checkpoints are unchanged.

The fixed diagnostic replays match registered loading means exactly. A loading
match can coexist with a sweep cap, unmatched activity and severe distortion
of reconstruction. Feature-specific residual precisions become extremely
small on affected genes; all posterior calculations can remain finite while
joint fitting is unstable. The fixed cases are diagnostic, not extra gate
observations.

| variance | recovered | false_active | converged | decreases | largest_mean_difference | first_decrease |
|---|---|---|---|---|---|---|
| 0.0025 | FALSE | 6 | TRUE | 1 | 0 | 2 |
| 0.01 | FALSE | 4 | FALSE | 17 | 0 | 2 |
| 0.05 | TRUE | 4 | FALSE | 84 | 0 | 2 |

| variance | expected_reconstruction_loss | min_tau | max_loading_energy_fraction |
|---|---|---|---|
| 0.0025 | 0.103143039382764 | 0.00554563962473137 | 0.0148814777294625 |
| 0.01 | 939.106205910707 | 4.16326340616329e-08 | 0.999074191437772 |
| 0.05 | 453803162951241 | 8.7483777770172e-20 | 0.999999999999997 |


The controlling gate requires a single grid arm to improve recovered-active
program counts by at least 3/10 in at least two large-p variance cells, with
mean unmatched survival-active count at most one and at least 90% convergence
in each qualifying cell. Full denominators and numerical errors are reported.

## Exact point of failure and current hypotheses

**The unresolved issue is statistical rather than only numerical.** A stable
nonzero posterior and a converged objective do not establish that the recovered
program is uniquely prognostic or that posterior activity is calibrated.

The main hypotheses supported by the diagnostics are:

1. **Tempered uncertainty.** Multiplying the Cox likelihood by γ increases its
   conditional coefficient precision. Away from a prior boundary, posterior
   SD scales roughly as γ^(-1/2). A nonvanishing slab prevents absorbing zeros,
   but cannot by itself stop weak or mixed directions from becoming apparently
   precise survival effects. This mechanism is plausible, not a completed
   frequentist coverage analysis. In the fixed-projection Normal approximation,
   z_γ ≈ √γ z_1. If the untempered null statistic is approximately standard
   Normal, P(|z_γ|>1.96) ≈ 2Φ(−1.96/√γ), which tends to one. Spike selection,
   projection uncertainty and learning F complicate this approximation, but
   ordinary credible-interval cutoffs are not automatically calibrated after
   increasing γ solely to balance reconstruction precision. A tiny coefficient
   can cross the registered z cutoff if its conditional variance is even
   smaller; changing to an absolute-effect threshold would change the
   pre-registered activity definition.
2. **Program mixing.** Factors need not align one-to-one with the planted
   program. Several mixtures can carry its survival signal and be flagged
   active, even when loading recovery is poor. The factorization lacks an
   identification condition forcing one uniquely prognostic column. Nearly
   collinear saved projections also expose a mean-field limitation: within a
   fixed included-coefficient pattern with joint Gaussian precision H,
   marginal variance is (H⁻¹)_kk, while an independent coordinate uses 1/H_kk.
   A structured coefficient posterior would require the corresponding cross-
   coefficient terms in the loading update; replacing SDs afterwards would
   not be a coherent correction to the registered approximation.
3. **Approximate objective.** Freezing and then recomputing the loading
   normalizer, and damping F while fully weighting β, does not guarantee
   ascent of one fixed normalized objective. Recorded decreases and sweep-cap
   fits remain meaningful limitations even after numerical posterior repair.

4. **Survival-only directions.** The model allows L_k to become exactly zero
   while F_k remains nonzero in the survival projection. Such a column is
   absent from genomic reconstruction but still predicts survival. Its scale
   is now identified numerically, but it is not a reconstructed expression
   program. Saved fixed-case factor tables identify these columns. Forcing
   joint support, discarding them, or redefining the recovery estimand would
   change the statistical model and cannot be justified by a numerical patch.

5. **Feature-noise adaptation.** In the fixed saved cases, very large
   reconstruction errors coexist with near-zero residual precisions and
   loading energy concentrated on a small number of genes. The prescribed
   τ_j=n/Q_j update downweights those features as their errors grow. This is
   finite arithmetic under the specified model, not a lost-variance bug.
   A noise prior, fixed precisions or an acceptance objective would change
   the method and needs a stated statistical choice.

The normalizer issue has a concrete mathematical form. If r_k=‖F_k‖ and
Z_k=Yc F_k/r_k, then

∂Z_k/∂f_jk = Yc_j/r_k − Z_k f_jk/r_k².

The registered coordinate update holds r_k fixed and does not include the
second term. Correcting this is not a scalar coding patch: it changes the
nonlinear optimization and the posterior approximation that supplies the
normal-means messages. A fixed normalized objective and its uncertainty
interpretation must be stated before claiming exact variational ascent.

Another unregistered search over weights, prior bounds or outcome-based gene
screens would not distinguish these hypotheses. No further numerical change
is supported by the repaired failing fixtures alone. Changing the coefficient
learning rate, reporting untempered conditional uncertainty, adding program
identification constraints, or switching to a latent predictor changes the
statistical method and needs a new stated decision and registration.

## Questions to resolve before another fitting study

- Should the next step derive a single normalized objective, including the
  derivative of the loading normalization, and an acceptance rule?
- Is γ a generalized-Bayes learning rate whose uncertainty needs calibration,
  or a discriminative loading penalty with a separate untempered coefficient
  fit? These imply different models and inferential claims.
- Should the planned fallback η=Lβ with fold-in prediction be evaluated next?
- Is recovery of one uniquely prognostic factor the right estimand if mixed
  factor representations are allowed, and what identification assumptions
  would make that claim defensible?

## Reproduction and evidence

Read the registered attempts log before running anything. Raw checkpoints,
complete-denominator summaries and gate tables are under
`results/benchmark_sim/outputs/yfb_supervision_2026-10-30/`.
Use `attempt7/gate_cells.csv`, `attempt7/gate.csv`, `attempt7/counts.csv`,
`attempt7/simulation_units.csv`,
`numerical_failure.csv`, `stable_damping_diagnosis.csv` and the fixed diagnostic
replays. The execution report identifies source revisions and test output.
The precision audit is conditional on a transformed saved training state;
it is not a newly fitted real-data model.

Focused mathematical checks can be run by sourcing `tests/test_helpers.R`,
`tests/test_fit_yfb_supervised.R` and `tests/test_yfb_stable_exponential.R`.
The exact registered engine is identified in `provenance.csv`; later code
adds diagnostics and config wiring without changing registered fitting values.
Do not overwrite original attempt directories. No patient-level inputs are
needed to reproduce these synthetic numerical fixtures or gate tables.
