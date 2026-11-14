# Fixed-prior survival supervision findings

**Paused draft.** Andrew requested a usage-limit pause on 2026-11-14. The real-data screen is complete. Synthetic Step 1 diagnostics and final audit/review are incomplete. Resume from the [pause handoff](../reviews/2026-11-14-bounded-risk-resume-handoff.md).
Date: 2026-11-14. Author: Andrew Walther. Registration: S-17.

## Result

No tested Normal prior scale passes the registered real-data convergence screen.
All four unit-weight Normal pilots converge stably, but every scale has a failed
higher-target fit. Normal priors with σ_b = 0.5 and 1 support stable target-0.05
fits within the coefficient limit. Neither supports the required targets 0.25 and
0.5. The Laplace comparison also has no passing scale.

The registration requires a scale to pass every target before training CV or the
ten-seed simulation and null screen. Steps 2 and 3 are therefore not run. No new
arm is selected. The prespecified full program-structure analysis and predictive evaluation
under these priors are not run. The working model is unchanged. External scoring is not run.

[The registration](../plans/2026-11-14-bounded-risk-supervision-plan.md) and
[execution prompt](../prompts/2026-11-14-codex-bounded-risk-prompt.md) fix the
families, scales, targets and stopping rule. The reviewing notes at the end of the
[S-16 handoff](../reviews/2026-11-13-share-targeted-weight-handoff.md) replace
S-16's risk-variance null gate. S-17 does not use that historical gate.

## Step 1: real-data convergence

Training uses 273 TCGA and CPTAC patients, 139 deaths and the prepared 2,218-gene
D4 matrix. The model has K = 5, point-Laplace L and F, pooled survival without
strata and the A1 survival block weight. Every same-prior pilot uses weight one.
The target loading vector is selected at that pilot and tracked through at most
three warm-started recalibrations. There is no ramp. Every fit has a maximum of
500 sweeps.

The survival predictor is η = (YF)β. The monitored criterion before
canonicalization is expression log likelihood plus w_S times the Cox term,
minus loading and score KL terms and w_S times the coefficient KL. Weighting
both coefficient terms leaves the conditional β update untempered at fixed F,
while increasing survival curvature in the F update. The target share is a
working quadratic precision ratio along the selected loading direction. It is
not a fraction of total Cox information or a measured change in orientation.

Normal coefficients use a fixed shared N(0, σ_b²) prior. Laplace coefficients use
fixed scale σ_b/√2, giving the same variance. Fixing these distributions imposes
shrinkage. They have unbounded support and are not hard bounds on coefficients.

The table gives final target passes. Target cells list **sweeps, stable
convergence, maximum absolute coefficient mean**. All pilots are stable, and their
cells give sweeps and maximum absolute coefficient mean. A scale passes only
when its pilot and every final target pass are stable and their maximum
coefficient means are at most 4σ_b.

| Prior | σ_b | 4σ_b limit | Pilot | Target 0.05 | Target 0.25 | Target 0.5 | Scale gate |
|---|---|---|---|---|---|---|---|---|
| Normal | 0.1 | 0.4 | 143, 0.219 | 500, no, 0.610 | 500, no, 0.556 | 500, no, 60.358 | fail |
| Normal | 0.25 | 1 | 66, 0.350 | 500, no, 1.556 | 500, no, 1.598 | 500, no, 24.977 | fail |
| Normal | 0.5 | 2 | 67, 0.383 | 64, yes, 0.421 | 500, no, 4.133 | 500, no, 8.844 | fail |
| Normal | 1 | 4 | 67, 0.392 | 64, yes, 0.434 | 500, no, 8.818 | 500, no, 16.320 | fail |
| Laplace | 0.1 | 0.4 | 145, 0.281 | 500, no, 21.319 | 500, no, 11.043 | 268, yes, 17.349 | fail |
| Laplace | 0.25 | 1 | 98, 0.358 | 500, no, 15.194 | 500, no, 14.939 | 500, no, 23.525 | fail |
| Laplace | 0.5 | 2 | 66, 0.373 | 500, no, 14.969 | 500, no, 18.341 | 500, no, 17.523 | fail |
| Laplace | 1 | 4 | 66, 0.384 | 65, yes, 0.438 | 500, no, 19.598 | 500, no, 20.141 | fail |

Of the 24 final target states, four are stable and three also meet the
coefficient limit. All eight pilots meet both conditions. No complete scale
passes all targets.

The Laplace σ_b = 0.1, target-0.5 arm uses two calibration fits. Both converge
stably, at 330 and 268 sweeps, but the final coefficient maximum is 17.349 against
a limit of 0.4. Target 0.05 also uses two fits for Normal σ_b = 0.25 and Laplace σ_b =
0.25 and 0.5. Those initial passes converge, but recalibration fails. The 32 final
real arm records therefore represent 36 actual fits. No real pilot triggers the
no-signal guard. The final Laplace σ_b = 0.1, target-0.5 pass reaches the 100,000
weight cap. Its achieved tracked share is 0.444248, with loading correlation
0.029019 to the pilot, so the tracked program identity is weak. No real fit emits
a warning or task error.

For Normal σ_b = 0.1, target 0.05 and target 0.25 meet the risk-stability
thresholds at the sweep cap, but fail the current ELBO criterion. Risk stability
alone is not counted as convergence. At target 0.5 the same prior ends with
maximum absolute β of 60.358 and maximum absolute η of 19.606. In its final sweep,
the program with that coefficient has β = 36.046 before canonicalization and a
score SD of 1.674, giving β = 60.358 after scaling. This is a finite-run diagnostic,
not proof that the objective is unbounded.

[Detailed fit records](../../results/benchmark_sim/outputs/bounded_risk_2026-11-14/screen_fits.csv)
retain the native ELBO flag separately from stability, all calibration passes,
weights, warnings and elapsed seconds. The per-program share table reports every
program, including inactive programs. A fixed continuous prior does not switch
off survival curvature through an estimated zero prior variance. Its activity
counts therefore have a different interpretation from the working model.

Three of the eight final target-0.05 arms meet both the stability and coefficient
conditions. They give a limited positive result. Fixed priors permit
moderate directional shares that failed under the S-16 estimated prior. Their
tracked loading remains close to the same-prior unit-weight pilot:

| Prior | σ_b | Final weight | Achieved directional share | Absolute loading correlation with pilot |
|---|---|---|---|---|
| Normal | 0.5 | 392.274 | 0.042223 | 0.997342 |
| Normal | 1 | 372.044 | 0.062140 | 0.997177 |
| Laplace | 1 | 391.387 | 0.044017 | 0.996892 |

These correlations describe the single direction already tracked for
calibration. They do not replace the unrun all-program comparison, annotation
or CV analysis in Step 2. Target attainment is approximate: the inherited rule
accepts an achieved directional odds ratio within a factor of two of the target
odds, or stops after three converged passes. Its tolerance is not an additional
convergence gate.

## Step 1: synthetic convergence diagnostics

Synthetic diagnostic results and the completeness audit are being assembled.
This section will be replaced after the registered Step 1 grid finishes.

## Steps 2 and 3: not reached

No Normal σ_b passes Step 1, so training CV, program-structure analysis, the
weak/null ten-seed screen, paired-weight nulls and the S-14 recovery comparison
are not run. Their criteria are not failed results. They are untested. No
adoption candidate is established, including for the individually stable
Normal target-0.05 arms.

## Implementation and verification

The opt-in fitter trace records the raw SD of each YF score before canonicalization, the scale actually applied, and coefficient means before and after rescaling. For a live program, beta_after = c_k beta_before. A collapsed score uses scale one, so the raw SD and applied scale are both retained. These are aggregate per-program quantities. The inputs are the current training matrix and posterior moments. The output is one row per program per sweep in history$s16.

The new stable_convergence_window argument defaults to zero. A positive value requires tracing and ELBO stopping. For S-17 its value is 50. The stopping block first checks the existing current-sweep relative ELBO criterion. It then requires all 50 relative L2 changes in eta to be below 0.001 and the endpoint change in max absolute eta to be below 1%. The denominator is the beginning-of-window maximum, floored at machine epsilon. The algorithm keeps iterating until both criteria hold or the 500-sweep cap is reached. history$native_elbo_converged separately records the current ELBO criterion. No default fit or update formula changes.

The shared S-16 fitting helper accepts prior_beta through fit_control. Normal remains its default. The fixed shared prior is passed unchanged to each same-prior pilot and warm-started calibration fit. The calibration helper can retain pass_fits, including the complete trace, for every pass. The default remains false. The dataset and real-data helpers also accept separate two-step reference controls. A separate unit-weight fit is still needed for the estimated per-program working-model reference. The fixed-prior pilot must not be mislabeled as that reference.

The new screen runner separates real and simulation stages and partitions each by family. Each process is serial. The Normal prior is created with ashr::normalmix(1, 0, sigma), and the equal-variance Laplace prior with ebnm::laplacemix(1, 0, sigma/sqrt(2)). These existing packages supply the fixed priors accepted by the fitter, without a new dependency. yaml reads the registered ranks and fitting tolerance from project configuration. Base R handles traces, checkpoint serialization and aggregate CSVs.

Each case is keyed by dataset, family and sigma. The private checkpoint contains source/input hashes, package versions, controls, fits, warnings, elapsed time, posterior second moments and final Cox weights. Atomic writes preserve completed arms. A stale signature is rejected. A nonconverged pilot remains a reported failure and leaves target fits unavailable because the weight calculator requires a converged pilot. An unexpected error is saved and rethrown. A no-signal guard reuses the pilot at weight one and is explicitly marked. Its requested target has not been attained.

The exporter independently recomputes the 50-sweep stability check from the saved trace, deduplicating its repeated sweep-level quantities across programs. It reports the final and peak risk/coefficients and the largest within-window deviation. Only endpoint risk drift enters the registered 1% gate. The real-data pair gate uses the pilot and final pass of all three targets and the final maximum coefficient mean limit of four prior SDs. Synthetic results do not determine which priors proceed.

Fixing a Normal or Laplace scale imposes shrinkage, not a hard coefficient bound. Canonicalization preserves LF' and eta, but rescaling q(beta) while holding the prior fixed generally changes its KL divergence. The native objective is computed before this rescaling. These fits test the registered algorithm and do not establish an exact CAVI derivation under a fixed prior expressed throughout in standardized coordinates. The added trace makes this limitation observable. The four-SD rule is an empirical acceptance threshold.


## Reproduction and evidence

Final commands, timing, audit totals and output ledger will be inserted after
all Step 1 cases finish. Patient-level matrices, risks and checkpoints stay in
ignored data/. Only aggregate outputs are committed.

## Decisions and limits

The study uses a new grid orchestrator around the existing `s16_fit()` and
`s16_calibrate()` helpers. The historical diagnostic script initializes its
fixed output folders when sourced and fits isolated weights without
recalibration. Extending that top-level workflow would mix the historical
S-16 namespace with the new case grid. The S-17 wrapper instead supplies
opt-in prior controls, tracing and retained calibration passes through the
shared fitting path. `s16_diagnostic_fit()` and its historical CLI remain
unchanged and are not used for S-17. This file-layout choice uses the prompt's
allowance for a new file where extending the old workflow is awkward.


The 50-sweep stability window must be observed before a fit can pass. An opt-in
stopping control continues after an early ELBO threshold until both the current
ELBO criterion and the registered window hold, or the fit reaches 500 sweeps.
This operationalization was independently reviewed and recorded before fitting.
It prevents an unobserved window from becoming either a pass or an automatic
failure. Default fitting behavior is unchanged.

The 1% maximum-risk check compares the beginning and end of the window. The
largest intermediate deviation is reported separately and does not tighten the
registered gate. The four-SD coefficient threshold applies to final means, while
peak coefficients are also reported. A no-signal guard is a weight-one fallback,
not an achieved stronger target.

The prepared D4 screen and cohort scaling are fixed historical preprocessing.
There is no nested preprocessing or new CV performance estimate here. The
500-sweep limit cannot distinguish every slow fit from a fit without a stable
fixed point. The synthetic strong-signal coefficient is about 6.7, outside all
registered fixed-prior scales. That cell assesses convergence, not recovery under
a correctly specified prior.

## Text for the shared status files

S-17 tested fixed shared Normal and equal-variance Laplace priors on survival
coefficients at four scales, with same-prior pilots and share-targeted weights.
No Normal scale passed the complete real-data stable-convergence and coefficient
screen. Moderate target-0.05 fits converged for Normal σ_b = 0.5 and 1, but the
required stronger targets failed. The registration therefore stopped before
training CV, program characterization and the ten-seed null/recovery study.
The unit-weight working model remains unchanged. No new external score was
computed and no fixed-prior arm is adopted.

## Draft paragraph for the progress book and Chapter 4

**Draft, for review.** We tested whether fixing the prior variance of the
standardized survival coefficients would permit stronger supervision of the
molecular loadings. Priors with standard deviations of 0.1, 0.25, 0.5 and 1 were
compared at unit survival weight and at weights calibrated to directional
curvature shares of 0.05, 0.25 and 0.5. The unit-weight Normal pilots converged
stably. A target share of 0.05 also converged under the two weaker Normal priors,
with standard deviations of 0.5 and 1. Neither prior supported stable fits at
both stronger targets, and no prior scale passed the complete registered
screen. We therefore did not evaluate these models by training cross-validation
or the registered null and recovery checks. The experiment does not establish
a new supervision weight or show that predictive or molecular-program structure
improved. Fixing a prior scale supplies shrinkage, while canonicalization still
rescales the coefficient posterior, so it does not itself impose a hard bound
on the standardized coefficients.
