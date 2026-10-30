# Targeted YFB mechanism screen

Twelve synthetic fits and a separate fixed-Gaussian-design coefficient study,
registered in the targeted diagnostic plan before execution. No patient data
or external-cohort scoring. All source-fitting defaults remain numerically
unchanged against the saved pre-change fixture.

- `design.csv`, `study_controls.csv`: registered cases and controls.
- `v*_*.rds`: complete fitted states, warnings, metrics and rejected proposals.
- `summary.csv`: all twelve cases; raw recovery can coexist with nonconvergence.
  `success` requires all pre-registered screening conditions. No case succeeds.
- `noise_comparison.csv`: matched adaptive/fixed-τ reconstruction comparison.
- `fit_checks.csv`: posterior finiteness and guarded accepted-objective checks.
- `beta_metrics.csv`, `beta_summary.csv`: fixed-design Gaussian working model,
  not actual Cox coverage. Independent-coordinate and exact finite-pattern
  posteriors use the same fixed prior and likelihood weight. Nonconvergence
  is explicit. Planted-signal runs rule out an always-zero apparent solution.
- `checkpoint_hashes.csv`: immutable synthetic fit hashes.
- `trajectories.png`: accepted-state trajectories; short lines may be stalled,
  not converged. Objectives from different criteria are not compared.
- Compressed logs and `test_results.csv` retain red/green and full verification.

The guard's normalized moment criterion uses a rigorous fixed-normalizer Cox
curvature bound but recomputes normalization from posterior means. It is not
claimed to solve the fixed-generative-model normalization/ELBO question.
The first completed feasible sweep supplies its initial comparison state.
The guard changes both working curvature and acceptance; gains cannot be
attributed to acceptance alone. Fixed τ is a diagnostic, not an adopted prior.

No result establishes the original complete-grid Gate S. Never pool these
fixed cases with that grid or select an arm from synthetic validation C.
