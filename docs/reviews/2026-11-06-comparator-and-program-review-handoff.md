# Review handoff for the real-data comparators and DeSurv-style analyses

Review the two registrations and the two results commits on
`fix/yfb-supervision`. The findings documents list all aggregate sources and
the commands that regenerate them. The registered recovery, false-activity
and convergence gate still governs adoption of survival supervision.

## Independent checks and limitations

- The fixed S-11 reference C values were checked against every saved S-12
  cohort value. The new pooled comparison reproduces the S-12 two-step
  contrast and donor interval. The remaining contrasts use the same
  estimator and match patients within each cohort.
- The cohort-stratified pooled HR, the equally weighted paired C contrast
  and the cohort-specific calibration slopes are different estimands.
  In particular, a larger pooled HR does not establish better discrimination.
  The HR intervals are model-based and do not use the C contrast's
  donor bootstrap.
- Elastic-net CV is conditional on the S-11 prepared feature matrix.
  Inspect the candidate path warnings before interpreting the alpha
  comparison. DeSurv uses the pinned author pipeline and published
  settings, with fitting restricted to the aligned training samples.
  The absence of a common clinical covariate limits comparison to the
  fitted molecular models.
- The two OOF cutpoints are model-specific. The risk-group HR estimates
  are descriptive for those thresholds. The local historical DeCAF labels
  distinguish permissive and restraining CAF calls. The subtype comparisons
  and enrichment tests are exploratory and their p-values are unadjusted.
- The program variance allocation is signed and can be negative. The S-11
  precision calculation gives an upper bound at the saved final state.
  S-13 and S-14 support a partial rule check only because program-level
  fold activity was not saved. S-14 lacks supervised per-program activity.
- The mixed-signal experiment distinguishes molecular signal from
  survival-effect strength. Its shared-loading construction gives the
  two directions equal standardized coefficients. A model can predict
  survival without recovering the planted near-threshold loading.
- These external data have been used in earlier project development.
  The locked scoring analysis did not choose settings from external
  outcomes, and its contrasts do not constitute a new independent
  validation of the overall research workflow.

The review priorities are the TCGA+CPTAC cohort alignment, the bootstrap
membership for repeated donors, the training-only cutpoint, the use of
precision bounds, and the denominators for each rule and simulation.
