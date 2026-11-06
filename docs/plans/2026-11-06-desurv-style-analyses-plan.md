# DeSurv-style analyses

Registration for [Prompt 7](../prompts/2026-11-06-codex-handoff-prompts.md).
Commit this registration with its index row and attempts-log entry before
running it. Start after Prompt 6 is complete. Use fixed S-11 fits, S-12
scoring, and the Prompt 6 cohort and donor mappings.

1. **Risk groups.** For each fixed score use one cutpoint for all external
   cohorts. Reconstruct S-11 out-of-fold predictions, scaling each with its
   fold-training mean and SD and estimating the w0 Cox coefficients on the
   fold's training rows. Maximize the absolute log-rank statistic over the
   midpoints between distinct standardized OOF scores, requiring at least
   20% of OOF patients on either side. Break ties by closest to zero, then
   by the smaller cutoff. Apply the cutpoint to the external risk using the
   full training score mean and SD. Report per-cohort KM and at-risk counts,
   and high-versus-low HRs by cohort and pooled with cohort strata. Check
   that all predictions were OOF and that the same cutoff was used in each
   external cohort.
2. **Programs, adjustment and enrichment.** For both fixed fits, report and
   plot marginal HRs for each active program per training SD, by cohort and
   pooled with cohort strata. Define active by absolute coefficient divided
   by its posterior SD or Cox SE of at least 1.96. Report continuous-risk
   HRs adjusted separately for PurIST subtype and DeCAF when calls are
   available, with coverage and a same-patient unadjusted estimate. Use
   two-sided Fisher tests for high-versus-low enrichment by cohort and
   report the proportions and counts. Record missing DeCAF sources explicitly.
3. **Rank comparison.** Fit EBMF on the same training D4 matrix with requested
   ranks 7, 15 and 40 and train Cox on the live projected program scores.
   Report achieved ranks and failures. Use the same S-12 external projection
   and Prompt 6 scoring and donor-resampling protocol.
4. **Program rule.** Compute standardized coefficient ratios using
   |βₖ| SD(YFₖ) relative to the largest training value. Use
   Cov(ηₖ, η)/Var(η) as the signed variance share, which allocates the
   covariance terms and sums to one for a nonconstant score. Use one-to-one
   absolute loading-correlation matching with minimum |r| = .70. Require
   posterior activity in at least four S-11 folds, ratio ≥ .25 and share
   ≥ .05. Compare the historical PVE < .01 label with survival dominance
   defined by median loading survival precision share > .50, using
   curvature from the precision audit and the fit's actual survival
   weight. Report raw values and both labels. Reconstruct S-13 and S-14
   checkpoint diagnostics without using held-out outcomes for classification.
   Report full-rule retention and false-active rejection where fold-level
   identity is available. Report the coefficient-and-share rule separately
   where folds were never fitted. S-14's saved coefficients and loadings
   do not identify posterior activity per program, so a full-fit exact
   replay would be required for an exact false-active denominator there.
5. **Mixed-signal simulation.** Extend the S-14 shared-loading DGP to two
   prognostic directions plus neutral programs at p = 2,000. Use shares
   .05 and .001 for the strong and near-threshold directions. Scale the
   combined log risk to training SD 1.1 and give the two standardized
   prognostic components equal coefficients so the weak direction refers to
   molecular detectability. Use seeds 1–10 and S-14's four arms. Keep
   S-14's training and validation sizes and censoring design. Match the
   two planted programs one-to-one by absolute loading correlation. Report
   overlap count and precision among the top 50 absolute fitted gene
   weights against the planted support, with deterministic ties. Track
   recovery, false actives, convergence and all failed or undefined cells.
   Use `coxph.control(timefix = FALSE)` for every Cox fit and
   `frozen_reverse_cindex()` for every C. Verify exact-time controls and
   matching on known synthetic examples before running the registered grid.
6. **Publish aggregate outputs.** Save findings, figures, CSVs, source
   inventory, regeneration commands, and fit and drop denominators. Run
   regression checks and the project test suite before the results commit.
   Keep individual external scores and group membership out of Git.
   Leave the prelim, writeup and progress book to the reviewer.
