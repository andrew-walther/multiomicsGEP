# Current project state chapter review

Reviewed the chapter from source against the fitter, prediction code, frozen
benchmark, registered supervision report, and the local DeSurv draft and
supplement. This review covers the chapter. The separate review required
before adopting or merging the experimental fitter remains open.

| Finding | Correction or evidence |
|---|---|
| The monitored expression used undefined prior and likelihood terms. | Define the Gaussian sampling model, stratified Breslow partial likelihood, empirical Bayes mixture priors, product-prior KL terms and point-estimated precision. Check against `code/fit_yfb_supervised.R` and `code/cox_working_logspace.R`. |
| The presentation could conflate historical YFB, the experimental projection and LB. | Label the historical loading feedback weight separately. Define the experimental centered normalized projection and distinguish the LB predictor. |
| The equation needed a scope limit for guarded experiments. | Identify the displayed approximation as the unguarded monitor and describe guarded fitting as a separate bound-based criterion. |
| The prior evidence summary had no completion inventory. | Separate historical and diagnostic deliverables from remaining acceptance evidence. Include methodology, simulations, real data, biology, sensitivity and conditional multimodal validation. |
| Risk of reading a progress fraction as scientific readiness. | State the counting rule and keep the failed supervision gate explicit. |

The frozen predictive result and selected diagnostic replays provide useful
anchors. The 0.627 figure agrees with the frozen benchmark README, and the
360 fits and 153 nonconverged fits agree with the October 30 execution report.
The completion checklist draws on DeSurv Figures 2–4 and supplement Figures
S1–S10. Its categories identify evidence to obtain for the final method,
including null controls, factor correspondence, cohort-level associations
and sensitivity. The checklist is a proposed inventory of deliverables
rather than a new registered experimental design.
