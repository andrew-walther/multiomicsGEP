# Attempt S-14: supervision near the detection threshold

Registered 2026-11-05, before any registered fit. Follows S-13 (running) and two pilot findings
from the same day, which are reported here because they motivate the design.

## Why

1. **On the S-9 test bed, supervision is not needed for recovery.** With the centered input the
   canonical fitter uses, the two-step fit (w_S = 0 loadings, then Cox) recovers the planted
   program in 10 of 10 datasets in every S-9 cell and in the S-13 weak cell
   (`yfb_supervision_strength_2026-11-05/twostep.csv`). The S-9 Gate S gain was measured
   against the S-7 reference, which was fitted to uncentered input and recovered 0 to 1 of 10.
   Against a matched two-step reference on the same input, the gain is 0.
2. **A masked program with too few factors is not a supervision problem in this model.** In a
   pilot with four high-variance neutral programs and one low-variance prognostic program,
   K = 4 misses the prognostic program with or without supervision (w_S up to 50, with or
   without the supervised start), while K = 5 finds it without supervision (|r| = 0.97). Swapping
   out a background program costs far more reconstruction likelihood than survival can supply
   short of overfitting.
3. **Near the detection threshold supervision helps (pilot, seeds 1–2).** With planted shares of
   0.05% and 0.1% at p = 2,000, K = 7, the two-step fit reaches |r| of 0.19 to 0.78 and validation
   C of 0.58 to 0.89. The w_S = 1 fit reaches |r| of 0.74 to 0.87 and C of 0.85 to 0.91, and
   w_S = 20 overfits (calibration slope about 0.1).

This is the regime the method is for: a program present in the expression data but too weak for
the factorization to resolve from noise alone, where survival supplies the missing information.

## Design

- **Generator:** `yfb_decision_data(seed, p = 2,000, share, truth_mode = "shared")`, z-scored per
  cohort as in S-9. Planted shares 0.05%, 0.1% and 0.2%. Seeds 1–10.
- **Survival strength:** "strong", the generator's own predictor (oracle validation C about
  0.94); "realistic", the true predictor rescaled so its training SD is 1.1 (oracle C about
  0.74, calibrated on seeds 101–106 from the generator alone), validation scaled by the same
  factor, outcomes redrawn (seed + 5000) with 30% censoring.
- **Arms (full training fits, K = 7, S-9 settings otherwise):**
  - two-step: canonical w_S = 0 loadings, then Cox on the live program scores;
  - A1 at w_S = 1 (the default and the S-9/S-11 choice);
  - A1 at w_S = 5;
  - A2 (joint q(β)) at w_S = 1.
- **Per fit:** best |r| with the planted direction, recovered (|r| ≥ 0.7, active at 1.96 SD for
  the supervised arms or Cox Wald for the two-step, correct sign), false-active programs,
  validation C, exact-times calibration slope, convergence.

## Registered criteria

1. **Primary (matched Gate S):** A1 at w_S = 1 recovers at least 3 more datasets than the
   two-step fit in at least 2 of the 3 strong-survival cells, with mean false-active at most 1
   and at least 9 of 10 converged.
2. **Secondary:** the mean paired validation C difference (A1 w_S = 1 minus two-step) is
   positive in at least 2 of the 3 strong cells.
3. Realistic-survival cells, A1 at w_S = 5 and A2 are reported without a threshold.

## Software

`results/benchmark_sim/run_yfb_near_threshold.R`, output
`results/benchmark_sim/outputs/yfb_near_threshold_2026-11-05/`. No fitter changes. Synthetic
data only.
