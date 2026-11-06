# Pruning after the flat β likelihood fix (KL-2)

[KL-2](../plans/2026-11-06-flat-beta-pruning-rerun-plan.md) repeated the
four-arm, 160-fit synthetic [KL-1 comparison](2026-11-06-kl-fix-pruning-rerun-findings.md)
at starting K = 7. The ordinary and protected pruning arms both failed the
registered recovery and rank criteria after the flat β fix.

| Arm | Noise | Recovered /20 | Mean false-active | Returned ELBO converged /20 | Median K |
|---|---:|---:|---:|---:|---:|
| Reference | 0.25 | 1 | 1.15 | 20 | 7 |
| Unpruned | 0.25 | 17 | 1.05 | 20 | 7 |
| Ordinary | 0.25 | 2 | 1.10 | 20 | 2 |
| Protected | 0.25 | 2 | 1.10 | 20 | 2 |
| Reference | 1 | 1 | 0.40 | 20 | 7 |
| Unpruned | 1 | 19 | 0.35 | 20 | 7 |
| Ordinary | 1 | 2 | 0.20 | 20 | 1 |
| Protected | 1 | 2 | 0.20 | 20 | 1 |

All 160 fits returned without errors or fit warnings. At noise 1, seed
20261017 hit one intermediate refit cap in each pruning arm. The final
fits converged. Both pruning arms passed the false-active criterion
and the returned-fit convergence criterion.

Every seed pruned in both arms. Both arms lost 15 of 17 initially recovered
programs at noise 0.25 and 17 of 19 at noise 1, with no gains among seeds that
started unrecovered. Both arms had the same per-seed recovery and retained
ranks. Some continuous diagnostics differed between the arms.

For ordinary pruning, the largest saved first-stage deletion gain was
positive in all 40 unpruned fits. Its candidate had a zero loading mean in
31 fits. These first candidates were classified from the rounded per-column
gains and squared loading norms in the paired unpruned rows. Of the 29
initially recovered fits in that group, 25 lost recovery. The other seven
initially recovered fits had a nonzero-mean first candidate and all seven
lost recovery. These summaries track the first candidate and final recovery.
They do not trace the planted program through each recursive refit.

## Implementation and review

For the fixed shared β prior, the coefficient likelihood is flat when its
unfloored quadratic and linear terms are both zero. The update then returns
the normal or Laplace mixture prior moments and exact zero negative KL. The
canonical fitter uses that value in the ELBO and omits the flat coordinates
from empirical Bayes prior estimation. If all coordinates are flat, it keeps
the prior from the preceding sweep. The deletion rule still compares the
full ELBO, which includes the molecular posterior terms. A zero loading mean
alone does not establish a zero deletion gain.

Review `code/update_beta.R`, `code/fit_ssbmf_yfb.R`, and
`tests/test_empty_program_kl.R`. The regression covers the prior second
moment, exact zero KL and a positive precision below the numerical floor.
`tests/test_collapse_diagnosis.R` also reconstructs the shared prior from
the informative entries in its trace.

The design, one row per fit, summary and run log are in
`results/benchmark_sim/outputs/flat_beta_pruning_rerun_2026-11-06/`.
The source snapshot used base `356347f` plus the two-file fitter patch
recorded in `provenance.txt`. This comparison uses current branch code,
which includes changes since KL-1. The result supports retaining the K = 7
unpruned working fit. Full-state refitting was outside this registration.
The historical Task1-5 evaluation tested that option and failed its
adoption gate.

Verification on the isolated source snapshot with the repository's existing
local test fixtures linked in was 781 passed and zero failed. The standalone
flat-likelihood regression passed 3/3. The full suite emitted a bounded β
optimizer profile-likelihood fallback warning in a separate diagnostic.
