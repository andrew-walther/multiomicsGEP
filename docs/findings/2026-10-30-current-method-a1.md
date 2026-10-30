# Does a larger survival weight make survival supervise the loadings?

Written 2026-10-30 for review at Gate A1 of
[`2026-10-26-parallel-work-sessions-plan.md`](../plans/2026-10-26-parallel-work-sessions-plan.md).
Branch `fix/current-method`, commits 6d1d674 (simulation) and 17ce0c5 (real data).
No external cohort was used.

## Question and answer

In the single-modality model, η = (YF)β, survival can enter the update of
the gene loadings F through a term multiplied by a survival weight w_S. With
w_S = 1, survival supplies only about 10⁻⁶ of the loading-update precision on
real data. The question is whether survival fails to shape F simply because
that weight is too small.

**No.** In the current code no value of w_S makes survival supervise F
without breaking something else:

- In simulation, w_S = 1 already gives survival almost all of the planted
  program's loading precision and recovers that program in 18–20 of 20
  datasets on the calibrated test bed. Larger weights recover it less often,
  flag nearly every program as survival-associated, and lower validation C.
- On the real training data, w_S = 10² raises the survival share only to
  about 10⁻⁴. The model then flags 6 of 7 programs as survival-associated, and
  cross-validated C falls slightly. At w_S = 10⁴ the share reaches about 0.5,
  but no fold converges, reconstruction fails and cross-validated C drops to
  0.555.

This matches what the rebuilt `ssbmf` package found with a survival weight
chosen by cross-validation. The extra survival-associated programs and the
mixed loadings are the same pattern, so the problem lies in the model's
balance between the two likelihoods, not in the current implementation.

## Does the update include the within-program coupling term?

It depends on the code path, and the experiments used the one that does.

- With `canonicalize = TRUE`, `fit_cox_on_yf()` updates F gene by gene
  (`multimodal_yfb_update_F_mk()`). Gene j's survival signal subtracts
  E[β_k²] Σᵢ wᵢ y_ij (Z_ik − y_ij E[f_jk]), where Z_k is program k's
  projection with the other genes' current loadings. The projection is
  refreshed after each gene, so the update includes the coupling between
  genes of the same program. Canonicalization also rescales every program to
  sd(Y E[F_k]) = 1, the role of the normalizer d_k in the 10/29 advisor
  review.
- With `canonicalize = FALSE`, `w_S` goes to the older vector update
  (`update_F_surv_YFB_k()`). It removes only the other programs' predictor
  and updates all genes at once, so it omits the coupling term. The 10/29
  review describes this omission for the `alpha_F` update, and the same
  omission applies here.

Every fit below uses `canonicalize = TRUE`.

## Simulation

**Design.** We used the existing test bed in
`results/benchmark_sim/run_supervision_diagnostics.R` (arm set `a1_ws`).
Each dataset has an expression block with p = 50 genes and three true
programs. Two programs carry most of the variance and have no effect on
survival. The third is the planted prognostic program, with at most 0.5% of
the variance. Each training and validation cohort is z-scored separately.

We fit K = 7 with point-Laplace priors on L and F, a per-program Normal β
prior and at most 300 sweeps. The grid crossed w_S ∈ {1, 10², 10⁴, 10⁶} with
the supervised start off and on (program K started from the null-model Cox
score direction). We ran 20 seeds at noise levels 0.25 and 1 on two beds:

- calibrated: n = 400 training patients, planted log hazard ratio 3.3, the
  smallest setting at which a survival-only lasso finds the program;
- original: n = 120, log hazard ratio 1.1, where even the lasso mostly fails.

That is 640 fits, with no errors.

**Metrics.**

- *Recovered:* the best-matching fitted program has |r| ≥ 0.7 with the planted
  loadings, a coefficient of the correct sign, and is survival-active
  (|E[β]| ≥ 1.96 SD).
- *False active:* the number of other programs that are survival-active.
- *Survival share:* w_S·A_surv / A_F, the survival share of the loading
  precision, as the median over genes at the last sweep, for the matched
  program. This is the quantity `loading_balance_terms()` computes. The
  trace records it inside the fit.
- *Reconstruction loss:* ‖Y − E[L]E[F]ᵀ‖² / ‖Y‖² on the training data. On this
  bed the unsupervised fit gives 0.83–0.97 and a rank-3 SVD gives 0.81–0.89,
  so values near 0.87–0.97 are normal. Values above 1 mean the fit is worse
  than predicting zero.
- *Validation C:* computed on a separate validation cohort of 120 patients.

### Calibrated bed (n = 400, log HR 3.3); 20 seeds per noise level

| w_S | Start | Recovered (0.25; 1) | Median \|r\| (0.25; 1) | False active per fit | Converged | Median recon. loss (max) | Validation C (0.25; 1) |
|---|---|---|---|---|---|---|---|
| 1 | off | 18; 19 | 0.88; 0.89 | 1.35; 0.40 | 19; 19 | 0.87 (0.91); 0.97 (0.99) | 0.771; 0.762 |
| 1 | on | **20; 20** | 0.87; 0.89 | 1.00; 0.25 | 19; 19 | 0.87; 0.97 | 0.770; 0.761 |
| 10² | off | 15; 20 | 0.75; 0.84 | 3.15; 0.90 | 20; 20 | 0.86; 0.99 | 0.756; 0.749 |
| 10² | on | 14; 20 | 0.77; 0.86 | 3.60; 1.65 | 20; 20 | 0.86; 0.99 | 0.756; 0.749 |
| 10⁴ | off | 16; 20 | 0.79; 0.84 | 6.05; 4.40 | 20; 20 | 0.90 (0.97); 0.95 | 0.756; 0.748 |
| 10⁴ | on | 16; 20 | 0.79; 0.84 | 6.10; 5.10 | 20; 20 | 0.88; 0.93 | 0.756; 0.748 |
| 10⁶ | off | 14; 20 | 0.80; 0.84 | 6.30; 5.65 | 20; 20 | **2.0 (1,220)**; 0.94 (4.6) | 0.756; 0.748 |
| 10⁶ | on | 14; 20 | 0.80; 0.84 | 6.30; 5.80 | 20; 20 | **4.3 (895)**; 0.92 (3.2) | 0.756; 0.748 |

The oracle validation C, using the true linear predictor, is 0.851 and 0.947.
A fit has K = 7 programs, so 6 false active programs means every program
except the planted one.

### Original bed (n = 120, log HR 1.1); 20 seeds per noise level

| w_S | Start | Recovered (0.25; 1) | Median \|r\| (0.25; 1) | False active per fit | Median recon. loss (max) | Validation C (0.25; 1) |
|---|---|---|---|---|---|---|
| 1 | off | 1; 6 | 0.17; 0.58 | 0.60; 0.55 | 0.87; 0.97 | 0.559; 0.651 |
| 1 | on | 1; 7 | 0.16; 0.59 | 0.60; 0.60 | 0.87; 0.97 | 0.557; 0.643 |
| 10² | off | 0; 1 | 0.31; 0.57 | 3.15; 1.60 | 0.87; 1.00 | 0.565; 0.640 |
| 10² | on | 0; 1 | 0.33; 0.58 | 3.00; 1.90 | 0.86; 0.98 | 0.564; 0.640 |
| 10⁴ | off | 0; 3 | 0.37; 0.59 | 6.95; 6.55 | 0.85; 0.90 | 0.564; 0.639 |
| 10⁴ | on | 0; 3 | 0.36; 0.59 | 6.90; 6.65 | 0.82; 0.87 | 0.564; 0.639 |
| 10⁶ | off | 0; 3 | 0.37; 0.59 | 7.00; 6.85 | **13.4 (366)**; **2.6 (68)** | 0.564; 0.639 |
| 10⁶ | on | 0; 3 | 0.37; 0.59 | 7.00; 6.85 | **17.3 (275)**; **2.3 (3,351)** | 0.564; 0.639 |

All 320 fits converged. Oracle validation C is 0.678 and 0.839. At noise 1
with w_S = 1, 2 fits (start off) and 1 fit (start on) gave a constant risk
score, so their validation C is undefined. The means above leave them out.

### What the simulation shows

- **The weight is not what limits recovery here.** At w_S = 1 the matched
  program's survival share is already about 1.0 on the calibrated bed and at
  noise 1 on the original bed. It is 0.007 only at noise 0.25 on the
  original bed, where the program is not found by any method, including the
  lasso. Raising w_S cannot add information the data do not have.
- **Larger weights overfit survival.** At noise 0.25 recovery falls from
  18–20 to 14–16 of 20 on the calibrated bed and from 1 to 0 on the original
  bed. Validation C falls by about 0.015 on the calibrated bed. The number of
  false survival-active programs per fit rises from about 1 to about 6, so
  nearly every program is flagged. On the original bed the median |r| rises from 0.17 to 0.37 while
  recovery falls to 0. The model then finds a mixture of the planted program
  and other directions, as the rebuild found (|r| 0.40).
- **Convergence is not the problem.** Every fit with w_S ≥ 10² converged.
  Sweeps that lowered the monitored objective also became rarer: 68–105 per
  arm of 20 fits at w_S = 10² and 0–28 at w_S ≥ 10⁴, against 189–250 at
  w_S = 1. At large weights the survival term dominates both the objective
  and the update. At w_S = 1, one calibrated
  dataset per noise level reached the 300-sweep cap, the same datasets
  reported on 10/15.
- **Reconstruction holds up to w_S = 10⁴ and fails at 10⁶.** The loss
  reached 1,220 times ‖Y‖² in one fit.
- **The supervised start makes almost no difference.** It adds 1–2 recovered
  datasets at w_S = 1 and none at larger weights.

## Real training data (TCGA + CPTAC), 5-fold CV

**Design.** `results/benchmark_sim/run_a1_real_ws_cv.R` uses the D4 training
matrix (273 patients, 2,064 genes, 139 events), with gene screening and
per-cohort z-scoring done once on all training patients. It also uses the
folds of the 10/15 training-CV check, stratified by cohort and event. The
configuration is that of the simulation.

Arms:

- the two best simulation settings: w_S = 1 with the start on, and with it
  off;
- w_S = 10² and 10⁴ with the start off. These two answer the weight question
  on data where the share at w_S = 1 is far below 1, because both w_S = 1
  arms use the same weight.

There were 20 fits in total.

**Metrics.**

- Held-out C: computed per fold, then averaged.
- Pooled out-of-fold C.
- CV partial log-likelihood: the Verweij–van Houwelingen form, summed over
  folds.
- Survival share: the median over survival-active programs of the share
  defined above.
- Training reconstruction loss.
- Held-out reconstruction loss: the fold's patients projected onto the
  fitted loadings by least squares.

| w_S | Start | Mean fold C (SE) | Pooled OOF C | CV partial log-lik. | Converged | Active programs (of 7) | Survival share | Train recon. loss | Held-out recon. loss | max \|β\| |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | on | 0.656 (0.024) | 0.652 | −779.2 | 4/5 | 2.6 | 1.2e-6 | 0.504 | 0.549 | 0.30–0.45 |
| 1 | off | 0.656 (0.023) | 0.652 | −778.8 | 4/5 | 2.6 | 1.2e-6 | 0.504 | 0.549 | 0.30–0.45 |
| 10² | off | 0.648 (0.025) | 0.644 | −781.1 | 5/5 | 6.2 | 1.1e-4 | 0.504 | 0.550 | 0.38–0.56 |
| 10⁴ | off | 0.555 (0.043) | 0.547 | **−56,241** | **0/5** | 7.0 | 0.53 | **19.6** | 0.952 | **13.7–16.9** |

The non-converged fold at w_S = 1 is fold 4, which reached the 300-sweep cap
in both arms. At w_S = 10⁴ all five folds reached the cap, and their training
reconstruction losses were 7.9–38.6 times ‖Y‖². Fold-level values are in
`outputs/a1_real_ws_cv/folds_*.csv`.

### What the real data show

- **The survival share grows about in proportion to w_S** (1e-6 → 1e-4 →
  0.5). Survival only reaches a share comparable to reconstruction at a
  weight of about 10⁴, and there the fit breaks down. The coefficients grow
  to about 15 log hazard ratio units per SD of score, and every program
  becomes survival-active. The loadings no longer reconstruct the data, and
  the held-out partial likelihood is 70 times worse.
- **At w_S = 10² the loadings barely move,** with a share of about 10⁻⁴ and
  the same reconstruction. Yet 6 of 7 programs become survival-active and C
  drops by 0.008. This is the same overconfidence seen in simulation.
- **The supervised start changes nothing on real data.** The fold C values
  differ by at most 0.006.
- **We did not run w_S = 10³.** By the proportional trend its share would be
  about 10⁻³, still far below reconstruction. We did not add arms to look for
  a value between 10² and 10⁴, because the plan forbids tuning this line.

## Why the simulation and real-data shares differ by six orders of magnitude

The share at w_S = 1 is about 1 in the 50-gene simulation and about 10⁻⁶ on
the 2,064-gene real data. Part of the gap comes from the number of genes.
Canonicalization fixes the SD of each program's score Y E[F_k] at 1. With
more genes, each loading is therefore smaller and each activity larger, so
the reconstruction precision τ_j Σᵢ E[l_ik²] grows with p. The survival
precision does not. The real programs also explain far more variance than
the planted program. We have not decomposed the gap further, so this
explanation is partial. It does mean that test beds with p = 50 do not
reproduce the balance problem seen on real data. A bed with p in the
thousands would, as the rebuild's 2,000-gene cells did.

## Comparison with the rebuilt model (`ssbmf`)

The rebuild reports a share of 10⁻⁵ to 4 × 10⁻⁴ in its 2,000-gene
simulation and about 5 × 10⁻⁷ on real expression at w_S = 1. With w_S chosen
by cross-validation from {1, 10, 100}, its 0.25%-variance program became
survival-active in 12 of 15 datasets, but with |r| of only 0.40. The 90%
β intervals then covered the truth only 27% of the time.

The current code behaves the same way:

| Rebuild | Current code |
|---|---|
| Larger w_S makes the program survival-active but gives mixed loadings | Same: \|r\| 0.17 → 0.37 on the original bed |
| Larger w_S gives overconfident β | Same: up to 7 of 7 programs flagged active |
| w_S ≥ 1,000 is slow (rejected message steps) | Not slow: w_S ≥ 10⁴ fits take 0.3–0.8 s in simulation against 2–6 s at w_S = 1 |

We did not measure β interval coverage. In this test bed the coefficient is
per SD of the fitted projection score, which has no single true value to
cover.

## How this bears on the 10/29 review and the 10/30 proposal

The 10/29 advisor review finds that survival does not supervise F in the
evaluated fit (α_F = 0) and that the α_F > 0 option omits d_k and the
coupling term. The canonical w_S path tested here has both, so the review's
objection to the update's form does not apply to these results.

The 10/30 proposal scales survival by γ = α/(1−α) · Σ_m v_m p_m / N_event.
For the real training data with α = 0.5, γ ≈ 2,064/139 ≈ 15. For α = 0.9 it
is about 130. These weights fall in the range tested here, between w_S = 1
and 10². On this evidence they would raise the survival share only to about
10⁻⁵ to 10⁻⁴ and make most programs survival-active. The proposal assumes
that γ makes the survival precision "O(1) relative to reconstruction". On
this data that holds only near γ = 10⁴, where the fit fails. These results
do not test the proposal's other parts, such as precision damping, the
non-collapsing β prior and point-exponential priors.

## Failures and limitations

- **The first real-data run failed completely.** All 20 fits used forked
  workers (`mclapply`) and segfaulted, which is consistent with forking after
  BLAS and the compiled sweep were loaded on macOS. The runner was switched
  to a socket cluster and rerun. The table reports the rerun only, and no
  result came from the failed run.
- **Non-converged fits:**
  - simulation: 4 of 80 at w_S = 1 on the calibrated bed (one dataset per
    noise level, in both start arms), each at the 300-sweep cap, and none of
    the other 560;
  - real data: 2 of 10 at w_S = 1 (fold 4 in both arms) and 5 of 5 at
    w_S = 10⁴.
- **Some C values are undefined.** In 3 of the 40 original-bed fits at
  w_S = 1 and noise 1, the risk score was constant, so validation C is undefined. These
  fits are left out of the C means and are counted in `c_validation_na`.
- **The real-data CV is conditional on screening.** Gene screening and
  scaling were done once on all training patients, so absolute CV C is
  optimistic. The comparison between arms shares the same folds.
- **Only one β prior was tested.** All arms use a per-program Normal β prior,
  not the shared point-Laplace prior chosen on 10/15. That choice changes how
  many programs are flagged active, so the false-active counts are specific
  to this prior.
- **Tests:** 710/710 pass. A first run in this worktree showed 18 failures
  because git-ignored inputs (`data/`, two D4 output files) and the isolated
  DeSurv library live only in the main checkout. With those linked in, all
  710 pass.

## Gate A1

**Question:** does any w_S make supervision of F work without breaking
convergence or reconstruction?

**Answer: no.**

- *In simulation:* w_S = 1 is already the best setting, and every larger
  weight lowers recovery and validation C.
- *On real data:* weights small enough to keep the fit stable leave survival
  at 10⁻⁶ to 10⁻⁴ of the loading precision. The weight that gives survival
  an equal share (about 10⁴) breaks convergence, reconstruction and
  prediction.

A single global multiplier on the Cox term does not supply what is missing.

For session B and the 10/30 proposal, the evidence points away from a
tuned weight and toward:

- a test bed with thousands of genes;
- a β prior that does not flag every program when survival is weighted up.

Decision for Andrew: keep, revise or drop this line. We recommend dropping
w_S > 1 as a fix for the current code and recording the result for the
rebuild.
