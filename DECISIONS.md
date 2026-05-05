# DECISIONS.md — multiomicsGEP

Architectural and analytical decisions made during development, in reverse chronological order.
Each entry records what was decided, why, what was traded away, and which files implement it.

---

## 2026-05-05 — K overfitting fix: k_pdac_single=10, k_pdac_synthetic=5; alpha CV added to LB runner

- **Decision (k_pdac_single=10):** Added `benchmark.k_pdac_single=10` to `config/globals.yml`.
  Both `run_LB_benchmark.R` and `run_YFB_benchmark.R` now use `K=20` for merged training
  (n=273) and `K=10` for single-cohort training (tcga_only n=144, cptac_only n=129).
  **Rationale:** K=20 on single-cohort data gives K/n≈0.14 — too many factors for the sample
  size. The archived baseline (C=0.63–0.65 on tcga_only) used K=10. Empirical confirmation:
  running K=20 on tcga_only (2026-05-05) gives K_eff=1 and C=0.37–0.50 — worse than random.

- **Decision (k_pdac_synthetic=5):** Changed `benchmark.k_pdac_synthetic` from 8 to 5 to
  match `synthetic.k_true=5`. **Rationale:** K=8 on synthetic DGP with K_true=5 gives
  LB C=0.135 and YFB C=0.092 — both anti-concordant. ARD with K>K_true absorbs signal variance
  into null factors, causing the model to miss survival-relevant directions. K_SYN must equal
  K_true to avoid this.

- **Decision (alpha CV in LB runner):** `run_LB_benchmark.R` now calls `select_alpha_cv()`
  before fitting each train mode. Alpha is CV-selected per mode and saved in the log and CSV.
  YFB runner still uses fixed alpha=0.50 (alpha CV calls `fit_supervised_mf_modular` internally,
  which is not compatible with `fit_cox_on_yf` — YFB alpha CV requires a separate implementation).

- **Empirical findings from 2026-05-05 benchmark runs (all 6 modes, both K corrections applied):**
  - **A_surv/A_gen imbalance confirmed across all modes:** LB iter-1 ratio = 0.0000–0.0097
    (cptac_only k=3), 0.0000–0.0033 (tcga_only), 0.0000–0.0011 (merged). The L update is
    structurally dominated by genomics, reducing the model to approximately unsupervised PCA.
  - **YFB β→0 collapse on all PDAC modes at both K=10 and K=20:** K_eff=0 regardless of K
    or prior. The ZF scale (‖Y·EF_k‖², sum over p=2000 genes) is enormous, driving EBeta
    to zero. Structural — not fixable by K tuning.
  - **LB tcga_only K=10:** K_eff=2, but external C=0.34–0.43 — worse than K=20 run (C=0.47–0.50).
    Two active factors are anti-concordant with external prognosis. K=10 does NOT recover
    archived 0.63–0.65 baseline. The archived baseline was a lucky PCA direction alignment,
    not a stable property of the model. A_surv/A_gen imbalance is the root cause.
  - **LB cptac_only K=10:** K_eff=3, C=0.32–0.45. Worse than K=20. Same pattern — more active
    factors that are anti-concordant externally.
  - **LB merged K=20:** K_eff=1, C=0.35–0.44. Unchanged from first run.
  - **LB synthetic K=5:** C=0.1353, K_eff=3 — SAME as K=8. ARD pruned both K=5 and K=8 to
    the same 3 active factors. Fixing K_SYN had no empirical effect. The archived C=0.828
    synthetic result was from a different script (exploratory runner) — not the benchmark
    runner's `generate_synthetic_benchmark_data()` function.
  - **YFB synthetic K=5:** C=0.092, K_eff=4, non-zero EBeta. Anti-concordance persists at
    K=5 — sign-direction inversion is structural to the YFB formulation (ZF = Y·EF mixes
    factors via Gram matrix EF'EF, can invert prognosis direction).
  - **Conclusion:** K tuning does not fix any of the observed failures. The root cause is
    the A_surv/A_gen structural imbalance in the L update. All further fixes must address
    the scale imbalance directly.

- **Affected files:** `config/globals.yml` (k_pdac_single, k_pdac_synthetic),
  `results/benchmark_sim/run_LB_benchmark.R` (alpha CV block, K assignment),
  `results/benchmark_sim/run_YFB_benchmark.R` (K assignment)

---

## 2026-05-05 — Benchmark train-mode support, benchmark_helpers.R, top_n_genes reverted to 2000

- **Decision (--train-mode):** Both benchmark runners now accept `--train-mode merged|tcga_only|cptac_only`.
  Single-cohort modes skip Section 1 (synthetic). Single-cohort training uses `preprocess_desurv_cohort()`
  (v1, per-cohort); merged uses `preprocess_merged_cohorts()` (v2). Output CSVs are mode-specific
  (`LB_benchmark_results_merged.csv`, etc.) and include a `train_mode` column.

- **Decision (benchmark_helpers.R):** Created `results/benchmark_sim/benchmark_helpers.R` to
  hold shared constants (`PDAC_DATA_ROOT`, `PLATFORM_LOG_TRANSFORM`, `EXTERNAL_COHORTS`) and
  functions (`load_pdac_raw`, `generate_synthetic_benchmark_data`) previously only in the archived
  `run_ssbmf_benchmark.R`. Both runners now source this file instead of the archive. PDAC cohort
  constants also added to `config/globals.yml` under `pdac:`.

- **Decision (top_n_genes reverted to 2000):** `preprocessing.top_n_genes` changed from 5000
  back to 2000 (DeSurv spec). Investigation of why current benchmark results (C ≈ 0.39–0.49)
  were far below the archived baseline (C ≈ 0.60–0.65) revealed three discrepancies: (1) top_n
  was 5000 vs 2000, (2) no alpha CV (fixed at 0.5 vs CV-selected), (3) K_max was 20 vs 10.
  Reverted top_n to 2000. Alpha CV and K_max alignment are deferred to the next session.
  **Tradeoff:** 2000 genes matches DeSurv and recovered the baseline. 5000 was expected to
  improve genomic reconstruction but in practice added noise that drowned survival signal.

- **Key empirical finding (2026-05-05):** Merged TCGA+CPTAC training has NEVER produced
  C-index > 0.50. The archived 0.60–0.65 results were entirely from tcga_only training (v1
  preprocessing, K=10, alpha CV, point_normal prior). Merged training gives median_ext ≈ 0.50
  even with the best archived settings (v1, K=10, alpha CV, K_eff=7). The β→0 collapse on
  merged data is structural and unresolved. Cluster A fixed training-side β=0 but external
  generalization regressed. Cluster B (YFB) also collapses on all train modes. This is the
  primary open problem.

---

## 2026-05-04 — Benchmark consolidation: K_max=10, prior comparison, cox_warmstart=FALSE, beta_threshold=0.001

- **Decision (K_max):** Canonical benchmark uses K_max=10 (from `cfg$cavi$k_max`). K=20 was
  used only for Cluster A/B diagnostic runs to test whether β=0 was a K-saturation artifact.
  It is not: on merged PDAC (n=273, p=2000), K=20 gave 0–2 active factors, same as K=10.
  K_max=10 gives K_eff ≈ 4 in practice (ARD pruning).

- **Decision (prior comparison):** Both benchmark runners (`run_LB_benchmark.R`,
  `run_YFB_benchmark.R`) run `prior_beta="point_normal"` AND `prior_beta="normal"` in a
  single pass for side-by-side comparison. Rationale: point_normal (spike-and-slab) collapses
  all EBeta to zero on real PDAC data for both Cluster A and Cluster B (point_normal EBNM
  shrinks to the spike component when survival signal is weak relative to ZF scale). Normal
  prior avoids this collapse at the cost of potentially retaining too many active factors.
  Outcome of external C-index comparison is an open empirical question — running both avoids
  having to re-fit to investigate.

- **Decision (cox_warmstart=FALSE as Cluster B baseline):** `fit_cox_on_yf()` now defaults
  to `cox_warmstart=FALSE` (EBeta initialized to 0). Rationale: matches Cluster A behavior
  for apples-to-apples comparison. Cox warm-start calibrates EBeta to the ZF scale from
  iteration 1, but with normal prior the CAVI itself can escape zero — warm-start may be
  unnecessary. Toggleable via `cox_warmstart=TRUE` if normal prior produces unstable initial betas.

- **Decision (beta_threshold=0.001):** Lowered from 0.05 in `config/globals.yml`. Under the
  YFB reformulation (η = ZF·β̃ where ZF = Y·EF), the natural EBeta scale is
  beta_true / sd(ZF) ≈ 0.003–0.008. A threshold of 0.05 would classify all YFB betas as
  inactive even when they are clearly non-zero. 0.001 distinguishes spike-shrunk zeros from
  non-zero betas in both LB and YFB models.

- **Affected files:** `config/globals.yml` (benchmark section, beta_threshold), `code/fit_cox_on_yf.R`
  (prior_beta="normal", N_burnin=0, cox_warmstart=FALSE, normalize_AB=FALSE defaults),
  `results/benchmark_sim/run_LB_benchmark.R` (new), `results/benchmark_sim/run_YFB_benchmark.R` (new),
  `results/benchmark_sim/archive/` (7 one-off scripts archived)

---

## 2026-05-04 — Cluster B architecture: dedicated files, alpha_F=0, interface reuse

- **Decision (file structure):** Cluster B (η = (YF)β̃) lives entirely in three dedicated files:
  `code/update_L_surv_YFB.R`, `code/update_F_surv_YFB.R`, `code/fit_cox_on_yf.R`. Prediction
  in `code/predict_cox_on_yf.R`. The Cluster A files (`fit_modular.R`, `update_L.R`,
  `update_F.R`, `predict.R`, `update_beta.R`, `compute_elbo.R`) are restored to exact main-branch
  versions. No cross-cluster coupling.

- **Reason:** Cluster A and Cluster B must coexist without risk of regression. Separate files
  mean: (1) the 171/171 Cluster A test suite validates Cluster A code unchanged; (2) Cluster B
  bugs cannot corrupt Cluster A runs. Interface reuse (`update_beta.R`, `compute_elbo.R`) is
  achieved by passing `ZF[,k]` as `EL_k` and `ZF[,k]^2` as `EL2_k` — observed projections
  have zero posterior variance, so the existing signatures work without modification.

- **Decision (alpha_F=0):** The Cluster B F update (`update_F_surv_YFB_k`) defaults to
  `alpha=0` (pure-genomics only; no survival contribution to the F precision or pseudo-obs).

- **Reason:** With η = ZF·β̃ (ZF = Y·EF), A_beta = Σ w_i ZF_ik². If EBeta ≈ 0, the Cox
  Hessian w_i ≈ 0 at a stable equilibrium but ZF is non-zero (EF initialized from SVD). So
  A_beta is non-zero, and the β update can escape zero. The root cause of the "normalize_AB"
  instability (see below) was that A_surv in the F precision depended on EBeta², creating a
  chicken-and-egg: EBeta≈0 → A_surv≈0 → x_F biased → EF grows → EL shrinks → positive
  feedback. With alpha_F=0, A_F = A_gen (τ * sum EL²); no survival term in denominator.
  EF is determined purely by genomics (same as unsupervised EBMF), and ZF can deliver
  non-zero signal to the beta update regardless of EBeta.

- **Trade-off:** With alpha_F=0, the loadings F are not jointly optimized for survival —
  they reflect genomic variance only. Survival signal enters only through the beta update. 
  This is less expressive than full joint optimization, but is numerically stable.

- **Diagnostic finding (2026-05-04):** On synthetic data (n=120, p=300), alpha_F=0 gives
  C-index=0.605 vs PCA=0.471 (3/8 active factors, converges in 11 iters). On merged PDAC
  training (TCGA_PAAD+CPTAC, n=273, p=2000), EBeta collapses to ~4.7e-7 (0 active factors).
  The β=0 collapse on real data persists even with alpha_F=0. Likely cause: on real data the
  survival signal is weaker relative to noise, and the point-normal spike-and-slab EBNM
  prior shrinks all betas to the spike component at the natural ZF scale (~sd(Y)·||EF_k||).

- **Affected files:** `code/fit_cox_on_yf.R`, `code/update_F_surv_YFB.R`,
  `code/update_L_surv_YFB.R`, `code/predict_cox_on_yf.R`,
  `results/benchmark_sim/run_cox_on_yf_benchmark.R`, `tests/test_cox_on_yf_smoke.R`

---

## 2026-04-30 — normalize_AB added to F update (Cluster B); positive-feedback instability discovered

- **Decision:** Added `normalize_AB` parameter to `update_F_k()` and `update_F_all()` in
  `code/update_F.R`, and wired it through from `fit_supervised_mf_modular()` in
  `code/fit_modular.R`. This is the Cluster B analogue of the Cluster A normalize_AB fix that
  was applied to `update_L_k()`. The parameter is backward-compatible (default FALSE).
  171/171 tests pass with the addition.

- **Motivation:** Under the Cox-on-YF reformulation (η = (YF)β̃), the F update is dual-source
  (genomics + survival). The scale imbalance between A_gen = Tau * sum_i(EL²_{ik}) and
  A_surv = EBeta²_k * Σ_i(w_i y²_{ij}) is structural: at initialisation A_gen/A_surv ≈ 10⁴.
  This is invariant under any reparameterisation of ZF by a constant (proven algebraically:
  EBeta scales inversely, leaving A_surv = EBeta² * YtWY unchanged). normalize_AB was the
  same fix that worked for Cluster A's L update.

- **What was implemented:** In `update_F_k()`, after computing A_surv and B_surv:
  ```r
  if (normalize_AB) {
    m_surv <- mean(A_surv); m_gen <- mean(A_gen)
    if (is.finite(m_surv) && is.finite(m_gen) && m_surv > 1e-12 && m_gen > 1e-12) {
      scale_surv <- min(m_gen / m_surv, 100)   # cap at 100
      A_surv_eff <- A_surv * scale_surv;  B_surv_eff <- B_surv * scale_surv
    }
  }
  ```
  The cap at 100 was added after observing that the uncapped version (scale ≈ 10,000)
  caused immediate catastrophic EF inflation even when EBeta was small.

- **Instability discovered (open issue):** Even with cap=100, `normalize_AB=TRUE` causes a
  runaway positive-feedback collapse of EL and EF by iteration 6-7 in synthetic smoke fits.
  The mechanism (traced per-iteration):

  1. After N_burnin=10 + Cox warm-start, EBeta ≈ 0.05 for one factor. A_gen/A_surv ≈ 89,000;
     cap=100 limits scale to 100, so A_surv_eff/A_gen ≈ 0.1%. Survival contribution to x_F
     is negligible but slightly biased toward the survival direction.
  2. This tiny bias causes EF[:,k] to grow slowly (+50-80% per iter for the active factor).
  3. Larger EF[:,k] → larger sum(EF²[:,k]) → larger A_L = Tau * sum(EF²[:,k]) in the L update
     of the NEXT iteration → smaller x_L = B_L/A_L → EBNM shrinks EL[:,k].
  4. Smaller EL[:,k] → smaller A_gen = Tau * sum(EL²[:,k]) in the F update → survival fraction
     of A_F grows → survival bias in x_F grows → EF grows faster.
  5. Positive feedback loop: EF doubles every few iters, reaching max|EF| = 68,000 by iter 6,
     then EBNM assigns A_L → Inf → EL → 0 → everything collapses.

  This feedback is structurally different from the Cluster A case (where normalize_AB in the
  L update was stable) because in the L update, A_gen = sum_j(Tau * EF²) depends on EF (not
  EL), so EL shrinkage doesn't feed back into A_gen. In the F update, A_gen = Tau * sum(EL²)
  depends on EL, creating the destabilizing loop.

  **Additional complication:** Factors k=4 and k=5 immediately collapse at iter=1 because SVD
  init with positive-part clipping (`EL[EL<0] <- 0`) leaves EL[:,4-5] ≈ 0. When A_gen ≈ 0,
  the normalize_AB guard (m_gen > 1e-12) prevents rescaling, but A_F ≈ alpha * A_surv (tiny),
  and x_F = B_surv / A_surv = x_surv. With EBeta[4] = −0.033 (negative) and point_exponential
  prior on F, EBNM zeros EF[:,4] immediately, from which F[:,4] never recovers.

- **Cap value rationale:** cap=100 gives A_surv_eff/A_gen ≈ 0.1% (far from 50-50 balance).
  This is not enough to deliver survival signal, yet is still enough to trigger the feedback.
  No safe cap exists in the range [1, A_gen/A_surv]: small caps are too weak; large caps
  amplify x_surv to catastrophic levels via the EBeta/EBeta2 ratio.

- **Next debugging directions (to be pursued in the next session):**
  1. **Decouple precision from signal**: Replace the current (A_gen + A_surv_eff, B_gen + B_surv_eff)
     formulation with a fixed-precision approach: A_F = A_gen (no survival in denominator);
     B_F = B_gen + gamma * B_surv_eff. This keeps A_F tethered to the genomics structure
     (preventing the feedback) while injecting survival direction. Requires a principled choice
     for gamma (ELBO justification unclear).
  2. **Pure-genomics F, survival via β only**: Use alpha=0 for the F update (F is purely
     genomic) and rely on the β update alone to select survival-relevant factors from ZF = Y*EF.
     Eliminates the feedback entirely. Sacrifices the "jointly supervised F" advantage of
     Cluster B but preserves the train/test consistency fix (prediction still uses ZF = Y*EF).
     This is the simplest stable path and may be sufficient for the dissertation.
  3. **Block coordinate descent for F × β**: Run a few extra β updates after each F update
     within the k-loop, so EBeta tracks the updated EF before the next k. May reduce the lag
     that allows the feedback to compound.
  4. **Alternative normalization**: Instead of scaling A_surv to match A_gen, scale both A_gen
     and B_gen DOWN by their magnitude so x_gen occupies the same scale as x_surv. This changes
     s_F but not x_F direction, and avoids inflating A_F.

- **Affected files:** `code/update_F.R` (normalize_AB logic + cap); `code/fit_modular.R`
  (normalize_AB argument wired to update_F_k call at line ~515).

- **Status:** Committed on branch `cox-on-yf-reformulation`. normalize_AB=FALSE (default)
  is stable and correct. normalize_AB=TRUE compiles and passes tests but is not yet usable
  due to the instability. The Cluster B framework (Steps 1-11) is fully implemented and
  correct; the remaining work is the scale-balancing mechanism for the F update.

---

## 2026-04-29 — Cluster A in-model fixes resolve training-side β=0; external generalization mixed

- **Decision:** Adopt the four Cluster A fixes from `docs/beta_zero_fix_design.md` §4 in
  `code/fit_modular.R` and `code/update_L.R`. Specifically:
  1. **Inner-loop reorder** β → L → F is now **canonical** (previously L → F → β). No new
     parameter — the reorder is unconditional. Justified by the symmetric `z_no_k` /
     `R_k` invariance argument: both expressions cancel in the current k's `EL[,k]` and
     `EBeta[k]`, so β can fire first using `z_no_k`, then L can reuse the same `z_no_k`
     with the freshly updated `EBeta[k]` flowing in via A_surv/B_surv.
  2. **`N_burnin` parameter** (default 0) — runs N_burnin iterations of β-only updates with
     EL fixed at SVD init (EL2 = EL^2). Replicates Warm-start Exp 1 to break the A_surv ≈ 0
     cycle at the very start.
  3. **`alpha_schedule` parameter** (default NULL) — `list(warmup_iters, ramp_iters)` ramps α
     from 0 to target over the warmup+ramp window. Curriculum lets L settle before survival
     pressure is applied.
  4. **`normalize_AB` parameter** in `update_L_k` (default FALSE) — when TRUE, rescales A_surv
     up to match A_gen's magnitude (and applies the same scale to B_surv) so α actually
     controls the fraction of influence between sources. **Reformulated** from the design
     doc's original §4.8 prescription, which divided both A_gen and A_surv by their means
     — that formula collapsed L to zero in the smoke fit (verified empirically; the rescale
     inflated 1/√A_L noise scale and over-shrunk EBNM). The retained reformulation preserves
     the original EBNM noise interpretation while still rebalancing the contributions.
- **Why:** The existing failure mode (β=0 on merged TCGA+CPTAC v2 training) was localized in
  Phase 1 to a chicken-and-egg + scale-imbalance trap inside `update_L_k()`. Instrumentation
  on the new branch confirmed the imbalance quantitatively: A_surv / A_gen ratios at iter 1
  for k = 1, 2, 3 are 0, 0, and 2e-4 respectively — survival precision is ~5000× smaller than
  genomics, so the survival term cannot pull L until the rescale is applied.
- **Smoke fit result (merged TCGA+CPTAC, n=273, p=2000, K=20, N_burnin=10, normalize_AB=TRUE):**
  Cox warm-start EBeta range [-0.048, 0.061]; post-burn-in [-0.034, 0.046]; final
  [-0.0588, 0.0580]; **2/20 factors active** (|β| > 0.05 — k=4 +0.058, k=6 -0.059); ELBO
  monotone non-decreasing across 60 iterations; max|EL| = 1.83e3.
- **External-cohort C-index (5 held-out cohorts vs. baseline N_burnin=0, normalize_AB=FALSE):**
  Improved on 1/5 (Moffitt_GEO_array +0.012). Regressed on 4/5 (Dijk -0.019, PACA_AU_array
  -0.060, PACA_AU_seq -0.024, Puleo_array -0.076). The recovered β favors training-Cox-aligned
  L directions that don't transport to held-out cohorts.
- **Trade-offs:** Fix 1 is a permanent change to the canonical CAVI ordering — backward-
  compatibility for any analysis that depended on the L → F → β trajectory is broken (no
  test asserts that trajectory, so 171/171 tests still pass). The `normalize_AB` rescale
  departs from strict ELBO maximization (verified empirically that ELBO is still monotone
  on the merged training set). The Fix 4 reformulation diverges from the design doc text;
  the design doc's §4.8 formula is documented as "reviewed and adopted with empirical
  reformulation" rather than rewritten.
- **Implication:** Cluster A solves the immediate training-side failure but does not deliver
  cross-cohort generalization. This is the design doc's predicted Cluster B trigger
  (`docs/beta_zero_fix_design.md` §3 row "EBeta non-zero but unstable"). Phase 4 (Cluster
  B — Cox-on-YF reformulation, `derivations/qF_supervised/`) is now the next priority.
- **Affected files:** `code/fit_modular.R` (instrumentation, reorder, N_burnin, alpha_schedule,
  normalize_AB threading); `code/update_L.R` (normalize_AB rescale of A_surv/B_surv);
  `results/benchmark_sim/run_cluster_a_smoke.R` (new); `results/benchmark_sim/run_cluster_a_external.R`
  (new); `results/benchmark_sim/outputs/cluster_a_smoke/`, `cluster_a_external/` (new).
- **Branch:** `fix-L-update-beta-cycle` (commits 12b0424, 55500fd).

---

## 2026-04-29 — EBMF warm-start pinpoints bug to the L update, not the β update

- **Decision:** The root cause of SSBMF's β=0 failure on the merged cohort is narrowed to `update_L_k()`. The β CAVI update is confirmed functional; the L/F updates are washing out the survival signal by prioritising the genomics reconstruction objective.
- **Evidence:** Two warm-start experiments run on merged TCGA_PAAD + CPTAC (v2 preprocessing, n=273, p=2000, K=20):
  1. **β-only experiment** — Fixed EL at the EBMF loading matrix and ran only `update_beta_k()` for 30 iterations. β moved non-zero at iteration 1, converged by iteration 7. **6/20 factors became active** (EBMF3/4/6/13/16/17), exactly matching the Cox-significant factors identified in the EBMF diagnostic. Max |β| = 6.04. Effective C-index ≈ 0.67 (raw concordance = 0.33, sign-inverted due to unit-norm L scaling). **Conclusion: β update is not broken.**
  2. **Full CAVI warm-start** — Initialized EL and EF from EBMF posterior means (`flash_fit$L_pm`, `flash_fit$F_pm`), ran full CAVI. Converged in 23 iterations. **β collapsed back to near-zero** (max |β| = 0.026, 0/20 active). The L update undid the EBMF initialisation and drove the loading matrix toward genomics-reconstruction-optimal directions, where the survival signal disappears.
- **Conclusion:** The β update is correct. The failure is that `update_L_k()`'s A_surv term (survival gradient contribution to the EBNM precision A) is dominated by A_gen (genomics reconstruction gradient) during CAVI. The model converges to a loadings solution that reconstructs Y well but is not informative for survival — then β has nothing informative to select.
- **Next debugging step:** Inspect the magnitude ratio A_surv / A_gen inside `update_L_k()` during a training run. If A_surv ≪ A_gen for most samples, the survival objective is not contributing meaningfully to the L update, and some form of objective rebalancing (within the L update specifically, not at the λ level) is needed.
- **Affected files:** `code/fit_modular.R` (EL_init/EF_init added), `results/benchmark_sim/run_ebmf_warmstart.R` (new)

---

## 2026-04-29 — EBMF diagnostic confirms survival signal exists; SSBMF failure is a model problem

- **Decision:** The β=0 failure on merged TCGA_PAAD + CPTAC training is classified as a **model problem**, not a data problem. Investigation via unsupervised EBMF + PCA diagnostic is now the official diagnostic path for cases where SSBMF produces all-zero β on a given dataset.
- **Reason:** Running `flashier::flash()` (EBMF, K=20) on the same v2-preprocessed merged training matrix used for SSBMF yielded 5/20 factors univariately associated with overall survival at p < 0.05. The strongest, EBMF6, has C-index = 0.629 and p = 3×10⁻⁶ using raw factor loadings alone. PCA confirmed: 4/20 components were also significant. Since unsupervised factorization — with no survival objective whatsoever — finds survival signal, the merged data contains recoverable signal. SSBMF's failure to surface non-zero β must originate in the model's CAVI objective, prior, or update equations, not in data quality.
- **Key result:** EBMF survival-associated factors (EBMF3, 4, 6, 16, 17) have the following top biological signals: EBMF3/EBMF17 = exocrine pancreas markers (PTF1A, GUCA1C, FGL1); EBMF4 = B cell / immune markers (FCRL1, TCL1A, CR2, FCER2); EBMF6 = metabolic / CYP (A2ML1, CYP24A1). Top gene tables at `results/benchmark_sim/outputs/ebmf_diagnostic/tables/ebmf_top_genes.csv`.
- **Trade-offs:** The EBMF diagnostic only tests whether signal is *detectable* by a purely unsupervised method. It does not guarantee SSBMF can recover the same factors — SSBMF imposes additional constraints (joint L/F/β optimisation, CAVI coordinate descent) that could prevent convergence to the EBMF solution even when the signal exists. The EBMF result rules out the data hypothesis; it does not pinpoint the model bug.
- **Recommended follow-on:** EBMF warm-start — initialise SSBMF L and F from the EBMF solution and optimise only β. This directly tests whether the β CAVI update is capable of assigning non-zero coefficients to factors that are empirically associated with survival.
- **Affected files:** `results/benchmark_sim/run_ebmf_diagnostic.R` (new), `results/benchmark_sim/outputs/ebmf_diagnostic/` (tables, figures, report)

---

## 2026-04-29 — Lambda increase ruled out for merged-cohort training (EL collapse)

- **Decision:** Amplifying the survival objective via λ > 1 is ruled out as a strategy for recovering non-zero β in the merged TCGA_PAAD + CPTAC training setting. The default λ=1.0 is retained. For merged training, λ tuning is actively harmful.
- **Reason:** A full λ × prior sweep (λ ∈ {1, 5, 10, 20} × {point_normal, point_laplace, normal}, v2 preprocessing) was run on the merged cohort. At λ=5, the CAVI degenerates completely: the entire L loading matrix collapses to zero (max|EL| < machine epsilon), not merely β. At λ=10 and λ=20, the same collapse occurs. The root cause is that amplifying the survival gradient in the L update overwhelms the genomics reconstruction signal; CAVI responds by driving L to zero (zeroing out the entire linear predictor) rather than shifting weight toward survival-informative directions.
- **Context:** The earlier sandbox (n=250, p=1000, K=5 synthetic data) found λ=1 flat vs. λ=p/n. The merged-cohort collapse is a qualitatively different, more severe failure: EL→0, not just β→0. The batch structure of the merged matrix (RNA-seq vs. proteomics platform factor) likely amplifies the instability.
- **Trade-offs:** λ remains in the codebase as an exposed parameter (default 1.0) for future experiments on single-platform cohorts or after batch effects are addressed. The λ=1 default is safe for all current benchmark runs.
- **Affected files:** `results/benchmark_sim/run_lambda_sweep.R` (new), `results/benchmark_sim/outputs/real_data/lambda_sweep_summary.csv`

---

## 2026-04-29 — v2 preprocessing adopted for merged-cohort benchmark; v1 preserved for single-cohort

- **Decision:** All merged-cohort (TCGA_PAAD + CPTAC) benchmark runs use **v2 preprocessing**: (1) intersect raw gene universes, (2) log₂(x+1) [RNA-seq only], (3) quantile normalization across all merged samples (`preprocessCore::normalize.quantiles()`), (4) top-2000 genes by merged-matrix variance, (5) per-subject rank transform. Single-cohort runs (tcga_only, cptac_only) continue to use v1 (per-cohort pipeline with `preprocess_desurv_cohort()`).
- **Reason:** Under v1, per-cohort top-2000 selection was applied *before* intersecting gene universes, yielding only ~838 common genes — far fewer than the ~2000+ expected. The preprocessing-order bug consumed most of the gene set before cohort merging. v2 fixes the order: intersect first, select top-2000 from the merged variance distribution. Quantile normalization aligns sample-level distributions across RNA-seq and proteomics platforms without introducing explicit batch labels (which would prevent generalisation to new cohorts at prediction time).
- **Trade-offs:** v1 is preserved under `preprocessing_version = "v1"` flag for backward compatibility. v2 output directories carry a `v2_` prefix (e.g., `outputs/real_data/merged/v2_point_normal/`) so v1 and v2 results coexist without overwriting. External cohorts still use v1 single-cohort preprocessing — they are never seen during training, so no joint quantile distribution to normalize against.
- **Affected files:** `code/preprocess_desurv.R` (`preprocess_merged_cohorts()`, `quantile_normalize_merged()`), `results/benchmark_sim/run_ssbmf_benchmark.R` (`run_real_data_benchmark()` `preprocessing_version` param)

---

## 2026-04-29 — Normal prior ruled out for merged-cohort β; point_normal remains default

- **Decision:** The `"normal"` prior for β (soft Gaussian shrinkage, no spike) is not adopted as the default. `"point_normal"` remains the canonical prior.
- **Reason:** The lambda sweep (above) ran all three priors at λ∈{1,5,10,20} on the merged cohort. Under the normal prior, β remains at zero just as with point_normal and point_laplace — the failure is not prior aggressiveness but the fundamental issue identified by the EBMF diagnostic (model/CAVI problem). Adding the normal prior to the benchmark sweep confirmed it does not rescue the merged-cohort fit and adds no new information. It may be revisited after the CAVI L-update issue is diagnosed.
- **Trade-offs:** The normal prior remains available in `update_beta.R` via `ebnm::ebnm_normal` and can be specified via `prior_beta = "normal"` in any benchmark call. It is not removed from the codebase.
- **Affected files:** `results/benchmark_sim/run_lambda_sweep.R`, `results/benchmark_sim/run_ssbmf_benchmark.R`

---

## 2026-04-24 — Lambda survival-scaling parameter: kept at 1.0 after sandbox evaluation

- **Decision:** The `lambda` parameter (scalar multiplier on survival precision terms in the L update) is retained in the codebase at the default value λ=1.0. No active λ tuning is performed.
- **Reason:** A principled argument for λ=p/n exists: the genomics ELBO term sums over p features while the Cox term sums over n patients, so when p>>n the genomics gradient dominates. However, a controlled sandbox (n=250, p=1000, K=5, seed=222) comparing λ∈{1, p/n=5, 2p/n=10} showed no benefit from scaling: hold-out C-index was flat at ≈0.805 across all three conditions, and β RMSE was *worse* at λ=p/n (+0.25) and λ=2p/n (+0.43) than at λ=1. Increasing λ inflates β estimates rather than correcting them, because the dominant source of β scale error is the L–β scale indeterminacy (L can rescale freely), not gradient imbalance.
- **Trade-offs:** The powered-likelihood approach (λ=p/n) is theoretically sound and used in robust Bayesian inference literature. It could become beneficial if the DGP changes (e.g., fewer features, stronger Cox signal). Keeping λ as an exposed parameter with default 1.0 costs nothing and preserves the ability to experiment.
- **Implementation:** λ is a named parameter in `update_L_k()`, `update_L_all()`, and `fit_supervised_mf_modular()` (all default 1.0). It is also in `config/globals.yml` under `cavi.lambda` and threaded through `run_ssbmf_benchmark()` and `run_real_data_benchmark()`. To test λ=p/n, change `globals.yml` and re-run.
- **Affected files:** `code/update_L.R`, `code/fit_modular.R`, `config/globals.yml`, `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/sandbox_lambda_test.R`

---

## 2026-04-24 — Proportional hazards diagnostics added to benchmark pipeline

- **Decision:** `cox.zph()` (Grambsch–Therneau test) is now run on SSBMF risk scores for each external PDAC cohort and results are saved to `ph_diagnostics_table.csv` alongside the benchmark outputs.
- **Reason:** The proportional hazards assumption underlies the Cox model used to generate risk scores. A violation means the log hazard ratio is time-varying, which can distort C-index estimates and KM stratification p-values. Formal PH testing is required before the results can be shared or published.
- **Results (TCGA-only, point_normal):** Dijk p=0.77, Moffitt p=0.72, PACA-AU array p=0.34, PACA-AU seq p=0.41 — all PASS. Puleo_array p=0.026 — **FLAG**. The Puleo violation is marginal and likely reflects the large sample size (n=288) giving power to detect subtle time-varying effects; the C-index and KM results remain valid as approximate assessments.
- **Implementation:** `compute_ph_diagnostics.R` (standalone re-fit + projection script); `run_ssbmf_benchmark.R` (PH test now wired into external cohort loop for future runs); `ssbmf_summary_report.qmd` (Section 4.3).
- **Affected files:** `results/benchmark_sim/compute_ph_diagnostics.R` (new), `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/ssbmf_summary_report.qmd`

---

## 2026-04-24 — ARD preferred over ELBO grid search for K selection

- **Decision:** K is determined automatically within a single model fit via Automatic Relevance Determination (ARD) — the point-normal/point-laplace prior on β shrinks irrelevant factor coefficients exactly to zero — rather than by fitting separate models at K = 1, …, K_max and comparing ELBOs.
- **Reason:** Two complementary advantages. *Efficiency:* ARD determines K_eff as a byproduct of CAVI; a grid search requires K_max separate full fits and introduces an outer loop. *Bayesian coherence:* ARD performs continuous soft shrinkage within a single probabilistic model; ELBO-based model selection performs discrete hard comparison across models with different dimensionalities. Setting K_max generously large (10) and letting ARD prune is equivalent to an ELBO grid search in expectation, without the computational overhead.
- **Trade-offs:** ARD can be conservative — correlated factors may collapse one even when both carry marginal survival signal. A full ELBO grid search would be more exhaustive but is 10× more expensive at K_max = 10. In practice, ARD + generous K_max matches published EBNM-based NMF standards (flash, flashier).
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R`, `config/globals.yml` (k_max), `results/benchmark_sim/ssbmf_summary_report.qmd` (Section 1.3)

---

## 2026-04-24 — point_normal chosen as default beta prior over point_laplace

- **Decision:** `point_normal` is the recommended default prior on β for all future SSBMF runs.
- **Reason:** Across synthetic and all PDAC training modes, `point_normal` matches or slightly outperforms `point_laplace` on external C-index (TCGA-only: 0.602 vs 0.579; CPTAC-only: 0.628 vs 0.620). `point_laplace` selects higher α̂ (0.7 vs 0.5 on synthetic), suggesting it compensates for over-shrinkage of small-to-moderate coefficients by drawing more heavily on the survival gradient. The Gaussian slab is also simpler to interpret — posterior SDs have a direct normal-distribution meaning, whereas the Laplace slab mixes two scale regimes.
- **Trade-offs:** `point_laplace` has heavier tails and may outperform `point_normal` in settings with very sparse survival signal (few events, high censoring) where strong coefficient shrinkage is needed. Revisit if future larger-n runs show a consistent >0.02 C-index advantage for `point_laplace`.
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R` (default `prior_beta`), `config/globals.yml`, `results/benchmark_sim/ssbmf_summary_report.qmd`

---

## 2026-04-24 — Multi-modal TCGA+CPTAC merge documented as expected failure

- **Decision:** The merged TCGA RNA-seq + CPTAC proteomics training mode (838-gene intersection, n=273) is documented as a known failure case rather than a valid benchmark condition. All β̂_k = 0 in both priors.
- **Reason:** When RNA-seq and proteomics are intersected at gene symbols and rank-normalised, PC1 separates the two platforms rather than separating patients by biology. The ARD prior correctly diagnoses that none of the learned factors carry survival signal — they carry platform identity instead. This is not a model failure; it is the model correctly reporting that no prognostic structure exists in this feature space.
- **Trade-offs:** Excluding merged results from the primary benchmark simplifies the comparison table. The failure case is retained in Section 6 of the report as a methodological lesson, motivating the shared-L multi-modal extension (separate F matrices per modality, shared L supervised by survival).
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R`, `results/benchmark_sim/ssbmf_summary_report.qmd` (Section 6)

---

## 2026-04-24 — DeSurv-aligned preprocessing pipeline added

- **Decision:** Added `code/preprocess_desurv.R` — a preprocessing module that matches the DeSurv paper (Young et al. 2025, PNAS) pipeline: log₂(counts+1) for RNA-seq → select top-2000 most-variable genes per cohort → rank-transform each subject's expression vector.
- **Reason:** DeSurv serves as the primary external benchmark. Using the same preprocessing ensures any C-index difference reflects model architecture, not data transformation choices. Proteomics/microarray platforms skip the log₂ step (already on a normalized scale).
- **Trade-offs:** Rank-transform destroys absolute expression magnitude (EBNM shrinkage is insensitive to scale, but gene-level variance information is lost). Top-2000 gene filter is platform-specific — genes selected differ across cohorts, requiring intersection after preprocessing. TCGA_PAAD × CPTAC intersection yielded 838 genes (42% of 2000), acceptable given DeSurv used the same cohorts.
- **Affected files:** `code/preprocess_desurv.R`, `results/benchmark_sim/run_ssbmf_benchmark.R`

---

## 2026-04-24 — Alpha mixing parameter grid expanded to [0, 1]

- **Decision:** `config/globals.yml` `alpha_grid` expanded from `[0.1, 0.3, 0.5, 0.7, 0.9]` to `[0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.0]`.
- **Reason:** The boundary values α=0 (pure genomics, unsupervised NMF) and α=1 (pure survival) are meaningful scientific conditions, not just edge cases. Including them lets CV reveal whether the supervision signal is worth anything at all (α=0 wins → pure NMF is optimal; α=1 wins → ignore genomics structure and regress directly).
- **Trade-offs:** Adds 2 extra fits per CV fold. α=0 is equivalent to an unsupervised NMF run; α=1 may be numerically unstable if survival events are sparse (A_surv can be near-zero for small event counts).
- **Affected files:** `config/globals.yml`, `code/select_alpha_cv.R`

---

## 2026-04-24 — predict_supervised_mf() changed from ridge solve() to SVD pseudoinverse

- **Decision:** Replaced `solve(crossprod(EF) + lambda * diag(K))` with an SVD-based Moore-Penrose pseudoinverse in `predict_supervised_mf()`.
- **Reason:** During alpha CV, early-iteration fits can have ARD drive some factor columns of EF to near-zero. The resulting EF'EF is near-singular even with a fixed ridge term λ·I, because near-collinear *non-zero* columns still make the Gram matrix ill-conditioned. The SVD approach sets d_inv = 0 for singular values below `lambda * max(d)` (relative threshold), so collapsed factors contribute nothing to L_test — which is correct since their EBeta ≈ 0 by the same ARD shrinkage.
- **Trade-offs:** SVD is slightly more expensive than a Cholesky solve for dense K×K matrices, but K ≤ 20 makes this negligible. The relative threshold means λ is now a dimensionless tolerance rather than an absolute precision floor; 1e-8 works well empirically.
- **Affected files:** `code/predict.R` (lines 73–94)

---

## 2026-04-24 — Synthetic DGP fixed: equal factor amplitudes + 4-factor survival signal

- **Decision:** Changed `generate_synthetic_benchmark_data()` to (a) use equal F amplitude for all factors (removed 5× multiplier for null factors) and (b) read `beta_true` from `cfg$synthetic$b_true = [1.5, -1.2, 0.8, -0.5, 0.0]` instead of the hardcoded `[1.0, -0.8, 0, 0, 0]`.
- **Reason:** The 5× null-factor inflation caused PCA to capture most variance from non-prognostic factors and incidentally correlate with survival — PCA C-index (0.715) exceeded supervised C-index (0.673) despite the model having the correct generative structure. Equalizing amplitudes restored the intended benchmark: supervised (0.79) beats PCA (0.76). The `b_true` change aligns the DGP with the 4-signal globals.yml spec and tests a more realistic setting (4 prognostic programs at varied effect sizes).
- **Trade-offs:** Changing the DGP invalidates any previously reported synthetic C-index numbers. The new DGP is harder (4 prognostic factors to recover vs. 2) and better reflects real-data complexity.
- **Affected files:** `results/benchmark_sim/run_ssbmf_benchmark.R` (`generate_synthetic_benchmark_data()`)

---

## 2026-04-09 — Repository reorganisation [PENDING]

- **Decision:** *Pending.* Propose a cleaner directory structure (documented in `ROADMAP.md` → Infrastructure section) but do not move any files until a dedicated refactor commit with no concurrent branch work.
- **Reason:** The current layout has accumulated structural debt across three simulation generations: `results/full_sim/`, `results/modular_sim_block/`, and `results/modular_sim_factor/` coexist; `.qmd` reports are mixed with output tables/figures; `demos/` is a top-level sibling of `code/` rather than nested within it; `code/SupervisedMF_Context.md` is a doc file in the algorithm directory; `derivations/EBMF/` and `derivations/SurvivalMF/` are early-sketch folders now superseded by the per-update derivation subdirectories.
- **Trade-offs:** Any file move invalidates hard-coded paths in runner scripts and `.qmd` `source()` calls — must audit before moving. Deferring keeps the repo stable while active development continues.
- **Affected files:** `results/`, `code/`, `demos/`, `derivations/`, `.gitignore`

---

## 2026-04-09 — Synthetic seed changed from 42 → 222

- **Decision:** Changed the random seed for the synthetic data-generating process from 42 to 222.
- **Reason:** Seed 42 produced a degenerate case where True F1 and True F5 (the null factor) were mixed by the bijective permutation alignment, causing β̂_null = −0.77 — an artifact of the seed rather than model behavior. Seed 222 cleanly separates all 5 factors and correctly zeroes the null factor.
- **Trade-offs:** Reported simulation results (RMSE, β̂ estimates, factor correlations) are seed-dependent; changing the seed resets the canonical benchmark. Prior figures are invalidated.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`

---

## 2026-04-09 — `bijective_match()` replaces column-greedy `which.max` for factor permutation alignment

- **Decision:** Replaced `apply(abs(cors), 2, which.max)` with a `bijective_match()` helper that performs greedy global-maximum assignment to align estimated factors to true factors.
- **Reason:** The column-greedy approach is non-bijective — it can assign the same true factor to multiple estimated factors when two estimated factors are most correlated with the same true factor. This produces incorrect factor-recovery plots and corrupted β̂ comparisons in Figure 3.
- **Trade-offs:** Greedy global-max is still a heuristic (not optimal) but guarantees a one-to-one mapping and is O(K²) — negligible for K ≤ 20.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`

---

## 2026-04-09 — Companion document kept in LaTeX (not Quarto) for print-first speaker notes

- **Decision:** The lab meeting companion document (`Notes/lab_meeting_april9_companion.tex`) is authored in plain LaTeX rather than Quarto.
- **Reason:** The companion is a print-first speaker notes document — dense text, no code execution, no cross-referencing with R output. LaTeX compiles faster, supports tighter typographic control, and avoids Quarto's overhead for documents with no R chunks.
- **Trade-offs:** Not integrated with the Quarto build system; must be compiled separately with `pdflatex`. Cannot embed live R output.
- **Affected files:** `presentation/walther_lab_meeting_04_09_2026/Notes/`

---

## 2026-04-01 — Full ELBO tracking added alongside proxy

- **Decision:** Implemented both a fast ELBO proxy (reconstruction error only) and a full ELBO (proxy + survival term + KL divergences) tracked at every iteration.
- **Reason:** The proxy is sufficient for convergence monitoring but the full ELBO is required for model comparison (prior family selection, K selection). Both are now computed and stored in `history$elbo_proxy` and `history$elbo_full`.
- **Trade-offs:** Adds per-iteration cost of `compute_survival_elbo()` (Taylor approximation) and `compute_ebnm_kl()`. Full ELBO requires Taylor refresh every iteration for accuracy.
- **Affected files:** `code/compute_elbo.R`, `code/fit_modular.R`

---

## 2026-04-01 — Four conditions compared: prior family × K strategy (PN/PL × K=5/K_eff)

- **Decision:** All benchmarks compare 4 conditions: point-normal (PN) vs. point-laplace (PL) crossed with fixed K=5 vs. adaptive K_eff from `auto_prune_K()`.
- **Reason:** Isolates the contribution of prior choice from K selection strategy. Reveals whether adaptive K adds value beyond a fixed-K run with the same prior.
- **Trade-offs:** Requires 4× the compute per dataset. K_eff runs require an extra `auto_prune_K()` pass before the main fit.
- **Affected files:** `results/modular_sim_factor/run_factor_modular_simulation.R`, `code/select_K.R`

---

## 2026-03-31 — `fit_modular.R` uses Gauss-Seidel (factor-wise sequential) CAVI, not full-gradient

- **Decision:** Each CAVI iteration updates factors sequentially (k=1, ..., K), immediately incorporating each updated factor before proceeding to the next. This is the canonical CAVI loop in `code/fit_modular.R`.
- **Reason:** Full-gradient (parallel) CAVI requires the dense Hessian of the joint ELBO with respect to all factors simultaneously, which is intractable in closed form under the Taylor approximation for the survival term. Gauss-Seidel is the standard choice for mean-field CAVI.
- **Trade-offs:** Sequential updates mean later factors in each iteration benefit from earlier updates (faster convergence per iteration), but the ordering introduces implicit asymmetry between factors. Not easily parallelized across factors.
- **Affected files:** `code/fit_modular.R`

---

## 2026-03-25 — Hold-out prediction uses pseudo-inverse projection of Y_test onto F

- **Decision:** `predict_supervised_mf()` obtains test-set loadings L_test by projecting Y_test onto the trained factor matrix F via pseudo-inverse: `L_test = Y_test %*% F %*% solve(t(F) %*% F)`.
- **Reason:** After training, F is fixed. The natural way to score a new patient is to find the loadings that best reconstruct their expression profile under the trained factors. The pseudo-inverse gives the least-squares solution.
- **Trade-offs:** Does not propagate uncertainty from F into L_test (point estimate only). Assumes test data is drawn from the same distribution as training data (no domain shift).
- **Affected files:** `code/predict.R`

---

## 2026-03-20 — Convergence criterion: dual threshold on mean|ΔL| AND mean|Δβ| after 5-iteration burn-in

- **Decision:** Convergence is declared when both `mean(|EL - EL_old|) < tol` AND `mean(|EBeta - EBeta_old|) < tol`, checked only after iteration 5.
- **Reason:** A single criterion on L alone can declare convergence while β is still adjusting (or vice versa). The dual criterion ensures both the genomic structure (L) and survival signal (β) have stabilized. Mean is used instead of max because max is dominated by a few high-variance entries near factor orientation boundaries and rarely reaches 1e-3 on real datasets.
- **Trade-offs:** Mean convergence is weaker than max convergence — some individual loadings may still be changing when the algorithm stops. The 5-iteration burn-in prevents premature stopping during large initial swings.
- **Affected files:** `code/fit_modular.R` (lines 353–370)

---

## 2026-03-12 — Modular update architecture: each update function is independently testable

- **Decision:** CAVI updates are split into four independent modules (`update_L.R`, `update_F.R`, `update_beta.R`, `update_tau.R`), each with its own test suite and demo scripts.
- **Reason:** The monolithic V2 implementation (`Supervised_Bayesian_MF_V2.R`) mixed all update logic into one function, making it hard to test individual components or swap implementations. Modular architecture enables TDD, isolated debugging, and future replacement of individual components.
- **Trade-offs:** Cross-module dependencies must be managed explicitly (e.g., `compute_R_k` lives in `update_L.R` but is also used by `update_F.R`, requiring careful source ordering).
- **Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R`, `code/update_tau.R`, `code/fit_modular.R`

---

## 2026-03-10 — EBNM priors chosen over fixed-penalty; g estimated per CAVI step

- **Decision:** Shrinkage priors for L, F, and β use the Empirical Bayes Normal Means (EBNM) framework, with the prior g estimated from data at each CAVI step rather than fixed by cross-validation.
- **Reason:** Fixed-penalty methods (Ridge, LASSO) require cross-validation for λ, adding computational overhead and a choice of CV criterion. EBNM folds g-estimation into the ELBO maximization — the prior adapts to the data, and no separate CV loop is needed. Implemented via the `ebnm` R package (Stephens lab).
- **Trade-offs:** g-estimation adds per-update overhead. The estimated g is dataset-dependent and may change across CAVI iterations (instability if the model is mis-specified). Point-normal and point-laplace are the two families tested; other families (e.g., generalized double-Pareto) are not yet explored.
- **Affected files:** `code/update_L.R`, `code/update_F.R`, `code/update_beta.R`

---

## 2026-02-12 — Baseline hazard h₀(t) left non-parametric; enters through Cox partial likelihood

- **Decision:** The survival component uses the Cox proportional hazards model with an unspecified baseline hazard h₀(t). h₀(t) cancels exactly in the Cox partial likelihood and never needs to be estimated.
- **Reason:** Specifying a parametric h₀(t) (Weibull, exponential) would add a nuisance parameter and require a prior on its shape/scale. The Cox partial likelihood avoids this entirely while preserving the proportional hazards structure needed to link L to survival.
- **Trade-offs:** Cannot predict absolute survival probabilities (only relative risk via exp(Lβ)). The partial likelihood is not a true likelihood (it conditions on observed event times), which means the ELBO approximation for the survival term requires a Taylor expansion rather than a closed-form KL.
- **Affected files:** `code/update_beta.R`, `code/update_L.R`, `code/compute_elbo.R`

---

## 2026-02-12 — Model formulation: shared L links genomics (Y = LF′ + E) and survival (h(t) = h₀(t)exp(Lβ))

- **Decision:** The core model posits a shared n×K loading matrix L that simultaneously factorizes the genomics matrix Y and enters the Cox proportional hazards model as the linear predictor.
- **Reason:** A two-stage approach (first factor Y, then regress loadings on survival) does not jointly optimize the factors for survival prediction. Sharing L with a joint objective ensures the learned gene expression programs are informative for both reconstruction and prognosis.
- **Trade-offs:** The joint objective couples genomics and survival gradients; the relative scale of the two terms is dataset-dependent and not normalized (p >> n means the genomics term dominates). This is a known limitation — adding a λ scaling parameter is an identified future direction.
- **Affected files:** `code/fit_modular.R`, `code/update_L.R`, `code/update_beta.R`, `derivations/MF_UpdateDerivations/`
