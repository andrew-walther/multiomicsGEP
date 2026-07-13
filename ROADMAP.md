# multiomicsGEP — Research Roadmap

> A living task list of potential extensions, methodological improvements, and dissemination
> goals for the multiomicsGEP project. Organized by theme. Add, edit, and check off items
> as the project evolves.
>
> **Status as of 2026-06-16.** Core model complete (modular CAVI, 246/246 tests passing).
> Two model variants fully implemented, benchmarked, and externally validated (5 held-out PDAC
> cohorts across RNA-seq, microarray, and proteomics platforms). Recommended configuration
> finalized: YFB + survival-ranked gene selection, K=7, mean external C=0.636. Multi-cohort simulation
> study complete: shared vs. study-specific factor recovery validated across 3 scenarios and
> 5 arms; proposal and key results in `docs/reports/multicohort_sim_proposal_06_14_26.pdf`.
> See `docs/reports/desurv_alignment_report_05_27_26.pdf` for the real-data benchmark.
> 6/18/2026 lab-meeting deck delivered and merged to main (`905279b`); unsupervised EBMF→Cox
> external baseline added (mean external C=0.564, K=20 — DECISIONS.md 2026-06-15). Adverse/protective
> program directions corrected to the marginal (YF)-projection convention (DECISIONS.md 2026-06-16).
>
> **Recommended configuration — D4:** YFB (η = (YF)β), DeSurv-aligned gene selection
> (combined mean+variance rank, top-3000 per cohort before per-platform z-standardization,
> 2064 genes after intersection), K=7 (CV-selected with biological floor K≥3), no cohort
> indicator. Mean external C=0.636 across 5 held-out PDAC cohorts. Identifies 2 active gene
> expression programs (K_eff=2). By marginal (YF)-projection survival association — the convention
> used in the talk and the biologically sensible one — **Program 7 is adverse** (MET/ITGA3/glycolytic)
> and **Program 3 protective** (epithelial). The joint-β signs (β̂₇=−0.041, β̂₃=+0.011) are *opposite*
> the marginal direction, a suppression effect among correlated programs (see DECISIONS.md 2026-06-16).
> New-patient scoring: η_new = Y_new F β̂ (exact, no approximation). Factor weight matrix F
> provides directly interpretable gene signatures for pathway enrichment.
>
> **Sensitivity — D3:** LB (η = Lβ), same DeSurv gene set, K=7, mean external C=0.622.
> Reproduces the same 2-program adverse/protective structure, confirming the finding is not
> specific to the YFB parameterization.
>
> **Key preprocessing finding:** Per-platform z-standardization applied before merging is
> required for mixed RNA-seq + proteomics training. 10 of 12 non-per-platform preprocessing
> configurations collapse to β=0 (K_eff=0). YFB point_normal K-CV still returns C=0.5 for
> all K (spike-small-n collapse; open item). Both models use SVD initialization with
> post-convergence sign check (global β negation if training C<0.5).

---

## How to Use This Document

Each item follows this format:

```
- [ ] **Title** `[Priority: High/Medium/Low]` `[Effort: Small/Medium/Large]`
  Description of what this involves and why it matters.
  *Notes: dependencies, caveats, or relevant files*
```

**Priority:** How much this would strengthen the work (High = strengthens core contribution; Medium = meaningful addition; Low = nice-to-have)

**Effort:** Rough implementation cost (Small = hours; Medium = days; Large = weeks or HPC run)

Move completed items to the [Completed](#-completed) section at the bottom.

---

## 🔥 Immediate Priorities

- [x] **Phase 3 (K-parsimony curve on real data)** *(Complete, 2026-07-13, branch `phase3-k-parsimony`)*
  Built `results/benchmark_sim/run_k_parsimony_curve.R`: refits YFB D4 (per-platform z-std, DeSurv
  combined-rank gene selection, top-3000 per cohort, no cohort_id) at K ∈ {2,3,4,5,7} and re-runs
  external validation against the same 5 held-out PDAC cohorts for each K — the curve the internal
  CV table doesn't give us. **Result: K=7 is not free to shrink, in this single-seed comparison.**
  Mean external C-index: K=2→0.541, K=3→0.594, K=4→0.541, K=5→0.596, K=7→**0.627** (SE=0.020) — no
  smaller K reaches within 1 SE of K=7 (margin 0.607). K=7 wins on every one of the 5 individual
  external cohorts, not just on average. The curve is non-monotonic (K=4 ties K=2 despite one more
  factor) and K_eff doesn't track external performance (K=5 has K_eff=3 but underperforms K=7's
  K_eff=2). **Caveat:** each K used a single seed/init; K=2 and K=4 converge suspiciously fast (7-9
  iterations) to a tiny beta_max, consistent with the Phase 2 CAVI factor-collapse failure mode —
  a multistart rerun at those K values is needed before treating them as a hard ceiling (K=5's
  result, which converges normally, is less suspect and already falls short of K=7 on its own).
  Conclusion: the K=7-vs-DeSurv's-K=3 gap is not simply "extra factors we could prune for free" —
  refitting at smaller K directly shows real performance loss in this run, though see the caveat
  before treating K=2/K=4 specifically as conclusive. Details: DECISIONS.md 2026-07-13.

- [ ] **Re-verify Phase 3's K=2/K=4 results with multistart before treating as conclusive**
  `run_k_parsimony_curve.R`'s K=2 and K=4 fits converge in 7-9 iterations with beta_max≈0.009 —
  consistent with the CAVI factor-collapse failure mode documented for Phase 2. A best-ELBO
  multistart rerun at K=2/K=4 (and ideally K=3/K=5 too, for consistency) would confirm whether these
  specific numbers reflect a genuine ceiling or a collapsed single-seed fit. Low urgency: K=5, which
  converges normally, already falls short of K=7 on its own, so the qualitative conclusion likely
  survives — this is about tightening the specific K=2/K=4 numbers, not re-opening the headline
  finding.

- [x] **Phase 2 (joint model vs. two-step value-add, survival-strength sweep)** *(Complete, 2026-07-12;
  comprehensively extended same day)*
  Post-lab-meeting action plan Phase 2, branch `validation-two-step`. Simulation sweep scaling the
  true prognostic effect from 0 to large, comparing YFB (joint) against new PCA+Cox and existing
  EBMF+Cox two-step baselines, across **10 seeds and 4 DGP scenarios** (`default` real EBMF
  templates, `sparse_synthetic`, `low_snr`, `high_K`; extended from an initial 5-seed/1-scenario
  pass). Honest results: equivalence at zero signal holds within ~1-2 SE; under `default`/`low_snr`/
  `high_K` the joint model's advantage is not immediate (tied with baselines at very weak signal) but
  emerges and grows reliably from moderate signal onward; the α=0 internal control is exact —
  max|F diff|=0 across all 240 scenario×strength×seed combinations. **`sparse_synthetic` reverses the
  ordering** (YFB ~5 SEs below both two-step baselines at strength=4) — root-caused with a persisted
  diagnostic script to a genuine, understood CAVI factor-collapse vulnerability shared by LB and YFB
  (see the follow-up item below), not a survival-coupling bug and not an amplitude-hierarchy effect
  (amplitude is identical across all shared factors in every scenario); the real differentiator is
  loading density/structure (real EBMF templates are dense and cross-correlated, `sparse_synthetic`
  is exactly sparse and disjoint). YFB has not collapsed in any real-template scenario tested
  (`default`/`low_snr`/`high_K`, 0/180 fits), but **LB has, occasionally, even with real templates**
  (`low_snr`/`high_K`, 10/180 fits) — so this is not purely a contrived-synthetic-data phenomenon.
  Full report: `docs/reports/joint_vs_twostep_sweep_07_12_2026.{qmd,pdf,html}`. New reusable code:
  `results/multi_cohort_sim/fit_pca_cox.R`, `results/multi_cohort_sim/plot_survival_strength_sweep.R`,
  `results/multi_cohort_sim/diagnose_factor_collapse.R`. Details: DECISIONS.md 2026-07-12.

- [ ] **Investigate CAVI factor-collapse vulnerability (discovered via Phase 2's comprehensive sweep)**
  Both LB and YFB's shared `update_L.R`/`update_F.R` joint-CAVI fitting can collapse multiple factors
  to exactly zero and get permanently stuck when factors have near-equal amplitude and disjoint
  (non-overlapping), low-density support — confirmed independent of the survival objective (identical
  dead-factor count at α=0 in every seed tested) and not fixable by switching to random init (makes
  it uniformly worse: 4/4 dead vs. 1-2/4 for SVD init). LB collapses more on average than YFB (mean
  2.4/4 vs. 1.8/4 dead over 5 seeds) but not in every individual seed. EBMF avoids this via greedy
  (residual-based) factor-at-a-time fitting; PCA avoids it by doing no sparsity-inducing shrinkage.
  Candidate fixes to evaluate: (1) multistart best-ELBO selection via the existing
  `code/fit_modular_multistart.R` (built but not yet tested against this specific failure mode — an
  attempted test hit an unrelated setup error, not yet resolved); (2) a greedy/sequential init variant
  analogous to EBMF's. Priority: elevated from "no evidence it affects real fits" — LB has shown this
  collapse occasionally on real-template scenarios (`low_snr`, `high_K`), not only on contrived
  synthetic data, so this is a real (if so-far-occasional) robustness gap worth investigating, not
  purely hypothetical. Details: DECISIONS.md 2026-07-12 (Phase 2
  same-day follow-up entry).

- [x] **Phase 1 (objective normalization, λ retirement, train/test preprocessing fix)**
  *(Complete, merged to `main` 2026-07-12)*
  Post-lab-meeting action plan Phase 1 (`docs/plans/ssbmf_post_lab_meeting_action_plan_07_08_2026.md`):
  normalized the unnormalized ~p-fold genomics/survival ELBO-term scale gap by boosting the
  survival contribution for the ELBO monitor only (not shrinking genomics — shrinking genomics
  collapses LB's L/F to exactly zero via a bilinear feedback loop; boosting β's own precision was
  tried too but found unjustified and reverted, see below); retired the redundant `lambda`
  survival-scale multiplier (`alpha` already plays that role); fixed a train/test preprocessing
  mismatch where external cohorts were rank-transformed while training was per-platform
  z-standardized (the reverse), across 3 benchmark scripts. **Honest net effect on D4 (recommended
  config): external mean C-index 0.636 → 0.627, K=7, K_eff=2 (unchanged)** — the small decline is
  attributable entirely to the preprocessing fix (a genuine bug fix, kept regardless); objective
  normalization itself has zero effect on YFB's fitted output (there was no real per-coordinate
  imbalance in YFB to fix — see the K_eff root-cause entry below). 267/267 tests passing. Full
  investigation (including two rejected normalization directions and the beta-boost correction) in
  DECISIONS.md 2026-07-12.

- [x] **Re-run benchmarks with top_n=2000 to restore single-cohort baseline** *(complete 2026-05-05)*
  `top_n_genes` was inadvertently raised to 5000. Reverted to 2000. K=20 on single-cohort
  data (tcga_only n=144) still gave K_eff=1 and C=0.37–0.50 — K overfitting, not top_n.
  Fixed by adding `k_pdac_single=10` and re-running. See DECISIONS.md 2026-05-05.

- [x] **Restore alpha CV to benchmark runners** *(complete 2026-05-05)*
  `select_alpha_cv()` added to `run_LB_benchmark.R`. The sign correction must not run inside
  CV folds — doing so causes fold-to-fold sign inconsistency that inflates apparent concordance
  and causes the 1-SE rule to select α=1.0 (degenerate pure-survival mode). Fix: set
  `sign_correction=FALSE` in CV fold fits. All three LB training modes select α=0.50 via 1-SE
  rule. YFB runner uses fixed α=0.50 (YFB-compatible alpha CV not yet implemented).

- [x] **Fix K overfitting and re-run benchmarks** *(complete 2026-05-05)*
  `k_pdac_single=10` added; all 6 modes re-run. **K tuning did not recover archived baseline.**
  LB tcga_only K=10: K_eff=2, external C=0.34–0.43 (worse than K=20 which gave 0.47–0.50).
  More active factors → more anti-concordant external predictions. The archived 0.63–0.65 was
  a lucky PCA direction alignment, not a stable property. Root cause is A_surv/A_gen structural
  imbalance — must be addressed directly. See DECISIONS.md 2026-05-05.

- [x] **Fix YFB β→0 collapse on merged TCGA+CPTAC** *(Complete — Phase 1, 2026-05-06)*
  Per-platform z-standardization (normalize TCGA_PAAD and CPTAC separately before merging)
  plus removal of the per-subject rank transform resolves the collapse. With K=3, merged YFB
  training achieves external C=0.53–0.67 (median 0.64), matching or exceeding the DeSurv
  benchmark range of 0.60–0.65. Both priors agree at K=3 (factors strong enough to escape
  the spike component). Implementation: `--per-platform-norm --no-rank` flags in
  `results/benchmark_sim/run_YFB_benchmark.R`; `per_platform_standardize=TRUE` and
  `rank_transform=FALSE` in `preprocess_desurv_cohort()`.
  *Files: `code/preprocess_desurv.R`, `results/benchmark_sim/run_YFB_benchmark.R`,
  `config/globals.yml` (k_pdac_yfb_merged: 3). Results: `outputs/YFB_benchmark_perplatform/`.*

- [x] **Cohort dummy-variable extension: shared L with cohort-specific F columns** *(Complete — 2026-05-22)*
  Corner-point encoded cohort indicator columns appended to L; F_cohort rows estimated via
  Normal conjugate update; cohort_id and sigma_F_cohort parameters in both models.
  Synthetic: offset absorption |cor|=0.995, factor recovery +88–109% vs baseline.
  Real PDAC: LB_cohort gains on RNA-seq external cohorts (PACA_AU_seq +0.019, PACA_AU_array +0.007)
  but loses on Dijk/Puleo (mean C 0.604 vs LB_base 0.618). YFB β→0 unchanged.
  Full evaluation: `docs/reports/cohort_lmm_benchmark_report.qmd` and DECISIONS.md 2026-05-22.
  *Files: `code/update_F_cohort.R`, `code/compute_elbo.R`, `code/fit_modular.R`, `code/fit_cox_on_yf.R`,
  `tests/test_update_F_cohort.R`, `tests/test_fit_modular_cohort.R`, `tests/test_fit_yf_cohort.R`.
  Validation: `results/cohort_lmm_sim/run_synthetic.R`, `results/benchmark_sim/run_cohort_lmm_benchmark.R`.*

- [x] **YFB β→0 on merged data: all CAVI-local fixes exhausted — formally closed** *(2026-05-25)*
  Frozen-F β pre-conditioning (N_frozen parameter, V9–V11) tested as the final CAVI-local
  strategy. All three frozen-F variants converge at iteration 9 — during the frozen phase,
  before EF unfreezes. Root cause: SVD initialization on merged TCGA+CPTAC is itself
  platform-dominated, so ZF_SVD has no survival correlation and β→0 during the frozen phase.
  Exhaustive record: V0 (baseline) through V11 (frozen-F N=30+warm-start) all fail.
  Practical solutions: (1) per-platform YFB (Phase 1, C_dijk≈0.573); (2) LB_cohort K=5
  (mean C=0.614); (3) LB_base K=20 (mean C=0.618). Possible future fix: initialize EF
  from a single-cohort TCGA-only YFB fit (not platform-dominated SVD) and freeze during
  merge training. Not implemented — deferred as low-priority.
  *Files: `code/fit_cox_on_yf.R` (N_frozen param, backward compat), `tests/test_fit_yf_frozen_f.R` (8 tests),
  `results/benchmark_sim/run_yfb_beta_fix_diagnostic.R` (V9–V11 added). DECISIONS.md 2026-05-25.*

- [x] **Evaluate cohort extension benefit on lower-K regimes** *(2026-05-25)*
  Run `run_cohort_lmm_benchmark.R --low-k` (K_LB=5, K_YFB=3, 300 iters).
  Key finding: at K=5, LB_base also collapses to β→0 (K_eff=0). LB_cohort at K=5 rescues
  β (K_eff=3, mean C=0.614), nearly matching LB_base K=20 (mean C=0.618) with 4× fewer
  biological factors. The cohort extension is *essential* in the interpretable low-K regime
  and *neutral to marginal* at K=20 where ARD handles platform effects implicitly.
  Practical recommendation: for merged multi-platform PDAC analysis, use K=5 + cohort_id.
  *Files: `results/benchmark_sim/run_cohort_lmm_benchmark.R` (--low-k flag added),
  `results/benchmark_sim/outputs/cohort_lmm_benchmark_low_k/`. DECISIONS.md 2026-05-25.*

- [x] **Merged-cohort comprehensive preprocessing benchmark (18 configurations)** *(Complete — 2026-05-25)*
  Extended 6-config comparison to 18 by adding joint quantile-no-rank, joint z-standardization,
  and log-only preprocessing. Biological K floor K_final = max(K_1se, 3) applied throughout.
  Winner: M5 (YFB × per-platform z-std × no cohort_id, K=3), mean external C=0.626, K_eff=2/3.
  Key finding: per-platform z-standardization is the *only* viable preprocessing for mixed
  RNA-seq + proteomics training data. 10 of 12 non-per-platform configs collapse to β=0 (K_eff=0).
  Rank transform is not the driver of YFB β→0 on joint-QN data (M13/M14 also collapse without it).
  Recommended for manuscript: M5 (primary), M4 (sensitivity). See DECISIONS.md 2026-05-25.
  *Files: `results/benchmark_sim/run_merged_kcv.R`, `results/benchmark_sim/run_merged_benchmark.R`,
  `code/preprocess_desurv.R` (normalize_method param), `docs/reports/merged_benchmark_report.qmd`.
  Results: `outputs/merged_benchmark/merged_benchmark_results_extended.csv`. DECISIONS.md 2026-05-25.*

- [x] **Align gene selection with DeSurv: combined mean+variance ranking, top-3000** *(Complete — 2026-05-27)*
  Implemented combined_rank method in select_top_variable_genes() and per-cohort
  selection (before normalization) in preprocess_merged_cohorts(). D4 (YFB DeSurv-aligned)
  mean external C=0.636 vs M5=0.624 (delta=+0.012, 5 cohorts). DeSurv gene selection
  adopted as new primary config. Also evaluated D5 (YFB DeSurv + cohort indicator, K=8):
  mean C=0.614 (−0.022 vs D4); per-platform z-std already absorbs the platform offset,
  making the cohort indicator counterproductive at K=7–8. See DECISIONS.md 2026-05-27.
  *Files: code/preprocess_desurv.R, results/benchmark_sim/run_desurv_comparison.R,
  docs/reports/desurv_alignment_report_05_27_26.qmd,
  docs/progress_report/SSBMF_Status_Update_5_28_26.qmd*

- [x] **Factor diagnostics report: convergence, survival associations, gene programs** *(Complete — 2026-05-27)*
  Full diagnostic evaluation of D1–D5 on merged TCGA+CPTAC training cohort. All five
  configurations converge to K_eff=2 (two active programs) regardless of K requested (3–8),
  confirming the ARD prior concentrates survival signal into at most 2 directions. D4 active
  factors: Factor 3 (β̂=+0.011) and Factor 7 (β̂=−0.041), both log-rank p<0.05.
  **[Corrected 2026-06-16: by the marginal (YF)-projection direction adopted in the 6/18 deck,
  Factor 7 is *adverse* (MET/ITGA3/glycolytic) and Factor 3 *protective* (epithelial) — opposite
  the joint-β-sign labels shown in this 05-27 report (suppression among correlated programs).
  See DECISIONS.md 2026-06-16.]**
  Signal replicates across all 5 external cohorts. Proportional hazards assumption not
  violated (Schoenfeld global p>0.05). Heatmaps and top-gene tables for all active factors
  provided as starting point for pathway enrichment.
  *Files: docs/reports/desurv_factor_diagnostics_05_27_26.{qmd,pdf}*

- [ ] **A/B comparison: SSBMF vs unsupervised EBMF on real PDAC data** `[Priority: High]` `[Effort: Medium]`
  Fit an unsupervised EBMF model (no survival objective) on the same merged TCGA+CPTAC training
  data and DeSurv gene set. Evaluate both on the same 5 held-out PDAC cohorts (η = YFβ for SSBMF;
  Cox PH on YF factor scores post-hoc for EBMF). Goal: quantify the performance and interpretability
  lift from the supervised component on real data. The multi-cohort simulation provides preliminary
  evidence: in the hybrid scenario (only some factors prognostic — the realistic case) the joint
  model achieves C=0.81 vs EBMF→Cox C=0.72 (+0.09); at equal signal the two-stage baseline nearly
  matches. Real-data confirmation would strengthen the manuscript's core claim.
  Key comparisons: (1) external C-index; (2) factor stability and biological coherence (do the
  same 2 active programs emerge unsupervised?); (3) pathway enrichment concordance.
  **Update 2026-06-16:** comparison (1) is done — `run_ebmf_cox_external.R` gives EBMF→Cox mean
  external C=0.564 (K=20) vs YFB 0.636 across the same 5 cohorts (DECISIONS.md 2026-06-15). The
  same-protocol external comparison is shown in the 6/18 deck. Parts (2) factor stability and
  (3) pathway concordance remain.
  *Files: YFB fit in `results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds`;
  EBMF via `flashier`; prediction via `code/predict.R`; baseline runner `run_ebmf_cox_external.R`.*

- [ ] **Pathway enrichment on D4 active factors** `[Priority: High]` `[Effort: Small]`
  Submit top-weighted genes from D4 Factor 7 (adverse: MET/ITGA3/glycolytic) and Factor 3
  (protective: epithelial) to MSigDB hallmark gene sets / KEGG / GSEA. Top-gene lists already
  exported to `presentation/walther_lab_meeting_06_18_2026/assets/active_factor_genes.csv`.
  Cross-reference with Moffitt basal/classical scores and
  Bailey et al. (2016) four-subtype annotations in TCGA_PAAD to assess concordance with
  established PDAC molecular axes. Patient loadings L̂_{i,3} and L̂_{i,7} serve as continuous
  molecular scores for subtype comparison. This is the primary interpretability task remaining
  before the manuscript methods section can be finalized.
  *Files: results/benchmark_sim/outputs/desurv_comparison/desurv_comparison_fits.rds (D4 fit)*

- [ ] **Prior comparison follow-up** `[Priority: Medium]` `[Effort: Small]`
  Current benchmarks test point_normal vs normal for both models. Key finding: for the YFB model
  on real PDAC data, the normal prior is strictly necessary — the spike-and-slab (point_normal)
  collapses all β to zero because the survival signal is too weak to exceed the spike threshold.
  For the LB model, both priors give identical results at the current signal level. A fuller
  comparison should include point_laplace and evaluate whether a β-only warm-up (fitting β while
  holding L fixed) helps the spike-and-slab prior find the active set.
  *Files: `results/benchmark_sim/run_LB_benchmark.R`, `results/benchmark_sim/run_YFB_benchmark.R`*

---

## 🧬 Model Specification

- [ ] **Incorporate F'Y term in survival model** `[Priority: Medium]` `[Effort: Large]`
  Replace $L_i$ in the Cox linear predictor with $(F^\top y_i)$ — the projection of patient $i$'s
  observed expression onto the factor space. This decouples factor learning from survival supervision:
  L is learned purely from genomics, and survival scoring is done post-hoc via the factor projection.
  Requires a new q(L) CAVI derivation (the survival term drops out of the L update) and a post-hoc
  scoring step. May improve identifiability and reduce the sign-flip issue observed in synthetic
  validation.
  *Notes: High effort — requires new derivation and update to `code/update_L.R`. See `derivations/qL/qL_update_derivation.pdf` for current L update derivation.*

- [x] **EBMF warm-start initialization for SSBMF** `[Priority: High]` `[Effort: Medium]` *(Complete — 2026-04-29)*
  Initialize SSBMF L and F from an EBMF solution and run CAVI from that starting point.
  Two experiments run: (1) β-only with EL fixed → β non-zero at iter 1, 6/20 factors active,
  confirms β update is functional. (2) Full CAVI warm-start → β collapses to zero in 23 iters,
  L/F updates wash out EBMF structure. Root cause: `update_L_k()` A_surv dominated by A_gen.
  See DECISIONS.md 2026-04-29 warm-start entry.
  *Notes: `fit_modular.R` extended with `EL_init`/`EF_init` params. Driver: `run_ebmf_warmstart.R`.*

- [x] **Investigate and resolve β=0 collapse in LB model** *(Complete — 2026-04-29)*
  Instrumentation confirmed a ~5000× scale imbalance between survival and genomics precision
  contributions to the L update at iteration 1 (A_surv/A_gen ~ 0–2e-4). Three fixes were
  implemented and evaluated: (1) Gauss-Seidel update ordering (β before L in the inner loop);
  (2) β-only warm-up iterations with L held fixed; (3) rescaling survival and genomics
  precision contributions to comparable magnitude. On merged TCGA+CPTAC: 2/20 factors active,
  ELBO monotone. Training-side β=0 resolved. External C-index was mixed (1/5 cohorts improved),
  motivating the Cox-on-YF reformulation below.
  *Files: `code/fit_modular.R`, `code/update_L.R`*

- [x] **Cox-on-YF (YFB) reformulation** *(Complete — 2026-05-04, merged to `main`)*
  See completed section below.

- [x] **Add λ scaling parameter to balance genomics vs. survival objectives** `[Priority: High]` `[Effort: Medium]` *(Implemented and evaluated — fixed at λ=1.0)*
  λ is implemented as an exposed parameter in `update_L_k()`, `update_L_all()`, and
  `fit_supervised_mf_modular()` (default 1.0), and registered in `config/globals.yml`.
  A controlled sandbox (n=250, p=1000, K=5) tested λ∈{1, p/n=5, 2p/n=10}: hold-out C-index
  was flat at ≈0.805 across all conditions; β RMSE was *worse* at λ=p/n (+0.25) and λ=2p/n (+0.43).
  The dominant β scale error is L–β scale indeterminacy, not gradient imbalance. λ=1 retained.
  *Notes: To experiment, change `cavi.lambda` in `config/globals.yml`. See `DECISIONS.md` entry 2026-04-24 for full analysis. Sandbox: `results/benchmark_sim/sandbox_lambda_test.R`.*

---

## 🌍 Validation & Generalization

- [x] **Cross-cohort validation: train on TCGA-only or CPTAC-only, evaluate on 5 external cohorts** `[Priority: High]` `[Effort: Medium]` *(Complete — see `results/benchmark_sim/outputs/real_data/`)*
  Merge CPTAC and TCGA_PAAD (batch-corrected) as a training set and hold out Dijk,
  Moffitt_GEO_array, PACA_AU_array, PACA_AU_seq, and Puleo_array as independent test cohorts.
  This tests generalisation of learned gene expression programs across platforms (RNA-seq,
  microarray, proteomics) and patient populations. More training data than the current 80/20
  within-cohort split. Requires R4 (inter-dataset normalisation) as a prerequisite.
  *Notes: `code/train_test_split.R` and `code/predict.R` are already in place for the projection step. Batch correction (ComBat or similar) not yet implemented — add to R4. See `results/modular_sim_factor/PDAC/` for current within-cohort results.*

- [x] **Inter-dataset normalisation before cohort merging** `[Priority: High]` `[Effort: Medium]` *(Complete — v2 preprocessing, 2026-04-29)*
  Before merging cohorts (required for R3), apply quantile normalisation or z-score
  standardisation across datasets beyond the current per-cohort column-centring. Addresses
  platform-specific mean/variance shifts between RNA-seq, microarray, and proteomics. Essential
  for any analysis that combines data across assay types.
  *Notes: Implemented as `preprocess_merged_cohorts()` in `code/preprocess_desurv.R`. Pipeline: intersect raw gene universes → log₂(x+1) [RNA-seq only] → `preprocessCore::normalize.quantiles()` across all merged samples → top-2000 by merged-matrix variance → per-subject rank transform. Gated by `preprocessing_version = "v2"` in `run_real_data_benchmark()`. See DECISIONS.md 2026-04-29 entry.*

---

## 📐 Model Selection

- [ ] **Joint (K, α) tuning via Bayesian optimization, to match DeSurv's search procedure** `[Priority: Low-Medium]` `[Effort: Medium]`
  Plan only, not implemented: `docs/plans/joint_k_alpha_bayesopt_plan_07_12_2026.md`. Motivation:
  DeSurv jointly tunes `k, α, λ` via Bayesian optimization; we tune only K via grid-search CV with
  α fixed at 0.5. This confirmed methodology difference is part of why K=7-vs-DeSurv's-k=3 isn't a
  fair apples-to-apples comparison (`DECISIONS.md` 2026-07-12). A joint `(K, α)` BO search (no
  direct analogue of DeSurv's λ exists for us — see the plan doc) would let that comparison stand
  on equal footing. Deferred behind Phase 2/3; revisit sooner if the manuscript (Phase 6) needs a
  defensible answer to "why does SBMF use more factors than DeSurv."

- [x] **K_eff 2→4 root-caused and resolved (not a `beta_threshold` calibration issue)** *(Complete — 2026-07-12)*
  D4's K_eff rose from 2 to 4 after the initial Phase 1 merge, alongside an apparent external
  C-index gain (0.636→0.642). Root cause, confirmed directly (not `beta_threshold` miscalibration
  as first suspected): an implementation detail in Phase 1a's objective normalization
  (`boost_beta=TRUE`) boosted β's own EBNM precision by a factor of p in *both* models, even
  though β's coordinate update has no genomics term competing with it in its own formula in
  either model — the genomics/survival imbalance Phase 1a targets does not structurally exist
  for β at all, so this boost corrected no real imbalance; it only reduced EBNM shrinkage,
  inflating K_eff without reflecting genuine new survival signal. Fix: `boost_beta` now defaults
  to `FALSE` in both `fit_supervised_mf_modular()` and `fit_cox_on_yf()`. Confirmed on real
  data: with `boost_beta=FALSE`, D4's β values are essentially identical to the pre-Phase-1
  baseline (+0.0115/−0.0404 vs. previously +0.011/−0.041) and K_eff is back to 2. The honest
  post-Phase-1 external C-index is **0.627** (down slightly from 0.636, attributable entirely to
  Phase 1c's preprocessing fix, not 1a) — see `DECISIONS.md` 2026-07-12 for the full comparison
  table and `CLAUDE.md` "Current model status" for the corrected headline numbers.

- [x] **Re-evaluate K=7's parsimony now that the objective/preprocessing fixes are settled** *(Complete — 2026-07-12)*
  Fresh `select_K_cv()` run (YFB, D4 preprocessing, K grid 2:10, no floor imposed a priori) under
  the corrected code (`boost_beta=FALSE`, fixed preprocessing): **K=7 is genuine, not an
  artifact.** A real, non-noise jump exists between K=6 (mean C=0.593) and K=7 (0.633) — larger
  than any fold's SE (~0.03) — with K=8 the actual peak (0.651) and K=7 the 1-SE-simplest choice
  tied with it. K=2 through K=6 are all meaningfully worse, not just noisier. Kept K=7. The
  K=7-vs-DeSurv's-k=3 gap is attributable to two confirmed methodology differences (DeSurv jointly
  tunes k/α/λ via Bayesian optimization with a fixed elastic-net penalty; we tune only K via CV
  with α fixed and no penalty at all) rather than evidence that the datasets themselves demand
  more latent structure — a fair comparison needs a matched-protocol re-run (Phase 2/6), not a
  different K-CV on our side. K_eff=2 (survival-active factors) already tracks DeSurv's own ~1
  reasonably well; the parsimony story that holds up is about *effective*, not total, factors.
  Full table and reasoning: `DECISIONS.md` 2026-07-12.
  *Files: `results/benchmark_sim/run_merged_kcv.R` (existing infra), a focused ad-hoc re-run for
  the D4 (YFB DeSurv-aligned) config specifically.*

- [x] **K selection via cross-validated C-index** *(Complete — Phase 4, 2026-05-06)*
  `select_K_cv()` implemented in `code/select_K.R` with 1-SE rule; accepts `model="LB"` or
  `model="YFB"`. Five-fold CV over K∈{2,…,10,15,20} on TCGA_PAAD (n=144):
  - LB normal prior: flat plateau (range=0.032 < per-K SE), K=3 selected (1-SE).
  - LB point_normal: K=8 selected (spike-small-n artefact; C=0.5 at K≤3). Practical default: K=5.
  - YFB point_normal: K=2 (β→0 collapse in all folds, C=0.5 everywhere).
  - YFB normal: K=5 selected (1-SE; C=0.56–0.61 across all K) after sign fix (2026-05-20).
  - YFB point_normal: C=0.5 for all K (β→0 collapse in folds; open item).
  193/193 tests passing (KCV-T15 added for YFB sign fix verification).
  *Files: `code/select_K.R`, `tests/test_select_K_cv.R`, `results/benchmark_sim/run_K_cv.R`.
  Results: `results/benchmark_sim/outputs/K_cv/`. Report: `docs/reports/ssbmf_summary_report_05_06_26.pdf`.*

- [x] **YFB sign correction in K-CV** *(Complete — 2026-05-20)*
  `sign_correction=FALSE` parameter added to `fit_cox_on_yf()`; passed in YFB branch of
  `select_K_cv()`. CV folds deliver raw SVD-oriented EBeta; concordance evaluated via
  `I(-pred$risk_scores)`. After fix: normal prior yields C=0.56–0.61, K=5 selected (1-SE).
  *Files: `code/fit_cox_on_yf.R`, `code/select_K.R`, `tests/test_select_K_cv.R` (KCV-T15).
  Results: `results/benchmark_sim/outputs/K_cv/K_cv_table_YFB_normal.csv`.*

- [ ] **YFB point_normal K-CV collapse** `[Priority: Medium]` `[Effort: Medium]`
  Under the spike-and-slab (point_normal) prior, YFB β→0 in all CV folds regardless of K
  (C=0.5 for all K ∈ {2,…,20}). The reduced training fold size (~115 subjects) is insufficient
  to escape the spike component. This is a separate issue from the sign fix. Possible
  approaches: β-only warm-up iterations before the spike prior activates; adaptive spike
  weight schedule; or accepting that point_normal is not viable for YFB K-CV on small-n data.
  *Files: `code/fit_cox_on_yf.R`, `code/select_K.R`*

- [ ] **YFB preprocessing alignment** `[Priority: Medium]` `[Effort: Medium]`
  Single-cohort and merged YFB training use different preprocessing contracts:
  single-cohort uses the per-cohort DeSurv pipeline (`preprocess_desurv_cohort()`);
  merged uses per-platform z-standardization (`per_platform_standardize=TRUE`,
  `rank_transform=FALSE`). These are not interchangeable — a model trained in single-cohort
  mode cannot be applied with merged-mode preprocessing at prediction time. A unified
  preprocessing API that makes the contract explicit and enforces train/test consistency
  would reduce downstream errors.
  *Files: `code/preprocess_desurv.R`, `results/benchmark_sim/run_YFB_benchmark.R`*

- [ ] **Per-platform noise variance (τ) investigation** `[Priority: Low]` `[Effort: Medium]`
  The current model uses a single shared noise precision τ (estimated per gene, shared across
  subjects and cohorts). In merged multi-platform training (TCGA RNA-seq + CPTAC proteomics),
  the true noise variance almost certainly differs between platforms. A per-platform τ would
  better reflect the generative structure and may reduce the genomics-term scale imbalance
  that drives β→0 in the L update. Requires a new q(τ) derivation for the grouped case.
  *Files: `code/update_tau.R`, `code/fit_modular.R`, `code/fit_cox_on_yf.R`*

- [ ] **Move `load_pdac_raw` and `plot_cohort_loading_heatmap` to `code/`** `[Priority: Low]` `[Effort: Small]`
  Both functions currently live in `results/benchmark_sim/run_ssbmf_benchmark.R` (the data-loading
  hub). They are general enough to belong in `code/` for reuse across runner scripts without
  sourcing the full hub. Requires updating all `source()` calls in runner scripts.
  *Notes: Do in a dedicated commit; coordinate with any active Longleaf paths.*

---

## ⚙️ Infrastructure

- [ ] **Investigate {targets} workflow for pipeline reproducibility** `[Priority: Low]` `[Effort: Medium]`
  Adopt the R `targets` package to manage computational dependencies across the pipeline (data
  loading → fitting → tables → figures → reports). Tracks which outputs are stale and re-runs only
  what is needed. Particularly valuable as the number of cohorts, conditions, and reports grows —
  currently all cohorts are re-run together even when only one has changed. Also facilitates
  Longleaf HPC submission of individual targets as SLURM jobs.
  *Notes: Low priority — investigate feasibility and overhead before committing to a refactor. The `{targets}` package works best when the pipeline DAG is well-defined; this is approximately true for the current runner structure. See `results/modular_sim_factor/run_factor_modular_simulation.R` for the current monolithic runner.*

- [x] **`results/` directory cleanup** `[Priority: High]` `[Effort: Small]` *(complete 2026-05-05)*
  The `results/` directory now mixes code scripts, CSV tables, PNG/PDF figures, and
  QMD/PDF/HTML reports within the same subdirectories, making it hard to navigate and
  prone to path errors. Proposed split (within `results/benchmark_sim/`):
  ```
  results/benchmark_sim/
    scripts/        ← MOVE run_ssbmf_benchmark.R, compute_ph_diagnostics.R
    reports/        ← MOVE ssbmf_summary_report.qmd/.pdf/.html
    outputs/        ← benchmark CSVs + figures (already structured by mode/prior)
  results/modular_sim_factor/
    scripts/        ← MOVE run_factor_modular_simulation.R, run_prior_k_comparison.R
    reports/        ← MOVE .qmd/.pdf/.html out of PDAC/ and synthetic/ subdirs
  ```
  Update all `source()` calls, `quarto render` paths, and README references when moving.
  *Notes: Dedicate a single refactor commit. Don't split across sessions. Update CLAUDE.md
  quick-reference table after moving. Coordinate with any active Longleaf HPC paths.*

- [x] **Repository reorganisation** `[Priority: Low]` `[Effort: Small]` *(results/ portion complete 2026-05-05; code/ and derivations/ portions deferred — low priority)*
  The current directory layout has accumulated structural debt across three simulation generations. 
  Need to clean up the results/ directory as multiple stages of tables/figures/reports are floating 
  around without clear structure Proposed clean structure (no file deletion — move and rename only):

  ```
  code/                     ← algorithm only (update_*.R, fit_modular.R, etc.)
    legacy/                 ← V1, V2 monolithic (already there); ADD Supervised_Bayesian_MF_V2.R
    demos/                  ← MOVE top-level demos/ here (executable examples per module)
  docs/                     ← companion .qmd/.pdf/.html (already clean)
    context/                ← MOVE code/SupervisedMF_Context.md, docs/PDAC_data_audit.qmd here
  derivations/              ← mathematical derivations (already clean per-update)
    archive/                ← MOVE derivations/EBMF/, derivations/SurvivalMF/ early sketches
                               MOVE dated .tex drafts from MF_UpdateDerivations/ (keep only REVISED + Companion)
  results/
    reports/                ← MOVE .qmd/.pdf/.html report files out of modular_sim_factor/
    figures/                ← figures by cohort (already structured)
    tables/                 ← tables by cohort (already structured)
    legacy/                 ← MOVE results/full_sim/, results/modular_sim_block/ here
  tests/                    ← test suite (already clean)
  presentation/             ← slide decks (already clean)
  longleaf_setup/           ← HPC scripts (already clean)
  paper/                    ← manuscript (already clean)
  config/                   ← globals.yml (added in Session 15)
  ```

  Additional housekeeping: add `.obsidian/`, `.claire/`, `.claude/worktrees/` to `.gitignore`
  (IDE/tool artifacts should not be tracked). Remove `.Rhistory` files scattered in subdirectories.
  *Notes: Propose only — do not move files until a dedicated refactor commit. Flag as pending in `DECISIONS.md`. Coordinate with any active branches before moving paths that appear in runner scripts or `.qmd` files.*

---

## ✅ Completed

- [x] **Multi-cohort simulation study: shared vs. study-specific factor recovery** *(Complete — 2026-06-14)*
  Validated that YFB natively recovers the shared/study-specific distinction across 3 scenarios
  (all-shared, hybrid, nothing-shared), 5 arms (YFB/LB ± cohort indicator + unsupervised EBMF→Cox),
  and 5 seeds. Key findings: specificity-classification accuracy 0.97–1.00 (hybrid, YFB); held-out
  C-index 0.81 vs EBMF→Cox 0.72 in the hybrid regime; β false-positive rate 0.03–0.07 in the
  null scenario. Signal-ratio sweep (1×–16× specific-to-shared variance) shows FP rate is bounded
  and does not grow monotonically. Recommendation from simulation: include cohort indicator as
  default and use |β̂| magnitude thresholding rather than binary nonzero classification.
  *Files: `results/multi_cohort_sim/`, `config/globals.yml` (`synthetic_multicohort` block),
  `docs/reports/multicohort_sim_proposal_06_14_26.{qmd,pdf}`. Branch: `multi-cohort-sim`.*

- [x] **YFB K-CV sign fix (2026-05-20)** — `sign_correction=FALSE` parameter added to
  `fit_cox_on_yf()` and passed in `select_K_cv(model="YFB")`. Before fix: C<0.5 for all K
  (Phase C oriented EBeta, then `I(-.)` re-flipped, producing anti-concordant predictions).
  After fix: YFB normal prior yields C=0.56–0.61 across K ∈ {2,…,20}; K=5 selected by 1-SE
  rule on TCGA_PAAD. YFB point_normal remains C=0.5 (spike-small-n collapse; separate open
  item). 193/193 tests passing.
  *Files: `code/fit_cox_on_yf.R`, `code/select_K.R`, `tests/test_select_K_cv.R`.*

- [x] **Training concordance sign correction + alpha CV fix (2026-05-05)** — Both LB and YFB
  were producing anti-concordant predictions (C < 0.5) due to sign ambiguity in SVD initialization:
  the positive-part truncation of the SVD can produce factor loadings whose inner product with β
  is inversely related to survival risk. Fix: after CAVI convergence, compute training concordance;
  if C_train < 0.5, negate β globally. Implemented as `sign_correction=TRUE` parameter in both
  `fit_supervised_mf_modular()` (LB) and `fit_cox_on_yf()` (YFB).
  A related bug in the alpha CV routine was fixed in the same pass: applying the sign correction
  inside CV folds caused fold-to-fold inconsistency, making the 1-SE rule select α=1.0 (pure-
  survival, degenerate). Fix: `sign_correction=FALSE` in fold fits.
  YFB synthetic result: C-index 0.12 → 0.91 (K_eff=4). LB external results (all modes, α=0.50):
  TCGA-only K_eff=2, C=0.55–0.66; CPTAC-only K_eff=3, C=0.55–0.67; merged K_eff=3, C=0.51–0.67.
  *Files: `code/fit_modular.R`, `code/select_alpha_cv.R`, `code/fit_cox_on_yf.R`,
  `code/predict_cox_on_yf.R`. Report: `docs/reports/ssbmf_summary_report_05_05_26.pdf`.*

- [x] **Benchmark pipeline consolidation (2026-05-04)** — Replaced scattered diagnostic scripts
  with two canonical benchmark runners: `run_LB_benchmark.R` (LB model, η = Lβ) and
  `run_YFB_benchmark.R` (YFB model, η = (YF)β). Both follow a 4-section structure (synthetic
  validation + PDAC training + external validation + CSV output) and test point_normal vs normal
  prior side-by-side. Archived 7 retired scripts to `results/benchmark_sim/archive/`. Updated
  `config/globals.yml` with benchmark defaults (K=10, α=0.5, λ=1.0, beta_threshold=0.001).

- [x] **Cox-on-YF (YFB) model implementation (2026-05-04)** — Reformulated linear predictor to
  η = (Y·F)·β, replacing η = L·β. Training and prediction both use the observed expression matrix
  projected onto the factor space, eliminating the train/test mismatch where training used
  EBNM-shrunk loadings but prediction used OLS-projected loadings. Implementation: `code/fit_cox_on_yf.R`,
  `code/predict_cox_on_yf.R`, `code/update_L_surv_YFB.R`, `code/update_F_surv_YFB.R`,
  derivations in `derivations/cox_on_YF/`. Synthetic: C-index 0.605 vs PCA 0.471 (3 active factors).
  Full PDAC benchmark comparison in `docs/reports/ssbmf_summary_report_05_05_26.pdf`.
- [x] **Core modular CAVI implementation** — `code/fit_modular.R` with four independently-tested update modules. 171/171 tests passing. *(Completed March 2026)*
- [x] **Hold-out prediction pipeline** — `code/predict.R` (`predict_supervised_mf()`), `code/train_test_split.R` (`stratified_split()`). 80/20 stratified hold-out. *(Completed March 2026)*
- [x] **Prior family comparison (PN vs PL × K=5 vs K_eff)** — Four-condition benchmark across 7 PDAC cohorts. Point-laplace preferred at fixed K; K selection dominates prior choice. *(Completed April 2026)*
- [x] **Full ELBO tracking** — Both proxy and full ELBO (genomics + survival + KL) tracked per iteration. `code/compute_elbo.R`. *(Completed April 2026)*
- [x] **DeSurv benchmark — synthetic + PDAC cross-cohort** — DeSurv-aligned preprocessing, alpha CV (1-SE rule), SVD pseudoinverse fix, prior sensitivity (point_normal vs point_laplace). Synthetic: supervised 0.79 > PCA 0.76. PDAC external median C-index 0.60 across 5 cohorts. *(Completed April 2026)*
- [x] **Prior sensitivity report** — point_normal vs point_laplace compared on synthetic and all PDAC training modes. point_normal recommended as default. *(Completed April 2026)*
- [x] **24-page DeSurv benchmark report** — `results/benchmark_sim/ssbmf_summary_report.pdf`. Includes ARD justification over ELBO K-grid-search, alpha gradient notation, per-figure takeaways, side-by-side prior comparison (Fig 19), multi-modal failure documented. *(Completed April 2026)*
- [x] **C-index honest reporting fix** — `get_cindex_comparison()` uses model's own `EL %*% EBeta`; corrected `concordance()` direction convention. *(Completed April 2026)*
