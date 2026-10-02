# Review: `paper/prelim/project3-ssbmf.qmd` (prelim proposal chapter, Project 3)

Reviewer: independent review, 2026-10-02. Branch `codex/multimodal-yfb`.
Scope: numbers against source files, claims against evidence, Methods against code, structure, citations, rendering, repo writing rules. No files other than this review were edited. The PDF was rendered from a scratchpad copy using RStudio's bundled Quarto and TinyTeX. It produced 22 pages with no LaTeX or citeproc warnings.

---

## Summary verdict

Nearly every number in the chapter matches its source file. The problems are in what the numbers are used to claim.

1. **The headline single-modality comparison uses a superseded baseline.** ΔC = 0.042 compares the K = 7 joint model with a two-step baseline at K = 20 that uses an unpenalized Cox regression. DECISIONS.md (2026-08-20) calls this an apples-to-oranges comparison and records two fairer ones: K-matched, +0.011 and not significant; independent K = 40 with LASSO, +0.026 and just significant. The chapter reports none of this.
2. **"Joint fitting" cannot explain that gain.** With α_F = 0, survival never enters the L or F updates in `fit_cox_on_yf.R`. The single-modality factorization is therefore unsupervised. The model differs from the two-step pipeline only in its rank, its priors and its empirical-Bayes-shrunk, tempered Cox stage. The Results heading in §3.1 and the Discussion ("fitting … together gave better external discrimination than factorizing first") go beyond what the α_F constraint allows.
3. **The multimodal real-data claims go beyond the evidence.** Pruning removed no factors on real data (K_final = K_init at 3, 7 and 10). The K = 7 and K = 10 fits hit the 500-sweep cap without converging. Concordance depends on K_init (0.557 to 0.685). Contribution 4's "chooses the number of factors without cross-validating K, and it substantially improves external discrimination" is supported only by the 2-seed simulation.
4. **Two internal contradictions and one dangling cross-reference.** These are the "regularized feature precisions" claim, the "avoided … a cross-validated search over K" claim, and the reference to "Section 4.4".

Structurally, the chapter has the expected Introduction, Methods, Results and Discussion. It is missing the elements a committee expects in a *proposal*: aims for the remaining work, a timeline, a simulation design table, and a DeSurv comparison that the Introduction promises. The prose is mostly bullet lists, which reads as notes rather than a chapter.

---

## Must-fix

**M1. The two-step baseline in Table 1 and the ΔC = 0.042 claim are superseded.**
- *Location:* §3.1 (lines 189–200): "against 0.581 for an unsupervised factorization …", "The joint model's concordance was higher in every cohort", "ΔC = 0.042 (95% CI 0.013 to 0.071)". Also Discussion §4.1 (line 288): "consistently across cohorts and platforms".
- *Problem:* The numbers match `external_cindex_ci.csv` and `external_paired_diff_ci.csv`. However, `run_ebmf_cox_external.R` uses its default `K_EBMF_MAX = k_pdac = 20`, so the results CSV has K = 20, and its second stage is an unpenalized `coxph()` on 20 scores with about 140 training events. DECISIONS.md 2026-08-20 records the following:
  - K matched to 7: two-step mean C = 0.623. ΔC = +0.011 (−0.013, 0.033), not significant. The two-step model is ahead in PACA-AU array and PACA-AU RNA-seq. Source: `results/benchmark_sim/outputs/ebmf_cox_external/yfb_vs_ebmf_k7_matched_ci.csv`.
  - Independent K = 40 with a LASSO Cox stage: mean C = 0.601. ΔC = +0.026 (0.0002, 0.0498). DECISIONS.md calls this the fairest comparison, and it was the 8/27 deck's headline.
  - The 2026-07-16 entry says its K = 20 result is "superseded" in its specifics.

  A committee member who asks how K was chosen for the baseline will find this.
- *Fix:* Use the LASSO, K = 40 baseline as the primary comparator, or report all three side by side: K = 20 unpenalized, K = 40 LASSO, and K = 7 matched. State plainly that the advantage shrinks to +0.011 and is not significant when the baseline is given the same K. Remove "higher in every cohort" and "consistently across cohorts and platforms".
- *Source:* DECISIONS.md lines ~1070–1150 and ~1696–1705; `ebmf_cox_external_results.csv` (K = 20); `yfb_vs_ebmf_k7_matched_ci.csv`; `ebmf_cox_regularized_results.csv`.

**M2. The gain is credited to "joint fitting", but with α_F = 0 the factorization is unsupervised.**
- *Location:* §3.1 heading, "Joint fitting improves external discrimination over a two-step pipeline in five held-out cohorts". Discussion line 288, "Fitting the factorization and the Cox model together gave better external discrimination than factorizing first". Contribution 1 also says "A joint model …".
- *Problem:* In `fit_cox_on_yf.R`, η = Ẑ_F β depends only on F, and the F update is pure genomics (line 662: "survival terms receive zero weight"). L enters only the reconstruction. So q(L) and q(F) are updated exactly as in unsupervised EBMF, and survival affects only q(β). The difference from the two-step pipeline lies in K, the point-exponential priors, the empirical-Bayes Normal prior on β, and α = 0.5. It does not lie in joint discovery. This is the α_F framing constraint in `Prelim_Proposal_Plan_10_1_26.md`. The Limitations section respects it, but the heading and the Implications bullet do not.
- *Fix:* Retitle §3.1 to describe the finding itself, for example "The supervised model's external concordance exceeds that of an independently tuned two-step pipeline". Say in §2.1 or §3.1 what actually differs between the two arms. Reword line 288 the same way.
- *Source:* `code/fit_cox_on_yf.R` lines 19–20, 170, 646–662; plan constraint 1.

**M3. Contribution 4 and §3.6 claim parsimony that the real data do not show.**
- *Location:* line 69, "The framework … chooses the number of factors without cross-validating K, and it substantially improves external discrimination". Also the §3.6 heading and the Table 2 caption.
- *Problem:*
  - In `real_sweep_summary.csv` the pruned real fits have K_final = K_init = 3, 7 and 10. The logs show `n_pruned = 0`, and the 10/2 chapter (line 36) says "Pruning removed no factors at these starting ranks".
  - The K = 7 and K = 10 fits have `converged = FALSE`, with 500 iterations reached.
  - External C depends on K_init: 0.557, 0.685 and 0.650.
  - Table 2 reports only K_init = 7, the best of the three, with no paired test against the original fit. The intervals overlap: (0.579, 0.779) against (0.319, 0.599).
  - The rank-independence result in §3.5 comes from simulation only: 2 seeds, one noise level, 50 + 50 features.
- *Fix:*
  - Limit Contribution 4 to what is shown: in simulation, pruning makes K_final independent of K_init. Say that on real data it has so far removed no factors and that the reason is under study (the nullcheck table in the 10/2 chapter, lines 271–287).
  - Report all completed K_init in §3.6. State that K_init = 7 was pre-specified to match the single-modality model, if that is true.
  - Replace "substantially improves" with the observed numbers plus the caveats: n = 50, 30 deaths, fits not converged.
  - Add a sentence to §3.6 that pruning removed no factors.
- *Source:* `results/multimodal_real/outputs/real_sweep_summary.csv`, `run_pruned_K7.log`, `docs/progress_book/chapters/2026-10-02.qmd` lines 36 and 271–287.

**M4. Single-modality concordance differs on the same donors between Table 1 and Table 2.**
- *Location:* Table 1, PACA-AU RNA-seq at 0.657, against Table 2, "Single-modality YFB fitter, expression only" at 0.491.
- *Problem:* DECISIONS.md 2026-10-01 says ICGC survival "agrees exactly with `PACA_AU_seq.survival_data.rds` for every RNA sample". The n = 50 ICGC primary set therefore overlaps almost entirely with the n = 52 PACA-AU RNA-seq cohort. Readers will see the same model family score 0.657 and then 0.491 on nearly the same patients, with no explanation. The differences are real:
  - training on TCGA alone (n = 144) rather than TCGA + CPTAC (n = 273);
  - gene screening on TCGA only;
  - log2 scale rather than per-platform z-standardization;
  - the expression-only fit hit its 100-iteration cap.

  DECISIONS.md 2026-09-04 already shows that TCGA-only training drops mean external C to 0.514.
- *Fix:* Add a sentence to §3.6 or the Table 2 caption that names these differences and cites the TCGA-only sub-analysis (0.514).
- *Source:* DECISIONS.md 2026-10-01 (matched cohorts) and 2026-09-04 (training-set sub-analysis); `run_pruned_expr_only_K7.log` (`converged FALSE`, 100 iterations).

**M5. "Regularized feature precisions" contradicts the Methods and the code.**
- *Location:* Discussion line 290, "a baseline intercept, regularized feature precisions, and evidence-based removal of factors".
- *Problem:* The framework keeps unregularized per-feature maximum-likelihood τ_mj (`tau_model = "feature"`). DECISIONS.md 2026-10-01: "Per-feature τ is kept. With the intercept in place, the degeneracy did not recur at moderate signal-to-noise". The empirical-Bayes Gamma alternative failed, with τ reaching about 10¹⁷. §2.2 item 3 says τ is "estimated per feature", which is correct. §3.5 lists unregularized τ as a cause but never says that the intercept, not regularization, resolved it.
- *Fix:* Delete "regularized feature precisions" from line 290. In §3.5, add that the intercept removed the degeneracy and that τ remains unregularized.
- *Source:* DECISIONS.md 2026-10-01, parsimony framework entry, `tau_model` bullet.

**M6. "Avoided … a cross-validated search over K" contradicts §2.1.**
- *Location:* line 289, "A single over-specified fit, followed by shrinkage, avoided both merged programs and a cross-validated search over K". Contribution 2 frames K selection the same way.
- *Problem:* §2.1 says K_init "is chosen by a consensus of the ELBO, BIC, bi-cross-validated genomic log-likelihood and cross-validated concordance, followed by a one-standard-error rule". That is a cross-validated sweep over K. DECISIONS.md 2026-08-20 also notes that K = 7 was originally chosen by cross-validating the survival outcome.
- *Fix:* Say what the procedure avoids: per-K model comparison for the *final* rank, since K_eff comes from one fit. Do not claim that it avoids cross-validation.

**M7. Broken cross-reference.**
- *Location:* line 102, "Removing these weights is part of the remaining work (Section 4.4)".
- *Problem:* There is no Section 4.4. Future work is §4.3.
- *Fix:* Change to Section 4.3, or better, use Quarto `@sec-` labels throughout so references cannot drift. Line 71 ("Section 4.3") is correct.

**M8. Contribution 3 promises a DeSurv comparison that the chapter does not contain.**
- *Location:* line 68, "with comparisons to a two-step factorization-then-Cox baseline and to DeSurv". Line 189: "Concordance in this range is comparable to that reported for DeSurv in its validation cohorts [@Young2026]".
- *Problem:* `results/benchmark_sim/outputs/desurv_comparison/` compares variants of the authors' own models under DeSurv-aligned preprocessing. It contains no fitted DeSurv model. The line 189 sentence gives no DeSurv numbers, and the cohorts and protocols differ.
- *Fix:* Either run DeSurv on the same training and validation split, or reword Contribution 3 to "using DeSurv's preprocessing and gene-selection rule" and move the head-to-head comparison to future work. If the line 189 sentence stays, give DeSurv's reported C values and say that they are not directly comparable.

---

## Should-fix

**S1. KM separation is overstated.** Line 202 says "separates patients by survival in each held-out cohort". Log-rank p is 0.075 for PACA-AU array and 0.17 for Moffitt, so only 3 of 5 cohorts reach p < 0.05 (`km_risk_tertile_stats.csv`). The text omits the PACA-AU array value. Suggested wording: "Risk tertiles were ordered as expected in all five cohorts. The log-rank test was significant in Puleo, Dijk and PACA-AU RNA-seq (p = 0.021)."

**S2. False survival-active factors are not reported.** §2.5 defines "the number of false survival-active factors" as a metric, but §3.5 does not report it. In `pruning_comparison_summary_quick.csv`, pruned fits average 0.5 to 1.0 false survival-active factors per fit, which is not better than the original fit in the informative scenarios. The 10/2 chapter also notes survival-only duplicates at large K_init and occasional over-pruning at K_init = 3. Report these alongside the recovery numbers.

**S3. Posterior uncertainty is not calibrated.** Contribution 1 says "Posterior uncertainty in the survival coefficients is estimated within the model". In the single-modality model this has two caveats:
- α = 0.5 multiplies the Cox precision A_k, which inflates the posterior variance of β by a factor of 1/α (`update_beta.R` lines 144–152).
- Ẑ_F is treated as fixed, with no uncertainty in F.

Add one clause saying the intervals are tempered and conditional on Ẑ_F. §2.1 states the second point but not the first.

**S4. Cohort-specific β prior.** Line 109 says "shrunk toward a common value". `update_beta_cohort.R` makes one vectorized `ebnm()` call per factor with `prior_family = "normal"`. As far as I can tell, ebnm's Normal family fixes the mode at 0 by default, so the C coefficients share a prior *variance* and shrink toward zero, not toward a common mean. Please check this and reword to "share one estimated prior variance (partial pooling of scale)". The code comment at lines 38–39 makes the same claim and may also need correcting.

**S5. The cohort-attribution result is vaguer than the source.** Line 234 says "attributed the effect to the correct cohort in 15–17 of 20 settings". The source reports 15/20 at K_init = 6 and 17/20 at K_init = 12, counted as seed-by-cohort cases. It also records one divergent seed whose estimates were 20 to 200 times the scale of the others, which is why medians were reported. State both counts and the divergent seed. The next sentence, "With larger cohorts, effects should be shared across cohorts", is an unsupported assertion. Replace it with the actual caution: training cohorts of 144 and 129, and no interval on the per-cohort β. DECISIONS.md 2026-09-04 also records that `beta_cohort_id` does not change external C (+0.022, CI includes 0), which is worth one clause.

**S6. The survival-active definition differs between the models and is not stated for the multimodal one.**
- Single modality: |E β_k| > 0.001 (`select_K.R`, `beta_thresh`). This cutoff depends on the units of Ẑ_F, so give the value and the reason for it.
- Multimodal: |E β_k| / SD(β_k) ≥ 1.96 (DECISIONS.md 2026-09-18). §3.5 and §3.6 use this definition, but §2.2 never gives it.

**S7. The "this explains" claim in Limitations is too strong.** Line 296 says the 10,000-to-1 term imbalance "explains why low-variance prognostic programs are not recovered". The simulation in which the program is lost has n = 120 and p = 50 + 50 (`run_multimodal_pruning_comparison.R` line 102), a ratio of 100 to 1. The imbalance is therefore a plausible contributor, not a demonstrated explanation. Write "is consistent with". The source for "< 3 nats" (factor 4: −2.9) and "10⁴–10⁵" is the 10/2 chapter's nullcheck table. Consider putting that table in the supplement.

**S8. The Methods do not match the reported runs.**
- §2.3 describes K_init = 3, …, 10, 15, two noise levels and three seeds.
- §3.5 reports K_init ∈ {3, 5, 7, 10, 15}, the moderate noise level only, and 2 seeds.

The [Preliminary] note covers seeds but not the K grid or noise level. §2.3 also omits n, p, censoring and the 4-SD baseline shift. The single-modality simulation never states the true K (6 in the hybrid scenario: 2 shared plus 2 + 2 specific) or the sample sizes, yet §3.3 relies on "fewer factors than truth".

**S9. Pooling methods are described inconsistently.** §2.5 says "a pooled, cohort-weighted bootstrap". §3.1 says "a fixed-effect, inverse-variance pooling". The CSV label is "POOLED (fixed-effect, weighted by 1/SE^2)". Use one description in both places and give the number of bootstrap replicates (2000 for the paired test in DECISIONS.md, 1000 for the multimodal fits).

**S10. Notation clash.** β denotes both the Cox coefficients and methylation β-values (line 175). Use "methylation beta-values (M)" or similar.

**S11. Convergence statement.** Line 281 says "The fits had stable risk scores but had not met the formal convergence criterion". The K = 3 pruned fit did converge, while K = 7 and K = 10 hit the 500-sweep cap. The evidence for "stable risk scores" is not in the outputs I checked. Either cite it or state only the iteration cap. Also note that the nullcheck ran on non-converged fits, because `prune` runs after the loop whether or not it converged.

**S12. Only one of the two single-modality programs appears in the multimodal fit.** §3.6 reports that the survival-active factor matches Program 7 (|r| = 0.676 on 1,781 genes ✓). The K = 7 fit had only one survival-active factor, so the protective Program 3 was not recovered as survival-active. Say so.

**S13. The DeSurv contrast is incomplete.** Line 60 lists only DeSurv's weaknesses. The literature review (§1.3.6) says DeSurv's survival term *does* enter the factorization objective and "directly addresses the discovery limitation". Given α_F = 0, the chapter should say that DeSurv supervises the loadings and that the current SSBMF does not. The distinctive points of SSBMF are then posterior uncertainty in β, shrinkage-based K_eff, and the multimodal extension. Stating this plainly will hold up better in an oral exam than leaving it implicit.

**S14. Minor accuracy points.**
- Line 175, "0.07% missing values within each cohort": the source gives 0.07% for TCGA only.
- Line 121, "Each modality–program pair has its own loading prior, chosen from …": the family is chosen per modality (`control$prior_F`), and the prior parameters are fit per modality–factor pair.
- Line 263 lists "signed point-Laplace loadings" as part of the framework, but the §2.2 parsimony list (items 1–3) omits it.

---

## Structure and completeness

What a prelim committee will expect that is missing or thin:

1. **Specific aims and success criteria for the remaining work.** §4.3 is a list. Turn it into two or three aims, each with a testable criterion. For example:
   - The α_F-free supervision test: recovery of the planted low-variance program, |r| ≥ 0.7 and survival-active in at least X% of seeds.
   - The full multimodal simulation.
   - The real-data multimodal fit, with a pre-specified K_init and a convergence requirement.
2. **Timeline** to the defense. Even a short table would do.
3. **Simulation design table.** A placeholder is already flagged at line 147. Give the scenarios, n, p, K_true, β, noise and seeds for both simulation studies.
4. **Abstract or chapter summary paragraph** at the top.
5. **A methods paragraph on the comparators:** what the two-step arm is (flashier, its K, the Cox stage) and why it is a fair comparison.
6. **Prediction and preprocessing of new cohorts:** per-platform standardization of the validation data, and centering on training means for the multimodal model. Contribution 1's point about predicting from measurements alone depends on this.
7. **Initialization and convergence rules** for both fitters: SVD initialization; relative ELBO change below 1e-5 with max_iter 100 for the single-modality model; the scale-invariant 1e-4 rule with 500 sweeps for the multimodal model.
8. **Results ordering.** The order is reasonable (main result, programs, K, cohorts, then multimodal). Consider adding a short "Summary of evidence" table that marks each claim as established or preliminary.
9. **Prose.** Most of Methods, Results and Discussion is bullet lists. Double spacing makes this look fragmented, with half-empty pages 6–7 and 16–18. Committee chapters are normally written as paragraphs. Keep bullets for parameter lists only.
10. **Discussion.** It does state that the project is in progress (line 302). Implications should be rewritten after M1, M2 and M6.
11. **Related multimodal work is not cited in this chapter.** Contribution 4 needs MOFA (`Argelaguet2018`) and iCluster (`Shen2009`), both already in the bib, so readers can see what is new. A brief cross-reference to the literature review is enough.

---

## Citations

- **Citations that hold up, judged from the bib titles:** Collisson2011, Moffitt2015, Bailey2016 and Raphael2017 for subtypes; Rashid2020 (PurIST) for single-sample classifiers; Brunet2004 and Kotliar2019 for factorization into programs; Prentice1982 and Carroll2006 for attenuation; Bair2004; Tibshirani1997 and Simon2011; Nygard2008; Wang2021; Willwerscheid2025 (ebnm); Cox1972 and Cox1975; Breslow1974 for ties; Blei2017; Harrell1982; BhattacharyaDunson2011 and GriffithsGhahramani2011 for factor-level shrinkage.
- **SEER2026** supports the five-year survival figure (13.7%). It does not support "outcomes vary widely among patients diagnosed at the same stage". Move the subtype citations there, or split the sentence.
- **Young2026** is a `@misc` entry with a `journal` field. It renders as "In *PNAS*." with no volume, pages or DOI. Make it `@article` with full details, checked against the literature review's verified bib.
- **Missing citations:**
  - flashier, for the greedy and backfit steps and the "flashier-style nullcheck";
  - the cohorts: CPTAC-PDAC, Puleo 2018, Dijk; PACA-AU via Bailey2016 is acceptable;
  - k-nearest-neighbor imputation (Troyanskaya 2001, or the `impute` package);
  - bi-cross-validation (Owen & Perry 2009);
  - BIC (Schwarz 1978);
  - the one-standard-error rule (Breiman et al. 1984, or Hastie et al.);
  - the percentile bootstrap (Efron & Tibshirani);
  - Weibull survival simulation (Bender et al. 2005);
  - MOFA and iCluster for Contribution 4.
- The literature-review section references (§1.2, §1.3, §1.3.3) resolve correctly. The review's third `##` section is "Supervised Bayesian Factorization for Prognostic Genomics", and its third subsection is "Survival Modeling and the Two-Step Pipeline", which contains the attenuation argument. Good.

---

## Rendering and format

- Renders cleanly: 22 pages, no unresolved citations, no missing glyphs, and the references section is present.
- **Placeholders visible in the PDF:** "[CROSS-REFERENCE: simulation design table]" (p. 7), "[CROSS-REFERENCE: appendix or supplement]" (p. 6), "[Preliminary: …]" (p. 15) and "[PENDING: remaining starting ranks.]" (p. 17).
- **Table 1 header:** "Two-step EBMF → Cox" wraps onto two lines, misaligned with the other headers. Shorten it to "Two-step" and define it in the caption.
- **Figure 2 (heatmap):** it fills almost all of p. 12 with small gene labels. Its column labels are "P1_Program 3 (Protective)", "P2_Genomics-only (A)", "P3_Genomics-only (B)" and "P4_Program 7 (Adverse)". The column headed "P3" is *not* Program 3, which will confuse readers. The embedded title, "K=7 fit's 4 kept factors … (3 fully-pruned factors omitted)", is internal phrasing. Relabel the columns "Program 3 (protective)", "Stromal", "Immune" and "Program 7 (adverse)", remove the title, and reduce the height. Note also that the "immune" column includes stromal genes (COL14A1, SLIT3, FBLN5, TNXB, SPARCL1, CHRDL1), so "immune/quiescent stroma" may be a more accurate label than "immune" in the text.
- **Figure 4 (ΔC by K_init):** the title "Joint YFB vs. unsupervised EBMF", the legend levels `all_shared`, `hybrid` and `nothing_shared`, and the axis label "K_init" are code-style. Relabel them and drop the embedded title.
- **Figure 5 (per-cohort β):** the title "(beta_cohort_id fit)", the caption text "ARD-kept" and the legend "TCGA_PAAD" are code names. The sign convention is also unstated. In the figure, Program 7 is positive and Program 3 negative, which is the risk = −η orientation. In DECISIONS.md the raw signs are the opposite (β₇ < 0). Add "positive = higher hazard" to the caption.
- **Figure 6:** the legend "Pruned (partial log-lik)" and "Original (no intercept, no pruning)" use internal wording, and the fonts are small.
- **Title block:** no author, date or affiliation. This is acceptable if the chapter is assembled into the prelim template, but the YAML comment marks it [TO CONFIRM].

---

## Writing rules (CLAUDE.md)

- **Internal labels in the text and tables:** "Joint YFB" (Table 1), "Single-modality YFB fitter" (Table 2), "New framework" and "Original multimodal fit" (Table 2 and §3.6), "the frozen model" (§3.1). "YFB" is never defined as a name. Use descriptive labels, for example "Supervised model (shared programs, η = (YF)β)" and "Multimodal, with intercept, signed priors and pruning".
- "EBMF" is used before it is defined. Write "empirical-Bayes matrix factorization (EBMF)" once in §1.3.
- There is no mannered prose. The writing is direct throughout.

---

## Verified claims (match source)

- **Table 1:** all ten cohort C values and CIs, the n values, and the means 0.627 and 0.581 (`external_cindex_ci.csv`). Per-cohort significance only in Puleo, and pooled 0.042 (0.013, 0.071) (`external_paired_diff_ci.csv`). Note M1: these are correct for the K = 20 baseline.
- **Orientation statement in §2.5:** consistent with DECISIONS.md 2026-10-02 (one global risk = −η, frozen orientations reproduce both CSVs).
- **KM:** Puleo p < 0.0001, Dijk 0.0005 (0.00046), Moffitt 0.17 (`km_risk_tertile_stats.csv`). 616 held-out patients (90 + 123 + 63 + 52 + 288). Program 3 pooled p = 0.004 (10/2 chapter).
- **Programs:** four retained from K = 7, two survival-active, three shrunk to zero. Program 7 genes (MET, ITGA3, glycolytic), Program 3 classical, stromal genes (collagens, FAP, SPARC, CTHRC1), and IKZF1, IL10RA and HCLS1 all appear in the heatmap. The classical and basal-like enrichment matches DECISIONS.md 2026-07-15.
- **Merges:** 7 of 27 at K_init = 2–4 and 0 of 114 at K_init = 6, 12 and 20 (DECISIONS.md 2026-09-04).
- **ΔC by K_init (figure):**
  - Hybrid: positive at every K_init ≥ 3, with a maximum of about 0.157 at K_init = 4 and about 0.095 from K_init = 6.
  - All shared: positive only from K_init = 6.
  - Nothing shared: within about 0.013.
  - 15 seeds.
  - The two-step arm in these simulations also used an unpenalized `coxph`, per DECISIONS.md 2026-08-20. Consider noting this.
- **Per-cohort β:** Program 7 is similar in both training cohorts, and Program 3 is concentrated in CPTAC (0.0304 against 0.0004).
- **Multimodal simulation (quick summary):**
  - Original K_final = K_init.
  - Pruned K_final between 2 and 3.5.
  - All three true programs recovered from K_init ≥ 5 in the informative scenarios.
  - Direction recall 0.75–1.
  - All pruned fits converged.
  - Low-variance program recall 0 in both arms.
  - Scenario β vectors and score scale 0.15 match `simulate_multimodal_yfb.R`.
  - The moderate noise level is labeled "factor PVE ≈ 10%" in the figure.
- **Table 2:** all five values and CIs (`real_sweep_summary.csv`). K_init 3/7/10 give 0.557, 0.685 and 0.650. |r| = 0.68 against Program 7 (0.676, 1,781 genes).
- **Data:** TCGA 144 patients with 75 deaths. ICGC 50 with 30 deaths, and 67 with 40. 305,352 CpGs. 3,000 genes and 10,000 CpGs screened on TCGA only (DECISIONS.md 2026-10-01).
- **Multimodal methods:**
  - E[Z²] includes Σ y² Var(f) in the β update: `projection$VZ` is passed to `multimodal_yfb_update_beta_k`.
  - No α weights.
  - Unit-SD projection rescaling (`multimodal_yfb_canonicalize_factors`).
  - Intercept re-estimated each sweep.
  - Pruning ELBO = reconstruction + Breslow log PL at E[η] − KL, then refit and repeat.
- **Discussion numbers:** "< 3 nats" (factor 4: −2.9), "10⁴–10⁵", "13,000 features", "> 10,000 to one" (about 1.9 million / 144) (10/2 chapter nullcheck table).
- **Single-modality methods:**
  - Point-exponential L and F priors with empirical-Bayes fits.
  - Normal empirical-Bayes β prior (`prior_beta = "normal"`).
  - Ẑ_F = Y E[F] D⁻¹ with unit-norm columns.
  - Ẑ_F treated as fixed, with no second-moment correction.
  - α = 0.5 and α_F = 0.
  - Breslow ties.
  - Classification: |β| threshold and PVE ≥ 1%.
- **α_F constraint:** satisfied in §2.1 (line 101) and Limitations (line 294). It is violated in spirit by the §3.1 heading and the line 288 Implications bullet (M2).
