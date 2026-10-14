# Independent review of the initial writeup

Reviewed 2026-10-14 on `writeup/draft`, clean at review start; draft `0a7cfc9`, monochrome revision `a9f4b47`. This review follows `review_prompt.md` and the newer author feedback in `implementation_plan.md`. Nothing below is marked fixed. Only this review file is changed; no report, code, result, inventory, or figure is edited.

## Summary

The draft consolidates substantial evidence and handles several difficult distinctions correctly: frozen risk orientation, the principal regularized comparator, the concatenated versus cohort-pooled estimands, the proteomics training input, overlapping validation cohorts, lack of direct survival feedback into the single-modality loading updates, and unconverged multimodal fits. It is a useful basis for revision, but is not ready for approval in its current form.

The main scientific corrections concern the Cox working curvature, the distinction between a random-loading joint model and the evaluated plug-in algorithm, normalization of the saved predictor, and what the implemented monitoring/removal objectives actually calculate. Both appendices give many correct final updates but skip the intermediate steps the author now requires. The numerical audit also identifies an overstatement of recovery and a rounding error. The manuscript-style changes requested by Andrew are substantive restructuring, not just proofreading.

Priority convention: **P1** blocks a defensible central conclusion; **P2** is a substantive scientific or methodological correction required before approval; **P3** is a smaller precision or presentation correction. No new fitting is needed to address the report corrections. Implementation defects identified here should be tracked for subsequent method work, without silently repairing historical fits.

## Part A — scientific and numerical findings

### A1 [P2] Define the joint model separately from the plug-in algorithm

**Location:** Methods §§1.1–1.3, `writeup.qmd:146–187`; Appendix A §§A.1, A.6, A.10 (`single_modality.qmd:10–16,92–105,149–163`).

Appendix A first defines the supposed joint survival model using posterior-mean loadings and fixed projections, then derives a survival term for random loadings. Under the first definition that loading term does not exist. A normalized random-loading predictor is also not exactly quadratic in one loading because its norm changes with that loading.

**Correction:** First define the intended partial-likelihood joint target using random (F\) and (\beta\):

\[
p(Y\mid L,F,\tau)\,\mathrm{PL}(t,\delta\mid Y,F,\beta)
\prod_{ik}g_{L,k}(L_{ik})\prod_{jk}g_{F,k}(F_{jk})\prod_kg_{\beta,k}(\beta_k).
\]

Explicitly specify whether the normalizer is a function (D(F)\) or an externally fixed scale. A Cox partial-likelihood target is not itself a fully specified survival-time likelihood. Then describe the evaluated algorithm separately: at sweep (r\), freeze (D^{(r)}\) and (\tilde Z^{(r)}=Y\bar F^{(r-1)}[D^{(r)}]^{-1}\), update coefficients on these covariates, and update activities/loadings from reconstruction only. A Gaussian survival-informed loading kernel requires the additional fixed-normalizer approximation. Preserve the distinction throughout both algorithms and Appendix A's comparison table.

### A2 [P2] The working curvature is not the exact Cox Hessian diagonal

**Location:** `writeup.qmd:181–185`; Appendix A §A.2; Appendix B §B.2 and §B.11.

The text says that only off-diagonal curvature is dropped. In fact the implemented (w_i=\theta_i\hat H_0(t_i)\) also omits a term from the exact diagonal. For distinct event times (a\), event counts (d_a\), risk sets (R_a\), define

\[
P_{ai}=\frac{1\{i\in R_a\}e^{\hat\eta_i}}{\sum_{j\in R_a}e^{\hat\eta_j}},\qquad
u_i=\delta_i-\sum_a d_aP_{ai}.
\]

The exact negative Hessian is

\[
H=\sum_a d_a\{\operatorname{diag}(P_a)-P_aP_a^\top\};\quad
H_{ii}=\sum_a d_aP_{ai}(1-P_{ai}).
\]

The implementation instead uses (a_i=\sum_a d_aP_{ai}\) as working diagonal curvature. `code/fit_cox_on_yf.R:103–119` further floors it: (w_i=\max(a_i,10^{-6})\), retaining the score (u_i=\delta_i-a_i\). Consequently writing (u_i=\delta_i-w_i\) is inaccurate at the floor. The multimodal implementation uses the division-free linear term (h=u+w\hat\eta\), without that same weight floor.

**Replacement:** “We retain the Breslow score and replace its coupled curvature matrix by a Poisson-style diagonal working approximation.” Show the exact Hessian, the replacement, the floor where applicable, and the square completion in B14 below. Do not claim this is the exact second-order Taylor expansion with only off-diagonals discarded.

### A3 [P2] The single-modality monitoring statistic does not contain the exact stated KL contributions

**Location:** Appendix A §§A.3, A.8; Methods stopping/objective descriptions; rank-sweep interpretation.

`code/compute_elbo.R:78–86` computes the expected pseudo-Gaussian log density with `0.5*sum(log(A/(2*pi)))`. In the activity update, `A` is scalar and the pseudo-observation/posterior vectors have length (n\) (`code/update_L_surv_YFB.R:78–96`; fitter `:707–708`). The quadratic term recycles the scalar precision correctly, but the log-normalization is included once rather than (n\) times. For (d\) pseudo-observations sharing precision (A\), the correct term is

\[
E\log p(x\mid\theta)=\frac d2\log\frac A{2\pi}
-\frac A2\sum_h(x_h^2-2x_h\bar\theta_h+\overline{\theta_h^2}).
\]

For the negative-KL contribution \(\kappa=\log p(x)-E\log p(x\mid\theta)\), this implies

\[
\kappa_{\rm code}=\kappa_{\rm correct}
+\frac{d-1}{2}\log\frac A{2\pi}.
\]

This discrepancy changes with precision and is not a harmless fixed offset. It enters the monitored stopping statistic; effects on stopping/rank preferences have not been assessed. It does not directly alter the coordinate posterior updates.

There is a second limitation: the survival monitor uses the partial log-likelihood at the sweep-start predictor while coefficient means subsequently change. Its recorded expression is bookkeeping, not the full expected quadratic at the updated mean. That expectation would also include
\(u^\top(E\eta-\hat\eta)-\tfrac12\sum_iw_i(E\eta_i-\hat\eta_i)^2\).

**Fix:** Retain the mathematically correct KL identity, but label the reported historical statistic as the implemented monitoring score. Do not identify its assembled terms as an exact current-state ELBO. Keep the numerical sweep as a description of the saved statistic, with this limitation in the technical appendix/limitations; do not rerun fits in the revision. The caption of Table 4 correctly calls it a monitor, but Results §3.3 still calls it “the evidence lower bound”; use consistent wording. Track the normalization defect for method work separately.

### A4 [P2] Use the stored prediction normalizer, and remove the feature-count invariance claim

**Location:** Table 1, Equation 2, Algorithm 1 output/prediction; `writeup.qmd:152–180,213–215`; Appendix A §A.9.

The displayed final predictor defines (D\) as the norm of the final loadings. The saved fit instead stores norms computed at the **start of the last sweep**, before its final loading updates. Appendix A.9 already acknowledges this. Evidence: `code/fit_cox_on_yf.R:643–646,734–737,926`; `code/predict_cox_on_yf.R:104–120`. Inspection of the saved fit confirms small nonzero discrepancies for all four retained columns; this is an actual difference, not merely a theoretical possibility.

**Fix:** Define a frozen stored normalizer (D_*\) for the saved predictor and use it consistently in Table 1, the equation, algorithm output and prediction. Distinguish this from the within-sweep definition (D^{(r)}\). State that numerical regularization is suppressed in displayed expressions, or include the (10^{-10}\) squared-norm ridge. Zero loading columns are not unit-length directions.

Also replace “keeps [the projection] on a scale that does not grow with (p\)” with “removes the arbitrary Euclidean magnitude of each loading column.” Indeed,
\(\operatorname{Var}(y^\top f/\|f\|)=f^\top\Sigma f/\|f\|^2\), which can grow with feature count. Identical standardized features and (f=\mathbf1_p\) give \(\sqrt p\,y_1\). The normalization/shrinkage development history is not a general mathematical guarantee.

### A5 [P2] Multimodal removal uses cached penalties from earlier priors

**Location:** Methods Equation 6 (`writeup.qmd:274–276`); Appendix B §§B.2, B.6, B.11.

`code/multimodal_yfb_updates.R:359–367,491–502` updates posteriors under the existing priors, estimates new priors, but returns KL values calculated using the old-prior marginal likelihood. This applies to both activities and loadings. The fitter and deletion checks sum those cached penalties (`code/fit_multimodal_yfb.R:528–531,577–581,734,742–754`), rather than recalculating KL against all current priors. B.6 acknowledges this for loadings only.

**Fix:** Describe “cached coordinate KL penalties calculated against the priors used for the latest coordinate updates,” and disclose the subsequent prior refit without penalty recomputation. Deletion improves the **implemented approximate removal score** before refitting. Replace the unqualified “conservative”/“evidence does not fall” language: neither exact marginal evidence nor an ELBO at the current posterior/prior is guaranteed to improve. No analytical result should be regenerated to conceal this distinction.

### A6 [P2] Qualify the variance correction and centering invariance

**Location:** `writeup.qmd:268,319`; Appendix B §§B.1–B.3, B.11 (`multimodal.qmd:21,25–37,147`).

Jensen's inequality establishes
\(E\ell_C(\eta)\le\ell_C(E\eta)\), but does not establish that the gap grows monotonically with the marginal variances. For example, adding the same random shift (U\) to every predictor increases all marginal variances while leaving the Cox likelihood and its Jensen gap unchanged. The full second-order expectation correction is

\[
-\tfrac12\operatorname{tr}\{H\operatorname{Cov}_q(\eta)\},
\]

not generally \(-\tfrac12\sum_iw_i\operatorname{Var}(\eta_i)\). Patients share uncertain loadings and coefficients, producing cross-patient covariance, and (w\) is itself the approximate curvature in A2.

**Replacement:** “Concavity makes the plug-in survival term an upper approximation to its posterior expectation. The alternative subtracts a variance penalty using the diagonal working curvature; the full second-order correction would also involve the coupled Cox Hessian and cross-patient covariances.”

Raw versus intercept-centered **frozen** predictions differ by a patient-independent constant, so exact Cox partial likelihood and concordance agree. That is correct. But the expected diagonal working surrogate need not preserve this invariance when the common shift is random under (q(F,\beta)\). Thus the reconstruction-only intercept update is an implemented approximation, not justified simply by saying survival never depends on the intercept. State that fitting uses centered projection moments and the intercept update maximizes reconstruction likelihood.

### A7 [P2] Weaken the “every shared program” simulation claim

**Location:** Results §3.1, `writeup.qmd:414`; compare Table 15 and the saved multicohort simulation output.

At starting rank at least 6, shared recovery is not always 100%: all-shared rank 12 is 0.97778 and hybrid rank 20 is 0.96667. The current statement contradicts the draft's own Table 15.

**Replacement:** “At starting ranks 6, 12 and 20, mean shared-program recovery was approximately 97–100%, with nearly complete recovery of cohort-specific programs.” Keep the scenario/rank-specific table and avoid implying every replicate recovered every program. The outcome-oriented concordance caveat and the absence of direct loading supervision are correctly stated.

### A8 [P2] Both appendices need the derivation additions in B14

**Location:** Appendix A §§A.1–A.8; Appendix B §§B.1–B.9.

Many final coordinate formulas agree with the implementation, but the route from likelihood and priors to the normal-means inputs is incomplete. Examples include A.6's “up to factors of (d_k\),” the absence of the complete variational objective, “same argument” substitutions, and omitted posterior-moment and variance-estimation derivations. This falls short of the author's explicit requirement that both appendices be reviewable step by step. B14 supplies the missing algebra and distinguishes exact identities from approximations.

### A13 [P2] Do not turn the absence of supervision into an unevaluated recovery failure

**Location:** Discussion `writeup.qmd:668` ("Neither model recovers ...") and Continuing work `:704` ("cannot").

The single-modality experiment described in Evaluation uses shared/cohort-specific programs; it does not test the multimodal planted low-variance scenario. The multimodal failure at lower noise is supported, but a failure of **both** models on this experiment is not demonstrated. Nor does α_F=0 imply that an unsupervised factorization can never recover low-variance structure; it means survival cannot directly help the loading updates select it.

**Replacement:** “The multimodal framework did not recover the planted low-variance prognostic program at lower noise. The evaluated single-modality loading update receives no survival information, so this configuration does not test whether survival supervision improves recovery of such programs.” In Continuing work, describe improving recovery as a goal rather than claiming that low-variance recovery is mathematically impossible without supervision.

### A9 [P3] Correct two numerical/algorithm details

- **Concatenated C:** `writeup.qmd:432` gives 0.629 for the single-modality model. Frozen-score recomputation gives **0.6284820**, which rounds to **0.628**, not 0.629, at three decimals. The comparator gives 0.6025621 → 0.603 and the difference 0.0259199 → 0.026. Use four decimals (0.6285 versus 0.6026) where helpful to avoid rounding confusion.
- **Multimodal stopping denominator:** Appendix B.10 uses only the new fitted quantity. The code uses \(\max(1,\|X_{\rm old}\|_\infty,\|X_{\rm new}\|_\infty)\) (`code/fit_multimodal_yfb.R:64–69,514–540`). Apply this definition to each reconstruction block and the mean predictor.

### A10 [P3] State the scale and feedback boundaries precisely

**Location:** Appendix B §§B.6, B.9; main scale discussion.

- Survival feedback depends on both \(\bar\beta_k\) and \(E\beta_k^2\). A near-zero mean with substantial variance does not switch it off. Say “small mean and second moment”; an exact point mass at zero turns off both terms.
- Complete scale transformations with \(\sigma_{\beta,k}^2\mapsto c_k^2\sigma_{\beta,k}^2\), \(EZ_k\mapsto EZ_k/c_k\), and \(\operatorname{Var}Z_k\mapsto\operatorname{Var}Z_k/c_k^2\). Explain pushforward distributions/Jacobians rather than invariant continuous density values. KL is invariant when both distributions are transformed consistently.
- The code substitutes (c_k=1\) for effectively constant projections (`fit_multimodal_yfb.R:197`); these columns do not acquire SD one.
- Say survival does not enter the single-modality reconstruction **updates**. The survival-dependent stopping statistic can still influence the final iterate, so the result need not equal a separately converged unsupervised fit.

### A11 [P3] Complete Revisit coverage while preserving the balanced markers

All 9 Revisit IDs occur in the main source. All START/END pairs are properly nested and balanced across all three sources. All existing figure definitions/captions and tables fall within markers. However, method-dependent unmarked passages remain:

- `writeup.qmd:132–134`: projection/activity notation and interpretation — add R2 (and R5 where cohort/data definitions change).
- `writeup.qmd:324`: evaluated configuration and future method statement — move to Continuing work under R9/R1, or mark before moving.
- `single_modality.qmd:3`: evaluated-method reconciliation introduction — R1/R2 if retained after process-language removal.
- `single_modality.qmd:163`: exactness/conclusion and future supervision — R1/R9.

The generic normal-means/KL identities in A.3 are method-independent and need not be wrapped solely to make every line marked. Appendix B is currently one enclosing R3 block with local R2 blocks; retain those but add local dependencies such as R4 for rank/removal and R9 for limitations if they become independently movable. Update R9's inventory of affected passages to Discussion/Continuing work, with the gap table removed. Recheck marker coverage after all moves and caption replacements.

### A12 [P3] Avoid implying uniformly ordered survival curves

**Location:** `writeup.qmd:454`, “low-risk tertile has the longest survival in every cohort.”

The plotted curves cross, notably in Moffitt; no saved tertile medians were found in the reviewed KM summary. A significant log-rank test is not proof of uniform ordering or a longer median for every group.

**Replacement:** “Survival curves generally favored the low-risk tertile; log-rank tests indicated differences among tertiles in three of five cohorts.” Preserve the exact cohort-specific p-values and descriptive-cutoff caveat.

### Independently checked drafter statements and numerical scope

| Statement | Review outcome |
|---|---|
| Principal +0.026 is concatenated, not inverse-variance pooled | Confirmed from saved scores and producer scripts. Both models are evaluated on 616 observations/414 events; cross-cohort comparable pairs enter this estimand. This does not establish cohort-stratified superiority. |
| Saved per-cohort C-indices equal frozen-orientation values | Confirmed for the established model and K40, K7 and K20 comparators. This supports these reported values, not outcome-based orientation as a generally valid validation procedure. |
| Established model C=0.685 on the matched 50 donors has no saved source | Qualified confirmation: the chapter/plot contains that number, but a matching saved computation/score source was not located in the reviewed sources. Omit it as a verified result. Do not assert exhaustive repository-wide or historical nonexistence. |
| Risk-scale Program 3/7 coefficients agree with marginal labels | Confirmed after the fixed risk = −η transformation: −0.01150643 and +0.04040517, respectively. Do not repeat the old sign disagreement. Conditional and marginal effects remain distinct. |
| Stored coefficient variances are zero | Confirmed exactly in the saved fit. The historical-code explanation is consistent with corrective comments, but its original executable was not reconstructed. |
| Program 3 enrichment includes both DeSurv lists | Confirmed; the classical label rests on concordance and association as well, rather than specific enrichment alone. |
| K=7 has the highest primary-ICGC C in the joint sweep | Confirmed among the nine saved starting ranks. State the selected favorable comparison and broad sensitivity, not an established optimum. |

Other checked material includes the saved fit dimensions/iteration count and retained-program moments; comparator cohort C values and the distinction among the three summary estimands; all Table 15 means and the simulation recovery exceptions above; all Table 5 aggregates from the 486 raw sweep rows (162/162 selected-variant fits converged versus 27/162 base fits); multimodal rank, convergence, classification and program-table values; and biological HR ranges/significance counts, subtype correlations, overlaps and shares. The final revision should retain exact source selections in the separate inventory even as visible IDs are removed. No fitting or new enrichment/correlation analyses were run. Full bootstrap distributions and the historical fitting environment were not recreated.

## Part B — requested manuscript revision, item by item

These are the author's requested changes, not claims that they have already been applied. Current exhibit numbers below refer to the reviewed draft and will change after restructuring.

### B1. Title and subtitle

**Locations:** YAML `writeup.qmd:2–3`; PDF title page, HTML title and PDF metadata.

Use exactly **Survival-Supervised Bayesian Matrix Factorization for Prognostic Molecular Programs in Pancreatic Cancer**. Remove the subtitle in both formats and metadata. The aspirational title must not lead Methods/Discussion to claim direct loading supervision for the evaluated single-modality algorithm.

### B2. Authored appearance

**Locations:** YAML running header `:90`; status box `:100–106`; Summary `:108–126`; process descriptions in §§1.3, 1.8, 2.4–2.5, 3.3–3.5, 4.1–4.2, 5 and all appendix introductions/reconciliation tables/captions.

Delete the status box and running draft header. Remove drafting/checking/reviewer language and author-list discussion from the scientific body. Keep verification/review status in repository documentation. Replace “no saved score file ... was found during drafting” with omission of the unverified number, not a paragraph describing the unsuccessful search. Replace “not used in this report” with a scientific distinction if needed, or cut it. The manuscript should state what was done and what the evidence supports.

### B3. Remove repository references and internal labels

**Affected locations:** Table 1; both algorithm headings/steps; §§1.3–1.8, 2.1–2.5, 3.1–3.5, 4.1–4.2, 5; every current caption containing sources; both appendix introductions and reconciliation tables; Appendix C. Search the entire included source, not just the main text.

Remove visible S01–S62 identifiers **entirely**, without moving them into hidden comments. Retain the separate source inventory. Retain only the required Revisit comments. Translate implementation labels to scientific descriptions:

| Current wording | Proposed wording |
|---|---|
| D4 | Single-modality model trained on pooled expression and proteomics with platform-specific standardization |
| D4 P3 / P7 | Single-modality Program 3 / Program 7 |
| flash k5 | Factor 5 of the unsupervised expression reference (define reference once) |
| expr-only | Expression-only model |
| `pruned_expression_only` | Expression-only fit of the multimodal framework |
| `original` / September model | Base multimodal model without feature intercepts or explicit pruning |
| `survival_elbo = plugin` | Survival term evaluated at the posterior-mean predictor |
| `prior_F = point_laplace` | Point-Laplace loading prior |
| `ebnm_point_laplace` | Empirical-Bayes normal-means estimation with a point-Laplace prior |
| `flip = FALSE` | Higher predicted scores denote higher risk, with direction fixed before validation |

Software package names can remain, but function names, paths, options, selectors and output-directory names cannot. Preserve the scientific prior, update order and realized settings rather than deleting those facts together with their code syntax. In Tables 8–9 distinguish the two models' program numbers explicitly; they are not a common numbering scheme.

### B4. Remove dated/process descriptions

**Locations:** §§1.3 and 1.3.1 (historical orientation), §1.8; §2.1; §§3.3–3.5; §4.1 direction discussion and Program 6; Figure 2/3/4/5 captions and Figure 2 panel title; §5; both appendix introductions, A.9/A.10 and B.13.

Describe base versus extended models by their components. Replace “September noise” by the numerical noise multiplier. Describe the historical frozen prediction convention directly: “For the evaluated single-modality fit, higher η denotes longer survival; reported risk is −η.” Keep zero stored variances as a fit limitation, without a date-stamped debugging history. Replace “the August 21 report” with the actual biological evidence. Drop meeting/progress/decision references and the unsupported per-program objective decomposition paragraph in §3.5; retain the high-dimensional objective imbalance as a limitation, not a measured quantitative decomposition.

### B5. Table 5 and Figure 2 noise labels

**Locations:** §2.1 `:334`; Summary; §3.4 `:495–523`; Table 5; Figure 2 caption/alt text and raster panel title; B.8/B.13 where relevant.

Use **Lower noise (0.25× baseline feature-noise SD)** and **Higher noise (1× baseline feature-noise SD)**, with “multiplier on each feature's noise standard deviation” in the caption. Define the baseline precision/noise distribution in Evaluation. Do not call a relative multiplier an absolute SD or imply one common numeric SD across heteroscedastic features. Preserve the same seeds/scenarios/means/ranges; redraw the saved-output figure labels only.

### B6. Self-contained captions and undefined abbreviations

The existing captions do not define the following terms locally. Expand or remove each term in each retained exhibit; a definition in Methods alone is insufficient. Gene symbols need not all be expanded as acronyms, but selected gene lists must be identified as such. Define nonstandard mathematical quantities/thresholds as well.

| Current exhibit | Undefined/insufficiently defined locally | Additional caption/note content needed |
|---|---|---|
| Table 1, model summary | TCGA-PAAD, CPTAC, RNA-seq, DNA, SD, Cox, KL where used; package `ebnm` needs a plain description | Dimensions, bar notation, projection scale, stored normalizer, evaluated priors; remove D4/source IDs |
| Table 2, cohorts | TCGA-PAAD, CPTAC, ICGC, PACA-AU, PDAC, RNA-seq, CpG, IPMN, KNN if named | Explicitly pooled training deaths; define assay/platform labels and overlap; clarify log transformation and screening scope |
| Table 3, external performance | C, CI, LASSO, PACA-AU, RNA-seq | Harrell concordance (higher is better), percentile bootstrap/B, concatenated versus cohort-pooled difference, paired patient sampling, frozen directions; don't round the positive lower bound to an apparently exact zero without a note |
| Table 4, rank sweep | ELBO, BIC, CV, PL, C | Define starting rank, classification thresholds, monitoring statistic rather than exact ELBO, model/dataset, CV folds, metric directions; clarify which criteria are maximized/minimized |
| Table 5, multimodal simulation | C; recovery/sign/classification measures are only partly specified | Training/validation size, true rank, noise multiplier, number of datasets versus fits, baseline parentheses, coefficient-active rule, no uncertainty interval supplied |
| Table 6, matched validation | TCGA, ICGC, PDAC, C, CI, Cox | Training sample/events, assay blocks, validation events, frozen direction, patient bootstrap/B, convergence criterion and cap; no selectors |
| Table 7, multimodal starting ranks | ICGC, PDAC, C | TCGA training n/events, both validation n/events, final/active rank definitions, convergence criterion/cap, absence of intervals; highlight this is one fit per starting rank |
| Table 8, single-modality biology | fgsea, RTK, PI3K, ERK, EMT, ECM, GO, HR, PurIST | Program classification rules, coefficient units/risk direction, association versus causal interpretation; keep pathway names comprehensible |
| Table 9, multimodal biology | expr, meth, HR, SD, TCGA, ICGC, PDAC; D4/P/flash k are internal labels | Signed-gene convention, named correspondence references, absolute r, non-bijective matching, rank-one-share denominator, class thresholds, 14 unadjusted tests |
| Table 10, completed/proposed | CpG; model/reference labels | Prefer prose rather than retaining this workflow table; if retained, define model/population and proposed status |
| Table 11 | EBMF, NMF, Cox, K, ELBO, C, PDAC and internal labels | Delete under B9; no caption repair required |
| Table 12, single-modality reconciliation | KL, EBNM, Cox; D4 is internal | Recast as mathematical target versus evaluated approximation; define weighted-monitor components and scales |
| Table 13, multimodal reconciliation | EM, SD, Cox, PL, KL where retained | Recast by component rather than dates; current priors, centered projections, cached penalties; optionally fold into appendix prose |
| Table 14 | fgsea, NES, padj, BH, ORA, GO-BP, TCGA, ICGC, PDAC, HR, SD, Cox, PurIST | Remove table/Appendix C; define these in the new Evaluation prose as needed |
| Table 15, single-modality simulation | C, Cox; “active,” “recovered,” “model/two-step” | Put model/two-step ordering explicitly in the note; number of datasets, split sizes, true shared/specific ranks by scenario, recovery and activity thresholds, outcome-oriented C; define SD/SE if displayed |
| Figure 1, risk tertiles | C; cohort/platform labels including RNA-seq/PACA-AU | Frozen risk direction, five cohort sample sizes, 616 total, tertile construction, survival/time units, meaning of censoring marks, log-rank test and whether any bands are shown |
| Figure 2, multimodal rank | No acronym requiring expansion if C/K labels are written out; internal scenario/noise names need replacement | n=120 training/120 validation, 50+50 features, true rank 3, three datasets per setting, seed ranges rather than confidence intervals |
| Figure 3, program survival | E(L), PDAC/Cox if used, program-number reference | Define programs as single-modality projections, pooled n=616, within-cohort grouping versus cohort-stratified test, time units, censoring marks, unadjusted p-values |
| Figure 4, modality shares | Program labels/asterisk, expr/meth if abbreviated | TCGA training n=144, K=7, rank-one squared-norm denominator, survival-active rule, no survival/PVE interpretation |
| Figure 5, simulation C difference | YFB, EBMF, C, K_init, scenario names in artwork; SE in caption | 15 paired seeds, 75 test patients/seed, true rank6, definition of scenarios, mean±one standard error, outcome-oriented concordance and unpenalized comparator |

Suggested expansions include empirical-Bayes matrix factorization (EBMF), nonnegative matrix factorization (NMF), evidence lower bound (ELBO; only where referring to a genuine bound), Bayesian information criterion (BIC), cross-validation (CV), partial likelihood (PL), proportion of variance explained (PVE), hazard ratio (HR), standard deviation (SD), confidence interval (CI), normalized enrichment score (NES), adjusted p-value (padj), Benjamini–Hochberg (BH), over-representation analysis (ORA), Kaplan–Meier (KM), epithelial–mesenchymal transition (EMT), extracellular matrix (ECM), and Gene Ontology biological process (GO-BP). Describe PurIST as a pancreatic-cancer subtype classifier; do not invent an expansion. Define consortium/cohort names in notes, and DeSurv as the comparator method attributed to Young et al., submitted. Use “the fgsea algorithm for ranked gene-set enrichment” rather than treating its software name as unexplained biological evidence.

### B7. Algorithms

**Locations:** Algorithm 1 spans PDF pp.4–5; Algorithm 2 p.7; source `:197–215,291–309`.

Replace callout containers with plain **Algorithm 1/2** headings followed by numbered, indented steps, at normal readable type size. Use “Inputs,” ordered updates and “Outputs/prediction.” Separate dataset-specific n/p/K settings into Evaluation or a compact note so the algorithm remains reusable. Remove implementation function names. In PDF keep a step with its continuation; if the full algorithm fits on one page, keep it together or start it on a fresh page. Do not use a framed block that continues onto a second page. Plain semantic ordered lists work in both formats and avoid a new LaTeX-only algorithm dependency. Incorporate A1–A6 corrections and preserve the relevant markers.

### B8. Figure-by-figure quality assessment

I rendered all 38 PDF pages to inspect layout, plus larger views of pp.12 and 19. `pdfimages -list` confirms every figure is raster. Their embedded pixel counts are adequate at their *current small sizes*, but font shrinkage and rasterized text remain visible. Increasing PNG DPI alone will not fix undersized labels.

| Figure | Pixels / effective PDF resolution | Current placement and issue | Recommendation |
|---|---|---|---|
| 1, external risk curves | 1350×900, 315 ppi; ~4.29 inches wide | Five panels and annotations compressed into 66% width. The producer uses 9×6 inches/base font11, so shrinking to ~4.29 inches reduces nominal text to ~5.2 pt. This explains poor readability despite adequate pixel density. | Redraw from saved risk scores, preserving groups/tests, at final physical size; use full available width and readable 8–10 pt plot text, perhaps a 2-column arrangement. Vector PDF plus ≥300-dpi PNG for HTML. |
| 2, multimodal rank | 1800×880, 346 ppi; ~5.20 inches wide | Two panels are reasonably sharp; base11 at 9-inch export shrinks to ~6.4 pt. Date-based panel label, long title/legend, non-manuscript scenario names. | Retain numerical values, ranges and color; replace labels, enlarge text at final size; vector PDF/high-resolution PNG. |
| 3, program survival | 1350×630, 315 ppi; ~4.29 inches wide | Two panels remain small with tiny test annotation and duplicated plot title. | Redraw at ~6–6.5 inches with larger text; vector PDF/PNG; scientific panel labels and caption instead of repeated title. |
| 4, modality shares | 1200×600, 298 ppi; ~4.03 inches wide | Seven category labels and legend are small; program identifiers need model context. | Enlarge and simplify labels, use common expression/methylation colors, define star in caption; vector PDF/PNG from saved shares. |
| 5, single-modality simulation | 1530×935, 267 ppi; ~5.73 inches wide | Historical title/legend use internal YFB/EBMF/scenario names; plot is noticeably less legible than body text. | If retained in §3.1, redraw from aggregate output with manuscript labels, same paired mean±SE; export vector PDF/PNG. |

Standardize plot typography, margins and legend positions. Figures may retain colors; text, captions and surrounding document remain black/white. Use Quarto format-conditional image references (or an equivalent existing mechanism) to select PDF versus PNG without fitting anything. Recheck output at actual reading size. No figure was regenerated during review.

### B9. Remove Table 11 and redistribute its six rows

**Locations:** §5.2 and Table 11, `:672–689`, Summary's `@sec-gaps` reference, any R9 mapping claiming a required table.

| Former claim row | Where its useful substance belongs |
|---|---|
| EBMF/NMF alone | Results/Biology: the actual loading/pathway correspondence, with restricted interpretation. Continuing work: matched comparisons under a defined prediction rule. Drop “the manuscript needs to prove” framing. |
| Two-step Cox | Results §3.2: exact principal difference/interval and comparator. Discussion: limited magnitude and cross-cohort estimand. Continuing work: final-method robustness comparisons. |
| DeSurv | Biology: overlap within a selected universe. Limitations: no matched predictive head-to-head evaluated. Continuing work: perform a matched benchmark. Remove directory-name audit narrative. |
| Rank selection | Methods/Results: actual CV rule, criterion disagreement, external evidence and starting-rank sensitivity. Limitations: dependence on initialization/validation and approximate objectives. Continuing work: survival-plus-reconstruction rank selection. |
| Multimodal added value | Results §3.5/Biology: same-donor performance and additional molecular structure. Limitations: convergence, precision and rank sensitivity. Continuing work: establish robustness and methylation-specific interpretation. Omit the unverified 0.685 benchmark. |
| Biology | Evaluation: characterization methods. Biology: supported labels/associations. Limitations: dependent inputs, unadjusted tests, provisional labels. Continuing work: characterize final programs and CpGs. |

The separate inventory remains the verification record; do not replace the deleted table with another source-audit table in the manuscript.

### B10. Discussion structure

**Locations:** current §5 and duplicated summary/limitations text elsewhere.

Proposed heading: **5 Discussion**, with:

1. **Interpretation and implications.** The single-modality model represents recognizable molecular programs and transfers a risk score across platforms; the tuned two-step comparison is small. The multimodal model captures additional molecular structure, with no demonstrated improvement in prediction. Explain why those are distinct contributions without claiming causal mechanism or successful survival-guided loading discovery.
2. **Limitations.** Frozen historical coefficient uncertainty; the approximate update/monitor/removal criteria; small multimodal sample/events and incomplete convergence; selection and rank sensitivity; feature screening; cohort heterogeneity/overlap; fixed-fit bootstrap intervals do not include model-training uncertainty; outcome-oriented simulation C; descriptive/dependent biology. Avoid reproducing every table number or development date.
3. **Continuing work.** One location for prospective supervision and the statement that α_F=0 is not the intended final method; joint survival/reconstruction rank selection, initialization/prior/feature robustness, matched comparators, final program characterization and CpG interpretation. Describe aims, not a sequence of code edits, dates or meeting decisions.

Cut §5.2 wholesale after redistributing its content. Move the unmarked future-method paragraph at `:324` here. Reduce the Summary's open-problems list and repeated model descriptions rather than restating Discussion verbatim. Keep the factual evaluated α_F=0 setting and its immediate meaning in Methods; that is necessary scientific disclosure, distinct from the prospective statement reserved for Continuing work.

### B11. Remove Appendix C and integrate its useful content

**Locations:** §2.5; §3.1; Appendix C `:730–775`; Figures 3 and Table 8 cross-references to Table 14.

Move biological characterization into **Evaluation §2.5**, with the following compact paragraphs/list. Each should name its input, population, statistic, multiplicity and interpretation:

- **Ranked enrichment and sensitivity:** nonnegative posterior-mean loadings across 2,064 selected genes; predefined collections; one-sided ranked enrichment, set sizes10–500; within-program/collection BH control at adjusted p<0.10; top50/100/150 over-representation against the selected universe. Interpret as gene-set correspondence. Define leading-edge and normalized enrichment score if used in Biology.
- **Subtype and overlap:** training latent activities versus subtype probability/calls (Spearman/Kruskal–Wallis, unadjusted); top270 loading genes and DeSurv lists restricted to the shared universe (Jaccard/hypergeometric, unadjusted). Training concordance and shared gene inputs are not independent validation.
- **External association and curves:** leading-edge signature mean after within-cohort gene standardization, univariable HR per unit; distinguish it from frozen program projections and overall risk. Explain the post hoc exclusion of the full custom PDAC collection for nonprognostic-program controls, grouping by within-cohort tertiles, and cohort-stratified log-rank tests for pooled program curves. Avoid rendering helper function names or internal C2/C4 labels.
- **Multimodal characterization:** manual signed-expression-gene interpretation, rank-one within-modality shares, absolute loading correspondence and 14 unadjusted per-SD joint-projection Cox tests; no completed enrichment/CpG annotation. Explain input/sign/matching limitations and define active classes.
- **Caveat paragraph:** dependent evidence and selected universe; within-training subtype check; post hoc controls; within-collection rather than global testing control; availability-driven choice of characterization methods and need for further assessment. State this as a limitation of the analyses, without citing an internal report.

Move Table 15 and Figure 5 into §3.1. Prefer keeping the figure for rank-dependent differences and a shortened table for recovery/false activity, removing duplicate concordance columns if they add no information. Explain true shared/specific counts per scenario. Delete Appendix C, its pagebreak and all referring text. Keep R6/R8 markers around the relocated paragraphs/exhibits. Replace “rows of Table 14” in Table 8 with named supporting analyses.

### B12. Remove the References section

**Locations:** YAML bibliography/link-citation settings, status box, §1.5 attribution, all `@Young2026` citations, References block `:715–718`, resulting PDF pp.23–24.

Use **“DeSurv (Young et al., submitted)”** at first substantive attribution and as needed, without a bibliography. Do not call its factor lists “published” (status box and Biology currently do). Preserve attribution but remove bibliography rendering from this document; leave the existing project bibliography unchanged.

### B13. Length and organization

The supplied PDF actually has **38 pages**: main scientific text begins on p.2 after title/contents and continues through p.23; References spills to p.24; Appendices A/B start on pp.25/29, and Appendix C starts on p.35. The “21-page main text” description does not fully capture the present layout. Main-text length should be counted again after restructuring.

Recommended sequence: **Methods → Evaluation → Results → Biological findings → Discussion**, with no Introduction for now, and only Appendices A/B. A short optional findings paragraph should not become a second Discussion.

Concrete cuts/reorganization:

- Remove the status box and long Summary (currently spans pp.1–2), or reduce the latter to one brief paragraph only if Andrew retains an overview.
- Keep notation dimensions once. Retain a shorter model-comparison table only for genuine single/multimodal differences; avoid repeating full prose and every implementation setting in Table 1.
- Delete the two-page evidence-gap table and its surrounding audit narrative.
- Replace the completed/proposed table by concise Biology/Continuing work prose if it duplicates Discussion.
- In Results, state the contrast and its interpretation, then refer to the table; do not repeat every cell. Keep the rank and convergence warnings.
- Shorten orientation/debugging history to the actual frozen convention plus the uncertainty limitation. Preserve the technical detail in Appendix A where needed to explain evaluated results, without dates/code names.
- Fold biological methods into Evaluation; place useful simulation exhibits in Results. Their new space is offset by the foregoing cuts, not by shrinking fonts.
- Expand derivations in Appendices A/B without treating appendix length as part of the 15–20-page main-text target.

There are also specific layout defects: Algorithm1 continues as a framed fragment on p.5; Table6 splits onto p.16; Table11 spans pp.21–22; References has one orphan entry line on p.24; Table13 has only its final rows at the top of p.34; Table14 spans pp.35–37; Table15 continues for three rows onto p.38. Deletion/restructuring removes several; use shorter tables, sensible page starts and unboxed algorithms for the rest. Do not solve them by reducing type further. Fix the malformed plain text “in -mm-sim, 0.04” at `writeup.qmd:500` to a proper table reference (currently present in HTML and PDF).

### B14. Complete, gap-free derivation additions

The following are the missing mathematical steps to insert, not permission to modify model code or silently change the evaluated method. Use the same notation across appendices, with (\bar\theta=E_q\theta\), \(\overline{\theta^2}=E_q\theta^2\), and conditional expectation (E_{-\theta}\) over all coordinates except the one being updated. Distinguish the random residual from its mean. State numerical floors and point-mass/prior-only boundaries separately from the algebra.

#### Common likelihood, working approximation and normal-means bridge

**A.1/B.1–B.2:** Write the Gaussian log likelihood explicitly:

\[
\ell_Y=\sum_{mij}\left\{\tfrac12\log(\tau_{mj}/2\pi)
-\tfrac12\tau_{mj}(y_{mij}-\mu_{mj}-\sum_kl_{ik}f_{mjk})^2\right\}.
\]

For the single block omit (m\) and set the intercept to zero. Define event/risk sets and Breslow ties:

\[
\ell_C(\eta)=\sum_a\left\{\sum_{i\in D_a}\eta_i
-d_a\log\sum_{i\in R_a}e^{\eta_i}\right\}.
\]

Show likelihoods/priors, the mean-field family, and the intended target
\(\mathcal L=E_q\ell_Y+E_q\ell_C-\mathrm{KL}_L-\mathrm{KL}_F-\mathrm{KL}_\beta\).
Specify the single-block normalization choice as in A1. Then give the exact score/Hessian from A2 and introduce the *replacement* curvature. Expand

\[
\ell_C(\hat\eta)+u^\top(\eta-\hat\eta)
-\tfrac12\sum_iw_i(\eta_i-\hat\eta_i)^2
=C+\sum_i\{(u_i+w_i\hat\eta_i)\eta_i-\tfrac12w_i\eta_i^2\}.
\]

For positive (w_i\), complete the square to obtain (z_i=\hat\eta_i+u_i/w_i\) and
\(C'-\tfrac12\sum_iw_i(z_i-\eta_i)^2\).
Explain the single-block weight floor versus the multimodal division-free (h_i=u_i+w_i\hat\eta_i\). Neither the working replacement nor the removal score should be identified with the exact target.

**A.3 and an explicit B cross-reference:** Derive the coordinate optimum and complete the square:

\[
\log q^*(\theta)=E_{-\theta}\ell+\log g(\theta)+C
=B\theta-\tfrac12A\theta^2+\log g(\theta)+C,
\]
\[
B\theta-\tfrac12A\theta^2
=-\tfrac12A(\theta-B/A)^2+\frac{B^2}{2A}.
\]

Hence
\(q^*(\theta)\propto N(x;\theta,s^2)g(\theta)\),
\(x=B/A\), \(s=A^{-1/2}\), for (A>0\).
Show shared-prior estimation explicitly:

\[
\hat g=\arg\max_g\sum_h\log\int N(x_h;\theta,s_h^2)g(d\theta).
\]

For the posterior corresponding to that same prior,
\(\mathrm{KL}(q\|g)=E_q\log N(x;\theta,s^2)-\log p_g(x)\).
For vectors sum the expectation over every observation, including the repeated Gaussian normalizer; then disclose A3/A5's implementation bookkeeping limitations. If (A=B=0\), specify the prior-only boundary rather than dividing by zero. Where floors are used, describe them as implementation rules.

#### Appendix A, section by section

**A.4, coefficient.** Do not say expectation is taken over the free coordinate itself. With a frozen projection (a_i=\tilde Z_{ik}\) and (H_i=\eta_{i,-k}\), expand

\[
E_{-\beta_k}(z_i-H_i-a_i\beta_k)^2
=C_i-2a_i\beta_k(z_i-EH_i)+a_i^2\beta_k^2.
\]

Multiply by \(-\alpha w_i/2\), sum, and collect
\(A=\alpha\sum_iw_ia_i^2\),
\(B=\alpha\sum_iw_ia_i(z_i-EH_i)\).
Write (x_\beta=B/A\), (s_\beta=A^{-1/2}\) before the prior/posterior calculation. Cancellation of α in x assumes α>0 and an inactive precision floor; code floors A at (10^{-10}\).

For the fitted normal prior, integrate it to get (x\sim N(0,s^2+\sigma^2)\), differentiate

\[
\ell(\sigma^2)=-\tfrac12\log(s^2+\sigma^2)
-\frac{x^2}{2(s^2+\sigma^2)}+C,
\quad
\ell'(\sigma^2)=\frac{x^2-(s^2+\sigma^2)}{2(s^2+\sigma^2)^2},
\]

and enforce the boundary:
\(\hat\sigma^2=\max(0,x^2-s^2)\).
Completing the posterior square gives mean
\(x\hat\sigma^2/(\hat\sigma^2+s^2)\) and variance
\(\hat\sigma^2s^2/(\hat\sigma^2+s^2)\). At equality \(|x|=s\), the fitted prior is also the point mass.

**A.5, activities.** Let
\(r_{ij}^{(-k)}=y_{ij}-\sum_{a\ne k}L_{ia}F_{ja}\) be random and (R_{ij}^{(-k)}=E r_{ij}^{(-k)}\). Expand the square and invoke mean-field independence:

\[
E_{-L_{ik}}(r_{ij}^{(-k)}-L_{ik}F_{jk})^2
=C_{ij}-2L_{ik}R_{ij}^{(-k)}\bar F_{jk}
+L_{ik}^2\overline{F_{jk}^2}.
\]

Multiplying by \(-\tau_j/2\) and summing (j\) produces
\(A_{L,k}=\sum_j\tau_j\overline{F_{jk}^2}\),
\(B_{L,ik}=\sum_j\tau_jR_{ij}^{(-k)}\bar F_{jk}\).
Give (x_{L,ik}=B_{L,ik}/A_{L,k}\), (s_{L,k}=A_{L,k}^{-1/2}\), the shared point-exponential prior and its EB optimization across patients.

**A.6, reconstruction loadings.** Show the counterpart explicitly rather than “same argument”:

\[
E_{-F_{jk}}(r_{ij}^{(-k)}-L_{ik}F_{jk})^2
=C_{ij}-2F_{jk}R_{ij}^{(-k)}\bar L_{ik}
+F_{jk}^2\overline{L_{ik}^2}.
\]

Collect (A^{\rm rec}_{F,jk}=\tau_j\sum_i\overline{L_{ik}^2}\),
(B^{\rm rec}_{F,jk}=\tau_j\sum_iR_{ij}^{(-k)}\bar L_{ik}\), then give x and s and the EB prior fit across genes.

**A.6, joint survival term.** Replace “up to factors of d” with a complete fixed-normalizer derivation. Set (a_{ij}=y_{ij}/d_k\), (C_i=\sum_{h\ne j}a_{ih}F_{hk}\), (H_i=\eta_{i,-k}\); then
\(\eta_i=H_i+\beta_k(C_i+a_{ij}F_{jk})\).
The coordinate-dependent expected square is

\[
a_{ij}^2E\beta_k^2 F_{jk}^2
-2a_{ij}F_{jk}\{\bar\beta_k(z_i-EH_i)-E\beta_k^2EC_i\}.
\]

Therefore

\[
A^{\rm surv}_{F,jk}=E\beta_k^2\sum_iw_i y_{ij}^2/d_k^2,
\]
\[
B^{\rm surv}_{F,jk}=\bar\beta_k\sum_iw_i(y_{ij}/d_k)(z_i-EH_i)
-E\beta_k^2\sum_iw_i(y_{ij}/d_k)EC_i.
\]

State that (d_k(F)\) varying with F invalidates this quadratic reduction. Compare this fixed-scale coupled expression with the implemented unnormalized/uncoupled survival branch, then combine its actual A/B with weights (1-\alpha_F,\alpha_F\), complete the square, and specialize to α_F=0. Equal half weights give the same x and half precision of the *same approximate unit-weight kernel*, not an exact posterior for the normalized joint model.

**A.7, precision.** Derive the expected residual using independence across factors:

\[
E(\sum_kL_kF_k)^2=(\sum_k\bar L_k\bar F_k)^2
+\sum_k(\overline{L_k^2}\,\overline{F_k^2}-\bar L_k^2\bar F_k^2),
\]
\[
Q_j=\sum_i\left[(y_{ij}-\sum_k\bar L_{ik}\bar F_{jk})^2
+\sum_k(\overline{L_{ik}^2}\,\overline{F_{jk}^2}-\bar L_{ik}^2\bar F_{jk}^2)\right].
\]

Differentiate (n\log\tau_j/2-\tau_jQ_j/2\) to get
\(n/(2\tau_j)-Q_j/2=0\), hence \(\hat\tau_j=n/Q_j\).
State the evaluated residual denominator floor (n\,10^{-8}\).

**A.8–A.10.** Explain the recorded monitor component versus a true expected working quadratic (A3), its weighted reconstruction/survival scales, and the negative-KL bookkeeping defect; distinguish this statistic from the intended target. Use the frozen (D_*\) for prediction and state the saved variance limitation without dated implementation-history prose. Keep the mathematical-versus-evaluated comparison accurate after these changes.

#### Appendix B, section by section

**B.3, intercept.** Show

\[
\frac{\partial E\ell_Y}{\partial\mu_{mj}}
=\tau_{mj}\sum_i(y_{mij}-\mu_{mj}-\sum_k\bar l_{ik}\bar f_{mjk})=0,
\]

and solve for the stated mean residual. Qualify the expected-surrogate centering issue from A6; exact Cox shift invariance is not a proof for that approximate objective.

**B.4, projection and predictor moments.** Explicitly expand (E Z_{ik}^2=(EZ_{ik})^2+\sum_{mj}\tilde y_{mij}^2\operatorname{Var}f_{mjk}\). Independence across factors and of β/F gives

\[
\operatorname{Var}\eta_i
=\sum_k\{E\beta_k^2 E Z_{ik}^2-\bar\beta_k^2(EZ_{ik})^2\}
=\sum_k\{E\beta_k^2\operatorname{Var}Z_{ik}+(EZ_{ik})^2\operatorname{Var}\beta_k\}.
\]

This derives the correct final formulas already present. Do not imply predictors for different patients are independent.

**B.5, activities.** Apply the random-residual square expansion from A.5 to both modalities; sum (m,j\) to obtain
\(A_{L,k}=\sum_{mj}\tau_{mj}\overline{f_{mjk}^2}\),
\(B_{L,ik}=\sum_{mj}\tau_{mj}R_{mij}^{(-k)}\bar f_{mjk}\).
State x/s, complete the square and show the prior/posterior construction below.

**B.6, reconstruction loading.** Repeat A.6's square expansion within modality m to obtain
\(A^{\rm rec}=\tau_{mj}\sum_i\overline{l_{ik}^2}\),
\(B^{\rm rec}=\tau_{mj}\sum_i\bar l_{ik}R_{mij}^{(-k)}\).

**B.6, survival loading.** Set (a_i=\tilde y_{mij}\), (U_i=Z_{ik}^{(-j)}\), (V_i=\eta_{i,-k}\), (f=f_{mjk}\). Expand

\[
\eta_i=\beta_k(a_if+U_i)+V_i,
\quad
\eta_i^2=a_i^2f^2\beta_k^2+2a_if(\beta_k^2U_i+\beta_kV_i)
+(\beta_kU_i+V_i)^2.
\]

Mean-field independence gives
(E(\beta_k^2U_i)=E\beta_k^2 EU_i\),
(E(\beta_kV_i)=\bar\beta_kEV_i\).
In (E[h_i\eta_i-w_i\eta_i^2/2]\), collect linear coefficient

\[
a_i[\bar\beta_k(h_i-w_iEV_i)-w_iE\beta_k^2EU_i]
\]

and quadratic precision (w_ia_i^2E\beta_k^2\). Sum i and add reconstruction terms to get B.6's A/B. Complete the square and give (x_{mjk}=B_{mjk}/A_{mjk}\), (s_{mjk}=A_{mjk}^{-1/2}\). Explain sequential updates and projection refresh after each feature, and old-prior posterior versus subsequently refitted prior.

**B.5/B.6, mixture weights and moments.** Show the half-line completion

\[
-\frac{(\theta-x)^2}{2s^2}-\lambda\theta
=-\frac{[\theta-(x-\lambda s^2)]^2}{2s^2}
-\lambda x+\frac{\lambda^2s^2}{2},
\]

and its negative-half-line counterpart after replacing θ by −θ. Integrate the normalized Gaussian on each half-line to obtain the displayed slab weights, with the spike weight proportional to \((1-\pi)N(x;0,s^2)\). Normalize all weights. For a positive truncated normal with location m, define (a=m/s\), (r=\phi(a)/\Phi(a)\), and show its moments

\[
M_+(m)=m+sr,\qquad S_+(m)=m^2+s^2+msr.
\]

These follow by integrating a standard normal and using
\(\int_{-a}^\infty z\phi(z)dz=\phi(a)\),
\(\int_{-a}^\infty z^2\phi(z)dz=\Phi(a)-a\phi(a)\).
For (m_+=x-\lambda s^2\), (m_-=-x-\lambda s^2\), use
\(Ef=w_+M_+(m_+)-w_-M_+(m_-)\),
\(Ef^2=w_+S_+(m_+)+w_-S_+(m_-)\).
The point-exponential case has only spike and positive components. Explain numerical tail evaluation after giving these mathematical moments. Show the shared-prior marginal-likelihood fit and the timing/cached-KL caveat for both L and F.

**B.7, coefficient.** Expand \((Z_{ik}\beta_k+V_i)^2\), average over F/other coefficients, and collect
\(A_\beta=\sum_iw_i E Z_{ik}^2\),
\(B_\beta=\sum_iEZ_{ik}(h_i-w_iEV_i)\).
Give x/s explicitly and reuse the normal-prior marginal likelihood derivative and posterior square completion from A.4. Show the normal KL expression by taking expectations of the log density ratio, including the point-mass boundary. Explain that \(1.96\) classification is an operational posterior mean/SD rule, not an adjusted multiple-test threshold.

**B.8, precision.** Sum the residual expansion from A.7 over patients within each modality-feature pair, differentiate the Gaussian log likelihood and solve for τ. State the zero-residual unbounded boundary and any implemented numerical protection rather than implying it is excluded by intercepts as a theorem.

**B.9–B.13.** Demonstrate reconstruction and predictor invariance under the scale map by substituting each transformed factor, then transform first/second moments, projection moments and all prior parameters (A10). Show Jacobian/pushforward treatment and collapsed-column boundary. Correct the stopping denominator, describe the actual cached-penalty deletion score, limit Jensen/centering statements, and distinguish exact raw-versus-centered frozen prediction from expected-surrogate behavior. Remove the dated base-versus-current audit table or rewrite it by model components, with no repository terminology.

## Checks performed, limits, and next handoff

- Confirmed `writeup/draft` and a clean tree before review. Reviewed the main source, both appendices, the 62-entry inventory and applicable author feedback. Separate read-only mathematical/numerical checks informed the findings above; the coordinating review checked the key curvature, scalar-KL and cached-prior code evidence directly.
- Inspected all 38 PDF pages as rendered page images/contact sheets and larger views of the external-KM/program-table pages. Identified five embedded raster figures, their effective resolution, small plot typography, split algorithm/tables and reference orphan. Mathematical glyphs appeared in the inspected PDF pages; this is not a proof of every line's typesetting correctness.
- Static HTML inspection found 651 math spans, embedded image/script assets, no duplicate IDs, no missing local fragment targets, no unresolved cross-reference tokens, and consistent column counts for all 15 tables. The malformed “-mm-sim” prose remains despite otherwise valid links. A live HTML browser/responsive-layout inspection was not completed in this review; static markup inspection does not establish MathJax runtime behavior or mobile readability.
- Numerical checks used saved aggregate outputs and local saved fit/score objects only. No models were fitted, no biological characterization was rerun, and no saved data/results were changed. Percentile bootstrap producer settings were inspected, but the complete bootstrap distributions/intervals were not independently regenerated. This limits claims of exhaustive computational reproduction. Cross-cohort donor overlap beyond the documented matched ICGC subset was not resolved; the draft appropriately discloses that uncertainty.
- The historical executable/environment that produced the evaluated fit was not reconstructed. Current code and cached settings establish the documented compatibility/limitations, not a byte-for-byte reproduction of its full fitting trajectory. Effects of the monitoring/KL defects on stopping and recommended rank remain untested.
- No external literature or new biological mechanisms were investigated. Gene-function labels remain interpretations of the existing characterization, with the stated uncertainty/dependence; this review does not validate novel methylation mechanism.

**Next:** Andrew authorizes the review commit if desired. Claude then checks each finding against its source, applies confirmed scientific corrections and all author-requested Part B revisions, records a disposition for every item, rerenders HTML/PDF, and requests Andrew's approval before merging. This review neither commits nor merges anything.
