# Verification of the revised writeup

**Reviewed:** 2026-10-14, branch `writeup/draft`, HEAD `17a8bb1`. The tree was clean at the start. Scope checked against `git diff 88cafda..17a8bb1 -- paper/writeup` and the complete changed-file list. Earlier findings refer to `review.md`; claimed fixes refer to `review_dispositions.md`.

## Verdict

**Substantial improvement, but further revision is needed before approval.** Most scientific corrections and all principal numerical corrections are present. The report now reads as an authored scientific account, has appropriately qualified comparisons, and uses sharp vector figures. The main performance and biological values checked below agree with their sources.

The dispositions overstate completion of the derivations, self-contained captions, and unsplit algorithms. There are also new or remaining scientific inaccuracies: an undefined zero-column normalized target, a contradictory orientation account, an incorrect statement about the actual multimodal precision floor, and an exact-target/working-quadratic conflation. These are report corrections; this verification does not authorize changes to model code or historical results.

Only this file was written. Nothing is marked fixed by this verification. No report, appendix, figure, source inventory, model code, default, data, or saved result was edited; no commit, push, or merge was made.

## Part A — finding-by-finding verification

“Verified” means the earlier requested correction is present and agrees with the inspected evidence. Additional wording refinements do not necessarily invalidate that correction. “Incomplete” means a material part of that finding remains unresolved. Severity and exact remaining corrections follow in the numbered problems below.

| Earlier finding | Verdict | Evidence in revision and implementation |
|---|---|---|
| **A1 Joint target versus evaluated algorithm** | **Incomplete** | `writeup.qmd:144–158` and Appendix A.1 now distinguish random loadings, the partial-likelihood target, and sweep-fixed posterior-mean projections. This is correct. However, A.1:29 divides by a random loading norm without defining the zero-column case under its spike prior; A.5:154 and Table 11/A.11:292 still attribute a second-moment quadratic update to the exact target. Fix V2/V4. |
| **A2 Cox working curvature/floor** | **Verified** | Appendix A.2:45–80 gives the Breslow score, full coupled negative Hessian, its diagonal, the replacement curvature, square completion, and separate raw score/floored weight. This matches `code/fit_cox_on_yf.R:103–119`. B.2:36–42 correctly uses unfloored weights and division-free h, matching `code/multimodal_yfb_helpers.R:51–86`. |
| **A3 Activity KL bookkeeping and monitoring** | **Verified as a disclosure; code defect remains** | A.3 gives the correct vector normalizer. A.9:251–261 gives the correct discrepancy sign and sweep-start survival limitation. The scalar activity precision from `code/update_L_surv_YFB.R:78–100` is passed into `compute_ebnm_kl()` at `fit_cox_on_yf.R:707–708`; `code/compute_elbo.R:78–86` counts its normalizer once. The issue was independently rechecked, not merely accepted from the earlier review. Results/Table 5 call the output a monitoring statistic. No corrective fitting has occurred. |
| **A4 Stored D_* and feature-count normalization claim** | **Verified, with residual wording precision** | Algorithm 1 output/prediction and A.10:271–275 use the saved preceding-sweep normalizer D_*; code `fit_cox_on_yf.R:643–646,734–737,926` and `predict_cox_on_yf.R:104–120` agree. The saved fit confirms the small norm discrepancies. The feature-count-invariance claim is gone. Table 1/Equation 2/Figure 4 still use “unit-length” too literally for final loadings divided by stored norms; see V9. |
| **A5 Cached old-prior KL in multimodal removal** | **Verified as a disclosure; code behavior remains** | Main text:242, B.6:132 and B.11:206 explicitly cover both activity and loading penalties, their earlier priors, and the absence of a current-prior ELBO guarantee. Matches `multimodal_yfb_updates.R:359–367,491–502` and fitter:528–531,577–581,734–754. The remaining non-strict-improvement wording is V10. |
| **A6 Jensen/variance correction and centering** | **Verified** | B.3:57 distinguishes exact frozen-score shift invariance from the expected surrogate. B.11:210 gives Jensen's direction, the common-shift counterexample, and the full Hessian/covariance correction before describing the diagonal alternative. The earlier unsupported monotonic variance-gap assertion is removed. |
| **A7 Shared-program recovery** | **Verified** | Results:383 and Table 3:402–416 report 97–100% and include rank 12. Raw S18 means at ranks 6/12/20 are 1/.97778/1 for all-shared and 1/1/.96667 for hybrid. All displayed recovery/activity cells match the saved output. |
| **A8 Complete derivations** | **Incomplete** | Most requested algebra has been added and checks correctly. However, the normal-prior convolution/posterior completion is still partly asserted, B.6 still skips its reconstruction expansion, B.8 skips differentiation, and multimodal A=0 boundaries are not separated explicitly from the single-block floor. See V5 and the section-by-section mathematical audit. B.8 also introduces the wrong precision-boundary description, V1. |
| **A9 Concordance rounding and stopping denominator** | **Verified** | Results:432 and Table 4 correctly use 0.6285 versus 0.6026 and difference 0.026. B.10:192–196 uses max(1, max abs old, max abs new), agreeing with fitter:64–69,514–540. A.9's stored stopping history also agrees with the saved sweep-20 stopping point. |
| **A10 Scale maps and feedback boundaries** | **Verified in the appendices; main-text qualification incomplete** | B.6:114 retains both coefficient moments. B.9:174–186 gives the moment/prior maps, pushforwards, KL invariance and constant-projection fallback; fitter:190–214 agrees. Main text/Table 1/Algorithm 2/B.13 still say SD becomes one without the exception. Main text:158 correctly allows stopping to affect the returned unsupervised-loading iterate. Add the compact qualifications in V9/V10. |
| **A11 Revisit coverage** | **Incomplete** | All pairs are balanced and nested; all nine IDs occur in the main source. Tables, figures, captions, revised notation and appendix introductions are covered. However, the single-model numerical-floor sentence at A.3:112 and the point-exponential prior-specific posterior section A.4:114–131 remain unmarked. The latter depends on the prior-family method work. Add R1 (and R2 where needed), while generic square-completion/KL identities may remain unmarked. |
| **A12 Survival-curve ordering** | **Verified** | Results:456 now says curves “generally favored” the low-risk tertile and reports three significant cohort tests. Recomputed log-rank p-values agree, and the new curves show the crossing/overlap for which the qualification is needed. |
| **A13 Unevaluated low-variance failure** | **Verified** | Discussion:675 confines the observed planted-program failure to multimodal simulation. Continuing work:695–696 no longer claims low-variance recovery is impossible for unsupervised estimation. The single-modality supervision question is explicitly untested. |

## Part B — all fourteen author items

The numbering here follows Part B of the earlier review and the review prompt (which separates noise labels from captions).

| Item | Verdict | Evidence and remaining action |
|---|---|---|
| **B1 Title/subtitle** | **Verified** | The requested title appears in source, PDF metadata/title page and live HTML. No subtitle remains. |
| **B2 Authored appearance** | **Verified** | No draft header/status box, reviewer/agent/inventory/untraced-source prose appears in the rendered scientific text. The short Summary is authored narrative. Technical discussion of implementation approximations is scientific content, not a drafting trace. |
| **B3 Repository references/internal labels** | **Verified** | Source-ID/path/function/argument/internal-label searches and live rendered HTML inspection found none of the prohibited material. D4, “flash k,” “expr-only,” visible S IDs and output selectors are gone. Defined SM/U/EO in Table 10 are scientific reference abbreviations, not unexplained old configuration labels. Required hidden Revisit comments remain. |
| **B4 Dated/process references** | **Verified** | Dated meeting/report/version references are removed. The ordinary title-page publication date remains; it is not a dated development reference. Words such as “meeting the stopping criterion” are not meeting-note references. |
| **B5 Noise levels** | **Verified** | Table 6, Figure 3/artwork, Evaluation and Results use lower 0.25× versus higher 1× baseline feature-noise SD, with Gamma(2,2) precisions stated. The labels describe heteroscedastic noise multipliers, not one common absolute SD. |
| **B6 Self-contained exhibits** | **Incomplete** | Captions are substantially better, but several abbreviations still lack definitions in their own exhibit, and a few symbol/population details remain. See the explicit caption checklist below. Definitions elsewhere in prose do not satisfy this requirement. |
| **B7 Algorithms** | **Incomplete** | Boxes are gone in HTML and PDF, and ordered steps are readable. **Algorithm 1 still spans PDF pp.3–4; Algorithm 2 spans pp.6–7, including an output sentence continued on p.7.** This fails the current explicit no-split requirement. See V6. |
| **B8 Figures** | **Verified for no model refitting, values, vector/PNG formats and printed readability; saved-input-only condition not met** | Both scripts avoid factor-model/Cox-model fitting. Five PDF figure files are vector; the report has zero raster image objects. Five HTML PNGs are 1950 px wide with 300-dpi metadata. Numerical values reproduce. However, `make_writeup_figures.R:95–109` reads local held-out expression to reconstruct frozen program projections; it does not read exclusively saved fits/scores/aggregate tables. S63 discloses this. See the figure audit and V11; no patient-level assets should be committed to hide the dependency. |
| **B9 Remove gap table/framing** | **Verified** | The gap table and claim-by-claim section are gone in both formats. Their useful content is distributed through Results, Biology, Discussion and Continuing work. |
| **B10 Discussion** | **Verified** | Section 5 has Interpretation and implications, Limitations, and Continuing work. The prospective α_F=0 statement is concentrated in Continuing work. Necessary factual disclosure of the evaluated setting remains in Methods/Results. |
| **B11 Remove Appendix C/integrate methods/simulation** | **Verified** | Only Appendices A/B remain. Evaluation §2.5 gives biological inputs, populations, statistics, multiplicity and caveats in prose/list form. Simulation Figure 1/Table 3 are now in Results §3.1. |
| **B12 No References** | **Verified** | No References section, bibliography rendering or formal citation block remains. DeSurv attribution is “Young et al., submitted.” Project bibliography is unchanged. |
| **B13 Organization/length/layout** | **Incomplete** | Organization and redundancy are improved, but the main report still occupies pp.1–27 (26 after excluding the title/contents/Summary page). The PDF is 42 pages overall; A begins p.28 and B p.36. Algorithms split and two tables have duplicate heads before captions. No table's data rows were found spread over multiple pages; the earlier apparent table split is a duplicate-header artifact on the same page. See V6 and ranked length recommendations. |
| **B14 Gap-free derivations** | **Incomplete** | Most final equations and new intermediate steps are correct, but the remaining steps in V5 and target/boundary corrections V1–V4 must be resolved. “Every step included” is too strong as a disposition. |

### Remaining caption checklist (both formats)

Use short local notes; expanding abbreviations does not require restoring source IDs or lengthy process explanations. Existing gene symbols need not all be expanded as acronyms.

| Current exhibit | Definitions/details still missing locally |
|---|---|
| **Table 1** | Define expr/meth in the data row, the dimensions/roles of L, F, n/p/K if the table is meant to stand alone; qualify D_* and collapsed projections as in V9. It presently sends the reader to §1.2 for part of its normalization definition. |
| **Table 2** | TCGA-PAAD, CPTAC, RNA-seq, PDAC and DeSurv are not expanded/described in this caption. ICGC/PACA-AU, IPMN, CpG and 450K are defined. Add the missing cohort/assay names rather than assuming their Methods definitions travel with the table. |
| **Table 3** | The simulation definition is largely self-contained. Define the tilde-beta symbol as the model's survival coefficient and explicitly name the model arm if extracted alone. No substantial missing acronym expansion. |
| **Table 4** | Define PDAC, CI (confidence interval), LASSO (least absolute shrinkage and selection operator), PACA-AU and RNA-seq. C is already defined. |
| **Table 5** | BIC and PL are defined, but CV is never explicitly expanded as cross-validation; define n and p used in the parameter count, and identify the 273-patient training population. |
| **Table 6** | Define SD locally as standard deviation. Identify the evaluated framework's intercept/signed-loading/prior/removal components, or explicitly name its algorithm as a supporting cross-reference while retaining the substantive description. C, fit counts, sample sizes and uncertainty absence are already stated. |
| **Table 7** | Define TCGA, ICGC, PDAC, CI and CpG locally. Expr./meth. and C are defined. The single-modality comparator's realized stopping setup should be given if retained (V12). |
| **Table 8** | Define TCGA, ICGC, PDAC and SD locally; state 75 training deaths and the 30/40 validation deaths if it is to stand independently of Table 7. C, classification and the cap are defined. |
| **Table 9** | Define PVE, BH, TCGA, PDAC, and “fgsea” as the ranked gene-set enrichment procedure. HR, PurIST, RTK, PI3K, ERK, EMT, ECM and GO-BP are defined. Avoid using “not prognostic” as a definitive null label where only a nonsignificant association was observed. |
| **Table 10** | Define TCGA, ICGC, PDAC and SD. SM/U/EO and R/S are already explained. Use the zero-SD classification convention from V10. |
| **Table 11 (Appendix A)** | Expand KL; clarify the intended target versus the *working-quadratic* second-moment statement in V4, and define its scale/weight symbols compactly if retained. |
| **Table 12 (Appendix B)** | Expand KL and sd, and state the collapsed-projection exception. |
| **Figure 1, simulation C difference** | Explicitly define C as Harrell's concordance and say higher C is better; outcome orientation and paired standard errors are already explained. |
| **Figure 2, external risk curves** | Define PDAC, PACA-AU and RNA-seq locally. Expand C to concordance index (the caption currently says “Harrell's C”). Grouping, censoring marks, direction and time units are clear. |
| **Figure 3, multimodal final rank** | Define SD locally. Means/ranges, true rank, noise multipliers and training/feature sizes are otherwise clear. |
| **Figure 4, program curves** | Expand PDAC and correct the literal unit-length definition (V9). Time, scores, grouping, censoring and stratified tests are explained. |
| **Figure 5, modality shares** | Define TCGA, CpG and SD locally. Share denominator and its interpretation are correctly specified. |

## New and residual problems, ranked

### V1 [P2, introduced] Multimodal precision boundary describes the wrong code path

**Location:** `appendices/multimodal.qmd:168`, rendered Appendix B.8, PDF p.40.

The new sentence says a zero expected residual makes the implementation stop with an error instead of imposing a floor. The evaluated fitter calls the chunked residual calculation and then `multimodal_yfb_tau_from_residuals()` (`code/fit_multimodal_yfb.R:505–509`). That helper floors residual sums at `.Machine$double.eps` before dividing (`:639–645`). The error belongs to the other update helper at `code/multimodal_yfb_updates.R:577–578`, not this fitting path.

**Exact replacement:** “Without numerical protection, Q_mj=0 gives an unbounded likelihood and no finite maximum-likelihood precision. The evaluated fitting path uses τ_mj=n/max(Q_mj, ε_machine), where ε_machine≈2.22×10^−16.” Distinguish the mathematical optimum from the protected computation. The removal score uses the same helper, so its description must agree.

### V2 [P2, remaining after A1] Define the random normalized target at a zero column

**Location:** `single_modality.qmd:27–35`; motivating target in `writeup.qmd:148`.

The spike-and-exponential prior gives a column positive probability of being identically zero, but the target divides by d_k(F)=||F_k||. The evaluated algorithm's ridge, described later for D_*, does not define that random-loading target. This is a genuine support/boundary issue, not a typographical one.

**Fix:** Define the normalized direction explicitly as F_k/||F_k|| when the norm is positive and the zero vector when the norm is zero. Alternatively define a strictly positive ridge normalizer as part of the target, then use that same definition throughout. State which target convention is adopted; do not imply that the historical evaluated fits are changed. Keep the fixed-normalizer approximation explicit when deriving Gaussian coordinates.

### V3 [P2, introduced] Orientation account is internally contradictory

**Location:** `single_modality.qmd:269`, PDF p.34.

It says coefficients are negated when η is discordant with risk, then says this check left the evaluated η increasing with longer survival. The present code uses risk-oriented concordance (`reverse=TRUE`) and flips below 0.5 (`fit_cox_on_yf.R:891–904`); the evaluated historical fit used a different convention. Removing dates must not erase that distinction.

**Exact replacement:** “For the evaluated single-modality fit, higher η denotes longer survival. All external evaluations therefore use risk=−η, with that direction fixed from training. The stored second moments equal the squared means, so this fit has no recoverable coefficient variance.” If current orientation behavior is discussed, put it in a separate clearly distinguished sentence. No dates, function names or debugging history need appear in the manuscript.

### V4 [P2, remaining] Second-moment coefficient statements still refer to the exact target

**Location:** `single_modality.qmd:154,292` (Table 11).

The exact E_q log partial likelihood is nonquadratic in a coefficient. Replacing a squared projection by E(Z²) is a property of the expected **working quadratic**, not an exact update of the joint target.

**Fix line 154:** “Under the working-quadratic approximation with random loadings, the quadratic coefficient involves E(Z_ik²); the evaluated plug-in update instead uses Z_ik².”

**Fix Table 11 cell:** “Nonconjugate under the exact target; the random-projection working quadratic uses E(Z_ik²).” Keep the evaluated-algorithm column as the plug-in approximation. This completes the target/algorithm distinction rather than suggesting a missing exact closed-form update.

### V5 [P2 for the requested completeness standard] Finish the few remaining derivation steps

The following are concise insertions, not a request for another wholesale rewrite.

1. **Normal-prior convolution and posterior square (A.5:156–167).** The marginal normal distribution is asserted rather than obtained by the requested integration. For v=σ²>0, insert

   $$
   \frac{(x-b)^2}{s^2}+\frac{b^2}{v}
   =\frac{(b-m)^2}{V}+\frac{x^2}{s^2+v},\qquad
   V=\frac{s^2v}{s^2+v},\quad m=\frac{xv}{s^2+v}.
   $$

   The normalized product is
   $$N(x;b,s^2)N(b;0,v)=N(x;0,s^2+v)N(b;m,V).$$
   Integrating b proves the marginal; retaining the normalized b factor proves the posterior. The already correct derivative, constrained maximum and v=0 point-mass boundary then follow. B.7 may reference this shared derivation explicitly.

2. **Reconstruction loading expansion (B.6:91–94).** Supply the currently omitted step, with r random and R its mean:

   $$
   E_{-f}(r_{mij}^{(-k)}-l_{ik}f)^2
   =c_{mij}-2fR_{mij}^{(-k)}\bar l_{ik}+f^2\overline{l_{ik}^2}.
   $$

   Multiply by −τ_mj/2, sum i, collect the displayed A_rec and B_rec, then add the survival coefficients and apply the normal-means square completion. Referring to A.7 explains why the formula is correct but does not meet the author's explicit requirement to show this step in B.

3. **Precision derivative (B.8:164–167).** Insert

   $$
   \partial_{\tau_{mj}}\{n\log\tau_{mj}/2-\tau_{mj}Q_{mj}/2\}
   =n/(2\tau_{mj})-Q_{mj}/2=0.
   $$

   Solve τ=n/Q for Q>0, then state the actual numerical floor from V1.

4. **Multimodal prior-only boundary (B.5–B.7).** State that x=B/A and s=A^−1/2 assume A>0; for A=B=0 the relevant update is prior-only, not division by zero. `multimodal_yfb_updates.R:39–49,181–187,524–527` implements these boundaries. Appendix A.3's immediately following single-modality precision floor must not be read as applying to B.

The remaining half-line completions, truncated moments, coefficient derivative, predictor moments, scale pushforwards and removal-score derivations are checked in the audit table below. No additional algebra errors were found in those checked expressions.

### V6 [P2 presentation requirement] Unsplit algorithms and clean table heads are not achieved

**PDF locations:** Algorithm 1 pp.3–4; Algorithm 2 pp.6–7; Table 9 p.22; Table 11 p.35; Appendix B.13 heading p.41/table p.42.

Both algorithms are unboxed, which is an improvement, but still cross page boundaries. Table 9 and Table 11 each have a stray column header **above** the caption followed by the normal column header below it on the same page. Their data rows do fit on a single page. Thus do not report these as multi-page data tables; the actual fault is the duplicated head/caption order. B.13's heading is orphaned before its table.

**Fix:** Keep each algorithm, including heading/inputs/outputs, together using a PDF-only keep-together/minipage block after confirming it fits at normal font size; start it at a suitable page position. Keep HTML's plain semantic lists. For tables that fit on one page, use a nonbreaking table layout with one caption/head, or correct the generated longtable head/caption handling. A blanket proportional `needspace` before a longtable is not proof of a clean result. Reinspect the actual PDF after rendering, including caption, first/last row and nearby headings. Do not shrink plot or body type to force a fit.

### V7 [P2 presentation requirement] Caption self-containment remains partial

Use the exhibit-by-exhibit checklist above. `review_dispositions.md` claims every caption defines its abbreviations, which is not true for TCGA/CPTAC/PDAC/LASSO/CI/CV/SD and other listed cases. Keep definitions concise, but local. Correct this disposition after the captions themselves are fixed.

### V8 [P3] Revisit comments miss two method-dependent blocks

**Locations:** `single_modality.qmd:112,114–131`.

Add a small R1 block around the single-block numerical precision-floor rule and an R1 block around the point-exponential/truncated-normal posterior section; add R2 if its parameterization changes with scaling. Generic identities in A.3 are not by themselves method-dependent. All current markers are balanced, so preserve nesting when adding the blocks. Everything else inspected after the restructuring—main equations, algorithms, figures, captions, tables, and appendix model/update sections—is covered.

### V9 [P3] Propagate normalization exceptions into summaries and captions

**Locations:** Table 1 `writeup.qmd:122–128`; Equation 2 description; main multimodal scale paragraph `:234`; Algorithm 2 `:267`; Figure 4 caption `:627`; B.13:233.

- For saved single-modality predictions use D_* consistently, or explicitly state that D in Table 1 is shorthand for D_*. The final F divided by the preceding sweep's stored norms is not exactly unit length, and inactive columns have zero direction. Replace Figure 4's “unit-length loading direction” with “loading column scaled by its fixed training normalizer.”
- For multimodal mean projections, say “unit training SD for nonconstant columns; effectively constant columns are left unchanged.” Appendix B.9 already states the c_k=1 fallback correctly. If an exact threshold is given, the code uses SD≤sqrt(ε_machine).

These are small refinements relative to the major A4/A10 fixes, not reasons to replace the saved normalizers or refit.

### V10 [P3] State classification/removal boundary conventions

**Locations:** main classification `writeup.qmd:246`, Algorithm 2 removal `:270`; B.7:154 and B.11:206; relevant captions.

- The coefficient classification uses |β̄|/max(SD(β),10^−12) in the actual fitter (`:553–556`). Bare 0/0 is undefined for a shrunk coefficient point mass. Include the protected denominator or explicitly declare a coefficient point mass at zero inactive.
- The deletion check requires K>1 (`fit_multimodal_yfb.R:591`); state that this implementation does not prune the final factor to rank zero.
- Accepted Δ≥0 means the implemented score **does not decrease**, not necessarily “improves.” Make B.11's sentence non-strict.

### V11 [P3 reproducibility qualification] Figure reconstruction is not aggregate-only

**Location:** `figures/make_writeup_figures.R:85–109`; source inventory S63.

The script reads the saved fit, training gene names and held-out saved risk scores, but also calls the local raw-expression loader to reconstruct frozen per-program projections. The script header and S63 disclose this correctly. The pruning-rank script reads the saved sweep CSV only. Neither refits a molecular factorization or a Cox model; Kaplan–Meier/log-rank calculations for the permitted redraw are separate descriptive computations.

The requested “only saved fits, scores and aggregates” assertion cannot be confirmed literally. Keep the dependency explicit in repository rendering/figure instructions; the already embedded HTML/PDF themselves do not need private data to open. If strict saved-only regeneration is required, use an existing saved program-curve/projection source if one exists; otherwise obtain a scope decision before creating an additional patient-level cache. Do not commit patient-level data to make the script self-contained. This is a qualification of provenance, not evidence that the redrawn curves are wrong.

### V12 [P3 completeness] State the isolated comparator's realized settings if retained

**Location:** Comparators `writeup.qmd:350` and Table 7 caption `:560`.

The identification as a single-block version of the model is correct: `code/run_multimodal_yfb_simulation.R:139–150` calls the same single-modality fitter with point-exponential loadings, a normal coefficient prior and its default survival weights. This should not be flagged as a different mathematical method. However, that wrapper uses an iteration cap of 100 by default, rather than the principal fit's 300, and receives unstandardized expression; the main runner does not pass a replacement control for this arm (`results/multimodal_real/run_multimodal_real_fit.R:143–145`). The text discloses the raw-expression scale but leaves its realized stopping setup unspecified.

**Fix if retaining this supporting row:** Call it a “single-block fit to unstandardized log2 expression” and give its 100-sweep cap and realized prediction convention in Evaluation or a table note. Do not imply it reproduces the principal pooled fit's complete configuration. Its saved numbers agree, and no methodological error is established merely by its current label. Dropping this supporting row is an optional length choice, not a required scientific correction.

### V13 [P3] Small remaining wording/rounding/reference corrections

- Summary `writeup.qmd:88`: “recovers the true rank ... at lower noise” is stronger than the detailed mean ranks 3.26/3.22/2.48. Say “returns mean ranks close to the true value at lower noise, with poorer recovery in the low-variance scenario.” Do not imply every fit selects exactly three factors.
- Main text `:179,190` calls the monitoring section “Appendix A.8”; its destination is now **A.9**. Links resolve, but their displayed labels are stale.
- Main text `:541` calls the matched difference 0.010. The saved difference is .6852071−.6757396=.0094675. Use “about 0.01” or “0.0095”; strict three-decimal rounding is 0.009. This imprecision predates the revision and does not change the conclusion.
- Table 9's “signature not prognostic”/“immune signature not prognostic” should be “no consistent association detected”/“no significant association detected,” retaining the distinction between nonsignificance and an established null.
- In Discussion, “no advantage over a rank-matched model” is clearest as “no statistically established advantage”; the estimated difference is positive but its interval includes zero.

## Mathematical audit of the complete appendices

This audit followed every displayed chain and the associated code, rather than checking only whether the earlier equations had been copied. Cross-references between A and B are acceptable for common identities, but the author's explicitly requested intermediate steps are itemized where still missing.

| Derivation | Verification |
|---|---|
| Gaussian reconstruction likelihoods, independent priors, mean-field family and partial-likelihood objective | A.1/B.1 are correct, subject to the random zero-norm target in V2. Cox partial likelihood is correctly distinguished from a full survival-time likelihood. |
| Breslow ties, score, coupled Hessian, replacement curvature and floor | A.2's derivatives and square completion are correct; the score is retained before the weight floor. B.2 correctly omits that floor and uses h. |
| Coordinate variational argument and normal-means bridge | A.3 correctly gives log q*=Bθ−Aθ²/2+log g+constant, completes the square and obtains x=B/A, s=A^−1/2. Shared-prior marginal likelihood and the KL identity use the same prior as the posterior. Add the distinct B zero-information boundary in V5. |
| Point-exponential slab and spike | A.4 correctly sets m=x−λs², completes the half-line square, integrates the positive half-line to λ exp(−λx+λ²s²/2)Φ(m/s), normalizes against the spike and obtains the posterior weights. |
| Truncated-normal moments | A.4 correctly uses a=m/s and r=φ(a)/Φ(a), conditional first/second standard-normal integrals, M_+=m+sr and S_+=m²+s²+msr. Spike-weighted unconditional moments are correct. Tail approximations are identified as numerical evaluation rules. |
| Fixed-projection normal coefficient update | A.5 correctly expands the residual over other coordinates, collects α-weighted A/B and restricts α cancellation to positive α with an inactive floor. The normal marginal derivative, constrained prior-variance MLE (x²−s²)_+, posterior moments and equality point-mass boundary are correct. Complete the explicit Gaussian integration/square in V5; qualify random projections as a working-quadratic result in V4. |
| Single-block activities | A.6 correctly distinguishes random residual from its mean and uses E(F²), not F̄², in A_L. Normal-means inputs and EB prior estimation agree with `update_L_surv_YFB.R`. |
| Single-block loadings | A.7's reconstruction expansion/A/B are correct. The fixed-normalizer survival expansion retains both d factors and feature coupling. The separately described implemented survival branch omits these terms as stated; α_F=0 recovers reconstruction-only updates. No claim of exact quadratic coordinates for varying d(F) should be added. |
| Single-block precisions and monitor | A.8 correctly expands expected squared residuals, cancels cross-factor terms and differentiates for τ; the n×10^−8 denominator floor agrees with `update_tau.R`. A.9's weighted J, negative-KL discrepancy and sweep-start displacement limitation match the fitter. A small no-fit arithmetic check of the KL helper reproduces the discrepancy. |
| Stored normalization, orientation and uncertainty | The saved D_* and zero coefficient variances are verified. Sweep 20 meets the stated stopping rule. Correct the orientation prose contradiction in V3; do not change the frozen scoring direction. |
| Multimodal intercept | B.3 correctly differentiates expected reconstruction likelihood and qualifies exact-score shift invariance versus the expected diagonal surrogate. |
| Projection/predictor moments | B.4 correctly derives Var(Z) by independent loading coordinates, E(Z²), and Var(η)=Σ_k{Eβ_k² Var(Z_ik)+Z̄_ik² Var(β_k)}. It explicitly acknowledges shared cross-patient uncertainty. |
| Shared multimodal activities | B.5's random-residual expansion, sums over both modalities, A/B and solver inputs are correct. Cached old-prior posterior versus new-prior timing is disclosed in B.6. |
| Multimodal loadings | B.6's survival expansion η=β_k(a_if+U_i)+V_i, its expected linear/quadratic coefficients, coupling term, sequential moment refresh and A/B are correct. The reconstruction half still skips the explicit expansion in V5. |
| Point-Laplace half-lines and moments | B.6 correctly handles the positive half-line using A.4 and the negative half by f=−f′, with m_−=−x−λs² and exponential constant +λx+λ²s²/2. Mixture normalization, signed first moment and unsigned second moment are correct. No sign error was found. |
| Multimodal normal coefficients and KL | B.7 correctly expands the working square, uses E(Z²), gives x/s, shares the constrained variance MLE and conjugate moments from A.5, and derives the normal KL by averaging log q−log g. Add the zero-A and zero-SD conventions. |
| Multimodal precisions | Expected residual formula is correct. Insert the derivative and correct the numerical floor as V1/V5. |
| Scale pushforwards | B.9 correctly verifies invariant reconstruction/predictor, transforms first/second moments, Z moments, exponential/Laplace rates and coefficient variance, and explains Jacobians and KL invariance. The c=1 fallback is correct in the appendix; propagate its exception to summary statements. |
| Stopping and removal | B.10's old/new denominator agrees with the fitter. B.11 correctly describes cached penalties, re-estimated intercepts/precisions, Jensen direction, the full Hessian/covariance correction and its diagonal replacement. Add K>1 and non-strict score language; apply the actual precision floor to the protected score calculation. |

### Numerical and figure verification

- **Simulation recovery/Table 3:** All displayed recovery/activity cells agree with the 540-row saved S18 output. At ranks 6/12/20 the all-shared recovery is 1.00000/.97778/1.00000; hybrid shared recovery is 1.00000/1.00000/.96667. Revised 97–100% is correct. Cohort-specific recovery and the false-activity qualifications agree.
- **Single-block external performance:** Saved-score recomputation gives concatenated C=.628481987 and .6025621; difference .0259 with the saved interval (.0002,.0498). Per-cohort table values and risk orientations agree with the earlier independently verified values. Principal concatenation and older inverse-variance pooling remain correctly distinguished.
- **Other changed/moved numerical material:** The new rank-12 simulation rows, cohort-specific-coefficient medians/counts, all six multimodal simulation rows/convergence counts (162/162 versus 27/162), matched-data C/CI rows and starting-rank classification/sensitivity agree with saved outputs. Table 10's modality shares and 14 HR/p-value cells, subtype correlations, Program 6 Kruskal–Wallis p, enrichment counts/p-values, signature ranges/significance counts and EBMF correspondence also agree. No changed numerical conclusion was found beyond the small remaining difference-rounding issue.
- **S63:** Independently replayed the figure computations without saving plots or creating loader temporary files. Reconstructed frozen η agrees with saved scores within 10^−8; saved/reconstructed times and statuses agree; totals are 616 patients/414 deaths. Risk-tertile log-rank p-values are .0004555969, .1705394, .07464884, .0208022 and .000007657 for Dijk, Moffitt, PACA-AU microarray, PACA-AU RNA-seq and Puleo. Pooled cohort-stratified program labels reproduce Program 7 p<.0001 and Program 3 p=.0041. All 18 paired simulation means/SEs reproduce (270 pairs). Shares are copied from the saved program table; its rounded methylation shares sum to .999, and the plot faithfully preserves them.
- **Figure formats and legibility:** All five PNGs are 1950 pixels wide with 11811 pixels/metre metadata (approximately 300 dpi). Their individual PDFs are 6.5 inches wide, one page each; `pdfimages -list` finds no raster objects in them or the 42-page report. The redrawn plots are sharper, have readable black labels at their intended size, retain color in graphical elements, and remove old internal/date-based titles. Captions need the local definitions listed above, but plot values and styles do not require another numerical redesign.
- **No refitting:** Script and call-path inspection found no factor-model or Cox-regression estimation in the two figure producers. Their Kaplan–Meier/log-rank calculations regenerate the permitted descriptive plots. The local raw-expression dependency is real and documented, as V11 explains.

## Length recommendations, least information lost first

**Opinion:** Keep the readable vector figures and the complete derivation appendices. The 27-page main report is too long relative to the agreed 15–20-page target, but the first remedy should be page composition and redundant material, not smaller fonts or cutting scientific qualifications. Excluding the title/contents/Summary still leaves 26 pages, so the discrepancy is not only counting convention. Estimates below are rough and overlap; only a rerender can establish the resulting length.

| Rank | Specific change | Information cost | Likely space effect |
|---|---|---|---|
| **1** | Replace blanket H figure pinning and oversized proportional `needspace` requests with controlled floats/keep-together blocks. Fit algorithms together and fix duplicate longtable heads. Let text fill space around sensibly placed exhibits. Review large blank areas on pp.11–12,14,17,19,21 and27. | None, if every caption/row remains intact and nearby text still identifies the exhibit. | Potentially several pages; verify rather than promising an exact saving. |
| **2** | Delete Table 9 as a redundant second account of the four single-modality programs, or replace it with a very compact program/class/coefficient summary. Keep the full biological interpretations in §4.1; the coefficients already appear in the risk-orientation/direction paragraphs. | Very low; most table evidence repeats the adjacent program paragraphs. | Roughly half to one page, including its layout overhead. |
| **3** | Remove Table 1's duplicate model summary or reduce it to four genuinely contrasting rows: data/reconstruction, loading prior, survival projection/uncertainty, rank handling. Notation and detailed Methods already define the same material. | Very low if the shared priors/precision facts stay in Methods. | Roughly half to one page. |
| **4** | Cut the separate cohort-specific-coefficient simulation paragraph in §3.1 to one sentence (or omit this supporting variant), and reduce threshold-grid sensitivity in §3.3 to “classification was stable over the examined practical thresholds,” retaining the threshold values in a concise note if needed. | Low; peripheral to the evaluated recommended configuration. Keep the uncertainty/convergence warning if the coefficient result is retained. | Several paragraphs, with additional savings from reflow. |
| **5** | Shorten Table 3 to ranks 6 and20 per scenario, retaining the rank12 exceptions in its note/text. Keep Figure1, because it shows concordance differences while Table3 shows distinct recovery/false-activity quantities. | Low; removes intermediate repetitions, not an entire endpoint. Do not drop recovery/false activity merely because the C-difference figure exists. | About a third to half a page. |
| **6** | Retain only the joint and expression-only fits of the same framework as the principal rows of Table7; summarize the less directly comparable base/isolated/two-step arms in a short supporting sentence or compact panel if still useful. | Moderate: less comparator detail, but better focus and clearer disclosure of the different supporting setups. | About a third to half a page, depending on caption length. |
| **7** | In §4.2 retain provisional labels/gene evidence in Table10; shorten prose repeating its HRs, correlations and shares to the two substantive findings (basal/classical survival axis and methylation-dominant structure). In §4.1 avoid repeating the direction/conditional-effect explanation twice. | Low to moderate; keep the no-enrichment/CpG and multiplicity caveats. | Roughly half a page. |
| **8** | Compress repeated algorithm/Methods exposition: retain scientific model/predictor and numbered update order in the main text; put long A/B formulas and detailed bookkeeping in A/B, with a brief approximation disclosure in Methods. Remove redundant appendix reconciliation tables if they repeat the now-complete derivation narrative. | Moderate; preserve the two algorithms, the evaluated α_F=0 meaning, frozen normalizer and approximate-score limits. | A further half to one page of main text; appendix-only cuts do not count toward the main target. |

Apply ranks 1–5 first, rerender and recount, then choose ranks6–8 only if needed. Keep 11-point body text, one-inch margins, readable figure labels, cohort/sample definitions, direct-versus-indirect supervision distinction, convergence/rank sensitivity, inferential limitations, and all required derivation steps. Do not delete missing-evidence qualifications to meet the page count. Caption abbreviation notes can be short, but must remain local under the author's requirement.

## Scope, rendered checks, and limits

- **Scope verified:** The complete `88cafda..17a8bb1` changed-file list contains only `paper/writeup/**`, `README.md`, `PROJECT_STATUS.qmd` and `CLAUDE.md`. No model code, defaults, raw data or saved analytical result outside writeup changed. `AGENTS.md` remains its symlink. The new review file in that range is the earlier authorized independent review, not an analytical change.
- **Rendered PDF:** Inspected all 42 pages via rendered images/contact sheets, with larger inspections of the risk-curve figure and the two duplicate-header tables. No boxed algorithms or References/Appendix C remain. Algorithms do split; table data rows fit per page, but Table9/Table11 have repeated heads before captions. Figures are vector and legible. No report was re-rendered or changed during verification.
- **Rendered HTML:** Opened the revised HTML in the in-app browser through a loopback preview exposing only that one self-contained file. Five embedded PNGs loaded, 848 MathJax containers rendered, and no math-error nodes were found. At the normal 1265-pixel document width, scroll width equaled viewport width. Static checks found no duplicate IDs or missing fragment targets. Computed colors for main prose, headings, cells, captions and links were black; no callout algorithm boxes remain. Searches of rendered main text found no prohibited source IDs, paths, function/argument names, internal labels, drafting/agent language or process references. A phone-width responsive inspection was not performed; the normal-width check does not establish mobile behavior.
- **Markers:** Main source has all R1–R9 IDs; all START/END pairs are balanced/nested across the three sources. Current figure/table definitions and captions are within markers. Only the method-dependent gaps in V8 were identified; generic identities are intentionally distinguished from update-specific statements.
- **Quantitative verification limits:** Bootstrap producer settings and saved intervals were checked, but full bootstrap distributions were not independently regenerated. No donor-overlap audit, new enrichment, new factor fitting or new scientific experiment was run. Numerical replay of figure computations checked existing results, not a new comparison or model.
- **Historical/code limits:** Saved fit moments, normalizers and stopping were inspected, but its historical executable/environment was not reconstructed. Effects of the acknowledged monitoring/KL issues on fitted stopping times or rank preferences remain untested, as the revision appropriately discloses. No user decision about changing those algorithms is implied by this document.

**Next action:** Andrew may authorize committing this verification. Claude should then verify and address V1–V13 and the incomplete A/B items, correct the dispositions to match the actual result, and rerender both formats. Choose the low-loss length/layout changes first. Recheck the resulting PDF, equations and captions before requesting approval to merge. This document does not grant merge or push authorization.
