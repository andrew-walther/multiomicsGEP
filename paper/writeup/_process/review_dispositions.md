# Review dispositions

Revision of the writeup on `writeup/draft` (2026-10-14) in response to `review.md` and the author feedback in `implementation_plan.md`. Each finding was checked against the code, saved outputs or the draft before it was applied. No model code, defaults, data or saved analytical results were changed, and no model was fitted. The figures were redrawn from the saved fit, saved scores and saved aggregates.

## Part A: scientific and numerical findings

| Finding | Disposition | Verification and change |
|---|---|---|
| A1 Joint model versus plug-in algorithm | **Fixed** | Methods §1.2 now states the joint target with random $F$ and $d_k(F)$, separately from the evaluated algorithm, which freezes $D^{(r)}$ and $\tilde Z^{(r)}$ within a sweep. Appendix A.1 defines the partial-likelihood target (@eq-a-target) and states that the Cox partial likelihood is not a full likelihood. A.7 derives the survival loading term only under a fixed normalizer, and A.11 compares target and algorithm. |
| A2 Working curvature is not the Hessian diagonal | **Fixed** (confirmed) | `code/fit_cox_on_yf.R:115–118`: `w = theta*H` floored at 1e-6, `u = status - theta*H`. A.2 derives the exact score and negative Hessian, shows $a_i = \sum_a d_aP_{ai}$ and $\mathcal H_{ii} = a_i - \sum_a d_aP_{ai}^2$, and calls the replacement a Poisson-style working curvature. It also states the floor and its effect on $u$. Methods §1.2 and B.2 (no floor, division-free $h$) match. |
| A3 Monitoring statistic KL bookkeeping | **Fixed in the text; implementation defect deferred to method work** | `code/compute_elbo.R:80–81` uses `0.5*sum(log(A/(2*pi)))`. For activities, `A` is a scalar (`update_L_surv_YFB.R:78–79`) and `x` has length $n$, so the normalizer is counted once rather than $n$ times. For loadings `A` is a $p$-vector, so the count is correct there. A.9 gives the discrepancy formula and the start-of-sweep survival limitation. Table 5 and §3.3 call it the "monitoring statistic", and "ELBO" no longer labels it anywhere. Effects on stopping and rank are listed as unassessed in Limitations. The code is unchanged. |
| A4 Stored prediction normalizer; feature-count claim | **Fixed** (confirmed) | Saved fit: stored `EF_norms` differ from the final-$\bar F$ norms by 2.9e-4, −1.6e-3, 5.7e-4 and 1.4e-4 for programs 3, 5, 6, 7 (computed from S58). $D_*$ is defined and used in §1.2, Algorithm 1 and A.10. The ridge is stated. The "does not grow with $p$" claim is replaced by the identifiability rationale. |
| A5 Cached KL penalties in removal | **Fixed in the text; implementation defect deferred to method work** | `multimodal_yfb_updates.R:359–367, 491–502`: the prior is refitted and the KL returned uses the old-prior `log_ml`, for both $L$ and $F$. §1.5, B.6 and B.11 describe cached penalties. "Conservative / evidence does not fall" is replaced by "improves the implemented removal score; no guarantee for exact evidence or the current-prior ELBO". |
| A6 Jensen variance claim; centering invariance | **Fixed** | §1.6 and B.11 now say the plug-in term is an upper approximation. B.11 adds the common-shift counterexample, the full $-\frac12\mathrm{tr}\{\mathcal H\,\mathrm{Cov}_q(\eta)\}$ correction, and how the diagonal penalty simplifies it. B.3 states that the expected surrogate need not be shift invariant and that the reconstruction-only intercept update is an implementation choice. |
| A7 "Every shared program" recovery | **Fixed** (confirmed) | Recomputed from S18: shared recovery at ranks 6/12/20 is 1.00/0.98/1.00 (all shared) and 1.00/1.00/0.97 (hybrid). The text says 97%–100%, and the rank-12 rows are added to Table 3. |
| A8 Appendix derivations incomplete | **Fixed** | Both appendices were rewritten from the likelihoods. They add Gaussian and Breslow likelihoods, priors, the mean-field target, random-residual expansions, the steps to each $A$, $B$, $x$ and $s$, normal-prior marginal likelihood and posterior, point-exponential and point-Laplace half-line completions, truncated-normal moments, precision MLE, predictor-moment identities, scale pushforward, and stopping and removal. See B14 below. |
| A9 Rounding and stopping denominator | **Fixed** (confirmed) | Concatenated C is now reported as 0.6285 versus 0.6026 (frozen recomputation in S11/S56). B.10 uses $\max(1, \max|X_{\rm old}|, \max|X_{\rm new}|)$ per `fit_multimodal_yfb.R:64–69`. |
| A10 Scale and feedback boundaries | **Fixed** (confirmed) | B.6: survival terms vanish only when both $\bar\beta_k$ and $\overline{\beta_k^2}$ are small. B.9: full moment, $\bar Z$, $\mathrm{Var}(Z)$ and $\sigma^2$ maps, pushforward/Jacobian explanation, KL invariance, and the $c_k = 1$ fallback (`fit_multimodal_yfb.R:197`). §1.2: survival does not enter the $L$/$F$ updates but can affect the returned iterate through the stopping rule. |
| A11 Revisit coverage | **Fixed** | Notation is wrapped in R2/R5. The future-method paragraph moved into Continuing work (R9/R1). The appendix introductions are wrapped in R1/R2 (A) and R3 (B), and B.11 removal in R4. Markers are balanced in all three sources (scripted check). |
| A12 Uniformly ordered survival curves | **Fixed** | §3.2 now reads "generally favored the low-risk tertile; log-rank tests indicated differences in three of five cohorts", with the exact p-values. |
| A13 Recovery failure attributed to both models | **Fixed** | The Discussion says only the multimodal framework was tested on the low-variance program. The single-modality configuration does not test whether supervision helps, and Continuing work states the goal without an impossibility claim. |

**Drafter statements checked by the review.** All seven were confirmed or qualified by the review, and the report follows them. The established-model C = 0.685 on the matched donors is omitted entirely. The risk-scale signs of Programs 3 and 7 are stated without the old disagreement narrative. The zero stored variances are stated as a fit limitation.

## Part B: author revision requirements

| Item | Disposition | Change |
|---|---|---|
| B1 Title/subtitle | **Fixed** | Title set to the working title. The subtitle and the bibliography/link-citation YAML are removed. |
| B2 Authored appearance | **Fixed** | Removed the running header, status box, long summary box, and all drafting, reviewer, inventory and "could not trace" language. |
| B3 Repository references and labels | **Fixed** | No source IDs, paths, function or argument names, configuration labels, D4/P3/flash k/expr-only labels or output-directory names remain (grep check). Software appears only as R packages (ebnm, flashier, glmnet, fgsea). Table 10 uses defined reference codes SM/U/EO. |
| B4 Dated and process references | **Fixed** | No month names appear except the title-page date. The base model is described by its components, and the orientation is stated directly. The per-program decomposition paragraph is dropped, and the objective imbalance is kept as a limitation. |
| B5 Noise labels | **Fixed** | Table 6 and Figure 3 now use "Lower (0.25×)" and "Higher (1×)" multipliers of the baseline noise SD. Baseline precisions ~ Gamma(2, 2) are defined in §2.2 and the Table 6 caption (`simulate_multimodal_yfb.R:106,181`). |
| B6 Self-contained exhibits | **Fixed** | Every caption defines its population, metric and direction, uncertainty method and abbreviations: ELBO, BIC, CV, PL, PVE, HR, SD, CI, C, NES, BH, ORA, EBMF, LASSO, PDAC, TCGA, ICGC, CPTAC, PurIST, EMT, ECM, GO-BP, RTK, PI3K, ERK, IPMN, CpG and 450K. Table 5 (rank sweep) now has preferred-direction arrows and definitions. |
| B7 Algorithms | **Fixed** | Callout boxes are replaced by plain numbered, indented lists. Dataset-specific settings are moved to the sentence after each algorithm, and function names are removed. |
| B8 Figures | **Fixed** | All five figures are redrawn at their printed size (6.5 in, 9-pt base text). PDF output uses vector graphics (`default-image-extension: pdf`; `pdfimages` lists no raster images) and HTML uses 300-dpi PNG. Figures are pinned in place (`fig-pos: H`). In-plot titles are removed, colours are kept, and text is black. |
| B9 Drop the evidence-gap table | **Fixed** | The table is removed. EBMF and two-step material is in Results and Discussion, DeSurv in Biology and Limitations, rank selection in Results and Continuing work, multimodal added value in Results, Discussion and Limitations, and biology in Evaluation §2.5 and Limitations. |
| B10 Discussion | **Fixed** | §5 now has three subsections: Interpretation and implications, Limitations, and Continuing work. The α_F = 0 statement is kept in Continuing work. |
| B11 Remove Appendix C | **Fixed** | The characterization table is now Evaluation §2.5 prose, with caveats. The simulation figure and table move to §3.1. The table drops its concordance columns, which the figure covers, and adds rank 12. The completed-versus-proposed table is folded into §4.2 and Continuing work. |
| B12 No References | **Fixed** | Changed to "DeSurv (Young et al., submitted)" in the text, with no bibliography. The project `.bib` is unchanged. |
| B13 Length/organization | **Partly addressed; see unresolved** | Order is Summary (one paragraph), Methods, Evaluation, Results, Biology, Discussion, Appendices A–B. Redundant summary, gap table, workflow table, Appendix C and process prose are cut. The malformed "-mm-sim" text is removed. No table or algorithm splits across pages. |
| B14 Gap-free derivations | **Fixed** | Every step listed in the review's B14 is included (A.1–A.11, B.1–B.13). |

## Unresolved or deferred

1. **Length (B13).** The PDF main text is 27 pages (title and contents on page 1; Discussion ends on page 27), against the 15–20-page target. The extra length comes from full-width figures sized for legibility, relocated Appendix C material, and white space reserved so tables do not split. Further cuts need Andrew's judgement about content. Candidates are Table 3 (the figure carries the C difference), the cohort-specific coefficient paragraph in §3.1, the threshold-sensitivity paragraph in §3.3, and the "Two-step, then Cox" rows of Table 7.
2. **Implementation defects (A3, A5).** The activity-KL normalizer in `compute_ebnm_kl` and the cached old-prior KL in the multimodal updates are documented in the report but not fixed. Fixing them changes the stopping and removal behaviour of saved fits, which belongs to the method work, not to this writing branch.
3. **Not independently regenerated:** full bootstrap distributions, the historical fitting environment of the evaluated single-modality fit, and the cross-cohort donor overlap beyond the matched ICGC subset. These are the same limits the review recorded.
4. **HTML runtime check.** HTML was inspected structurally (5 embedded figures, resolved cross-references). It was not opened in a browser at mobile width.

## Verification round (`verification.md`)

Codex verified the revision at `17a8bb1` and found that several dispositions above overstated completion: A8/B14 (derivations), B6 (captions), B7 (algorithm pagination) and A11 (markers). Each verification finding was checked against code or the rendered PDF before it was applied. Required corrections are kept separate from optional wording and length suggestions.

### Required corrections

| Finding | Disposition | Evidence and change |
|---|---|---|
| V1 Multimodal precision floor | **Fixed** (confirmed) | `fit_multimodal_yfb.R:508` calls `multimodal_yfb_tau_from_residuals()`, which uses `pmax(r2, .Machine$double.eps)` (`:640`). The error branch is in the unused helper `multimodal_yfb_updates.R:577`. B.8 now separates the mathematical boundary ($Q = 0$ makes the likelihood unbounded) from the protected computation $n/\max(Q, \epsilon_{\rm mach})$, used in fitting and in the removal score. |
| V2 Zero column in the normalized target | **Fixed** | A.1 defines $u_k(F) = F_{\cdot k}/d_k(F)$ for $d_k > 0$ and $u_k = 0$ for a zero column. It states that this convention defines the target, does not change the evaluated fits (whose stored norms include a ridge), and that the Gaussian updates additionally hold $d_k$ fixed. Methods §1.2 carries the same convention. |
| V3 Contradictory orientation account | **Fixed** (confirmed) | Current code flips when the risk-oriented `reverse = TRUE` training C is below 0.5 (`fit_cox_on_yf.R:891–904`). The saved fit was produced under a different convention. A.10 now uses Codex's replacement and states only the fixed convention of the evaluated fit ($\mathrm{risk} = -\eta$) and the absence of recoverable coefficient variance. |
| V4 Second-moment statement attributed to the exact target | **Fixed** | A.5 says the $E(\tilde Z^2)$ coefficient belongs to the random-projection working quadratic, and that the exact target is nonquadratic and nonconjugate. The Table 10 (A.11) cell is reworded the same way. |
| V5 Remaining derivation steps | **Fixed** | (1) A.5 now has the full Gaussian convolution and square completion, giving $N(x;b,s^2)N(b;0,v) = N(x;0,s^2+v)N(b;m,V)$, with each coefficient matched. B.7 references it. (2) B.6 shows the reconstruction-loading expansion explicitly. (3) B.8 shows the precision derivative. (4) The prior-only boundary ($A = 0$, hence $B = 0$) is stated for B.5–B.7. It matches `multimodal_yfb_updates.R:39–49` (prior moments returned) and `:524–527` (prior-only $\beta$), including the error raised when $A = 0$ and $B \ne 0$. A.3's floor is now labelled as single-modality only. |
| V6 Algorithm and table pagination | **Fixed** (checked in the PDF) | Each algorithm is wrapped in a PDF-only `minipage`, so it is unbreakable; the HTML keeps plain lists. Algorithm 1 is wholly on p.4 and Algorithm 2 on p.7, each with its inputs and outputs. The duplicate heads had one cause: a longtable started too low on a page, broke before its first row, and repeated `\endhead` on the next page. Table heights were measured from the rendered PDF (Table 7 ≈ 0.54, Table 8 ≈ 0.49, Table 9 ≈ 0.74 of the text height). Sized `\needspace` reserves are now used only where needed. The three figures that could float above those tables (Figures 3–5) are pinned in place, and the orphaned A.11 heading gets a reserve too. A scripted check confirms every table's caption and last row are on the same page, with no repeated head. |
| V7 Caption self-containment | **Fixed** | All items on the checklist are covered locally. Table 1 defines expr/meth, $n$, $p$, $K$, $L$, $F$, $D_*$ and SD. Tables 2, 4, 6, 7, 8 and 9 and Figures 1–5 define TCGA(-PAAD), CPTAC, RNA-seq, PDAC, ICGC, PACA-AU, CI, LASSO, CV, PL, SD, CpG and DeSurv as they occur. Table 5 gives $n = 273$ and $p = 2{,}064$. Table 3 defines $\bar{\tilde\beta}_k$ and the model arm. Tables 10–11 define KL, sd and their symbols. The earlier claim that every caption defined its abbreviations was not true at `17a8bb1`; it is true now. |
| V8 Revisit gaps | **Fixed** | Added an R1 block around A.3's single-modality floor and an R1 block around A.4 (point-exponential posterior). All three sources are balanced (scripted check). |

### Smaller corrections (P3)

| Finding | Disposition | Change |
|---|---|---|
| V9 Normalization exceptions | **Fixed** | Table 1 uses $D_*$ and states the constant-column exception. The predictor paragraph now says "loading columns, each divided by a norm". Figure 4 now says "loading column divided by its stored training norm". The multimodal scale paragraph, Algorithm 2 step iv, B.9 (with threshold $\sqrt{\epsilon_{\rm mach}}$, `fit_multimodal_yfb.R:197`) and B.13 state that numerically constant columns are not rescaled. |
| V10 Classification and removal boundaries | **Fixed** (confirmed) | Methods, B.7 and Table 9 use $|\bar\beta_k|/\max\{\mathrm{sd}(\beta_k), 10^{-12}\}$ (`fit_multimodal_yfb.R:553–556`): a point mass is inactive. Removal runs only while $K > 1$ (`:591`), as stated in Methods, Algorithm 2, B.11 and B.13. B.11 now says a deletion "does not decrease" the score. |
| V11 Figure provenance | **Not changed; disposition corrected** | Confirmed: `make_writeup_figures.R` reads local held-out expression to recompute frozen program projections, and no saved projection file exists. The earlier wording "saved fits, scores and aggregates only" is corrected here. The dependency is stated in the script header, in the README command (`PDAC_DATA_ROOT=...`) and in S63. No patient-level cache was created or committed. No model is fitted. |
| V12 Single-block comparator settings | **Fixed** (confirmed) | `run_multimodal_yfb_simulation.R:139–149` caps the fit at 100 sweeps and uses unstandardized input. It is now labelled "single-block fit … to unstandardized log2 TCGA expression (at most 100 sweeps)" in Comparators and in Table 7, where the row is "Single-block model". The row is kept. |
| V13 Wording and rounding | **Fixed** | Summary: "mean ranks close to the true value at lower noise, with poorer recovery when the prognostic program has low variance". Appendix links now point to A.9. The matched-donor difference is "about 0.01" (saved 0.0095). Discussion: "no statistically established advantage". The "not prognostic" table wording went with the deleted table. |

### Additional correction found during this round

The Discussion said the held-out cohorts were "measured on platforms not used in training". Two of the five are RNA-seq, which is a training platform. It now says "three of them on a microarray platform not used in training".

### Length and layout

Following Codex's ranking, least loss first:

1. **Layout (rank 1).** Removed blanket `fig-pos: H` and the oversized blanket `needspace` reserves. Relaxed the float fractions. Gave only the tables that need one a measured reserve. Moved explanatory paragraphs ahead of large tables so text fills the space they leave.
2. **Rank 2.** Deleted the single-modality program table; its content is in the §4.1 paragraphs.
3. **Rank 3.** Reduced Table 1 to four contrasting rows.
4. **Rank 4.** Cut the cohort-specific coefficient paragraph and the threshold-sensitivity paragraph to one sentence each.
5. **Rank 5.** Table 3 now shows ranks 6 and 20, with the rank-4 and rank-12 values in its note.
6. **Rank 7.** Removed repetition of Table 9's numbers in §4.2 and the duplicated direction explanation in §4.1.
7. **Other.** Shortened Table 8's row labels. Tightened Limitations and Continuing work without dropping a limitation or an aim. Reduced two figure heights slightly (program curves 3.1 → 2.9 in, modality shares 3.0 → 2.7 in) at the same 9-pt plot text.

Not applied: ranks 6 (dropping comparator rows from Table 7) and 8 (moving Methods exposition into the appendices). Both lose moderate information and need Andrew's decision.

**Result:** the PDF is 39 pages. The main text is pages 1–24 (title, contents and Summary on page 1; Methods through Discussion on pages 2–24), down from 27. Appendix A begins on page 25 and Appendix B on page 33. The 15–20-page target is still not met.

### Still unresolved

1. Main-text length is 24 pages against a 15–20 target. The remaining options (Codex ranks 6 and 8, or a shorter Methods) need author judgement.
2. Implementation defects A3 and A5 remain documented for the method work and are not fixed.
3. Figure regeneration needs the local held-out expression data (V11).
4. The HTML has not been checked at phone width. Codex confirmed MathJax renders at normal width.

## Second verification (`verification_round2.md`)

Codex verified the scientific corrections, added algebra, local caption definitions
and PDF pagination at `fbf8317`. Andrew authorized the remaining minor updates.

| Finding | Disposition | Change |
|---|---|---|
| N1 First sweep-normalizer formulas omit the ridge | **Applied** | Methods §1.2 and Appendix A.1 now explicitly use the square root of the preceding posterior-mean loading norm squared plus `10^-10`, matching `code/fit_cox_on_yf.R:643–646` and A.10. The random target's separate zero-direction convention is unchanged. |
| N2 Wrong disposition-file pointer | **Applied** | `PROJECT_STATUS.qmd` now points to this file instead of attributing the dispositions to `verification.md`. |

The 24-page main text and V11's local-expression dependency remain disclosed.
No further length cuts, figure regeneration, model changes or new analyses are
part of these corrections. Final approval remains pending.

Verification after correction: HTML and PDF renders both exited successfully.
Both corrected formulas appear in HTML; PDF pages 3 and 25 were inspected as
images. Only those two pages' extracted text changed. The PDF remains 39 pages,
with the main text ending on p.24, algorithms on pp.4/7 and every table on its
previous single page. Revisit markers remain balanced and nested.

## Final editorial check (`final_editorial_check.md`)

Andrew asked for the required corrections E1–E5 and Methods cuts M1–M4. Suggestions O1–O4 and cut M5 were not applied.

| Item | Disposition | Change |
|---|---|---|
| E1 Table 7 caption versus row label | **Fixed** | The caption now says "the single-block model is a fit of the model of @sec-single-model…". |
| E2 Base model missing from Comparators | **Fixed** | §2.3 now lists four further arms, starting with the base multimodal model (no feature intercepts or factor removal; nonnegative loadings; moment-based prior updates). |
| E3 Unsupported claim in Continuing work | **Fixed** | Item 1 now rests on the reported multimodal result (survival enters the loading update, but the planted low-variance program was not recovered) and gives the mechanisms as plausible reasons, not findings. |
| E4 Uncited exhibits | **Fixed** | Table 1 is cited in Notation, Table 10 in A.11, and Table 11 in B.13. A scripted check finds no uncited table or figure and no dangling reference. |
| E5 Abbreviations at first prose use | **Fixed** | §2.1 defines ICGC and PACA-AU in its first sentence and CpG at its first use. |
| M1 Algorithms out of the main text | **Applied** | Algorithms 1 and 2, with their Revisit markers, now open Appendices A and B under "Fitting algorithm". Each still fits on one PDF page (p.23 and p.32). Methods keeps a short "Fitting" paragraph per model with the update order, the evaluated settings and a link. Simulation and Table 6 references read "Algorithm 2, Appendix B". |
| M2 Single-modality predictor and Cox approximation | **Applied** | The sweep-fixed normalizer, ridge and zero-column details are left to A.1/A.10. The `@eq-cox-working` display is removed; Algorithm 1 now cites @eq-a-working. Methods keeps @eq-zf, the fixed-covariate and stored-$D_*$ statements, and a statement that the Poisson-style diagonal curvature is not the second-order Taylor expansion. |
| M3 Multimodal removal | **Applied** | The @eq-mm-elbo display moved to B.11, under a new anchor `sec-appendix-b-removal`; references from §3.5 and B.11 still resolve. Methods describes the score in words and keeps every qualification: plug-in survival term, cached KL penalties against earlier priors, does-not-decrease rule, unconverged fits, $K > 1$, and no evidence or ELBO guarantee. |
| M4 Approximations subsection | **Applied** | §1.6 is removed. Limitations now says neither algorithm is exact coordinate ascent on one ELBO or guaranteed to increase its monitored quantity. The Jensen and variance-penalty discussion remains in B.11. |

**Layout.**
- Figure 5 now follows the paragraph that cites it, and the "evidence therefore supports…" paragraph moved below Table 9. Both moves remove page-sized gaps around the large multimodal programs table.
- A.11 starts on a fresh page, because the earlier space reserve printed a stray table header.
- No table splits across pages or repeats its header.
- To fit the page, Continuing work items 1, 3 and 4 were shortened. Each still states its aim, and item 3 keeps "evaluated in replicated simulations".

**Result:** the PDF is 39 pages. The main text runs from p.1 to p.22 (from p.24 before these cuts). Appendix A starts on p.23 and Appendix B on p.32.

**Still open:** the main text is still above the 15–20-page target. Further reductions would need M5, condensing §2.2 Simulation designs, or O1–O4, all at Andrew's discretion.

## Third verification correction (`verification_round3.md`)

Andrew authorized correction R3-1. Continuing work item 1 now limits the failure to recover the planted low-variance prognostic program to the lower-noise multimodal simulation with strong reconstruction programs. This agrees with Section 3.4 and Table 6 and does not deny the higher-noise recovery result. Both renders were refreshed successfully; only PDF p.22's text changed, its layout was visually inspected, and the 39-total/22-main-text pagination is unchanged. No model code, defaults, data, figure assets or saved results changed. Required corrections are resolved; interim approval and acceptance of the length remain pending.
