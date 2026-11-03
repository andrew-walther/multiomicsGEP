# Talk outline: Survival-supervised Bayesian matrix factorization (prelim Part III)

Part III of the prelim oral exam (2026-12-07). The exam runs about 90 minutes with
questions, and each of the three parts gets 20–25 minutes. The committee is four
statisticians and one non-statistician. Draft 1 (2026-11-03), written before the deck
so the content can be reviewed first.

Part III is the one unfinished project, so the talk is a research talk about a proposed
project with its progress inside it. Each section gives what exists first, then what is
planned. Every slide is tagged **prelim** (main deck), **technical** (backup) or
**talk-only** (dropped when the deck joins the full prelim deck).

The slide entries below were written by hand for review. Once the deck
(`ssbmf_talk.qmd`) exists, they are regenerated from it with
`python3 presentation/prelim-slides/tools/outline_from_deck.py`, and `outline.html` with
`python3 presentation/prelim-slides/tools/build_outline_html.py`.

## Structure and timing

Timings assume about 140 spoken words per minute.

| Section | Slides | Minutes |
|---|---|---|
| Title and outline | 1–2 | about 0.5 |
| 1. Introduction | 3–6 | about 3.5 |
| 2. Methods | 7–13 | about 7.5 |
| 3. Results (simulation and PDAC application) | 14–23 | about 9.25 |
| 4. Discussion | 24–26 | about 2.75 |
| Main body | 1–26 | about 23.5 |
| Backup | B1–B14 | as needed |

## Author decisions (2026-11-03)

- Conventional order: Introduction, Methods, Results with the PDAC application,
  Discussion. Planned work sits inside each section after what exists, and the wording
  marks the difference, without "current" or "planned" labels.
- The main-deck results all come from the frozen October 14 fit, in which survival sets
  the coefficients but does not shape the gene loadings. Its programs are described as
  expression programs that survival weights, never as survival-guided discoveries.
- Every result on fits where survival shapes the loadings (October 15, October 30,
  November 3) stays in backup until Andrew has verified it. Those slides move to the main
  body once the supervised model is validated.
- The extension to several data types is presented as proposed work only: a model slide,
  its expected value on the aims slide, and future work. No multimodal results appear
  anywhere, since the model has not yet converged on real data.
- Multi-omic supervised factorization is presented as an open gap in the field, without
  singling out DeSurv.
- The proposed model uses a point-normal (spike-and-slab) prior on the survival
  coefficients. The fitted model used a normal prior, and the slide says so.
- The fitted model standardized each gene to mean 0 and SD 1 while its loadings are
  nonnegative. The methods slide explains this and the Limitations slide lists it.
- Plain language throughout, including backup: no internal labels or project shorthand,
  no question-and-answer phrasing. Titles state what the slide shows.
- Remaining work has no dates. It is work between the prelim and the final exam.
- Content cutoff 2026-11-03. Slide content freezes on 2026-11-30 for practice.

## Open items

- The four diagrams marked "draft" are progress-book diagrams converted for review. They
  carry in-figure titles and internal labels, so they will be redrawn as plain PNGs by a
  deck script before the deck is built.
- The iBAG description on slide 5 (it relates methylation-driven and other expression
  components of each gene to survival, gene by gene) should be checked against Wang et
  al. (2013) in the literature review.
- The number of unique donors among the 616 external records (568 in the tracker) should
  be confirmed before it goes on the Limitations slide.
- The title slide lists Andrew only. Dr. Liu, Dr. Rashid and the Rashid Lab are
  credited on the shared prelim thank-you slide.

## Sources

- Model and fitted configuration: `paper/writeup/writeup.qmd` (single-modality model
  section and Appendix A), `docs/progress_book/current-state.qmd`,
  `derivations/multimodal_YFB/multimodal_YFB_derivation.qmd` (objective with all weights
  equal to one), `paper/prelim/project3-ssbmf.qmd` (multimodal model).
- Simulation: `writeup.qmd` simulation designs and Figure `single_sim_delta`;
  Chapter 4 results on starting rank and cohort attribution.
- PDAC application: `results/benchmark_snapshot_2026-10-14/`,
  `results/benchmark_sim/outputs/ebmf_cox_external/`, `writeup.qmd` external validation
  table and biological findings.
- Backup on supervised fits: `docs/findings/2026-11-03-joint-yfb-verification-findings.md`,
  `results/benchmark_sim/outputs/yfb_supervision_2026-10-30/attempt7/summary.csv`,
  `docs/findings/2026-10-30-current-method-a1.md`,
  `docs/findings/2026-10-30-yfb-targeted-diagnostics-findings.md`, `DECISIONS.md`
  (October 15 entries).
- Remaining work: `ROADMAP.md`, `docs/plans/2026-11-02-package-and-paper-pathway-plan.md`.
- Figures: `presentation/prelim-slides/tools/make_figures.R` copies scripted figures and
  draws the concordance forest plot.

<!-- SLIDES (generated) -->

## 1. Introduction

### 1. Survival-supervised Bayesian matrix factorization identifies prognostic molecular programs
- Figure: none
- Slide text: Survival-supervised Bayesian matrix factorization identifies prognostic molecular programs · Andrew J Walther · Department of Biostatistics, University of North Carolina at Chapel Hill · Prelim Part III
- Notes: My third project is a model that finds gene expression programs in pancreatic tumors and relates them to patient survival in one fit. The method is still being finished, so I'll show what is built and validated, and then what remains.
- Time: 0.25
- Carry-over: talk-only

### 2. Outline
- Figure: none
- Slide text: Introduction: pancreatic cancer, two-step analyses and the gap · Methods: the joint model, its fitting, the finalized model and the extension to several data types · Results: simulation and five external PDAC cohorts · Discussion: findings, limitations and remaining work
- Notes: I'll start with the clinical problem and how it is usually analyzed, then the model, then results in simulation and in pancreatic cancer data, and I'll end with the limitations and the work left before the final exam.
- Time: 0.25
- Carry-over: talk-only

### 3. Pancreatic cancer is lethal, and tumors differ in ways that predict survival
- Figure: none
- Slide text: Five-year relative survival is 13.7% (SEER, 2016–2022 diagnoses) · Tumor expression separates basal-like and classical subtypes, and basal-like tumors have shorter survival (Collisson 2011, Moffitt 2015, Bailey 2016) · Bulk expression also carries signal from the stroma and immune cells around the tumor · Goal: find expression programs that predict survival and can score a new patient from their expression alone
- Notes: Pancreatic ductal adenocarcinoma, PDAC, has a five-year survival of about fourteen percent. Earlier work found that tumor expression splits into a basal-like and a classical subtype, and patients with basal-like tumors die sooner. A bulk tumor sample also contains stroma and immune cells, so its expression mixes several sources. The goal of this project is to find the expression programs that matter for survival and to score a new patient with them using nothing but that patient's expression.
- Time: 0.75
- Carry-over: prelim

### 4. The usual analysis finds programs first and tests survival afterward
- Figure: draft_two-step-vs-joint.png Draft diagram, to be redrawn without the panel titles. Top: two-step fit. Bottom: joint fit.
- Slide text: Two-step: factorize expression Y ≈ LFᵀ, then fit a Cox model to the program scores · A program that explains little expression variance but predicts survival can be missed, because cell-composition variation dominates an expression-only fit · The Cox step treats estimated scores as known, so its standard errors leave out the uncertainty in the scores · Joint fit: estimate the programs and the Cox model together
- Notes: The usual pipeline has two steps. First a matrix factorization finds programs in the expression matrix, and then a Cox model tests the program scores against survival. Survival never touches the first step. That has two consequences. A program that explains only a little expression variance can be lost, even if it predicts survival well, because large sources of variation such as cell composition dominate the factorization. And the Cox step treats the estimated scores as if they were measured exactly, so its standard errors are too small. A joint fit estimates both pieces together.
- Time: 1.0
- Carry-over: prelim

### 5. Existing methods each cover part of the problem
- Figure: none
- Slide text: Survival-supervised NMF (DeSurv, Young et al.): survival guides program discovery, with penalized point estimates and rank and supervision strength tuned by cross-validation · iBAG (Wang et al. 2013): relates methylation-driven and other parts of each gene's expression to survival, one gene at a time · MOFA (Argelaguet 2018) and iCluster (Shen 2009): shared factors across data types, without survival · Gap: no method learns programs across data types with survival guiding them, posterior uncertainty, and a number of programs estimated from the data
- Notes: Several methods cover pieces of this. Survival-supervised nonnegative matrix factorization from our group lets survival guide which programs are found, but it reports point estimates and tunes its rank and the strength of supervision by cross-validation. iBAG links methylation and expression to survival, but gene by gene rather than through programs. MOFA and iCluster find factors shared across data types, but they ignore survival. What none of them does is learn programs across data types with survival guiding them, carry uncertainty through, and estimate how many programs there are.
- Time: 1.0
- Carry-over: prelim

### 6. Project aims
- Figure: none
- Slide text: Aim 1: a joint Bayesian model of expression and survival that can find low-variance prognostic programs, reports uncertainty, and shrinks unneeded programs · Aim 2: validate the model in external PDAC cohorts and characterize its programs biologically · Aim 3: extend the model to several data types, starting with expression and methylation, to link methylation changes to expression programs
- Notes: The project has three aims. The first is the joint model itself: it should find programs that predict survival even when they explain little expression variance, report uncertainty, and decide how many programs to keep. The second is to validate it in external pancreatic cancer cohorts and interpret its programs. The third is to extend it to several data types, so that a program can include both genes and methylation sites.
- Time: 0.5
- Carry-over: prelim

## 2. Methods

### 7. Each patient contributes an expression profile and a survival outcome
- Figure: none (notation glossary and a small matrix sketch to draw)
- Slide text: Y: n × p expression, patients by genes · (tᵢ, δᵢ): observed time and event indicator · Y ≈ LFᵀ with K programs · L (n × K): program activity of each patient · F (p × K): gene loadings, the definition of each program · Zᵢₖ: patient i's expression projected onto program k's normalized loading · Risk score ηᵢ = Σₖ Zᵢₖ βₖ, with βₖ the program's log hazard ratio
- Notes: Here is the notation. Y has a row for each patient and a column for each gene, and each patient also has a survival time and an event indicator. We approximate Y as L times F transpose with K programs. A column of F is a program, its weight on each gene. A column of L gives how active that program is in each patient. The risk score projects a patient's expression onto each program's loading column, scaled to unit length, and multiplies by a coefficient beta, which is a log hazard ratio. Because the score uses expression and loadings only, a new patient can be scored without fitting anything.
- Time: 1.0
- Carry-over: prelim

### 8. The joint model adds a Cox partial likelihood to a matrix factorization
- Figure: none (equation slide)
- Slide text: yᵢⱼ = Σₖ Lᵢₖ Fⱼₖ + eᵢⱼ, eᵢⱼ ~ N(0, τⱼ⁻¹) · hᵢ(t) = h₀(t) exp(ηᵢ), η = (Y꜀ F D) β, D = diag(1/‖F·ₖ‖) · Objective: 𝒥(q) = λ_G E_q[log p(Y | L, F, τ)] + γ E_q[ℓ_Cox(η)] − KL(q_L ‖ g_L) − KL(q_F ‖ g_F) − KL(q_β ‖ g_β) · λ_G and γ weight the expression and survival terms · F appears in both terms, so survival can reach the loadings
- Notes: The model has two parts. The expression part is a Gaussian factorization with a separate noise precision for each gene. The survival part is a Cox model whose linear predictor is the projection score times beta. We fit a variational approximation q by maximizing this objective: the expected expression log-likelihood, plus the expected Cox log partial likelihood, minus a Kullback-Leibler penalty for each block of parameters. The weights lambda G and gamma set how much the expression data and the survival data count. The loadings F appear in both terms, and that shared F is how survival can shape the programs.
- Time: 1.25
- Carry-over: prelim

### 9. Empirical Bayes priors let the data decide how much to shrink
- Figure: none (prior formulas with small density sketches)
- Slide text: Lᵢₖ ~ (1 − π_Lk) δ₀ + π_Lk Exp(λ_Lk) · Fⱼₖ ~ (1 − π_Fk) δ₀ + π_Fk Exp(λ_Fk) · βₖ ~ (1 − π_β) δ₀ + π_β N(0, σ²_β) (proposed, while the fitted model used βₖ ~ N(0, σ²ₖ)) · Prior parameters are estimated from the data for each program (ebnm R package) · The point mass at zero lets a whole program, or its survival effect, shrink to exactly zero
- Notes: The activities and loadings get point-exponential priors: a point mass at zero mixed with an exponential, so they are nonnegative and sparse. The survival coefficients get a point-normal prior, which lets a program have no survival effect at all. The model fitted so far used a plain normal prior on beta instead. In each case the prior's parameters are estimated from the data, program by program, which is what makes this empirical Bayes. The point masses are also how the model decides how many programs it needs.
- Time: 1.0
- Carry-over: prelim

### 10. The model is fitted by coordinate ascent, starting with too many programs
- Figure: draft_variational-updates.png Draft diagram, to be redrawn with plain labels.
- Slide text: Each sweep: approximate the Cox term by weighted least squares around the current risk score · Update βₖ, then Lₖ, then Fₖ, one program at a time · Update the noise precisions τⱼ · Start with more programs than expected and let unneeded ones shrink to zero · Fix the direction of the risk score from training data only
- Notes: We fit the model by coordinate ascent variational inference. The Cox term is not quadratic, so at the start of each sweep it is replaced by a weighted least squares approximation around the current risk score. Then each program is updated in turn, its coefficient, its activities and its loadings, and the noise precisions are updated at the end of the sweep. We start with more programs than we expect and let the priors remove the ones that are not needed. The sign of the risk score is fixed once from the training data.
- Time: 1.0
- Carry-over: prelim

### 11. A new patient's risk score needs only their expression
- Figure: draft_held-out-prediction.png Draft diagram, to be redrawn with plain labels.
- Slide text: Keep the training genes and standardize within the new cohort · Project onto the fitted loadings: Z_new = Y_new F̂ D̂ · Risk score η_new = Z_new β̂, with the training direction · No patient-level parameters are refitted, so the score is the same formula in training and validation
- Notes: To score a patient in a new cohort, we restrict to the training genes, standardize within that cohort, project onto the fitted loadings and multiply by the fitted coefficients. Nothing about the new patient is estimated. That keeps validation clean, because the same formula is used for training and test patients.
- Time: 0.75
- Carry-over: prelim

### 12. So far survival sets the coefficients but not the loadings
- Figure: none (two columns: fitted model, finalized model)
- Slide text: Fitted model: survival sets β but does not update F · Genes standardized to mean 0, SD 1 within cohort, although loadings are nonnegative · Normal prior on β, fixed-score coefficient update · Finalized model: the survival term enters the loading update · Obstacle: with about 2,000 genes, the expression term outweighs the survival term in each loading update · Plan: rebalance λ_G and γ, fit on a nonnegative scale, and pass a simulation test before use on PDAC · Fallback: a risk score built on patient activities L, with new patients scored by fitting their activities to the fixed loadings
- Notes: The model fitted so far is a step short of the full proposal. Survival sets the coefficients, but the loadings are updated from expression only, so its programs are expression programs that survival then weights. It also standardized each gene, which gives negative values, while the loadings are nonnegative. That mismatch will be removed. The finalized model lets the survival term into the loading update. The obstacle is scale. With about two thousand genes, the expression part of each loading update outweighs the survival part, so survival barely moves the loadings. The plan is to rebalance the two weights and to accept the result only if it passes the simulation test I'll show in the results. If that fails, the fallback builds the risk score on patient activities instead.
- Time: 1.5
- Carry-over: prelim

### 13. The extension to several data types shares patient activities across them
- Figure: none (block diagram of shared L and stacked F⁽¹⁾, F⁽²⁾ to draw)
- Slide text: Y⁽ᵐ⁾ ≈ 1μ⁽ᵐ⁾ᵀ + L F⁽ᵐ⁾ᵀ for data types m = 1, …, M, with one shared L · Each data type has its own loading prior family and noise precisions, for example a sparser prior for many methylation sites · η = (Σₘ Y⁽ᵐ⁾ F⁽ᵐ⁾ D⁽ᵐ⁾) β · M = 1 gives the single data type model · One program can load on genes and methylation sites together
- Notes: The extension keeps one matrix of patient activities and gives each data type its own loadings. Each data type also gets its own prior family and noise levels, which matters because there are many more methylation sites than genes and most of them should have zero weight. The risk score sums the projections over data types. With one data type it is exactly the model I just described. The point of the extension is that one program can include genes and methylation sites at once. A derivation exists, and implementing it is part of the remaining work.
- Time: 1.0
- Carry-over: prelim

## 3. Results

### 14. Simulation: two cohorts with shared and cohort-specific programs
- Figure: none (design table)
- Slide text: Two cohorts of 150 patients, 1,000 genes, true rank 6 · Loadings taken from a factorization of the real training data · Weibull survival from the shared programs (β = 1.5, −1.2, 0.8, −0.5), 30% censoring · Scenarios: all shared, shared plus cohort-specific, nothing shared · Starting ranks 2, 3, 4, 6, 12, 20 with 15 datasets each · Compared with a two-step fit: factorization, then an unpenalized Cox model
- Notes: The first simulation has two cohorts of 150 patients and a thousand genes. The true loadings come from a factorization of the real data, so the programs look like real expression programs. Survival depends on the shared programs. Three scenarios change how many programs the cohorts share, and each is fitted at six starting ranks with fifteen datasets per setting. The comparison is a two-step fit, the factorization followed by an ordinary Cox model.
- Time: 0.75
- Carry-over: prelim

### 15. The joint fit predicts better than the two-step fit when cohorts have their own programs
- Figure: single_sim_delta.png
- Slide text: Shared plus cohort-specific: C higher by 0.10–0.16 from a starting rank of 3 · All shared: about 0.02 higher from a starting rank of 6, slightly lower at 3–4 · Nothing shared: within about 0.01 · The loadings were not survival-guided in either arm, so the gain comes from the priors and the shrunken survival step
- Notes: Each bar is the joint fit's held-out concordance minus the two-step fit's. When each cohort has programs of its own, the joint fit is ahead by 0.10 to 0.16 from a starting rank of three. When every program is shared it is ahead by about 0.02 once the starting rank reaches six. With no survival signal the two are the same, as they should be. In these fits survival did not shape the loadings either, so the gain comes from the priors and the shrinkage of the survival coefficients.
- Time: 1.0
- Carry-over: prelim

### 16. Starting with too few programs merges true programs, and starting with too many does not
- Figure: none (two large numbers)
- Slide text: Starting ranks 2–4: 7 of 27 estimated programs were merges of two true programs · Starting ranks 6, 12, 20: 0 of 114 · Cohort attribution of a survival effect was correct in 15 of 20 cases at a starting rank of 6 and 17 of 20 at 12 · Start high and let the priors shrink
- Notes: Starting below the true rank forces true programs together: seven of 27 estimated programs matched two true programs at once. At or above the true rank there were no merges in 114 programs. The model also attributed a survival effect to the right cohort in most cases. So the working rule is to start with too many programs and let the priors shrink the extras.
- Time: 0.75
- Carry-over: prelim

### 17. The finalized model will be judged on a planted low-variance program
- Figure: none (design and pass criteria in boxes)
- Slide text: One planted program explains under 1% of expression variance but has a strong survival effect · 400 training patients, noise at several levels, 50 and 2,000 genes · Pass: finds the planted program (loading correlation ≥ 0.7, correct direction, survival effect distinguishable from zero) · Pass: flags no spurious survival programs · Pass: the fit converges · The same test is run on a fit where survival does not shape the loadings
- Notes: The first simulation could not test supervision, because survival never shaped the loadings. This test is built for that. One planted program explains less than one percent of expression variance but strongly predicts survival, exactly the kind of program a two-step analysis misses. The finalized model passes only if it finds that program in the right direction, does not call other programs prognostic, and converges, at both fifty and two thousand genes. The same data are fitted without survival feedback for comparison. Runs are in progress, and the current numbers are in the backup slides.
- Time: 0.75
- Carry-over: prelim

### 18. PDAC application: train on two cohorts, validate on five
- Figure: real_data_workflow_single.png
- Slide text: Training: TCGA (n = 144) and CPTAC (n = 129), RNA-seq, 2,064 genes after screening each cohort for mean and variance · Validation: Dijk, Moffitt, PACA-AU microarray, PACA-AU RNA-seq, Puleo (616 records) · Starting rank 7 · Every choice, including the risk direction, made from training data
- Notes: For the application we trained on two RNA-seq cohorts, TCGA and CPTAC, 273 patients in total, using 2,064 genes that passed a mean and variance screen in both. We then scored five external cohorts across microarray and RNA-seq platforms. Everything, including the direction of the risk score, was fixed from training before any validation outcome was used.
- Time: 0.75
- Carry-over: prelim

### 19. External concordance matches a rank-matched two-step model
- Figure: external_c_forest.png
- Slide text: Mean external C: 0.627 for the joint model, 0.623 for a two-step model with 7 programs, 0.601 for a two-step model with 40 programs and a LASSO Cox step · Difference from the 7-program model: 0.011 (95% CI −0.013 to 0.033) · Difference from the 40-program model: 0.026 (0.0002 to 0.050) · Lines are 95% bootstrap intervals for the joint model
- Notes: Across the five external cohorts the mean concordance was 0.627. A two-step model with the same seven programs reached 0.623, and the difference is well inside its confidence interval. Against an independently tuned two-step model with forty programs and a lasso Cox step, the joint model was ahead by 0.026, with an interval that just excludes zero. In short, prediction is comparable to a well-chosen two-step analysis. That is expected here, because survival did not shape these loadings.
- Time: 1.0
- Carry-over: prelim

### 20. Predicted risk separates survival in every external cohort
- Figure: km_risk_external.png
- Slide text: Patients split into thirds of predicted risk within each cohort · The lowest third had the best survival in all five cohorts · Log-rank p < 0.0001 (Puleo), 0.0005 (Dijk), 0.021 (PACA-AU RNA-seq), 0.075 (PACA-AU microarray), 0.17 (Moffitt)
- Notes: Splitting each external cohort into thirds of predicted risk, the low-risk third had the best survival every time. The separation is significant in three of the five cohorts. Moffitt, the weakest, also had the lowest concordance.
- Time: 0.75
- Carry-over: prelim

### 21. Two programs match the basal-like and classical subtypes
- Figure: km_programs_external.png
- Slide text: Program 7, shorter survival: glycolysis genes (GAPDH, ENO1, LDHA) and integrin and growth-factor receptor genes (ITGA3, MET) · enriched for the basal-like gene list (adjusted p = 0.0006) · correlation 0.58 with the PurIST basal-like probability · Program 3, longer survival: epithelial and secretory genes · correlation −0.58 with PurIST · classical signature hazard ratio below 1 in all five cohorts · External tertiles: p < 0.0001 (Program 7), p = 0.004 (Program 3)
- Notes: Two of the retained programs carry the survival signal. Program 7 is associated with shorter survival. Its top genes are glycolysis genes and integrin and growth-factor receptor genes, it is enriched for the published basal-like gene list, and it correlates with the PurIST basal-like classifier. Program 3 is associated with longer survival and lines up with the classical subtype. In the external cohorts, patients in the top third of Program 7 had the shortest survival and those in the top third of Program 3 the longest. These programs were found from expression, so they reproduce the known subtype axis rather than discover a new one.
- Time: 1.25
- Carry-over: prelim

### 22. The other programs describe stroma and immune cells
- Figure: programs_heatmap.png
- Slide text: Program 5: fibroblast and collagen genes (COL1A2, COL3A1, SPARC), led by epithelial-to-mesenchymal transition (adjusted p = 0.0005) · Program 6: lymphoid genes (IKZF1, IL10RA) with vascular and basement-membrane genes (COL14A1, TNXB) · Neither predicts survival in the external cohorts · Together the four programs correspond to the published tumor, stromal and immune factors, with stroma and immune cells split apart · Program 3's survival association comes mainly from CPTAC
- Notes: The other two programs explain expression without a survival association. Program 5 is the stroma, dominated by collagen and fibroblast genes. Program 6 combines lymphoid immune genes with vascular genes. Neither is prognostic in the external cohorts. Together the four programs match the tumor, stromal and immune factors reported for pancreatic cancer, with the stroma and immune cells separated. One caution: with cohort-specific coefficients, Program 3's protective effect comes almost entirely from CPTAC.
- Time: 1.0
- Carry-over: prelim

### 23. What the finalized and multi-data-type models are designed to show
- Figure: none (designed aims slide: four colored boxes, each with its test)
- Slide text: Survival-guided programs → test: the planted-program simulation, then training cross-validation in PDAC against the two-step model · Uncertainty carried into survival effects → test: interval coverage in simulation · Programs that link methylation to expression → test: programs loading on genes and CpG sites, CpG sites mapped to genes and promoters, validation in matched ICGC data · Number of programs from the data alone → test: the same rank across starting ranks, without validation cohorts
- Notes: This slide is about results we expect rather than results we have. First, programs shaped by survival, tested by the planted-program simulation and then by cross-validation in the training cohorts. Second, uncertainty in the programs carried into the survival coefficients, checked by interval coverage in simulation. Third, programs that include both genes and methylation sites. Mapping those sites to genes and promoters would show which program genes are regulated through methylation, for example hypomethylated promoters of basal-like genes alongside higher expression. That is useful even if concordance barely changes. Fourth, a number of programs chosen without looking at the validation cohorts.
- Time: 1.25
- Carry-over: prelim

## 4. Discussion

### 24. Summary
- Figure: none
- Slide text: A joint Bayesian model of expression and survival scores new patients from expression alone · The fitted version, with expression-only loadings, reached external C = 0.627 in five cohorts, similar to a rank-matched two-step model · It recovers the basal-like and classical axis and separate stromal and immune programs · In simulation, starting with too many programs avoids merging true programs
- Notes: To summarize, the model fits expression and survival together and scores new patients from expression alone. The version fitted so far validates in five external cohorts at about the level of a well-matched two-step analysis. It recovers the main subtype axis and separates stroma from immune cells. The simulations support starting with many programs and letting the priors shrink them.
- Time: 0.75
- Carry-over: prelim

### 25. Limitations
- Figure: none
- Slide text: In the fitted model, survival does not shape the programs · Prediction is comparable to a two-step model, not better · Genes were standardized to mean 0 while loadings are nonnegative · The starting rank of 7 was chosen partly from validation performance · The two PACA-AU series share donors (616 records, about 568 people) · Not yet compared with DeSurv, PurIST, elastic-net Cox or clinical covariates
- Notes: The main limitation is the one I've repeated: in the fitted model survival does not shape the programs, so the benefit of supervision is not shown yet on real data. Prediction is comparable to two-step analysis, not better. The data scaling is inconsistent with nonnegative loadings. The starting rank used validation performance. Two external series share donors, so pooled intervals are optimistic. And several comparators are still missing.
- Time: 1.0
- Carry-over: prelim

### 26. Remaining work
- Figure: none
- Slide text: Survival-guided loadings: rebalance the expression and survival terms, pass the planted-program test, check training cross-validation in PDAC · Refit on a nonnegative data scale and choose the number of programs without validation data · Implement the multi-data-type model and fit matched expression and methylation (TCGA training, ICGC validation) · Complete the biological characterization: enrichment, subtype agreement, CpG-to-gene annotation · Wider comparison: DeSurv, PurIST, elastic-net Cox, clinical covariates, donor-aware intervals · An R package with one fitting function, tests and a vignette · A manuscript whose results are regenerated from the package
- Notes: Here is what remains before the final exam. First and most important, get survival to shape the loadings and confirm it with the simulation test and with cross-validation in the training data. Then refit on a nonnegative scale and choose the number of programs without the validation cohorts. Next, implement the extension to expression and methylation and fit it to matched TCGA and ICGC data. Alongside that, finish the biological characterization and the wider comparison with other methods. The model then goes into an R package, and the manuscript will regenerate its results from that package.
- Time: 1.0
- Carry-over: prelim

## Backup

### B1. With survival shaping the loadings, the model finds a planted low-variance program that two-step methods miss
- Figure: planted_matched_comparators.png
- Slide text: 400 training patients, 50 genes, planted log hazard ratio 3.3, 20 datasets per noise level · Planted program found in 17 of 20 and 19 of 20 datasets · With survival feedback off, 1 of 20 each · Four comparator methods found it in 0 of 20 at both noise levels · Strong effect and few genes; not yet a general result
- Notes: Use when asked whether supervision can work at all. In this calibrated simulation the supervised fit found the planted program in 17 and 19 of 20 datasets, against one for the unsupervised fit and none for four comparators. The effect is strong and there are only 50 genes, so this shows the mechanism, not performance at real-data size.
- Carry-over: technical

### B2. Recovery and held-out concordance in the planted-program simulation
- Figure: planted_yfb_performance.png
- Slide text: Held-out C: 0.763 and 0.765 with survival feedback, 0.641 and 0.540 without · 19 of 20 supervised fits converged at each noise level · The planted program explains at most 0.49% and 0.04% of expression variance
- Notes: Use when asked about prediction in the planted simulation. Held-out concordance was about 0.76 with survival feedback against 0.54 to 0.64 without.
- Carry-over: technical

### B3. At 2,000 genes the supervised model finds the planted program, with slightly too many spurious survival programs
- Figure: none (table from the October 30 runs)
- Slide text: Planted program found in 10 of 10 datasets at each of three planted variance shares (0.25%, 1%, 5%) · 0 of 10 without survival feedback · All supervised fits converged · Mean spurious survival programs 1.0, 1.5 and 1.1, against a pass limit of 1; the unsupervised fit had 1.3 to 1.6 · Not yet checked on PDAC at this setting
- Notes: Use when asked about the larger setting. At two thousand genes, the supervised fit found the planted program every time, while the unsupervised fit never did. It narrowly fails the spurious-program limit in two of three settings.
- Carry-over: technical

### B4. With survival feedback on, PDAC concordance matches the two-step model
- Figure: joint_external_concordance.png
- Slide text: Mean external C 0.621 (survival feedback on) and 0.623 (two-step, 7 programs) · Paired difference −0.002 (95% bootstrap interval −0.012 to 0.006)
- Notes: Use when asked what happens on real data with survival feedback switched on. Concordance is the same as the two-step model.
- Carry-over: technical

### B5. With survival feedback on, PDAC programs are almost identical to unsupervised programs
- Figure: joint_program_correspondence.png
- Slide text: The two survival-associated programs correlate 0.993 and 0.975 with unsupervised programs · Gene-loading correlations 0.991 and 0.976 · No additional low-variance prognostic direction found
- Notes: Use with the previous slide. On real data the supervised fit found the same programs as the unsupervised fit, which is what we expect when survival has almost no weight in the loading update.
- Carry-over: technical

### B6. Survival carries almost no weight in the loading update at 2,000 genes
- Figure: none (two numbers and a short table)
- Slide text: Share of the loading update's precision from the survival term: about 1 at 50 genes, about 1 in a million at 2,064 genes · Raising the survival weight 10,000-fold brings the share to one half, but no cross-validation fold converges and C falls from 0.656 to 0.555 · No single survival weight works at this dimension, so the plan rebalances the expression term instead
- Notes: Use when asked why survival does not shape the loadings on PDAC. The survival term is about a millionth of the information in each loading update, and simply turning up its weight makes the fit unstable.
- Carry-over: technical

### B7. An earlier real-data check could not tell the supervised and fitted models apart
- Figure: none
- Slide text: Supervised fit (October 15) minus fitted model: ΔC = −0.003 (95% CI −0.029 to 0.024) · The check required the lower limit to be above −0.02 · The interval was too wide to show the supervised fit was no worse
- Notes: Use when asked about the October supervised fit on real data. Its concordance was essentially the same as the fitted model's, but the interval was too wide to meet the margin we had set in advance.
- Carry-over: technical

### B8. Fixing the noise level stops divergence but does not make supervision work
- Figure: none
- Slide text: With free noise precisions, large survival weights made reconstruction error explode (to 4.5 × 10¹⁴) · Fixed noise kept it at 0.079 · Of 12 fits, 1 converged, 8 hit the iteration limit and 3 stalled · Null programs were still flagged as prognostic
- Notes: Use when asked about stability. Holding the noise level fixed prevents the blow-up, but the fits still do not converge reliably and still flag null programs.
- Carry-over: technical

### B9. Pass criteria for the planted-program simulation
- Figure: none
- Slide text: Recovery: loading correlation ≥ 0.7 with the planted program, correct direction, survival effect at least 1.96 posterior SDs from zero · Gain over the unsupervised fit of at least 3 of 10 datasets in at least two settings · At most 1 spurious survival program on average · At least 90% of fits converged
- Notes: Use when asked how the simulation test is scored. These criteria were set before the runs.
- Carry-over: technical

### B10. Choosing the starting rank
- Figure: none
- Slide text: Five-fold cross-validated training C: best at K = 8 (0.651, SE 0.026) · K = 7 the smallest within one SE (0.633) · Over starting ranks 2–20, external C reached a plateau at 7 · The final choice therefore also used validation performance
- Notes: Use when asked how seven was chosen. Cross-validation pointed to seven or eight, and the external plateau was also at seven, so the choice was not made from training data alone.
- Carry-over: technical

### B11. Cohort-specific survival coefficients
- Figure: percohort_beta.png
- Slide text: Program 7's association appears in both training cohorts · Program 3's comes mainly from CPTAC (0.030 against 0.0004 in TCGA) · Cohort-specific coefficients did not change external concordance
- Notes: Use when asked whether the programs are prognostic in both training cohorts.
- Carry-over: technical

### B12. Preprocessing
- Figure: none
- Slide text: Top 3,000 genes per cohort by combined mean and variance rank, intersection 2,064 · Each gene standardized to mean 0 and SD 1 within cohort · External cohorts standardized the same way and scored on the training genes they measure · Standardized data contain negative values, which nonnegative loadings cannot reconstruct
- Notes: Use when asked about scaling. The standardization helped stability in earlier fits but conflicts with nonnegative loadings, and the finalized model will work on a nonnegative scale.
- Carry-over: technical

### B13. Enrichment of the four retained programs
- Figure: none (table)
- Slide text: Program 7: basal-like list (adj. p = 0.0006), growth-factor receptor, PI3K and RAS–ERK pathways (0.001), Bailey squamous list (0.03) · Program 3: classical and basal-like lists (0.017 each), Moffitt classical list · Program 5: EMT (0.0005), stromal-immune list, extracellular matrix · Program 6: mostly immune GO terms · Normalized enrichment scores 1.05–1.29
- Notes: Use when asked for the enrichment detail. The enrichment scores are modest, and Program 3's enrichment is less specific than its subtype correlation.
- Carry-over: technical

### B14. A risk score built on patient activities
- Figure: none
- Slide text: Alternative: η = Lβ, risk from each patient's program activities · A new patient's activities are estimated against the fixed loadings · Mean external C 0.611, against 0.627 for the projection score · The fallback if survival cannot shape the loadings through the projection score
- Notes: Use when asked about the alternative model. Basing risk on activities lets survival reach the activities directly, at the cost of estimating activities for each new patient.
- Carry-over: technical
