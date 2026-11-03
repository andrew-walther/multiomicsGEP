# Talk outline: Survival-supervised Bayesian matrix factorization (prelim Part III)

Part III of the prelim oral exam (2026-12-07). The exam runs about 90 minutes with
questions, and each of the three parts gets 20–25 minutes. The committee is four
statisticians and one non-statistician. Draft 1 (2026-11-03).

Part III is the one unfinished project, so the talk is a research talk about a proposed
project with its progress inside it. Each section gives what exists first, then what is
planned. Every slide is tagged **prelim** (main deck), **technical** (backup) or
**talk-only** (dropped when the deck joins the full prelim deck).

The deck is `ssbmf_talk.qmd` (rendered to `ssbmf_talk.html`; press S for the speaker
view with the read-aloud script, one sentence per bullet). The slide entries below are
generated from the deck by `python3 presentation/prelim-slides/tools/outline_from_deck.py`,
and `outline.html` by `python3 presentation/prelim-slides/tools/build_outline_html.py`.
Figures come from `tools/make_figures.R` and `tools/make_diagrams.py`.

## Structure and timing

Slide timings come from each script's word count at 140 words per minute, rounded to the
nearest quarter minute.

| Section | Slides | Minutes |
|---|---|---|
| Title and outline | 1–2 | 0.5 |
| 1. Introduction | 3–6 | 3.25 |
| 2. Methods | 7–13 | 7.25 |
| 3. Results (simulation and PDAC application) | 14–23 | 8.25 |
| 4. Discussion | 24–26 | 2.25 |
| Main body | 1–26 | 21.5 (2,988 words) |
| Backup | B1–B23 | as needed |

The target is about 20 minutes, running a little long if needed (author, 2026-11-03).
Pauses on the equations and figures will add to the read-aloud time.

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

## Checked on 2026-11-03

- iBAG (Wang et al. 2013, cited in the literature review): splits each gene's expression
  into a part explained by its promoter methylation and the rest, and relates both to
  survival in an accelerated failure time model. Slide 5 now says this.
- External donors: the PACA-AU series share 47 donors and the RNA-seq series repeats one,
  so the 616 records represent 568 people (writeup cohort table and limitations).
- The title slide lists Andrew only. Dr. Liu, Dr. Rashid and the Rashid Lab are thanked
  on the shared prelim thank-you slide.
- The diagrams were redrawn with plain labels by `tools/make_diagrams.py`.

## Open items

- Slides 12 and 23 leave space at the bottom and could take larger text.

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

### 1. Title
- Figure: none
- Slide text: Survival-supervised Bayesian matrix factorization identifies prognostic molecular programs Andrew J Walther Department of Biostatistics, University of North Carolina at Chapel Hill Prelim Part III · December 7, 2026
- Notes: My third project is a model that finds gene expression programs in pancreatic tumors and relates them to patient survival in one fit. The method is still being finished, so I'll show what is built and validated, and then what remains.
- Time: 0.25
- Carry-over: talk-only

### 2. Outline
- Figure: none
- Slide text: • Introduction Pancreatic cancer, two-step analyses of expression and survival, and the gap this project addresses • Methods The joint model, its priors and fitting, the finalized model, and the extension to several data types • Results Simulation, and validation in five external pancreatic cancer cohorts with biological interpretation • Discussion Summary, limitations, and remaining work
- Notes: I'll start with the clinical problem and how it is usually analyzed, then the model, then results in simulation and in pancreatic cancer data. I'll end with the limitations and the work left before the final exam.
- Time: 0.25
- Carry-over: talk-only

### 3. Pancreatic cancer is lethal, and tumors differ in ways that predict survival
- Figure: none
- Slide text: A lethal cancer 13.7% Five-year relative survival for pancreatic cancer (SEER, diagnoses 2016–2022) Tumors differ Expression separates basal-like and classical tumors, and basal-like tumors have shorter survival A bulk sample also carries signal from the stroma and immune cells around the tumor Goal: find expression programs that predict survival and can score a new patient from their expression alone SEER Cancer …
- Notes: Pancreatic ductal adenocarcinoma, PDAC, has a five-year survival of about fourteen percent. Earlier work found that tumor expression splits into a basal-like and a classical subtype, and patients with basal-like tumors die sooner. A bulk tumor sample also contains stroma and immune cells, so its expression mixes several sources. The goal of this project is to find the expression programs that matter for survival and to score a new patient with them using nothing but that patient's expression. Moffitt and colleagues found these subtypes by factorizing bulk tumor expression, separating tumor programs from stromal programs. So any method that finds programs in bulk expression has to handle tumor, stroma and immune signal at the same time.
- Time: 0.75
- Carry-over: prelim

### 4. The usual analysis finds programs first and tests survival afterward
- Figure: diagram_two_step_vs_joint.png
- Slide text: A program that explains little expression variance but predicts survival can be missed, because cell-composition variation dominates an expression-only fit The Cox step treats estimated scores as known, so its standard errors leave out the uncertainty in the scores
- Notes: The usual pipeline has two steps. First a matrix factorization finds programs in the expression matrix, and then a Cox model tests the program scores against survival. Survival never touches the first step, and that has two consequences. A program that explains only a little expression variance can be lost, even if it predicts survival well, because large sources of variation such as cell composition dominate the factorization. And the Cox step treats the estimated scores as if they were measured exactly, so its standard errors are too small. A joint fit estimates both pieces together. Picture a program that explains one percent of the expression variance but doubles the hazard. An expression-only factorization has little reason to keep it, and the Cox step can only test the programs it is given.
- Time: 1
- Carry-over: prelim

### 5. Existing methods each cover part of the problem
- Figure: none
- Slide text: Survival-supervised NMF DeSurv (Young et al.): survival guides which programs are found. Penalized point estimates, with the rank and the strength of supervision tuned by cross-validation Methylation and survival, gene by gene iBAG (Wang et al. 2013): splits each gene's expression into a part explained by its promoter methylation and the rest, and relates both to survival in an accelerated failure time model …
- Notes: Several methods cover pieces of this. Survival-supervised nonnegative matrix factorization from our group lets survival guide which programs are found, but it reports point estimates and tunes its rank and the strength of supervision by cross-validation. iBAG links methylation and expression to survival, but one gene at a time rather than through programs. MOFA and iCluster find factors shared across data types, but they ignore survival. None of them learns programs across data types with survival guiding them, carries uncertainty through, and estimates how many programs there are. DeSurv is the closest to our approach, and we use its gene screening rule so the two can be compared on the same genes. What we add is a Bayesian version: posterior distributions instead of point estimates, and shrinkage instead of cross-validation to set the number of programs.
- Time: 1
- Carry-over: prelim

### 6. Project aims
- Figure: none
- Slide text: Aim 1. Joint model A Bayesian model of expression and survival that can find low-variance prognostic programs, reports uncertainty, and shrinks programs it does not need Aim 2. PDAC validation Validate the model in external pancreatic cancer cohorts and characterize its programs biologically Aim 3. Several data types Extend the model to expression and methylation together, so that methylation changes can be …
- Notes: The project has three aims. The first is the joint model itself: it should find programs that predict survival even when they explain little expression variance, report uncertainty, and decide how many programs to keep. The second is to validate it in external pancreatic cancer cohorts and interpret its programs. The third is to extend it to several data types, so that a program can include both genes and methylation sites.
- Time: 0.5
- Carry-over: prelim


## 2. Methods

### 7. Each patient contributes an expression profile and a survival outcome
- Figure: diagram_notation.png
- Slide text: Y n × p expression, patients by genes (tᵢ, δᵢ) observed time and event indicator L n × K program activity of each patient F p × K gene loadings, one column per program Zᵢₖ patient i's expression projected onto program k's normalized loading βₖ log hazard ratio of program k $$\eta_i = \sum_{k=1}^K Z_{ik}\,\beta_k \qquad\text{(risk score)}$$
- Notes: Here is the notation. Y has a row for each patient and a column for each gene, and each patient also has a survival time and an event indicator. We approximate Y as L times F transpose with K programs. A column of F is a program, its weight on each gene, and a column of L gives how active that program is in each patient. The risk score projects a patient's expression onto each program's loading column, scaled to unit length, and multiplies by a coefficient beta, which is a log hazard ratio. Because the score uses expression and loadings only, a new patient can be scored without fitting anything. One detail matters for later. The score uses the projection of expression onto F, not the activities L, and because F sits inside the risk score, survival information can flow back into the loadings.
- Time: 1
- Carry-over: prelim

### 8. The joint model adds a Cox partial likelihood to a matrix factorization
- Figure: none
- Slide text: $$y_{ij}=\sum_{k=1}^K L_{ik}F_{jk}+e_{ij},\quad e_{ij}\sim N(0,\tau_j^{-1}) \qquad\quad h_i(t)=h_0(t)\exp(\eta_i),\quad \eta = Y_c F D\,\beta$$ $$\begin{aligned}\mathcal J(q)={}&\lambda_G\,E_q\big[\log p(Y\mid L,F,\tau)\big]+\gamma\,E_q\big[\ell_{\mathrm{Cox}}(\eta)\big]\\&-\mathrm{KL}(q_L\Vert g_L)-\mathrm{KL}(q_F\Vert g_F)-\mathrm{KL}(q_\beta\Vert g_\beta)\end{aligned}$$ $q = q_L\,q_F\,q_\beta$ is a mean-field …
- Notes: The model has two parts. The expression part is a Gaussian factorization with a separate noise precision for each gene. The survival part is a Cox model whose linear predictor is the projection score times beta. We fit a variational approximation q by maximizing this objective: the expected expression log-likelihood, plus the expected Cox log partial likelihood, minus a Kullback-Leibler penalty for each block of parameters. The weights lambda G and gamma set how much the expression data and the survival data count. The loadings F appear in both terms, and that shared F is how survival can shape the programs. Two special cases help to read it. With gamma equal to zero this is an ordinary empirical Bayes matrix factorization, and with F held fixed it is a Bayesian Cox model on projected scores. The joint model sits between the two, and the weights decide where.
- Time: 1
- Carry-over: prelim

### 9. Empirical Bayes priors let the data decide how much to shrink
- Figure: none
- Slide text: $$L_{ik}\sim(1-\pi_{Lk})\,\delta_0+\pi_{Lk}\,\mathrm{Exp}(\lambda_{Lk}) \qquad F_{jk}\sim(1-\pi_{Fk})\,\delta_0+\pi_{Fk}\,\mathrm{Exp}(\lambda_{Fk})$$ $$\beta_k\sim(1-\pi_\beta)\,\delta_0+\pi_\beta\,N(0,\sigma^2_\beta)$$ Nonnegative, sparse programs Point-exponential priors on activities and loadings: a point mass at zero mixed with an exponential Survival effects can be zero Point-normal prior on β in the …
- Notes: The activities and loadings get point-exponential priors: a point mass at zero mixed with an exponential, so they are nonnegative and sparse. The survival coefficients get a point-normal prior, which lets a program have no survival effect at all. The model fitted so far used a plain normal prior on beta instead. In each case the prior's parameters are estimated from the data, program by program, which is what makes this empirical Bayes. The point masses are also how the model decides how many programs it needs. Empirical Bayes matters here because there are thousands of loadings and only a few hundred patients. Estimating the prior from all the genes in a program borrows strength across genes, so each loading is shrunk by an amount the data choose rather than a penalty we tune.
- Time: 1
- Carry-over: prelim

### 10. The model is fitted by coordinate ascent, starting with too many programs
- Figure: diagram_fitting.png
- Slide text: Each update is an empirical Bayes normal-means problem: the data give an estimate and its precision, and the prior shrinks it The direction of the risk score is fixed once, from the training data
- Notes: We fit the model by coordinate ascent variational inference. The Cox term is not quadratic, so at the start of each sweep it is replaced by a weighted least squares approximation around the current risk score. Then each program is updated in turn, its coefficient, its activities and its loadings, and the noise precisions are updated at the end of the sweep. We start with more programs than we expect and let the priors remove the ones that are not needed. The sign of the risk score is fixed once from the training data. Each update has a closed form: the data give an estimate and a standard error for one parameter, and the fitted prior shrinks it. The working approximation keeps each patient's Cox score exact but replaces the full curvature, which couples patients through their risk sets, with a diagonal one. That is what allows one coordinate at a time, and the equations are in a backup slide.
- Time: 1.25
- Carry-over: prelim

### 11. A new patient's risk score needs only their expression
- Figure: diagram_prediction.png
- Slide text: No patient-level parameters are refitted, so training and validation patients are scored with the same formula The direction of the score is set in training, before any validation outcome is seen
- Notes: To score a patient in a new cohort, we restrict to the training genes, standardize within that cohort, project onto the fitted loadings and multiply by the fitted coefficients. Nothing about the new patient is estimated. That keeps validation clean, because the same formula is used for training and test patients. Two-step methods also score new patients by projection, so differences in validation come from the fitted loadings and coefficients, not from the scoring rule.
- Time: 0.5
- Carry-over: prelim

### 12. So far survival sets the coefficients but not the loadings
- Figure: none
- Slide text: Model fitted so far Survival sets β, but F is updated from expression only Genes standardized to mean 0 and SD 1 within cohort, although the loadings are nonnegative Normal prior on β, with the projection scores held fixed in the β update Finalized model The survival term enters the loading update Obstacle: with about 2,000 genes, the expression term outweighs the survival term in each loading update Plan: …
- Notes: The model fitted so far is a step short of the full proposal. Survival sets the coefficients, but the loadings are updated from expression only, so its programs are expression programs that survival then weights. It also standardized each gene, which gives negative values, while the loadings are nonnegative, and that mismatch will be removed. The finalized model lets the survival term into the loading update. The obstacle is scale: with about two thousand genes, the expression part of each loading update outweighs the survival part, so survival barely moves the loadings. The plan is to rebalance the two weights and to accept the result only if it passes the simulation test I'll show in the results. If that fails, the fallback builds the risk score on patient activities instead. To make the scale problem concrete, the precision of each loading update has an expression part and a survival part. At fifty genes the survival part is most of it, and at about two thousand genes it is roughly a millionth of the total. Multiplying the survival weight by ten thousand makes the fit unstable, so the plan rebalances from the expression side instead.
- Time: 1.5
- Carry-over: prelim

### 13. The extension to several data types shares patient activities across them
- Figure: diagram_multi_data.png
- Slide text: $$Y^{(m)}\approx \mathbf 1\mu^{(m)\top}+L\,F^{(m)\top},\quad m=1,\dots,M$$ $$\eta=\Big(\textstyle\sum_m Y^{(m)}F^{(m)}D^{(m)}\Big)\beta$$ Each data type has its own loading prior family and noise precisions, for example a sparser prior for many methylation sites $M=1$ gives the single data type model One program can load on genes and methylation sites together
- Notes: The extension keeps one matrix of patient activities and gives each data type its own loadings. Each data type also gets its own prior family and noise levels, which matters because there are many more methylation sites than genes and most of them should have zero weight. The risk score sums the projections over data types, and with one data type it is exactly the model I just described. The point of the extension is that one program can include genes and methylation sites at once. A derivation exists, and implementing it is part of the remaining work. Cohorts and data types enter differently: more patients add rows of L, and more data types add blocks of F. Each feature also gets an intercept, which keeps spare programs from absorbing feature means, a problem we saw in an earlier version.
- Time: 1
- Carry-over: prelim


## 3. Results

### 14. Simulation: two cohorts with shared and cohort-specific programs
- Figure: none
- Slide text: Data Two cohorts of 150 patients, 1,000 genes, true rank 6 True loadings taken from a factorization of the real training data Survival Weibull proportional hazards from the shared programs, β = 1.5, −1.2, 0.8, −0.5 30% censoring, 75 test patients per dataset Settings All shared, shared plus cohort-specific, nothing shared Starting ranks 2, 3, 4, 6, 12, 20, with 15 datasets each Comparison: a two-step fit, the …
- Notes: The first simulation has two cohorts of 150 patients and a thousand genes. The true loadings come from a factorization of the real data, so the programs look like real expression programs. Survival depends on the shared programs. Three scenarios change how many programs the cohorts share, and each is fitted at six starting ranks with fifteen datasets per setting. The comparison is a two-step fit, the factorization followed by an ordinary Cox model. Starting rank means the number of programs a fit begins with, before shrinkage, and the true rank is six.
- Time: 0.75
- Carry-over: prelim

### 15. The joint fit predicts better than the two-step fit when cohorts have their own programs
- Figure: single_sim_delta.png
- Slide text: Shared plus cohort-specific: C higher by 0.10–0.16 from a starting rank of 3 All shared: about 0.02 higher from a starting rank of 6, slightly lower at 3–4 Nothing shared: within about 0.01 Survival did not shape the loadings in either arm, so the gain comes from the priors and the shrunken survival step
- Notes: Each bar is the joint fit's held-out concordance minus the two-step fit's. When each cohort has programs of its own, the joint fit is ahead by 0.10 to 0.16 from a starting rank of three. When every program is shared it is ahead by about 0.02 once the starting rank reaches six. With no survival signal the two are the same, as they should be. In these fits survival did not shape the loadings either, so the gain comes from the priors and the shrinkage of the survival coefficients. The two-step arm here uses an unpenalized Cox step, so part of the gain is shrinkage that a penalized two-step model would also get. That is why the PDAC comparison uses a penalized two-step model as well.
- Time: 1
- Carry-over: prelim

### 16. Starting with too few programs merges true programs, and starting with too many does not
- Figure: none
- Slide text: Starting ranks 2–4 7 of 27 estimated programs were merges of two true programs Starting ranks 6, 12, 20 0 of 114 estimated programs were merges With cohort-specific coefficients, the cohort carrying a survival effect was identified in 15 of 20 cases at a starting rank of 6 and 17 of 20 at 12 Working rule: start high and let the priors shrink the extra programs
- Notes: Starting below the true rank forces true programs together: seven of 27 estimated programs matched two true programs at once. At or above the true rank there were no merges in 114 programs. The model also attributed a survival effect to the right cohort in most cases. So the working rule is to start with too many programs and let the priors shrink the extras. On the real data we followed this rule: we started at seven programs, kept four, and three shrank to zero.
- Time: 0.5
- Carry-over: prelim

### 17. The finalized model will be judged on a planted low-variance program
- Figure: none
- Slide text: Design One planted program explains under 1% of expression variance but has a strong survival effect 400 training patients, several noise levels, 50 and 2,000 genes The same data fitted without survival feedback for comparison The model passes only if it finds the planted program: loading correlation ≥ 0.7, correct direction, survival effect distinguishable from zero flags no spurious survival programs converges …
- Notes: The first simulation could not test supervision, because survival never shaped the loadings, and this test is built for that. One planted program explains less than one percent of expression variance but strongly predicts survival, exactly the kind of program a two-step analysis misses. The finalized model passes only if it finds that program in the right direction, does not call other programs prognostic, and converges, at both fifty and two thousand genes. The same data are fitted without survival feedback for comparison. Runs are in progress, and the current numbers are in the backup slides. The variance share is set by rescaling the planted program's activities, so the same survival effect sits on a weaker and weaker expression signal across settings. Recovery is judged from the posterior of the coefficient, so a program has to be both found and called prognostic.
- Time: 1
- Carry-over: prelim

### 18. PDAC application: train on two cohorts, validate on five
- Figure: real_data_workflow_single.png
- Slide text: Training: TCGA (n = 144) and CPTAC (n = 129), RNA-seq, 2,064 genes after screening each cohort for mean and variance, starting rank 7 Validation: Dijk, Moffitt, PACA-AU microarray, PACA-AU RNA-seq and Puleo (616 records), with every choice fixed from training
- Notes: For the application we trained on two RNA-seq cohorts, TCGA and CPTAC, 273 patients in total, using 2,064 genes that passed a mean and variance screen in both. We then scored five external cohorts across microarray and RNA-seq platforms. Everything, including the direction of the risk score, was fixed from training before any validation outcome was used. The external cohorts differ in platform, size and death rate, from 52 to 288 patients, which is a fairly hard test for a model trained on RNA-seq.
- Time: 0.75
- Carry-over: prelim

### 19. External concordance matches a rank-matched two-step model
- Figure: external_c_forest.png
- Slide text: | Model | Mean external C | |---|---| | Joint model | 0.627 | | Two-step, 7 programs | 0.623 | | Two-step, 40 programs, LASSO Cox | 0.601 | Difference from the 7-program model: 0.011 (95% CI −0.013 to 0.033) Difference from the 40-program model: 0.026 (0.0002 to 0.050) Lines are 95% bootstrap intervals for the joint model
- Notes: Across the five external cohorts the mean concordance was 0.627. A two-step model with the same seven programs reached 0.623, and the difference is well inside its confidence interval. Against an independently tuned two-step model with forty programs and a lasso Cox step, the joint model was ahead by 0.026, with an interval that just excludes zero. In short, prediction is comparable to a well-chosen two-step analysis. That is expected here, because survival did not shape these loadings. The rank-matched model was ahead in both PACA-AU series and behind in Dijk, so no single cohort decides the comparison. The bootstrap resamples records rather than donors, so with the shared PACA-AU donors these intervals are somewhat too narrow.
- Time: 1
- Carry-over: prelim

### 20. Predicted risk separates survival in every external cohort
- Figure: km_risk_external.png
- Slide text: The lowest third of predicted risk had the best survival in all five cohorts Log-rank p < 0.0001 (Puleo), 0.0005 (Dijk), 0.021 (PACA-AU RNA-seq), 0.075 (PACA-AU microarray), 0.17 (Moffitt)
- Notes: Splitting each external cohort into thirds of predicted risk, the low-risk third had the best survival every time. The separation is significant in three of the five cohorts. Moffitt, the weakest, also had the lowest concordance. The cut points are thirds within each cohort, and no validation outcome was used to set them.
- Time: 0.5
- Carry-over: prelim

### 21. Two programs match the basal-like and classical subtypes
- Figure: km_programs_external.png
- Slide text: Program 7: shorter survival Glycolysis (GAPDH, ENO1, LDHA), integrin and growth-factor receptor genes (ITGA3, MET) Basal-like gene list enriched (adjusted p = 0.0006), correlation 0.58 with PurIST basal-like probability Program 3: longer survival Epithelial and secretory genes Correlation −0.58 with PurIST, classical signature hazard ratio below 1 in all five cohorts External tertiles: p
- Notes: Two of the retained programs carry the survival signal. Program 7 is associated with shorter survival. Its top genes are glycolysis genes and integrin and growth-factor receptor genes, it is enriched for the published basal-like gene list, and it correlates with the PurIST basal-like classifier. Program 3 is associated with longer survival and lines up with the classical subtype. In the external cohorts, patients in the top third of Program 7 had the shortest survival and those in the top third of Program 3 the longest. These programs were found from expression, so they reproduce the known subtype axis rather than discover a new one. PurIST is a published single-sample subtype classifier from our group, so its agreement is an outside check on the subtype label. Joint coefficients can change sign when correlated programs are in the model together, so the adverse and protective labels come from each program's own survival association.
- Time: 1
- Carry-over: prelim

### 22. The other programs describe stroma and immune cells
- Figure: programs_heatmap.png
- Slide text: Program 5: fibroblast and collagen genes (COL1A2, COL3A1, SPARC), led by epithelial-to-mesenchymal transition (adjusted p = 0.0005) Program 6: lymphoid genes (IKZF1, IL10RA) with vascular and basement-membrane genes (COL14A1, TNXB) Neither predicts survival in the external cohorts Together the four programs correspond to the published tumor, stromal and immune factors, with stroma and immune cells split apart …
- Notes: The other two programs explain expression without a survival association. Program 5 is the stroma, dominated by collagen and fibroblast genes. Program 6 combines lymphoid immune genes with vascular genes. Neither is prognostic in the external cohorts. Together the four programs match the tumor, stromal and immune factors reported for pancreatic cancer, with the stroma and immune cells separated. One caution: with cohort-specific coefficients, Program 3's protective effect comes almost entirely from CPTAC. The stromal program leans weakly toward basal-like, with a correlation of 0.30, and its mesenchymal signature is not consistently prognostic in the external cohorts. The immune program is not correlated with the subtype axis at all.
- Time: 0.75
- Carry-over: prelim

### 23. What the finalized and multi-data-type models are designed to show
- Figure: none
- Slide text: Survival-guided programs Programs that predict survival even when they explain little expression variance Test Planted-program simulation, then training cross-validation in PDAC against the two-step model Uncertainty carried through Uncertainty in the programs reflected in the survival coefficients Test Interval coverage of β in simulation Methylation linked to expression Programs that load on genes and CpG sites …
- Notes: This slide is about results we expect rather than results we have. First, programs shaped by survival, tested by the planted-program simulation and then by cross-validation in the training cohorts. Second, uncertainty in the programs carried into the survival coefficients, checked by interval coverage in simulation. Third, programs that include both genes and methylation sites. Mapping those sites to genes and promoters would show which program genes are regulated through methylation, which is useful even if concordance barely changes. Fourth, a number of programs chosen without looking at the validation cohorts. For methylation, array probes can be mapped to genes and promoter regions with the standard annotation, which is how a methylation loading becomes a statement about a gene. Each of these has a test that could fail, and that is how the remaining work is set up.
- Time: 1
- Carry-over: prelim


## 4. Discussion

### 24. Summary
- Figure: none
- Slide text: A joint Bayesian model of expression and survival scores new patients from expression alone The fitted version, with expression-only loadings, reached a mean external C of 0.627 in five cohorts, similar to a rank-matched two-step model It recovers the basal-like and classical axis and separate stromal and immune programs In simulation, starting with too many programs avoids merging true programs
- Notes: To summarize, the model fits expression and survival together and scores new patients from expression alone. The version fitted so far validates in five external cohorts at about the level of a well-matched two-step analysis. It recovers the main subtype axis and separates stroma from immune cells. The simulations support starting with many programs and letting the priors shrink them. What is not shown yet is the central claim of the proposal, that survival can shape the programs themselves.
- Time: 0.5
- Carry-over: prelim

### 25. Limitations
- Figure: none
- Slide text: In the fitted model, survival does not shape the programs Prediction is comparable to a two-step model, not better Genes were standardized to mean 0 while the loadings are nonnegative The starting rank of 7 was chosen partly from validation performance The two PACA-AU series share 47 donors, so the 616 external records represent 568 people Not yet compared with DeSurv, PurIST, elastic-net Cox or clinical covariates
- Notes: The main limitation is the one I've repeated: in the fitted model survival does not shape the programs, so the benefit of supervision is not shown yet on real data. Prediction is comparable to two-step analysis, not better. The data scaling is inconsistent with nonnegative loadings. The starting rank used validation performance. Two external series share donors, so pooled intervals are optimistic. And several comparators are still missing. The rank and donor issues have direct fixes: choose the rank inside training cross-validation only, and keep one record per donor or resample donors.
- Time: 0.75
- Carry-over: prelim

### 26. Remaining work
- Figure: none
- Slide text: Survival-guided loadings: rebalance the expression and survival terms, pass the planted-program test, and check training cross-validation in PDAC Model details: refit on a nonnegative data scale, and choose the number of programs without validation data Several data types: implement the model and fit matched expression and methylation (TCGA training, ICGC validation) Biology: complete the characterization, with …
- Notes: Here is what remains before the final exam. First and most important, get survival to shape the loadings and confirm it with the simulation test and with cross-validation in the training data. Then refit on a nonnegative scale and choose the number of programs without the validation cohorts. Next, implement the extension to expression and methylation and fit it to matched TCGA and ICGC data. Alongside that, finish the biological characterization and the wider comparison with other methods. The model then goes into an R package, and the manuscript will regenerate its results from that package. The survival work comes first because the package, the extension to several data types and the manuscript all depend on a single-modality model that passes its test.
- Time: 1
- Carry-over: prelim


## Backup

### B1. With survival shaping the loadings, the model finds a planted low-variance program that two-step methods miss
- Figure: planted_matched_comparators.png
- Slide text: 400 training patients, 50 genes, planted log hazard ratio 3.3, 20 datasets per noise level Found in 17 of 20 and 19 of 20 datasets. With survival feedback off, 1 of 20 each. Four comparator methods: 0 of 20 Strong effect and few genes, so this shows the mechanism, not performance at real-data size
- Notes: Use when asked whether supervision can work at all. The supervised fit found the planted program in 17 and 19 of 20 datasets, against one for the unsupervised fit and none for four comparators.
- Carry-over: technical

### B2. Recovery and held-out concordance in the planted-program simulation
- Figure: planted_yfb_performance.png
- Slide text: Held-out C 0.763 and 0.765 with survival feedback, 0.641 and 0.540 without 19 of 20 supervised fits converged at each noise level The planted program explains at most 0.49% and 0.04% of expression variance
- Notes: Use when asked about prediction in the planted simulation.
- Carry-over: technical

### B3. At 2,000 genes the supervised model finds the planted program, with slightly too many spurious survival programs
- Figure: none
- Slide text: | Planted variance share | Found, survival feedback | Found, no feedback | Spurious survival programs, feedback (no feedback) | Converged | |---|---|---|---|---| | 0.25% | 10 of 10 | 0 of 10 | 1.0 (1.6) | 10 of 10 | | 1% | 10 of 10 | 0 of 10 | 1.5 (1.3) | 10 of 10 | | 5% | 10 of 10 | 1 of 10 | 1.1 (1.1) | 10 of 10 | The pass limit is at most 1 spurious survival program on average, so this narrowly fails in two of …
- Notes: Use when asked about the larger setting. At two thousand genes the supervised fit found the planted program every time, while the unsupervised fit never did.
- Carry-over: technical

### B4. With survival feedback on, PDAC concordance matches the two-step model
- Figure: joint_external_concordance.png
- Slide text: Mean external C 0.621 (survival feedback on) and 0.623 (two-step, 7 programs) Paired difference −0.002 (95% bootstrap interval −0.012 to 0.006)
- Notes: Use when asked what happens on real data with survival feedback switched on.
- Carry-over: technical

### B5. With survival feedback on, PDAC programs are almost identical to unsupervised programs
- Figure: joint_program_correspondence.png
- Slide text: The two survival-associated programs correlate 0.993 and 0.975 with unsupervised programs, and their gene loadings 0.991 and 0.976 No additional low-variance prognostic program was found
- Notes: Use with the previous slide. On real data the supervised fit found the same programs as the unsupervised fit, as expected when survival has almost no weight in the loading update.
- Carry-over: technical

### B6. Survival carries almost no weight in the loading update at 2,000 genes
- Figure: none
- Slide text: $$A_{F_{jk}}=\underbrace{\tau_j\sum_i E[L_{ik}^2]}_{\text{expression}}+\underbrace{E[\beta_k^2]\sum_i w_i\,y_{ij}^2}_{\text{survival}}\qquad\text{(precision of the update for } F_{jk})$$ Share of this precision from the survival term: about 1 at 50 genes, about 1 in a million at 2,064 genes Raising the survival weight 10,000-fold brings the share to one half, but no cross-validation fold converges and C falls …
- Notes: Use when asked why survival does not shape the loadings on PDAC. The survival term is about a millionth of the information in each loading update, and simply turning up its weight makes the fit unstable.
- Carry-over: technical

### B7. An earlier real-data check could not tell the supervised and fitted models apart
- Figure: none
- Slide text: Supervised fit (October 15) minus the fitted model: ΔC = −0.003 (95% CI −0.029 to 0.024) The check, set in advance, required the lower limit to be above −0.02 The interval was too wide to show the supervised fit was no worse
- Notes: Use when asked about the October supervised fit on real data. Its concordance was essentially the same as the fitted model's, but the interval was too wide to meet the margin set in advance.
- Carry-over: technical

### B8. Fixing the noise level stops divergence but does not make supervision work
- Figure: none
- Slide text: With free noise precisions, large survival weights made the reconstruction error explode (to 4.5 × 10¹⁴) Fixed noise kept it at 0.079 Of 12 fits, 1 converged, 8 hit the iteration limit and 3 stalled Null programs were still flagged as prognostic
- Notes: Use when asked about stability. Holding the noise level fixed prevents the blow-up, but the fits still do not converge reliably and still flag null programs.
- Carry-over: technical

### B9. Pass criteria for the planted-program simulation
- Figure: none
- Slide text: Recovery: loading correlation ≥ 0.7 with the planted program, correct direction, survival effect at least 1.96 posterior SDs from zero Gain: at least 3 of 10 more datasets recovered than the unsupervised fit, in at least two settings Specificity: at most 1 spurious survival program on average Convergence: at least 90% of fits
- Notes: Use when asked how the simulation test is scored. These criteria were set before the runs.
- Carry-over: technical

### B10. Training and validation cohorts
- Figure: none
- Slide text: | Cohort | Platform | Patients or records | Deaths | Role | |---|---|---|---|---| | TCGA-PAAD | RNA-seq | 144 | 75 | Training | | CPTAC | RNA-seq | 129 | 64 | Training | | Dijk | RNA-seq | 90 | 81 | Validation | | Moffitt | Microarray | 123 | 83 | Validation | | PACA-AU | Microarray | 63 | 38 | Validation | | PACA-AU | RNA-seq | 52 | 31 | Validation | | Puleo | Microarray | 288 | 181 | Validation | 616 validation …
- Notes: Use when asked about the cohorts or the event counts.
- Carry-over: technical

### B11. Preprocessing
- Figure: none
- Slide text: Top 3,000 genes per training cohort by combined mean and variance rank, intersection 2,064 Each gene standardized to mean 0 and SD 1 within cohort External cohorts standardized the same way and scored on the training genes they measure Standardized data contain negative values, which nonnegative loadings cannot reconstruct. The finalized model will work on a nonnegative scale
- Notes: Use when asked about scaling. The standardization helped stability in earlier fits but conflicts with nonnegative loadings.
- Carry-over: technical

### B12. The Cox working approximation
- Figure: none
- Slide text: $$P_{ri}=\frac{\mathbf 1\{i\in R_r\}\,e^{\eta_i^0}}{\sum_{a\in R_r}e^{\eta_a^0}},\qquad u_i=\delta_i-\sum_r d_rP_{ri},\qquad w_i=\sum_r d_rP_{ri}$$ $$\ell_{\mathrm{Cox}}(\eta)\approx \ell_{\mathrm{Cox}}(\eta^0)+u^\top(\eta-\eta^0)-\tfrac12(\eta-\eta^0)^\top W(\eta-\eta^0),\qquad W=\operatorname{diag}(w)$$ Breslow handling of ties, with the score kept exact The exact negative Hessian couples patients through …
- Notes: Use when asked how the Cox term fits into coordinate ascent.
- Carry-over: technical

### B13. Each coordinate update is an empirical Bayes normal-means problem
- Figure: none
- Slide text: $$\log q^*(\theta)=-\tfrac{A_\theta}{2}\theta^2+B_\theta\,\theta+\log g(\theta)+C\quad\Rightarrow\quad x_\theta=B_\theta/A_\theta,\quad s_\theta^2=A_\theta^{-1}$$ $$A_{\beta_k}=\sum_i w_i\big(\bar Z_{ik}^2+v_{Z_{ik}}\big),\qquad B_{\beta_k}=\sum_i \bar Z_{ik}\,h_i^{-k},\qquad h_i=u_i+w_i\eta_i^0$$ The prior is fitted to the observations x with standard errors s (ebnm), and the posterior moments come back The β …
- Notes: Use when asked for the update equations.
- Carry-over: technical

### B14. How the fitted model was run
- Figure: none
- Slide text: Starting rank 7, initialized from absolute singular vectors of Y Survival weight 0.5 on the coefficient update, and no survival term in the loading update Stopped when a monitoring statistic changed by less than 10⁻⁵ (met at sweep 20). The statistic is not the evidence lower bound Coefficient posterior variances were not stored, so no coefficient intervals are reported for this fit
- Notes: Use when asked about convergence or settings of the fitted model.
- Carry-over: technical

### B15. Why the risk score projects expression onto the loadings
- Figure: none
- Slide text: Lₖ Fₖᵀ is unchanged when Lₖ is multiplied and Fₖ divided by the same constant, so the score uses each loading column scaled to unit length β is then a log hazard ratio per unit of normalized projection A new patient needs no fitted activity, so training and validation use the same formula The alternative, η = Lβ, lets survival act on the activities directly but needs each new patient's activities estimated (mean …
- Notes: Use when asked why the model does not use L in the Cox model.
- Carry-over: technical

### B16. Retained, survival-associated and inactive programs
- Figure: none
- Slide text: Survival-active: absolute posterior mean coefficient above 0.001 Expression-only: not survival-active, and explains more than 1% of the variance of Y Inactive: neither From a starting rank of 7: two survival-active (7, 3), two expression-only (5, 6), three inactive Both thresholds are operational and are not posterior probabilities
- Notes: Use when asked how programs were labeled.
- Carry-over: technical

### B17. Choosing the starting rank
- Figure: none
- Slide text: Five-fold cross-validated training C was best at K = 8 (0.651, SE 0.026), and K = 7 was the smallest within one SE (0.633) Over starting ranks 2–20, external C reached a plateau at 7 The final choice therefore also used validation performance
- Notes: Use when asked how seven was chosen. Cross-validation pointed to seven or eight, and the external plateau was also at seven, so the choice was not made from training data alone.
- Carry-over: technical

### B18. The two-step comparators
- Figure: none
- Slide text: Factorization: unsupervised empirical Bayes matrix factorization (flashier), point-normal priors, gene-specific variances, greedy fitting with backfitting 40 programs, LASSO Cox: penalty chosen by cross-validation, 15 of 40 programs kept 7 programs, Cox: rank matched to the joint model 20 programs, Cox: mean external C 0.581 Held-out patients scored with the fixed loadings and coefficients, and the risk direction …
- Notes: Use when asked how the comparators were built.
- Carry-over: technical

### B19. External concordance by cohort
- Figure: none
- Slide text: | Cohort | Records | Joint model (95% CI) | Two-step, 40, LASSO | Two-step, 7 | Two-step, 20 | |---|---|---|---|---|---| | Dijk | 90 | 0.634 (0.555, 0.708) | 0.600 | 0.577 | 0.592 | | Moffitt | 123 | 0.549 (0.467, 0.625) | 0.560 | 0.541 | 0.541 | | PACA-AU microarray | 63 | 0.648 (0.553, 0.747) | 0.623 | 0.665 | 0.597 | | PACA-AU RNA-seq | 52 | 0.657 (0.544, 0.772) | 0.603 | 0.690 | 0.573 | | Puleo | 288 | 0.645 …
- Notes: Use when asked for cohort-level numbers. The rank-matched two-step model was ahead in both PACA-AU series.
- Carry-over: technical

### B20. Cohort-specific survival coefficients
- Figure: percohort_beta.png
- Slide text: Program 7's association appears in both training cohorts Program 3's comes mainly from CPTAC (0.030 against 0.0004 in TCGA) Cohort-specific coefficients did not change external concordance
- Notes: Use when asked whether the programs are prognostic in both training cohorts.
- Carry-over: technical

### B21. Enrichment of the four retained programs
- Figure: none
- Slide text: | Program | Leading enriched sets (adjusted p) | |---|---| | 7 (shorter survival) | Basal-like list (0.0006), growth-factor receptor, PI3K and RAS–ERK pathways (0.001), Bailey squamous list (0.03) | | 3 (longer survival) | Classical and basal-like lists (0.017 each), Moffitt classical list | | 5 (stroma) | Epithelial-to-mesenchymal transition (0.0005), stromal-immune list, extracellular matrix | | 6 (immune, …
- Notes: Use when asked for the enrichment detail. Program 3's enrichment is less specific than its subtype correlation.
- Carry-over: technical

### B22. How this model differs from survival-supervised NMF
- Figure: none
- Slide text: | | Survival-supervised NMF (DeSurv) | This model | |---|---|---| | Factorization | NMF, penalized point estimates | Empirical Bayes factorization with posterior distributions | | Number of programs | Chosen by cross-validation | Shrinkage from an over-specified start | | Survival strength | A tuned weight | Likelihood weights, with survival-shaped loadings still to be shown | | Uncertainty | Not reported for …
- Notes: Use when asked how the method differs from DeSurv. Projection-based prediction is common to both, so it is not claimed as new here.
- Carry-over: technical

### B23. The model for several data types in more detail
- Figure: none
- Slide text: A feature intercept μ⁽ᵐ⁾ for each data type, re-estimated every iteration, so spare programs do not absorb feature means The loading prior family is chosen per data type (point-exponential, point-Laplace or normal) and fitted for each data type and program Reconstruction and survival terms enter each loading update as a plain sum, with no tempering weight Each program is rescaled after its update so its …
- Notes: Use when asked how the extension is specified. The derivation exists, and the implementation is part of the remaining work.
- Carry-over: technical
