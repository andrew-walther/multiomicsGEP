# Abstract — Working Draft
# Title: Supervised Bayesian Matrix Factorization Identifies Prognostic Gene Expression Programs
# Status: Draft for review — results to be confirmed against final rendered reports
# Note: "Prognostic" preferred over "Survival-Associated" for now; revisit at submission

---

## Abstract

Gene expression programs (GEPs) — coordinated patterns of transcriptional activity
reflecting underlying biological processes — are a natural unit of analysis in cancer
genomics, yet most matrix factorization methods identify GEPs without reference to
clinical outcomes, leaving their prognostic relevance to be assessed post hoc.
Frequentist approaches that incorporate survival supervision exist, but typically
require manual specification of regularization hyperparameters — such as sparsity
penalties and shrinkage scales — that are difficult to tune in high-dimensional
genomic settings and may not adapt appropriately across datasets or platforms. We
present Supervised Bayesian Matrix Factorization (SBMF), a probabilistic framework
that jointly decomposes high-dimensional gene expression data and a time-to-event
survival outcome within a single inference procedure. SBMF models gene expression as
a low-rank matrix factorization and links patient factor loadings to survival through
a Cox proportional hazards model, placing sparse adaptive priors on all latent
quantities via the Empirical Bayes Normal Means framework. Critically, the degree of
sparsity and the scale of each factor's prior are estimated directly from the data via
marginal likelihood maximization, rather than fixed in advance — allowing the model to
adapt to the signal present in each dataset without manual hyperparameter tuning. This
enables simultaneous discovery of latent transcriptional structure and automatic
identification of which programs carry prognostic signal, with non-prognostic factors
shrunk toward zero by the data-informed prior. Inference is performed via Coordinate
Ascent Variational Inference with a factor-wise update schedule that ensures tight
coupling between genomic structure and survival signal at each iteration.

We validate SBMF on synthetic benchmarks, confirming accurate factor recovery and
correct shrinkage of null factors, then apply the method to seven pancreatic ductal
adenocarcinoma (PDAC) cohorts spanning RNA-seq and microarray platforms. All cohorts
converge reliably, and the supervised model achieves higher discrimination than a PCA
baseline in the majority of cohorts, with held-out prediction confirming generalizability
of the learned factor structure. A heavier-tailed sparse prior consistently improves
held-out discrimination across cohorts, suggesting that asymmetric gene contributions
are a characteristic feature of PDAC expression programs.

SBMF provides a unified framework for the simultaneous discovery and prognostic
evaluation of gene expression programs, with direct applicability to any cancer type
with paired transcriptomic and survival data. By embedding clinical supervision
directly into the factorization, rather than treating it as a downstream analysis
step, the method recovers biologically meaningful and clinically relevant programs
that purely unsupervised approaches may miss.

---

## Notes for revision

- Confirm C-index results and cohort counts against final rendered reports before submission
- Add a sentence on pathway enrichment once fgsea/Enrichr analysis is complete
- "Prognostic" vs. "Survival-Associated": revisit at submission — "survival-associated"
  is more conservative if C-index improvements are modest across cohorts
- Word count: ~310 words (most journals target 250–350)
