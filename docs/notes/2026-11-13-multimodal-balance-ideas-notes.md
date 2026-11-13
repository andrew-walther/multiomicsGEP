# Ideas for balancing expression and methylation in the joint model

Date: 2026-11-13. Author: Andrew Walther.

Status: ideas only. Nothing here has been run. The multimodal work (task D) was paused on
2026-11-13 so the time could go to the advisor meeting and the prelim submission.

## The problem

The joint model shares the patient scores L across data types. The precision of L sums
over every feature, so with 10,000 CpG sites against 2,218 genes, methylation supplies most
of that precision. On TCGA at 10k CpGs it supplied 83–97% in three of five programs.

Raising the methylation reconstruction likelihood to a power w < 1 did not fix this. The
weight lowers methylation's influence on L, but it also lowers each CpG's evidence for its
own loading by the same factor. Empirical Bayes then refits the methylation prior to a
point mass, and the methylation loadings disappear before the data types are balanced. At
w = 0.7, 42 to 46 of 50 simulated methylation blocks were empty. Lower weights were worse.
See [the weight-path findings](../findings/2026-11-12-weight-path-findings.md).

Separate loading priors do help when the feature counts are balanced. At 50:50 they found
expression-only, methylation-only and joint programs more often than a pooled prior.

## Ideas, cheapest first

1. **Weight methylation only in the score update.** Use Σ_m w_m Σ_j τ_mj E[F²_mjk] for the
   precision and mean of L, but keep the F_m and τ updates at full weight. Each CpG keeps
   its full evidence for its own loading, so the collapse mechanism above does not apply.
   - *Cost:* this is no longer coordinate ascent on one objective, so the ELBO can no longer
     serve as the stopping rule. Convergence would be judged by parameter change, and the
     fixed point needs its own derivation.
   - *Effort:* a small change to the fitter. This is the first one to try.
2. **Protect the methylation prior from the weight.** Keep the power weight, but fit each
   methylation prior by empirical Bayes on the unweighted pseudo-observations, or put a
   floor on the slab weight. This keeps a coherent weighted update for the posterior means
   while stopping the prior from collapsing to a point mass.
   - *Cost:* the prior and posterior then answer slightly different likelihoods, which needs
     a short derivation.
3. **Summarize CpG sites before fitting.** Average β values within a region (promoter, CpG
   island or gene body), or within clusters of strongly correlated sites. Every site still
   contributes, the feature counts become closer, and loadings read at the region level,
   which matches the gene–CpG question.
   - *Cost:* site-level detail inside a region is lost. Region definitions need a fixed
     annotation and must be chosen without outcomes.
4. **Calibrate a weight to a target share, combined with idea 1 or 2.** Choose w so that
   methylation's share of the L precision matches a target, such as its share under equal
   feature counts. This only works once the collapse is prevented.
5. **Add programs specific to one data type.** Use shared scores L plus scores L_m that
   belong to one data type, in the style of JIVE or of MOFA's factor-by-view sparsity.
   Methylation-specific variation then has its own place and does not have to take over the
   shared programs.
   - *Cost:* a new model component with its own derivation, rank choice and tests. This is
     the largest change.

**Ruled out.** Rescaling the methylation data, or fixing its noise precision at a smaller
value, acts like the power weight on both L and F, so it would collapse in the same way.
Cutting the CpGs to the size of the gene list is also out. Task C found that the
methylation-specific immune program forms only from about 50k CpGs.

## Suggested order when the work resumes

1. Idea 1 in the existing weight-path simulation, Blocks A and B, with the same rules.
2. If idea 1 keeps methylation loadings non-empty and balances the shares, fit it on TCGA at
   10k and 100k CpGs.
3. Otherwise, idea 3 with promoter and island summaries.
4. Idea 5 only if the shared-only model keeps forcing methylation-specific structure into
   joint programs.
