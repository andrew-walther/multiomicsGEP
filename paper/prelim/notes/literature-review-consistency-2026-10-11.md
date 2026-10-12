# Chapter 4 consistency with the revised literature review

Author-approved wording and reference corrections on 2026-10-11. The canonical
chapter uses the same L/F convention as the literature review. Section 1.3.3
distinguishes possible coefficient bias from uncertainty omitted by treating
estimated scores as fixed. The multimodal derivation carries loading variance
into the quadratic term of the coefficient update. This does not establish an
attenuation correction. Integration methods are reviewed in Section 1.3.7.

Only prose and a section pointer changed. Model equations, fitting procedures
and results are unchanged. Replaced material is preserved below.

## Introduction: distinguish possible bias from omitted score uncertainty

Date: 2026-10-11. Reason: Introduction: distinguish possible bias from omitted score uncertainty.

```markdown
Treating the estimated scores as fixed Cox covariates can attenuate their associations (Section 1.3.3).
```

## Integration pointer: the reviewed methods are in Section 1.3.7

Date: 2026-10-11. Reason: Integration pointer: the reviewed methods are in Section 1.3.7.

```markdown
Existing integration methods, such as MOFA [@Argelaguet2018] and iCluster [@Shen2009], combine cohorts or data types without survival supervision (Section 1.3.6).
```

## Notation: both chapters already use L for patient scores and F for feature loadings

Date: 2026-10-11. Reason: Notation: both chapters already use L for patient scores and F for feature loadings.

```markdown
^[Section 1.3.2 of the literature review writes $X\approx LF$ with $L$ the gene loadings; this chapter transposes that convention to match the EBMF software, so here $L$ holds patient scores.]
```

## Second moment: describe loading-uncertainty propagation without an attenuation-correction claim

Date: 2026-10-11. Reason: Second moment: describe loading-uncertainty propagation without an attenuation-correction claim.

```markdown
This propagates loading uncertainty into the survival coefficients, the correction for attenuation motivated in literature review Section 1.3.3.
```
