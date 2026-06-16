# Revision Notes: `presentation/walther_lab_meeting_06_18_2026`

These notes direct the revision process for the 6/18/2026 lab meeting slide deck.

---

## General Comments

- Include section headings to break up the sections of the slide deck — specifically for "Method" & "Results"
- For all figures, a proper label & caption are necessary. Most of these figures are pulled from other reports where we already have strong informative captions (we don't need the entire thing but that should serve as a basis for specifying the figure in the slide deck)
- In a few spots the "point-normal" prior is specified. This isn't correct. We'd settled on different priors (verify this). I believe we use exponential for L & F and normal for B (very very important)
- Should we specify the number of genes retained from the "top 3000"? I believe it was 2064 after merging the PDAC training cohorts
- The process for how we select "K" factors isn't specified anywhere. This is deserving of a slide to specify that we settled on cross-validation to select K. Further specify how we arrive at "K_eff" from "K" for the joint model (very very important)
- There is a borderline excessive amount of "call out" boxes used in the deck. The content can be included in-line with slides if necessary but we don't need to overly edit these in this pass
- Call a reviewer subagent to review the presentation of the report with attention to the following items:
  - Text & figures are legible, no excessive white space on slides where content can be shifted or expanded
  - Each figure has a label & descriptive caption
  - The takeaways of the presentation are clear and obvious to a technical audience
  - All terms, equations are clearly defined
  - No extraneous artifacts (like the `:::{.notes}` markdown artifact) are included in the rendered slide deck

---

## Slide-by-Slide Notes

### Slide 1
- Title looks good
- The subtitle can just be "Project Update"
- Change the date to "2026-06-18"

### Slide 2
- Heading "Agenda" is good
- Change "Validation" to "Results"
- Drop "(complete)" from this line

### Slide 3
- Write "gene expression programs (GEPs)"
- Drop "We will demonstrate this with results - not assert it"
- Shift the content in this slide further down the slide; currently the entire bottom half of the slide is just blank white space so we can expand the content on the slide to fill it a bit more

### Slide 4
- Looks like the feedback about specifying h0(t) (and how it cancels out) is clarified here
- Are the models & parameters fully defined here?
- I like the note about how a program is prognostic only if B isn't 0

### Slide 5
- The two linear predictors look good; specify that YFB is useful for exact risk scoring on new samples and that the LB predictor requires an approximation
- I believe the priors we use are the normal prior for beta, and exponential for L & F (non-negative). This was some prior feedback so make sure we get this right!
- Also make sure to clearly specify what the q() terms actually are
- There is a lot of white space at the bottom of this slide; the content can be expanded or shifted downward for better symmetry/balance

### Slide 6
- This slide is a bit unclear; need some text to specify what the q() equation is
- The boxes for the combined likelihood could be removed and just present the information for the likelihood terms on the slide (the additional callout boxes get excessive!)
- Clearly specify what each term of the ELBO is (and what the ELBO is for)
- You have an artifact `:::{.notes}` that shouldn't be in the rendered output
- Make sure we are specifying the likelihood/ELBO for the YFB model! (not the LB model)
- In general this slide presents 2 big equations without much in terms of definitions or context to explain the terms:
  - Why do we have the "q()" expression — this needs to be well defined and explained
  - The combined objective needs to be defined and its utility (ELBO is the objective function for the iterative updates) needs to be clear as well

### Slide 7
- Clearly specify what the convergence criteria we use is and explain it in plain English as well
- I like the note about the sign correction (matching it with what was found in training)
- Is the convergence criteria saying that the mean of all values in L and mean of all values in B must have a shift under the threshold? Is this the proper way to evaluate convergence?

### Slide 8
- The point about the within-platform survival precision being dominated is good! Make sure to define what the components of this expression "A_k = sum_i w_i(YF)^2_{ik}" are for clarity. This is important!

### Slide 9
- There is a lot of white space at the bottom of this slide; the content can be expanded or shifted downward for better symmetry/balance
- Maybe make the title "L Extension: cohort-membership indicator" so it is clear that this is about the loadings matrix
- Address how the addition of this indicator requires an additional degree of freedom so we select one more factor than we would without it
- Provide the "Result on real data" on the body of the slide without a callout box. This will prevent the text from being really small and drive home the recommendation. In fact, this could have a small header like "Recommendation" and then state the finding below it.

### Slide 10
- There is a lot of white space at the bottom of this slide; the content can be expanded or shifted downward for better symmetry/balance
- Drop "(completed)" from multi-cohort simulation
- Change title to "Data - Synthetic Validation & PDAC Cohorts"
- Change "3 scenarios x 5 model arms x 5 seeds" to two lines about the "Scenarios" and "Model Arms" (drop the "5 seeds"):
  - **Scenarios** (list out the 3 scenarios below)
    - "All Shared Factors"
    - "All Cohort-Specific Factors"
    - "Shared & Cohort-Specific Factors"
  - **Model Arms** (list out the model arms/characteristics below)
    - EBMF
    - YFB_base
    - YFB_cohort
    - LB_base
    - LB_cohort
- Change "Real PDAC - external validation design" to "PDAC cohorts - External Validation"

### Slide 11
- Change title to "Synthetic Results - Shared & Cohort-Specific Factor Recovery"
- Need to specify what these figures actually show! It looks like they show the average correlation of the true vs. recovered factors across 5 different model types for the 3 scenarios?
- Each figure needs a caption below it with a figure number like "Figure 1: figure description"
- What is the distinction between factor recovery (cor) & specificity-classification accuracy? Do we need both plots? The descriptions on the slide need to motivate the need for both plots!!!
  - Is the specificity plot for the performance of identifying the correct factor type?
- Expand the content to remove white space at the bottom of the slide

### Slide 12
- Shift slide content to remove excess white space
- Change title to "Synthetic Results - SSBMF outperforms EBMF+Cox"
- The EBMF+Cox method results should be included in this figure if we are stating a comparative result. Fix this with an update!
- Change "Hybrid (realistic) scenario, mean C-index:" to "Hybrid (Shared & Cohort-Specific), mean C-index:"
  - Is the "mean c-index" referring to the average over 5 seeds?
- Move the call-out box to another line across the slide that has "Takeaway" as a heading and then "Joint supervision beats unsupervised two-step approach on data with known ground truth"
- I think it would be really useful to show the by-cohort factor recovery heatmap as well. This really illustrates which recovered factors are shared and specific to each cohort (this figure might fit better context-wise on slide 11 but there is more space on slide 12 — you can add a slide between 11 and 12 as well to include this figure!):
  - The cohort recovery heatmap is in the appendix — a version is also included on page 8 of `multicohort_sim_proposal_06_14_26` that has great presentation too (this could be a standalone slide with title "Synthetic Results - Factor Recovery Heatmap")

### Slide 13
- Change title to "PDAC Cohorts - Prognostic Gene Programs"
- The text says "recommended model" — instead specify what the model actually is
- Is this figure showing that 7 factors are selected with 3 active being factor 3, 5, and 7?
  - Does the 3rd program not have any survival relationship? (it looks like 3 are active but you mentioned just 2 — this needs to be clarified)
  - And then the rows are specific genes that are active in each factor?
- The explanation of this figure needs to be more clear as I think it has a potential to be a strong result!!!
- Use a note/overlay/comment to specify which of the active factors are adverse & protective (this helps create continuity with the next slide)

### Slide 14
- Address excess white space at the bottom of the slide (expand/shift content)
- Expand the K-M plots if possible and make sure the p-values are more readable (larger in size)
- Change title: "PDAC Cohorts - Prognostic Association"
- Under the "Adverse program" and "Protective program" headings, include some of the genes included in those programs
- I don't understand what the "Always reported together:" text is for. This needs to be clarified or dropped (but "Always reported together:" definitely doesn't belong in the slides)
- Give the K-M plots a proper figure caption and make the figure labels/legend larger if possible

### Slide 15
- Change title to "PDAC Cohorts - External C-index"
- Give the figure a proper figure caption
- The figure looks a bit stretched out width-wise and smushed top to bottom. There is plenty of white space on the slide to plot the figure at its proper aspect ratio (and it can be enlarged)
- Change "Projection model mean" to "Projection (YFB) model mean" and "loadings model" to "loadings (LB) model"
- The following line needs to be reworked: "Real-data analog of the synthetic result: the supervised joint model is competitive-to-better than the unsupervised two-step". Maybe something like "The projection-based supervised joint model outperforms the unsupervised 2-step approach across all 5 external PDAC cohorts"

### Slide 16
- Make the title "Model Comparison Summary" and place it at the top of the slide
- Expand the content to fill the slide
- You specify "DeSurv Selection" — this needs to be defined (it is the preprocessing method we use):
  - It's better to specify the gene selection process instead of referencing "DeSurv" — leave out references to "DeSurv" for now and just state the work we did
  - You could make a small note to specify that the selection process was motivated by the selection process for the DeSurv method in a call out box to the side
- We have to explain how we get from "K" to "K_eff". Maybe this belongs earlier in the method section. Yes — we definitely need to explain how "K" is selected given that we now specified it is chosen by cross-validation

### Slide 17
- The findings don't need to be included in a call out box. They can just be numbered as lines on the slide (or bullet points)
- When referencing the DeSurv approach, reference some results to compare it quantitatively if possible (otherwise omit the comparison comment — probably better to omit specifically stating the comparison to DeSurv for now)

### Slide 18
- Pathway enrichment is good
- I'm not sure what the "Finalize the EBMF / DeSurv head-to-head" line is for
- The manuscript preparation line is good. Maybe just say "Manuscript preparation" though
- Here is a good spot where you could mention that we addressed items from the "Next Steps" of the last presentation. Those were:
  - CV-based K selection (complete)
  - Multi-dataset shared factors (complete)
  - Full ELBO as model selection (complete)
    - Right? We do our updates based on the ELBO now?
- **Summary:**
  - Verify what priors we actually use (I don't think it's point-normal)
  - I don't like the wording of "Per-platform normalization is load-bearing - the p >> n fix". This should be reworked to specify that the pre-processing steps we implemented are crucial for model performance (and avoiding the precision term domination if I recall right)
  - "Validated on a completed multi-cohort simulation & 5 external PDAC cohorts" could include some C-index values to support the claim
  - "Recommended model" could be broken out onto multiple lines like this:
    - **Recommended model**
      - Linear Predictor: projection (YF)B
      - L Cohort Indicator: No
      - Pre-processing: per-platform z-std
      - Selected K: 7 (K_eff = 2), mean external C = 0.636

### Slide 19
- This slide is fine

---

## Appendix Notes

- Determine if these slides are truly necessary
- It looks like the appendix is largely composed of the derivation of the L/F/B/T updates as well as additional results figures:
  - I don't think the derivations are necessary
  - The prior/hyperparameter table isn't necessary. Anything that is necessary in this slide needs to already be specified & defined in the methods section (check this!)
  - "which preprocessing pipelines collapse" can be dropped (this isn't really relevant to the talk)
  - The results slides (last 3) can stay but each figure needs a label & caption with enough supporting context to make the presentation of the figures valuable. Additional relevant figures for the simulation study and/or PDAC application results can be included too
  - In the last appendix slide, the loading heatmap should be moved to the body of the slide deck (this has already been made clear)
- Review the appendix slides of the 4/9/26 presentation deck and determine if any similar or other relevant figures are present for the most up-to-date results with the multi-cohort simulation study & the PDAC data application
