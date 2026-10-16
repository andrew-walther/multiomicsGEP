# Corrected complete-cycle diagnostic

Task2-7 ran once after preregistrationc94c182 and full suite694/0.
Exact entire-fit/history identity, unchanged input/source hashes and complete
last50-cycle density coverage pass. Maximum ledger/native moment gap is
2.22e-16; zero warnings; runtime42.867s. The original fit still caps300.
Task2-6 remains invalid and is preserved separately.

The corrected Gaussian expectation retains between-program mean products.
Finite last50 reference-objective comparisons at the registered1e-8
descriptive tolerance give:

| Transition | Decreases / finite comparisons |
|---|---:|
| F posterior, fixed incoming prior/working terms |0/350|
| L posterior, fixed newly fitted prior/working terms |0/350|
| beta posterior, fixed newly fitted prior/working terms |0/350|
| F prior refit, posterior fixed |50/350|
| Cox working refresh, posterior fixed |18/50|
| Canonicalization with transformed reference priors |0/50|
| Noise-precision update |0/50|

The smallest late F-prior change is−.0253224; smallest refresh change
is−.00206559. This does not support blaming a defective fixed-prior F
coordinate step in this simulation, or treating refresh as the only source
of non-ascent. It supports checking the loading-prior hyperparameter step.

The2100 F-prior records match the earlier exact Task2-5 audit by iteration
and factor. On2074 finite comparisons, prior-density changes agree within
3.64e-12; that earlier audit's native/analytic slab gap was6.66e-16. This
reuses the earlier probability validation on the identical saved trajectory;
it is not a newly measured native-probability bound for every L/beta density.
Moment identity alone does not establish density/entropy accuracy in all
extreme states. These limits matter for interpreting the joint L/beta losses.

L/beta marginal-EB fits change prior and posterior jointly. Their intermediate
fixed-q prior losses are not alone evidence of an incorrect joint update.
The additional joint-block summary records17/350 negative L differences
(minimum−1.09106) and16/350 beta differences (minimum−1.13e-8, near the
descriptive tolerance). Optimizer accuracy, reconstructed density accuracy
and reference-prior interpretation need scrutiny before attributing these
to a new model/update defect. Do not omit these observations.

L/beta priors transformed for bookkeeping are mathematical references, not
implemented carried priors. Infinite support transitions remain explicit;
they cannot give finite ascent claims. The expected local quadratic is not
exact E_q logPL and is not a global lower bound. No common implemented ELBO,
convergent remedy or causal explanation of multimodal/M3 instability is
established. No default, threshold, scientific gate or external score changed.

Recommendation for Andrew's next decision: derive/review an opt-in fixed-q
point-Laplace F-prior EM step (pi=S/p, rate=S/T) with its corresponding KL,
rather than select a damping value or survival weight. A new implementation
or panel needs its own approval and prospective registration. Stop here.
