# Task2-9: valid bounded case, unresolved output stability

Executed once after preregistration commit02a658a. Exit0,2.659s, no warnings;
reference and all recorded source hashes unchanged. Detailed fit/telemetry
remain ignored; only compact aggregates and provenance are versioned.

The candidate uses fixed_q_em with original controls and reaches the native
ELBO rule at10 iterations (relative change9.997123e-6 versus1e-5). Its existing
output diagnostic remains false: max relative reconstruction/predictor change
.01914876 versus1e-4. Reference:300 cap, ELBO change.000172121, output change
.105693. Neither result supports saying that fitted outputs have settled.
Changing KL can change where an existing stopping rule triggers; no stopping
rule was changed or reclassified to make the result pass.

Both recover the planted direction (|r|.89463 candidate/.89564 reference),
with2 credible programs and1 false-active program. Simulated validation C is
.756846/.753722. One seed cannot establish a gain, recovery rate, absence of
runaway or adoption. Terminal loading norms.948/.948 and activity norms
88.87/85.39 are contextual, not a trajectory guarantee.

All67 actual F-update helper calls have finite KL and nonnegative prior gains
at the registered1e-8 descriptive tolerance. The fitter's existing dead-program
shortcut is not a helper call; no omitted helper result is silently discarded.
Four predefined informative scalar cases have max absolute mean error5.55e-17,
second-moment error1.39e-17 and normalizer error7.06e-14 versus independent
integration. These are selected checks, not certification of all states or
resolution of the independently demonstrated tail-kernel approximation.

The candidate monitor has3 negative increments. A local prior repair therefore
does not establish global ascent; L/beta prior bookkeeping and the changing
Cox surrogate remain distinct uncertainties. No real/external data or new
criterion was used. The benchmark is unchanged and no configuration adopted.

Recommendation: audit the stopping monitor against output movement and the
intended objective before a larger remedy panel. A new observation/remedy
requires its own prospective design; no replay/longer cap or new stopping
rule is authorized automatically. This registered case is terminal.
