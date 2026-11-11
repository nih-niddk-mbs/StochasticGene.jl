# First-Passage Numerical Validation

Numerical corrections to first-passage propagation are separate from changes
to the observation model. Propagation starts with the entrance distribution at
physical time zero, even when the first reported bin has a positive left edge.
Both solver paths return plain time-by-state matrices. Eigen propagation uses
matrix exponentials when the eigenvector basis is ill-conditioned. Default ODE
tolerances are `reltol=1e-8`, `abstol=1e-11`; low-level `time_evolve_diff`
callers can override them.

The previous ODE solution-adjoint wrapper used unusual slice orientation:
although its size was time-by-state, selecting states could return a
state-by-time matrix. Its `dims=1` reduction was therefore not automatically
wrong. The new plain-matrix representation and `dims=2` must be considered
together. Existing successful zero-left-edge simulation comparisons remain
meaningful; they do not test a positive cutoff with a nonexponential dwell law.

`pdf_from_cdf` retains the historical window-conditioned PMF. It also handles
decreasing survival probabilities from a killed generator. The entrance flux
is normalized after subsetting. Histogram likelihoods reject invalid
probabilities, preserve logarithms of tiny positive values, and give zero
contribution to an empty observed bin with zero model probability.
Small negative roundoff down to `-1e-10` is tolerated only in empty bins,
consistent with the ODE absolute tolerance. An occupied bin with a nonpositive
prediction still makes the likelihood invalid; it is not epsilon-clipped.

Positive bin right edges are required. Rejecting a zero-start, offset, or
nonuniform grid in the simulator is an intentional compatibility change:
the existing simulator bins durations using `ceil(duration/dt)` and does not
implement those other grids. Do not relabel experimental data automatically.

Run targeted checks with:

```julia
using StochasticGene
include("test/first_passage_tests.jl")
```

These tests include analytic exponential/Erlang laws, offset/nonuniform
prediction bins, dense/sparse matrices, zero-time propagation, entrance
normalization, and histogram likelihood/gradient checks. They are not a claim
that historical posterior estimates are unchanged. Compare saved predictions
and refit representative datasets before deciding on a broader migration.

Truncated, complete, and censored data require separate observation semantics.
These corrections do not infer censoring or a sample size from a normalized curve.
CDF subtraction near one can still lose rare-tail precision; use survival
differences when extreme tail accuracy is required.
