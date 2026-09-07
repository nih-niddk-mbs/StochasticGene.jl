# Version 2.0

## Version 2.0.5

Version 2.0.5 restores one coordinate convention across the legacy
`hierarchical=(2, fittedindividual, fixedeffects)` trace stack. The scientific
likelihood now receives only physical rates and noise parameters. Its first
hierarchical block is the vector of population means; its second block is the
vector of positive population coefficients of variation.

For a positive individual parameter, the hierarchy conditional is log-normal
and is parameterized by its arithmetic mean and CV. For a signed parameter,
such as a zero-median Gaussian noise mean, the conditional is normal. Every CV
is positive even when its corresponding individual parameter is signed. The
sampler still operates in unconstrained coordinates, and the individual
inverse-transform Jacobian is included exactly once in the posterior density.

The HMM layer now calculates only observation likelihoods. The trace likelihood
adds the hierarchy density above that layer, including for grid models, and
continues to return data-only pointwise terms for WAIC.

This is a correctness change relative to versions 1.5.0 through 2.0.4, which
could compare physical hyperparameters with transformed individual parameters.
Old rate and info files remain readable, but continuing one targets the
corrected posterior and may move away from the previous result.

## Version 2.0.4

Version 2.0.4 exposes the ON/OFF dwell-time CDFs that were already calculated
internally before deriving the normalized PDFs. `write_ONOFFhistograms` and
`write_ONOFFhistograms_key` now write `ON_CDF` and `OFF_CDF` alongside the
existing `ON` and `OFF` columns. Simulation output also includes `SimON_CDF`
and `SimOFF_CDF`.

Key-based and legacy folder generators accept `bin_size` and `maxtime` in
minutes. Use `bin_size=100/60` to match data acquired every 100 seconds. An
explicit `bins` vector remains supported and takes precedence.

## Version 2.0.3

Version 2.0.3 brings coupled simulation into parity with the coupled transition
matrix used by fitting and likelihood evaluation. The simulator now handles
connections to different target transitions independently, reciprocal
coupling, hidden three-unit models, and repeated entries in `unit_model`.
Stationary-state indexing works for any number of units, and coupled trace
simulation preserves `a_grid`.

The release also makes the two reporter-state conventions explicit:

- `Rany` uses one sentinel source state and multiplies the target transition by
  `1 + gamma` whenever at least one reporter position is occupied.
- `Rsum` uses one connection per reporter position. If `n` positions are
  occupied, a tied strength gives the factor `1 + n*gamma`.

The historical shorthand `make_coupling("R5", G, R)` expands to `Rsum`, not
`Rany`. In direct simulator calls, the rate vector contains one strength per
expanded connection; repeat a tied fitted gamma for each R-position connection.

Malformed unit-model maps, rate/noise layouts, and coupling sign metadata now
raise `ArgumentError` rather than being silently misinterpreted. These changes
can alter coupled simulated trajectories but do not change fitted rate-file
formats.

## Version 2.0.2

Version 2.0.2 corrects the parameter domain for `Rsum` coupled models. An
`Rsum` coupling shared by `m` R positions contributes through
`1 + γR₁ + ... + γRₘ`. Because those R positions can be occupied
simultaneously, the physical lower bound is `γ > -1/m`, not `γ > -1` for each
expanded connection independently.

The real-line parameter transform now detects tied R-position connections and
uses `(-1/m, 0)` for inhibitory `Rsum` couplings and `(-1/m, Inf)` for free
`Rsum` couplings. This applies to ordinary and hierarchical fits. Default CSV
priors are also kept strictly inside the appropriate interval. `Rany` models
and tied couplings between mutually exclusive source states retain the usual
lower bound of `-1`.

Historical `Rsum` fits are not modified on disk. Before using an older rate
file as an initial condition, check every tied group and migrate values for
which `1 + mγ <= 0` into the open valid interval. For new CSV workflows, use
explicit `Rsumk` or `Ranyk` tokens when the distinction matters.

## Version 2.0.1

Version 2.0.1 is a correctness patch for coupled-model simulation and
correlation analysis.

- Centered theoretical correlations now subtract the correct stationary means
  and approach zero at long lag.
- Coupled simulation now handles multiple active coupling connections to the
  same target transition correctly.
- Physical-time equilibration follows any requested event-count warmup and
  ends at a time cutoff rather than a reaction epoch.
- Sparse simulated traces retain the requested recording endpoint.
- `simulate_trials(...; nexperiments=N)` runs independent experiments in
  parallel and reports ON-ON correlation intervals across experiments.

These corrections can change coupled simulated traces and theoretical
correlation outputs relative to version 2.0.0. Fitted rate-file formats and the
public fitting interface are unchanged.

StochasticGene 2.0 is the first stable release of the API tested on the 1.11
beta line. It requires Julia 1.11.

## Highlights

- `fit` supports Metropolis-Hastings, NUTS, and ADVI through a common keyword
  interface and run-spec representation.
- Recursive/shared-parameter fits keep fitted parameters separate from the
  complete rate vectors consumed by likelihoods, with transition and emission
  cache groups for repeated HMM computations.
- Key-based `info_<key>.jld2` workflows support reproducible continuation,
  analysis, and scheduler command generation.
- RNA workflows support explicit legacy and non-truncated likelihood modes,
  metadata-aware decay and allele handling, genome-scale `makeswarm_genes`,
  and dataframe assembly.
- Trace analysis includes key-aware writers, generalized state-observable
  correlations with reusable HMM contexts, and posterior burst prediction.
- Coupled and shared analyses preserve complete per-unit or per-group rate
  outputs while fitting only the selected parameter representation.

## Migration Notes

- Use `stage_write_run_specs` plus `makeswarm`, or the CSV command-file helpers,
  for key-based batch jobs. `makeswarm_genes` remains the high-level RNA gene
  batch entry point.
- Prefer `trace_specs` and `dwell_specs` over legacy `traceinfo` and `dttype`
  metadata in new scripts.
- Prefer tuple/vector `datatype` and modality-keyed `datapath` for multimodal
  fits; legacy combined datatype strings remain available for compatibility.
- The terminal summary now labels the likelihood evaluated at coordinate-wise
  median rates explicitly. It is not the median of sampled likelihood values.
- Development-only notebook and documentation packages are no longer runtime
  dependencies. Documentation dependencies remain isolated in `docs/`.

## Release Checks

The release workflow checks that the package precompiles on Julia 1.11, all
exports resolve, the standard test suite passes, public exports have static
docstring coverage, and the Documenter manual builds without warnings.
