# Version 2.0

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
