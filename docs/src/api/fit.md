# fit Function

!!! note "Authoritative reference"
    The full keyword list, coupling description, `key`-based loading, and return values are maintained in the **in-source docstring** for `fit`. In the Julia REPL, use `?fit` after `using StochasticGene`. This page is a **short** summary; details may lag the code.

`fit` fits steady-state or transient GM/GRSM models to data for a single gene (or coupled units), writes results through `finalize`, and returns fit results and diagnostics.

For coupled transcribing units, arguments `transitions`, `G`, `R`, `S`, `insertstep`, and trace-related settings become **tuples** of the single-unit type (e.g. two units with `G = (2, 3)`).

## Syntax

```julia
fits, stats, measures, data, model, options = fit(; kwargs...)
```

The keyword form is the usual entry point. Positional overloads exist for advanced / legacy callers; optional inference keywords on the positional path are forwarded into `make_structures` (see the `_MAKE_STRUCTURES_OPTION_KW` constant in `fit.jl`).

## Inference methods

Posterior / variational inference is selected with **`inference_method`** (also accepted: plain symbols `:mh`, `:nuts`, `:advi`, or the aliases `INFERENCE_MH`, `INFERENCE_NUTS`, `INFERENCE_ADVI`):

| Method | Meaning |
|--------|---------|
| **`:mh` / `INFERENCE_MH`** | Metropolis–Hastings MCMC (`run_mh` / `metropolis_hastings`); default. |
| **`:nuts` / `INFERENCE_NUTS`** | NUTS HMC on the same transformed parameter space as MH (`run_nuts_fit`, AdvancedHMC). |
| **`:advi` / `INFERENCE_ADVI`** | Mean-field ADVI (`run_advi_fit`); returns a variational approximation (see notes below). |

Shared **budget** keywords are harmonized in `load_options` when building option structs:

- **`samplesteps`**: MH — posterior samples to collect; NUTS — `n_samples`; ADVI — `maxiter` unless you set **`maxiter`** explicitly in the run dict.
- **`warmupsteps`**: MH — discarded warmup (shares wall time with sampling via **`maxtime`**); NUTS — `n_adapts` unless **`n_adapts`** is set explicitly; ADVI — not the same semantic object (use **`n_mc`** for ELBO Monte Carlo draws).

Other cross-method options stored on `MHOptions`, `NUTSOptions`, and `ADVIOptions` in `common.jl`:

- **`device`**: `:cpu` or `:gpu` (GPU paths may error if unsupported for a method).
- **`parallel`** (alias **`parallelism`**): `:single`, `:threaded`, or `:distributed` — used with **`nchains`** for multi-chain NUTS/ADVI dispatch (`run_inference`); MH multi-chain still uses the existing distributed MH chain runner when `nchains > 1`.
- **`gradient`**: method-specific; e.g. `:finite`, `:ForwardDiff`, `:Zygote` for NUTS/ADVI; `:none` default for MH. String values in TOML/JSON-style dicts are coerced when possible.

`make_structures` merges **`_current_run_spec[]`** (from `fit(; key=...)`) with explicit `samplesteps` / `warmupsteps` / `maxtime` / `temp` and any extra `kwargs...`, then calls `load_options` on that dict. You can also call `load_options` directly in scripts or tests.

## Arguments (high level)

### Basic model parameters

- **`G::Int`**, **`R::Int`**, **`S::Int`**, **`insertstep::Int`**, **`transitions`**: Model topology (see in-source `fit` docstring).
- **`coupling`**, **`grid`**, **`hierarchical`**: Coupled / grid / hierarchical layouts.

### Data parameters

- **`datatype`**: String or symbol for a single legacy datatype, e.g. `"rna"`, `:trace`, `"rnadwelltime"`, `"tracejoint"`. A tuple or vector such as `(:rna, :dwelltime)` requests the v1.10 [`CombinedData`](combined_data.md) path.
- **`datapath`**: Path, vector/tuple of paths, or for `CombinedData` preferably a modality-keyed `NamedTuple`, e.g. `(rna = "smFISH", dwelltime = ["ON.csv", "OFF.csv"])`.
- **`datacond`**, **`cell`**, **`gene`**, **`nalleles`**, **`trace_specs`**, **`dwell_specs`**, …

### Fitting / inference parameters

- **`nchains::Int`**: Number of parallel chains for **`fit(nchains, data, model, …)`** dispatch (`run_inference`); for MH this matches the existing pooled-chain behavior; for NUTS/ADVI, multi-chain runs are merged when `nchains > 1` (see `src/inference_common.jl` in the package source).
- **`maxtime`**: Primary wall budget for **MH** (warmup + sampling); numeric seconds or strings like `"90m"`, `"2h"` (see `maxtime_seconds`). Interpretation for NUTS/ADVI is method-specific (many NUTS/ADVI controls are in the option structs from `load_options`).
- **`samplesteps`**, **`warmupsteps`**, **`temp`**: As above; **`temp`** is MH temperature; NUTS/ADVI paths use a neutral value for finalize when not applicable.
- **`propcv`**: MH proposal CV / covariance loading (see [Package overview](../package_overview.md#MCMC-proposal-covariance-and-warmup)).

### Priors, indices, outputs

- **`priormean`**, **`priorcv`**, **`noisepriors`**, **`fittedparam`**, **`fixedeffects`**, **`onstates`**, **`decayrate`**, …
- **`resultfolder`**, **`label`**, **`writesamples`**, **`burst`**, **`optimize`**, …

### v1.10 API changes

- **`CombinedData`**: New multimodal fits use `datatype = (:rna, :trace)` or `datatype = (:rna, :dwelltime)`. The order is canonicalized, likelihoods are evaluated per modality, scalar likelihoods are summed, and WAIC pointwise predictions are concatenated in canonical modality order.
- **Retired legacy input keywords**: `infolder` and `inlabel` are no longer part of the public API. Use `root`, `datapath`, `label`, and `resultfolder`. Older key-based run specs may still contain `infolder`; `fit(; key=...)` ignores it during migration.
- **Trace and dwell metadata**: Prefer `trace_specs` and `dwell_specs`. Legacy `traceinfo` and `dttype` entries may be consumed from old run specs, then are dropped from newly written run specs.

### Run specification and key-based naming

- **`key = nothing`**: When a string, `fit` loads `info_<key>.jld2` (companion to the marker TOML) and merges with explicit keywords (keywords win). Outputs use that stem. See [Run specification (info TOML)](../run_spec_toml.md).

## Returns

The top-level **`fit(; …)`** return tuple is:

```julia
fits, stats, measures, data, model, options = fit(; ...)
```

- **`fits`**: `Fit` — posterior samples (MH/NUTS) or variational mean column (ADVI); **`ll`**, WAIC-related fields, acceptance summaries depend on the method.
- **`stats`**, **`measures`**: `Stats`, `Measures` — comparable structs across methods where meaningful (ADVI may use single-point proxies for some diagnostics; see the `run_advi_fit` docstring in the source).
- **`data`**, **`model`**, **`options`**: The concrete `Options` subtype is **`MHOptions`**, **`NUTSOptions`**, or **`ADVIOptions`** depending on `inference_method`.

## Examples

### Basic RNA histogram (MH default)

```julia
fits, stats, measures, data, model, options = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/HCT116_testdata/",
    gene = "MYC",
    datacond = "MOCK",
)
```

### NUTS (same keyword surface; different `options` type)

```julia
fits, stats, measures, data, model, options = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/HCT116_testdata/",
    gene = "MYC",
    datacond = "MOCK",
    inference_method = :nuts,
    samplesteps = 500,
    warmupsteps = 250,
    parallel = :single,
    gradient = :ForwardDiff,
)
```

### Combined RNA + dwell-time data

```julia
fits, stats, measures, data, model, options = fit(
    datatype = (:rna, :dwelltime),
    datapath = (
        rna = "HBEC_smFISH",
        dwelltime = [
            "dwelltime/CANX_ON.csv",
            "dwelltime/CANX_OFF.csv",
            "dwelltime/CANX_ONG.csv",
        ],
    ),
    gene = "CANX",
    cell = "HBEC",
    datacond = "",
    resultfolder = "FISH",
    dwell_specs = [
        (
            unit = 1,
            onstates = [Int[], Int[], [2, 3]],
            dttype = ["ON", "OFF", "ONG"],
        ),
    ],
)
```

For new scripts, the keyed `datapath` form is preferred over the legacy positional `"rnadwelltime"` layout. See [v1.10 CombinedData API](combined_data.md).

### Trace fit with multiple MH chains

```julia
fits, stats, measures, data, model, options = fit(
    G = 3,
    R = 2,
    S = 2,
    insertstep = 1,
    transitions = ([1,2], [2,1], [2,3], [3,1]),
    datatype = "trace",
    datapath = "data/testtraces",
    cell = "TEST",
    gene = "test",
    datacond = "testtrace",
    trace_specs = [(unit = 1, interval = 1.0, start = 1.0, t_end = -1.0, zeromedian = true, active_fraction = 1.0, background = 0.0)],
    noisepriors = [40., 20., 200., 10.],
    nchains = 4,
)
```

## Notes

1. **Rate order** — G transitions, R transitions, S if present, decay, noise parameters (see package overview).

2. **Convergence** — R-hat, ESS, etc. apply to sample-based chains; interpret ADVI diagnostics separately.

3. **MH proposal covariance and warmup** — See [Package overview](../package_overview.md#MCMC-proposal-covariance-and-warmup).

4. **Key-based workflows** — Use `key="..."` for reproducible cluster runs and `write_run_spec_preset` / `makeswarmfiles`; include `inference_method` / `parallel` / `gradient` in the saved dict when needed.

5. **Cluster scripts** — `makeswarm` and related helpers can emit `fit(; key=..., inference_method=:nuts, …)` overrides; positional gene/coupled scripts append `; kw=...` suffixes for the same keywords. See [Cluster and batch workflows](../cluster_batch_workflows.md).

6. **Migration from older inputs** — Replace `infolder` with `datapath`; replace `inlabel` usage with `label`; keep output routing in `resultfolder`. Prefer `trace_specs` / `dwell_specs` over `traceinfo` / `dttype`.
