# fit Function

!!! note "Authoritative reference"
    The full keyword list, coupling description, `key`-based loading, and return values are maintained in the **in-source docstring** for [`fit`](@ref). In the Julia REPL, use `?fit` after `using StochasticGene`. This page is a **short** summary; details may lag the code.

Fit steady state or transient GM/GRSM model to RNA data for a single gene, write the result (through function finalize), and return fit results and diagnostics.

For coupled transcribing units, arguments transitions, G, R, S, insertstep, and trace become tuples of the single unit type, e.g. If two types of transcription models are desired with G=2 and G=3 then G = (2,3).

## Syntax

```julia
fits = fit(; kwargs...)
```

## Arguments

### Basic Model Parameters

- `G::Int = 2`: Number of gene states
- `R::Int = 0`: Number of pre-RNA steps
- `S::Int = 0`: Number of splice sites (must be ≤ R - insertstep + 1)
- `insertstep::Int = 1`: Reporter insertion step (must be ≥ 1; ignored when R = 0)
- `transitions::Tuple`: Tuple of vectors specifying state transitions

### Data Parameters
- `datatype::String = "rna"`: Data type:
  - "rna": mRNA count distributions
  - "trace": Intensity traces
  - "rnadwelltime": Combined RNA and dwell time data
  - "tracejoint": Simultaneous recorded traces
- `datapath::String = ""`: Path to data file or folder
- `datacond::String = ""`: Data condition identifier
- `cell::String = ""`: Cell type
- `gene::String = "MYC"`: Gene name
- `nalleles::Int = 1`: Number of alleles

### Fitting Parameters

- `nchains::Int = 2`: Number of parallel chains for **Metropolis–Hastings** (`inference=:mh`). For **NUTS** (`inference=:nuts`), must be `1` (single chain).
- `inference = :mh`: Posterior algorithm — [`INFERENCE_MH`](@ref) (MH, default), [`INFERENCE_NUTS`](@ref) (NUTS via AdvancedHMC; use `samplesteps` / `warmupsteps` as `n_samples` / `n_adapts`), or [`INFERENCE_ADVI`](@ref) (not wired through `fit`; call [`run_advi`](@ref)).
- `steady_state_solver = :augmented`: Passed to the likelihood when using NUTS.
- `ad_likelihood = nothing`: Passed to NUTS (`nothing` selects AD likelihood for RNA count data when appropriate).
- `maxtime = 60`: Maximum **total** wall time for the **MH** phase, including **warmup and sampling** together. A **numeric** value is **seconds**; you may also pass a **string** with suffix `m` (minutes) or `h` (hours), e.g. `"90m"`, `"2h"`. Set `samplesteps` large and use `maxtime` as the primary stop (e.g. cluster time limits). See [`maxtime_seconds`](@ref). (NUTS stopping is controlled by `samplesteps` / `warmupsteps`; `maxtime` is not applied to NUTS in the current `fit` wiring.)
- `samplesteps::Int = 1000000`: MH: max sampling steps (may stop earlier at `maxtime`). NUTS: `n_samples`.
- `warmupsteps::Int = 0`: MH: discarded warmup (shared `maxtime` with sampling). NUTS: `n_adapts`.
- `propcv::Float64 = 0.01`: Proposal distribution coefficient of variation
- `temp::Float64 = 1.0`: MCMC temperature

### Prior Parameters

- `priormean::Vector{Float64} = Float64[]`: Mean rates for prior distribution
- `priorcv::Float64 = 10.0`: Prior distribution coefficient of variation
- `noisepriors::Vector = []`: Observation noise priors
- `fittedparams::Vector{Int} = Int[]`: Indices of rates to fit
- `fixedeffects::Tuple = tuple()`: Fixed effects specification

### Trace Parameters (for trace data)

- `traceinfo::Tuple = (1.0, 1., -1, 1.)`: Trace parameters:
  - Frame interval (minutes)
  - Starting frame time (minutes)
  - Ending frame time (-1 for last)
  - Fraction of active traces
- `datacol::Int = 3`: Data column index
- `probfn::Function = prob_Gaussian`: Observation probability function
- `noiseparams::Int = 4`: Number of noise parameters
- `zeromedian::Bool = true`: Subtract median from traces

### Output Parameters

- `resultfolder::String = "test"`: Results output folder
- `label::String = ""`: Output file label
- `infolder::String = ""`: Folder for initial parameters
- `inlabel::String = ""`: Label of initial parameter files
- `writesamples::Bool = false`: Write MCMC samples

### Run specification and key-based naming

- `key = nothing`: When nothing, fit uses the keyword arguments you pass (and defaults). When a string (e.g. `key = "33il"`), fit looks for `info_<key>.toml` in the results folder; if found, loads that spec (kwargs override spec); if not found, uses kwargs and defaults. All outputs use that stem (e.g. `rates_<key>.txt`, `info_<key>.toml`). See [Run specification (info TOML)](@ref).

## Returns

- `fits`: MCMC fit results containing:
  - Posterior samples
  - Log-likelihoods
  - Acceptance rates

## Examples

### Basic RNA Histogram Fit

```julia
fits = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/HCT116_testdata/",
    gene = "MYC",
    datacond = "MOCK"
)
```

### Trace Data Fit

```julia
fits = fit(
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
    traceinfo = (1.0, 1., -1, 1.),
    noisepriors = [40., 20., 200., 10.],
    nchains = 4
)
```

## Notes

1. **Rate Order**
   - G transitions
   - R transitions
   - S transitions
   - Decay
   - Noise parameters

2. **MCMC Convergence**
   - R-hat should be close to 1 (ideally < 1.05)
   - Increase `maxtime` or `samplesteps` if R-hat is high
   - Use `warmupsteps` to improve proposal distribution
