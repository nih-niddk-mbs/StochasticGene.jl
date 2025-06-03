# fit Function

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
- `insertstep::Int = 1`: Reporter insertion step
- `transitions::Tuple`: Tuple of vectors specifying state transitions

### Data Parameters
- `datatype`: Data type (supported: :rna, :rnacount, :rnaonoff, :rnadwelltime, :dwelltime, :trace, :tracerna, :tracejoint, :tracegrid)
- `datapath`: Path to data files
- `root`: Root directory containing data and results folders
- `gene`: Gene name(s)
- `cell`: Cell type
- `decayrate`: mRNA decay rate (if < 0, will be looked up)

### Coupled Models
- `coupling`: Coupling structure between units
- `nalleles`: Number of alleles (automatically set to 1 if coupling is nonempty)

### Transition Types
- `TransitionType`: String specifying transition pattern (e.g., "KP", "cyclic")

## Returns

- Model parameters
- Posterior distributions
- Predicted data
- Goodness of fit metrics
- Model diagnostics
- Analysis results

## Notes

- For coupled models, many parameters become tuples (e.g., G = (2,3) for two units with 2 and 3 gene states)
- The number of splice sites S cannot exceed R - insertstep + 1
- If coupling is used, nalleles is automatically set to 1
- Default transitions can be generated using the `get_transitions` function

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

- `nchains::Int = 2`: Number of MCMC chains
- `maxtime::Float64 = 60`: Maximum wall time (minutes)
- `samplesteps::Int = 1000000`: Number of MCMC sampling steps
- `warmupsteps::Int = 0`: Number of warmup steps
- `propcv::Float64 = 0.01`: Proposal distribution coefficient of variation
- `annealsteps::Int = 0`: Number of annealing steps
- `temp::Float64 = 1.0`: MCMC temperature

### Prior Parameters

- `priormean::Vector{Float64} = Float64[]`: Mean rates for prior distribution
- `priorcv::Float64 = 10.0`: Prior distribution coefficient of variation
- `noisepriors::Vector = []`: Observation noise priors
- `fittedparams::Vector{Int} = Int[]`: Indices of rates to fit
- `fixedeffects::Tuple = tuple()`: Fixed effects specification
- `hierarchical::Tuple = tuple()`: Hierarchical model specification

### Trace Parameters (for trace data)

- `traceinfo::Tuple = (1.0, 1., -1, 1.)`: Trace parameters:
  - Frame interval (minutes)
  - Starting frame time (minutes)
  - Ending frame time (-1 for last)
  - Fraction of active traces
- `datacol::Int = 3`: Data column index
- `probfn::Function = prob_Gaussian`: Observation probability function
- `noiseparams::Int = 4`: Number of noise parameters
- `zeromedian::Bool = false`: Subtract median from traces

### Output Parameters

- `resultfolder::String = "test"`: Results output folder
- `label::String = ""`: Output file label
- `infolder::String = ""`: Folder for initial parameters
- `inlabel::String = ""`: Label of initial parameter files
- `writesamples::Bool = false`: Write MCMC samples

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
