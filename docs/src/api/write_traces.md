# write_traces Function

Generate model-predicted intensity traces for GRSM models.

## Syntax

```julia
write_traces(; kwargs...)
```

## Arguments

### Model Parameters

- `G::Int = 2`: Number of gene states
- `R::Int = 0`: Number of pre-RNA steps
- `S::Int = 0`: Number of splice sites
- `insertstep::Int = 1`: Reporter insertion step
- `transitions::Tuple = ()`: State transitions
- `rates::Vector{Float64}`: Model rates
- `nalleles::Int = 1`: Number of alleles

### Trace Parameters

- `ntraces::Int = 100`: Number of traces to generate
- `tspan::Tuple{Float64, Float64} = (0., 100.)`: Time span for traces
- `dt::Float64 = 1.0`: Time step
- `noise::Vector{Float64} = Float64[]`: Observation noise parameters
- `probfn::Function = prob_Gaussian`: Observation probability function

### Output Parameters

- `outfolder::String = "traces"`: Output folder
- `label::String = ""`: Output file label
- `write::Bool = true`: Write traces to file
- `returntraces::Bool = false`: Return traces array

## Returns

- `traces`: Array of generated traces (if `returntraces = true`)

## Examples

```julia
# Generate traces for a simple G model
write_traces(
    G = 2,
    R = 0,
    rates = [0.1, 0.2],  # G1->G2, G2->G1
    ntraces = 100,
    tspan = (0., 100.),
    dt = 1.0
)

# Generate traces for a GR model
write_traces(
    G = 2,
    R = 3,
    S = 2,
    insertstep = 1,
    rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    ntraces = 50,
    tspan = (0., 200.),
    dt = 0.5
)
```

## Notes

1. **Rate Order**
   - G transitions
   - R transitions
   - S transitions
   - Decay
   - Noise parameters

2. **Trace Generation**
   - Traces are generated using Gillespie algorithm
   - Noise is added according to specified `probfn`
   - Output format: time, intensity, G state, R state, S state
