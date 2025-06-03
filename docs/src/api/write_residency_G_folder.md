# write_residency_G_folder Function

Generate G state residency probabilities for GRSM models.

## Syntax

```julia
write_residency_G_folder(; kwargs...)
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

### Simulation Parameters

- `ntraces::Int = 1000`: Number of traces
- `tspan::Tuple{Float64, Float64} = (0., 1000.)`: Time span for traces
- `dt::Float64 = 1.0`: Time step
- `burnin::Float64 = 0.0`: Burn-in time
- `sampletime::Float64 = 1.0`: Sampling interval

### Output Parameters

- `outfolder::String = "residency"`: Output folder
- `label::String = ""`: Output file label
- `write::Bool = true`: Write residency data to file
- `returnresidency::Bool = false`: Return residency data

## Returns

- `residency`: Array of residency probabilities for each G state

## Examples

```julia
# Generate residency for a simple G model
write_residency_G_folder(
    G = 2,
    R = 0,
    rates = [0.1, 0.2],  # G1->G2, G2->G1
    ntraces = 1000,
    tspan = (0., 1000.),
    dt = 1.0
)

# Generate residency for a GR model
write_residency_G_folder(
    G = 2,
    R = 3,
    S = 2,
    insertstep = 1,
    rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    ntraces = 5000,
    tspan = (0., 2000.),
    dt = 0.5,
    burnin = 100.0,
    sampletime = 5.0
)
```

## Notes

1. **Residency Calculation**
   - Traces are generated using Gillespie algorithm
   - Residency is calculated as fraction of time in each G state
   - Results are averaged over all traces

2. **Rate Order**
   - G transitions
   - R transitions
   - S transitions
   - Decay
