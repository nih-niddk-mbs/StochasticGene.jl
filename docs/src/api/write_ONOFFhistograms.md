# write_ONOFFhistograms Function

Generate ON/OFF dwell time histograms for GRSM models.

## Syntax

```julia
write_ONOFFhistograms(; kwargs...)
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

### Histogram Parameters

- `ntraces::Int = 1000`: Number of traces for histogram
- `tspan::Tuple{Float64, Float64} = (0., 1000.)`: Time span for traces
- `dt::Float64 = 1.0`: Time step
- `bins::Int = 100`: Number of histogram bins
- `maxtime::Float64 = 100.0`: Maximum dwell time

### Output Parameters

- `outfolder::String = "histograms"`: Output folder
- `label::String = ""`: Output file label
- `write::Bool = true`: Write histograms to file
- `returnhist::Bool = false`: Return histogram data

## Returns

- `hist`: Tuple containing:
  - ON dwell time histogram
  - OFF dwell time histogram
  - Bin edges

## Examples

```julia
# Generate histograms for a simple G model
write_ONOFFhistograms(
    G = 2,
    R = 0,
    rates = [0.1, 0.2],  # G1->G2, G2->G1
    ntraces = 1000,
    tspan = (0., 1000.),
    bins = 100,
    maxtime = 100.0
)

# Generate histograms for a GR model
write_ONOFFhistograms(
    G = 2,
    R = 3,
    S = 2,
    insertstep = 1,
    rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
    ntraces = 5000,
    tspan = (0., 2000.),
    bins = 150,
    maxtime = 200.0
)
```

## Notes

1. **Histogram Generation**
   - Traces are generated using Gillespie algorithm
   - ON/OFF states are determined by gene state
   - Histograms are normalized to total number of transitions

2. **Rate Order**
   - G transitions
   - R transitions
   - S transitions
   - Decay
