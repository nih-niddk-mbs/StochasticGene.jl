# simulator Function

Simulate stochastic gene expression models using the Gillespie algorithm.

## Syntax

```julia
simulator(rin, transitions, G, R, S, insertstep; kwargs...)
```

## Arguments

### Required Arguments

- `rin::Vector{Float64}`: Initial transition rates
- `transitions::Tuple`: Tuple of vectors specifying state transitions
- `G::Int`: Number of gene states
- `R::Int`: Number of pre-RNA steps
- `S::Int`: Number of splice sites (must be ≤ R - insertstep + 1)
- `insertstep::Int`: Step where reporter becomes visible

### Optional Keyword Arguments

#### Simulation Parameters
- `warmupsteps::Int = 0`: Number of warmup steps before recording
- `nalleles::Int = 1`: Number of alleles
- `nhist::Int = 20`: Number of histogram bins
- `bins::Vector{Float64} = Float64[]`: Custom histogram bins
- `traceinterval::Float64 = 0.0`: Time interval for trace recording (0 = no traces)

#### Model Configuration
- `coupling::Tuple = tuple()`: Coupling parameters for multi-unit models
- `onstates::Vector{Int} = Int[]`: States where transcription is active
- `splicetype::String = ""`: Splicing configuration ("", "offeject", "offdecay")

#### Observation Model
- `probfn::Function = prob_Gaussian`: Probability function for observations
- `noise::Vector{Float64} = Float64[]`: Noise parameters
- `reportersteps::Vector{Int} = Int[]`: Steps where reporter is visible

#### Output Control
- `tspan::Tuple{Float64, Float64} = (0., 1000.)`: Time span for simulation
- `ntrials::Int = 1`: Number of simulation trials
- `resultfolder::String = ""`: Output folder for results

## Returns

- `histogram::Vector{Float64}`: Steady-state RNA count distribution
- `traces::Vector{Vector{Float64}}`: Intensity traces (if `traceinterval > 0`)
- `dtimes::Vector{Float64}`: Dwell times (if applicable)

## Examples

### Basic Two-State Model

```julia
using StochasticGene

# Simple two-state telegraph model
rates = [0.1, 0.2]  # G1->G2, G2->G1
transitions = ([1,2], [2,1])
G, R, S = 2, 0, 0
insertstep = 1

# Simulate RNA histogram
histogram = simulator(
    rates, transitions, G, R, S, insertstep,
    nhist = 50,
    ntrials = 1000
)
```

### GRS Model with Traces

```julia
# Gene-Reporter-Splice model
rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
transitions = ([1,2], [2,1])
G, R, S = 2, 3, 2
insertstep = 1

# Simulate with intensity traces
histogram, traces = simulator(
    rates, transitions, G, R, S, insertstep,
    traceinterval = 1.0,    # 1 minute intervals
    tspan = (0., 500.),     # 500 minute simulation
    noise = [10.0, 5.0],    # Gaussian noise parameters
    probfn = prob_Gaussian
)
```

### Coupled Model

```julia
# Two coupled transcriptional units
rates1 = [0.1, 0.2, 0.3]
rates2 = [0.15, 0.25, 0.35]
transitions1 = ([1,2], [2,1])
transitions2 = ([1,2], [2,1])

# Coupling: units share some transition rates
coupling = (2, 2, [1,2], [3,4], 2)  # Coupling specification

histogram = simulator(
    [rates1; rates2], 
    (transitions1, transitions2),
    (2, 2),        # G states for each unit
    (1, 1),        # R steps for each unit
    (0, 0),        # S sites for each unit
    (1, 1),        # Insert steps
    coupling = coupling,
    nalleles = 2
)
```

### Hierarchical Model

```julia
# Simulate data for hierarchical fitting
rates = [0.1, 0.2, 0.3]
transitions = ([1,2], [2,1])
G, R, S = 2, 1, 0
insertstep = 1

# Generate multiple traces for hierarchical analysis
ntraces = 10
all_traces = Vector{Vector{Float64}}()

for i in 1:ntraces
    # Add noise to rates for each trace
    noisy_rates = rates .* (1 .+ 0.1 * randn(length(rates)))
    
    histogram, traces = simulator(
        noisy_rates, transitions, G, R, S, insertstep,
        traceinterval = 0.5,
        tspan = (0., 200.),
        noise = [20.0, 10.0]
    )
    
    append!(all_traces, traces)
end
```

### Custom Observation Model

```julia
# Define custom observation function
function prob_Poisson(y, μ, σ)
    return pdf(Poisson(μ), round(Int, y))
end

# Simulate with Poisson observation noise
histogram = simulator(
    rates, transitions, G, R, S, insertstep,
    probfn = prob_Poisson,
    noise = [15.0],  # Poisson rate parameter
    nhist = 30
)
```

## Rate Ordering

The rate vector `rin` must follow this specific ordering:

1. **G transitions**: Rates between gene states
2. **R transitions**: Rates between pre-RNA steps
3. **S transitions**: Splicing rates
4. **Decay rates**: mRNA decay rates
5. **Noise parameters**: Observation noise parameters

### Example Rate Ordering

For a model with G=2, R=3, S=2:
```julia
rates = [
    # G transitions
    0.1,    # G1 -> G2
    0.2,    # G2 -> G1
    
    # R transitions
    0.3,    # R1 -> R2
    0.4,    # R2 -> R3
    0.5,    # R3 -> eject
    
    # S transitions
    0.6,    # Splice site 1
    0.7,    # Splice site 2
    
    # Decay
    0.05    # mRNA decay
]
```

## Performance Notes

1. **Memory Usage**: Large `nhist` values require more memory
2. **Simulation Time**: Longer `tspan` and smaller rates increase runtime
3. **Trace Recording**: `traceinterval > 0` significantly increases memory usage
4. **Parallel Processing**: Use multiple calls for embarrassingly parallel simulations

## Common Use Cases

### Parameter Estimation Validation
```julia
# Generate synthetic data for testing parameter estimation
true_rates = [0.1, 0.2, 0.3]
synthetic_data = simulator(
    true_rates, transitions, G, R, S, insertstep,
    nhist = 100,
    ntrials = 5000
)

# Use synthetic_data to validate fitting algorithms
```

### Model Comparison
```julia
# Compare different model structures
models = [
    (2, 0, 0),  # Simple telegraph
    (2, 1, 0),  # With pre-RNA
    (3, 1, 0),  # Three gene states
]

for (G, R, S) in models
    histogram = simulator(rates, transitions, G, R, S, insertstep)
    # Analyze and compare histograms
end
```

### Burst Analysis
```julia
# Simulate for burst size analysis
histogram, traces = simulator(
    rates, transitions, G, R, S, insertstep,
    traceinterval = 0.1,    # High temporal resolution
    tspan = (0., 1000.),
    onstates = [2]          # State 2 is transcriptionally active
)

# Analyze burst properties from traces
```

## Error Handling

The function includes validation for:
- Rate vector length consistency
- Valid state transitions
- Proper model dimensions (G, R, S relationships)
- Positive rate values
- Valid time spans

## See Also

- [`simulate_trace`](@ref): Generate traces only
- [`simulate_trials`](@ref): Multiple simulation runs
- [`fit`](@ref): Fit models to data
- [`prob_Gaussian`](@ref): Gaussian observation model