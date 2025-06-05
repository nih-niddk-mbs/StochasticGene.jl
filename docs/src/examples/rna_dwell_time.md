# RNA Dwell Time Analysis

This example demonstrates how to analyze RNA dwell time data to understand transcriptional dynamics.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("dwell_time_example")
cd("dwell_time_example")

# Generate example data using test_fit_rnadwelltime
predicted_hist, data_hist = test_fit_rnadwelltime(
    rtarget=[0.33, 0.19, 2.5, 1.0],  # Target rates for simulation
    transitions=([1, 2], [2, 1]),    # Simple two-state model
    G=2,                             # 2 gene states
    R=2,                             # 2 RNA states
    S=2,                             # 2 splicing states
    nRNA=100,                        # Number of RNA molecules
    nalleles=2,                      # Number of alleles
    onstates=[1],                    # States considered "on"
    dttype="on",                     # Analyze "on" dwell times
    fittedparam=[1, 2, 3],          # Parameters to fit
    nchains=1                        # Single MCMC chain for example
)

# Print results
println("Predicted histogram: ", predicted_hist)
println("Data histogram: ", data_hist)
```

## Data Preparation

Place your dwell time data in the `data/` directory. The data should be organized as follows:

```
data/
├── gene_name/
│   ├── condition1/
│   │   ├── dwell_times.csv
│   │   └── metadata.csv
│   └── condition2/
│       ├── dwell_times.csv
│       └── metadata.csv
```

Each dwell time CSV file should contain:
- Time intervals
- State labels (ON/OFF)
- Optional metadata

## Model Definition

We'll fit a model with:
- 2 gene states (G=2) - ON and OFF
- No pre-RNA steps (R=0)
- Simple transitions between states

```julia
# Define model parameters
G = 2  # Number of gene states (ON and OFF)
R = 0  # Number of pre-RNA steps

# Define state transitions
# Format: (from_states, to_states)
transitions = (
    [1, 2],  # From states
    [2, 1]   # To states
)

# Define transcription rates for each state
transcription_rates = [0.0, 1.0]  # No transcription in OFF state
```

## Fitting the Model

Now we can fit the model to our dwell time data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    transcription_rates = transcription_rates,
    datatype = "dwell_time",
    datapath = "data/",
    gene = "MYC",
    datacond = "CONTROL"
)
```

## Analyzing Results

### Basic Analysis

```julia
# Print basic statistics
println(stats)

# Plot the results
using Plots
plot(fits)

# Save results
save_results(fits, "results/")
```

### Dwell Time Analysis

```julia
# Analyze dwell time distributions
dwell_times = analyze_dwell_times(fits)
plot_dwell_times(dwell_times, "results/dwell_times/")

# Calculate state residence times
residence_times = analyze_residence_times(fits)
plot_residence_times(residence_times, "results/residence_times/")
```

### Burst Analysis

```julia
# Analyze transcriptional bursts
burst_stats = analyze_bursts(fits)
println(burst_stats)

# Plot burst statistics
plot_bursts(fits, "results/bursts/")
```

## Advanced Analysis

### State Transition Analysis

```julia
# Analyze state transitions
transitions = analyze_transitions(fits)

# Plot transition probabilities
plot_transitions(transitions, "results/transitions/")

# Calculate transition rates
rates = calculate_transition_rates(fits)
println(rates)
```

### Model Comparison

```julia
# Compare different model configurations
configurations = [
    (G=2, rates=[0.0, 1.0]),           # Basic ON-OFF
    (G=3, rates=[0.0, 0.5, 1.0]),      # Three-state
    (G=4, rates=[0.0, 0.3, 0.7, 1.0])  # Four-state
]

config_fits = []
for (G, rates) in configurations
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = ([1:G...], [2:G..., 1]),
        transcription_rates = rates,
        datatype = "dwell_time",
        datapath = "data/",
        gene = "MYC",
        datacond = "CONTROL"
    )
    push!(config_fits, (fits, stats))
end

# Compare configurations
compare_configurations(config_fits, configurations, "results/config_comparison/")
```

## Best Practices

1. **Model Selection**
   - Start with basic ON-OFF model
   - Use model selection criteria
   - Validate state assumptions

2. **Parameter Estimation**
   - Check parameter identifiability
   - Verify convergence
   - Consider parameter correlations

3. **Interpretation**
   - Relate states to biological mechanisms
   - Consider experimental validation
   - Document assumptions

## Common Issues and Solutions

### Parameter Identifiability
```julia
# Check parameter identifiability
identifiability = check_identifiability(fits)
plot_identifiability(identifiability, "results/identifiability/")
```

### State Validation
```julia
# Validate state assignments
validation = validate_states(fits)
plot_validation(validation, "results/validation/")
```

## Next Steps

- Try different state configurations
- Experiment with transition patterns
- Compare results across different genes

For more advanced examples, see:
- [Pre-RNA Steps](@ref)
- [Coupled Models](@ref)
- [Hierarchical Models](@ref) 