# Multi-State RNA Analysis

This example demonstrates how to analyze RNA count data using models with multiple gene states.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("multi_state_example")
cd("multi_state_example")

# Generate example data using test_fit_simrna with a multi-state model
fitted_rates, target_rates = test_fit_simrna(
    rtarget=[0.33, 0.19, 0.15, 2.5, 1.0],  # Target rates for simulation
    transitions=([1, 2], [2, 1], [2, 3], [3, 2]),  # Three-state model
    G=3,                             # 3 gene states
    nRNA=100,                        # Number of RNA molecules
    nalleles=2,                      # Number of alleles
    fittedparam=[1, 2, 3, 4],       # Parameters to fit
    nchains=1                        # Single MCMC chain for example
)

# Print results
println("Fitted rates: ", fitted_rates)
println("Target rates: ", target_rates)
```

## Data Preparation

Place your RNA count data in the `data/` directory. The data should be organized as follows:

```
data/
├── gene_name/
│   ├── condition1/
│   │   ├── counts.csv
│   │   └── metadata.csv
│   └── condition2/
│       ├── counts.csv
│       └── metadata.csv
```

## Model Definition

We'll fit a model with:
- 3 gene states (G=3)
- No pre-RNA steps (R=0)
- Complex transitions between states

```julia
# Define model parameters
G = 3  # Number of gene states
R = 0  # Number of pre-RNA steps

# Define state transitions
# Format: (from_states, to_states)
transitions = (
    [1, 2, 3],  # From states
    [2, 3, 1]   # To states
)

# Define transcription rates for each state
transcription_rates = [0.0, 1.0, 2.0]  # No transcription in state 1
```

## Fitting the Model

Now we can fit the model to our RNA count data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    transcription_rates = transcription_rates,
    datatype = "rna",
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

### State Analysis

```julia
# Analyze state probabilities
state_probs = analyze_states(fits)
plot_states(state_probs, "results/states/")

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
# Compare different state configurations
configurations = [
    (G=2, rates=[0.0, 1.0]),           # Two-state
    (G=3, rates=[0.0, 1.0, 2.0]),      # Three-state
    (G=4, rates=[0.0, 1.0, 2.0, 3.0])  # Four-state
]

config_fits = []
for (G, rates) in configurations
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = ([1:G...], [2:G..., 1]),
        transcription_rates = rates,
        datatype = "rna",
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
   - Start with simpler models
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