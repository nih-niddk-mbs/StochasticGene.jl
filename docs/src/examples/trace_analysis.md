# Trace Analysis

This example demonstrates how to analyze trace data from live cell imaging experiments.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("trace_example")
cd("trace_example")

# Generate example data using test_fit_trace
fitted_rates, target_rates = test_fit_trace(
    traceinfo=Dict("cell" => "example", "gene" => "test"),  # Trace metadata
    G=2,                             # 2 gene states
    R=2,                             # 2 RNA states
    S=2,                             # 2 splicing states
    transitions=([1, 2], [2, 1]),    # Simple two-state model
    rtarget=[0.33, 0.19, 2.5, 1.0],  # Target rates for simulation
    totaltime=1000,                  # Total simulation time
    ntrials=10,                      # Number of simulation trials
    fittedparam=[1, 2, 3],          # Parameters to fit
    nchains=1                        # Single MCMC chain for example
)

# Print results
println("Fitted rates: ", fitted_rates)
println("Target rates: ", target_rates)
```

## Data Preparation

Place your trace data in the `data/` directory. The data should be organized as follows:

```
data/
├── gene_name/
│   ├── condition1/
│   │   ├── traces.csv
│   │   └── metadata.csv
│   └── condition2/
│       ├── traces.csv
│       └── metadata.csv
```

Each trace CSV file should contain:
- Time points
- Fluorescence intensity values
- Optional background values

## Model Definition

We'll fit a model with:
- 2 gene states (G=2) - ON and OFF
- 1 pre-RNA step (R=1)
- Simple transitions between states

```julia
# Define model parameters
G = 2  # Number of gene states (ON and OFF)
R = 1  # Number of pre-RNA steps

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

Now we can fit the model to our trace data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    transcription_rates = transcription_rates,
    datatype = "trace",
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

### Trace Visualization

```julia
# Plot individual traces
plot_traces(fits, "results/traces/")

# Plot average trace
plot_average_trace(fits, "results/average_trace/")
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

### Time Series Analysis

```julia
# Analyze time series
time_series = analyze_time_series(fits)

# Plot time series
plot_time_series(time_series, "results/time_series/")

# Calculate autocorrelation
autocorr = calculate_autocorrelation(fits)
plot_autocorrelation(autocorr, "results/autocorrelation/")
```

### State Transitions

```julia
# Analyze state transitions
transitions = analyze_transitions(fits)

# Plot transition probabilities
plot_transitions(transitions, "results/transitions/")

# Calculate transition rates
rates = calculate_transition_rates(fits)
println(rates)
```

### Multiple Conditions

```julia
# Compare different conditions
conditions = ["CONTROL", "TREATMENT"]
condition_fits = []

for cond in conditions
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = transitions,
        transcription_rates = transcription_rates,
        datatype = "trace",
        datapath = "data/",
        gene = "MYC",
        datacond = cond
    )
    push!(condition_fits, (fits, stats))
end

# Compare conditions
compare_conditions(condition_fits, conditions, "results/condition_comparison/")
```

## Best Practices

1. **Data Quality**
   - Check for photobleaching
   - Verify signal-to-noise ratio
   - Consider background correction

2. **Model Selection**
   - Start with simple models
   - Use model selection criteria
   - Validate assumptions

3. **Analysis Pipeline**
   - Document preprocessing steps
   - Save intermediate results
   - Version control your analysis

## Common Issues and Solutions

### Photobleaching
```julia
# Check for photobleaching
bleaching = check_photobleaching(fits)
plot_photobleaching(bleaching, "results/photobleaching/")

# Correct for photobleaching
corrected_fits = correct_photobleaching(fits)
```

### Background Noise
```julia
# Analyze background noise
noise = analyze_background_noise(fits)
plot_background_noise(noise, "results/background_noise/")

# Apply background correction
corrected_fits = correct_background(fits)
```

## Next Steps

- Try different model configurations
- Experiment with preprocessing steps
- Compare results across different genes

For more advanced examples, see:
- [Hierarchical Models](@ref)
- [Joint Analysis](@ref)
- [Coupled Models](@ref) 