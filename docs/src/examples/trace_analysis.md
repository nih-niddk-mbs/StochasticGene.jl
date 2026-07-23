# Trace Analysis

!!! warning "Legacy conceptual example"
    This page is retained as a conceptual sketch. Some snippets use older helper
    names or simplified trace arguments and should not be copied verbatim. For
    current trace batch/key workflows, see [Cluster and batch workflows](../cluster_batch_workflows.md).

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

# Fit writes rates, measures, and parameter statistics under `resultfolder`.
```

### Trace Visualization

```julia
# For trace predictions, use `write_traces_key` or `write_traces` after fitting.
```

### Burst Analysis

```julia
# For burst summaries, run fits with `burst=true` and inspect the burst output files written by `fit`.
```

## Advanced Analysis

### Time Series Analysis

```julia
# For fitted key-based trace folders, use:
write_traces_key("results/my-trace-run")
write_correlation_functions_key("results/my-coupled-trace-run")
```

### State Transitions

```julia
# Inspect fitted transition rates in `rates_*.txt` or assembled summary CSVs.

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
- Hierarchical Models
- Joint Analysis
- Coupled Models
