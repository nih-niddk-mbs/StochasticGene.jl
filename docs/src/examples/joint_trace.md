# Joint Trace Analysis

This example demonstrates how to analyze joint trace data from multiple reporters (e.g., MS2 and PP7) in live cell imaging experiments.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("joint_trace_example")
cd("joint_trace_example")

# Generate example data using test_fit_tracejoint
fitted_rates, target_rates = test_fit_tracejoint(
    coupling=Dict("gene1" => "gene2"),  # Coupling between genes
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

Place your joint trace data in the `data/` directory. The data should be organized as follows:

```
data/
├── gene_name/
│   ├── condition1/
│   │   ├── ms2_traces.csv
│   │   ├── pp7_traces.csv
│   │   └── metadata.csv
│   └── condition2/
│       ├── ms2_traces.csv
│       ├── pp7_traces.csv
│       └── metadata.csv
```

Each trace CSV file should contain:
- Time points
- Fluorescence intensity values
- Optional background values

## Model Definition

We'll fit a joint model with:
- 2 gene states (G=2) - ON and OFF
- 2 pre-RNA steps (R=2)
- Simple transitions between states
- Reporter positions for MS2 and PP7

```julia
# Define model parameters
G = 2  # Number of gene states (ON and OFF)
R = 2  # Number of pre-RNA steps

# Define state transitions
# Format: (from_states, to_states)
transitions = (
    [1, 2],  # From states
    [2, 1]   # To states
)

# Define transcription rates for each state
transcription_rates = [0.0, 1.0]  # No transcription in OFF state

# Define reporter positions
reporter_positions = (
    ms2 = 1,  # MS2 reporter at first step
    pp7 = 2   # PP7 reporter at second step
)
```

## Fitting the Model

Now we can fit the joint model to our trace data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    transcription_rates = transcription_rates,
    datatype = "joint_trace",
    datapath = "data/",
    gene = "MYC",
    datacond = "CONTROL",
    reporter_positions = reporter_positions
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

### Reporter-Specific Analysis

```julia
# Analyze MS2 reporter
ms2_analysis = analyze_reporter(fits, "ms2")
plot_reporter(ms2_analysis, "results/ms2/")

# Analyze PP7 reporter
pp7_analysis = analyze_reporter(fits, "pp7")
plot_reporter(pp7_analysis, "results/pp7/")
```

### Cross-Correlation Analysis

```julia
# Calculate cross-correlation
cross_corr = calculate_cross_correlation(fits)
plot_cross_correlation(cross_corr, "results/cross_correlation/")

# Analyze time lag
time_lag = analyze_time_lag(fits)
plot_time_lag(time_lag, "results/time_lag/")
```

## Advanced Analysis

### Burst Analysis

```julia
# Analyze transcriptional bursts
burst_stats = analyze_bursts(fits)
println(burst_stats)

# Plot burst statistics
plot_bursts(fits, "results/bursts/")
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
        datatype = "joint_trace",
        datapath = "data/",
        gene = "MYC",
        datacond = cond,
        reporter_positions = reporter_positions
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

- Try different reporter configurations
- Experiment with preprocessing steps
- Compare results across different genes

For more advanced examples, see:
- [Hierarchical Models](@ref)
- [Coupled Models](@ref)
- [Time Series Analysis](@ref) 