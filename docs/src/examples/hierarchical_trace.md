# Hierarchical Trace Analysis

This example demonstrates how to analyze trace data using hierarchical models to account for cell-to-cell variability.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("hierarchical_example")
cd("hierarchical_example")

# Generate example data using test_fit_trace_hierarchical
fitted_rates, target_rates = test_fit_trace_hierarchical(
    traceinfo=Dict("cell" => "example", "gene" => "test"),  # Trace metadata
    G=2,                             # 2 gene states
    R=2,                             # 2 RNA states
    S=2,                             # 2 splicing states
    transitions=([1, 2], [2, 1]),    # Simple two-state model
    rtarget=[0.33, 0.19, 2.5, 1.0],  # Target rates for simulation
    totaltime=1000,                  # Total simulation time
    ntrials=10,                      # Number of simulation trials
    fittedparam=[1, 2, 3],          # Parameters to fit
    hierarchical=true,               # Use hierarchical model
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
│   │   ├── cell1.csv
│   │   ├── cell2.csv
│   │   └── metadata.csv
│   └── condition2/
│       ├── cell1.csv
│       ├── cell2.csv
│       └── metadata.csv
```

Each cell CSV file should contain:
- Time points
- Fluorescence intensity values
- Optional background values

## Model Definition

We'll fit a hierarchical model with:
- 2 gene states (G=2) - ON and OFF
- 1 pre-RNA step (R=1)
- Simple transitions between states
- Cell-specific parameters

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

# Define hierarchical structure
hierarchical = true
```

## Fitting the Model

Now we can fit the hierarchical model to our trace data:

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
    datacond = "CONTROL",
    hierarchical = true
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

### Cell-Specific Analysis

```julia
# Analyze cell-specific parameters
cell_params = analyze_cell_parameters(fits)
plot_cell_parameters(cell_params, "results/cell_parameters/")

# Plot individual cell traces
plot_cell_traces(fits, "results/cell_traces/")
```

### Population Analysis

```julia
# Analyze population statistics
pop_stats = analyze_population(fits)
plot_population_stats(pop_stats, "results/population/")

# Calculate population averages
pop_avg = calculate_population_average(fits)
plot_population_average(pop_avg, "results/population_average/")
```

## Advanced Analysis

### Variability Analysis

```julia
# Analyze cell-to-cell variability
variability = analyze_variability(fits)
plot_variability(variability, "results/variability/")

# Calculate correlation structure
correlations = calculate_correlations(fits)
plot_correlations(correlations, "results/correlations/")
```

### Model Comparison

```julia
# Compare hierarchical and non-hierarchical models
models = [
    (hierarchical=true,  name="Hierarchical"),
    (hierarchical=false, name="Non-hierarchical")
]

model_fits = []
for (hier, name) in models
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = transitions,
        transcription_rates = transcription_rates,
        datatype = "trace",
        datapath = "data/",
        gene = "MYC",
        datacond = "CONTROL",
        hierarchical = hier
    )
    push!(model_fits, (fits, stats, name))
end

# Compare models
compare_models(model_fits, "results/model_comparison/")
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

### Parameter Identifiability
```julia
# Check parameter identifiability
identifiability = check_identifiability(fits)
plot_identifiability(identifiability, "results/identifiability/")
```

### Convergence
```julia
# Check model convergence
convergence = check_convergence(fits)
plot_convergence(convergence, "results/convergence/")

# Analyze parameter correlations
correlations = analyze_parameter_correlations(fits)
plot_parameter_correlations(correlations, "results/parameter_correlations/")
```

## Next Steps

- Try different hierarchical structures
- Experiment with parameter priors
- Compare results across different genes

For more advanced examples, see:
- [Joint Analysis](@ref)
- [Coupled Models](@ref)
- [Time Series Analysis](@ref) 