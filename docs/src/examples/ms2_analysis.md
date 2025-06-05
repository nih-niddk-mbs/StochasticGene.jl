# MS2 Reporter Analysis

This example demonstrates how to analyze live cell imaging data using MS2 reporters.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
rna_setup("ms2_example")
cd("ms2_example")
```

## Data Preparation

Place your MS2 reporter image data in the `data/` directory. The data should be in a format compatible with the package (see [Data Types](@ref) for details).

## Model Definition

We'll fit a model with:
- 2 gene states (G=2)
- 1 pre-RNA step (R=1)
- MS2 reporter insertion at the first step

```julia
# Define model parameters
G = 2  # Number of gene states
R = 1  # Number of pre-RNA steps
transitions = ([1,2], [2,1])  # Gene state transitions
```

## Fitting the Model

Now we can fit the model to our MS2 data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    datatype = "trace",
    datapath = "data/",
    gene = "MS2",
    datacond = "CONTROL"
)
```

## Analyzing Results

Let's examine the fitting results:

```julia
# Print basic statistics
println(stats)

# Plot the results
using Plots
plot(fits)

# Save results
save_results(fits, "results/")
```

## Model Interpretation

The fitted model provides:
- Gene state transition rates
- Transcription rates
- Reporter maturation rates
- Model likelihood and fit statistics

## Advanced Analysis

### Time Series Analysis
```julia
# Analyze time series data
time_series = analyze_time_series(fits)
plot(time_series)
```

### Burst Analysis
```julia
# Analyze transcriptional bursts
burst_stats = analyze_bursts(fits)
println(burst_stats)
```

## Next Steps

- Try fitting models with multiple reporters
- Experiment with different reporter configurations
- Compare results across different genes

For more advanced examples, see:
- [Dual Reporter System](@ref)
- [Time Series Analysis](@ref)
- [Coupled Models](@ref) 