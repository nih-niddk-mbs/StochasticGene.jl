# Dwell Time Analysis

!!! warning "Legacy conceptual example"
    This page is retained as a conceptual sketch. New dwell-time workflows
    should prefer `dwell_specs` and the current `fit` keyword interface.

This example demonstrates how to analyze dwell time distributions from live cell imaging data.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
rna_setup("dwell_example")
cd("dwell_example")
```

## Data Preparation

Place your dwell time data in the `data/` directory. The data should be in a format compatible with the package (see Data Types for details).

## Model Definition

We'll fit a model with:
- 2 gene states (G=2)
- No pre-RNA steps (R=0)
- Simple transitions between states

```julia
# Define model parameters
G = 2  # Number of gene states
R = 0  # Number of pre-RNA steps
transitions = ([1,2], [2,1])  # Gene state transitions
```

## Fitting the Model

Now we can fit the model to our dwell time data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    datatype = "dwell",
    datapath = "data/",
    gene = "PP7",
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

# Fit writes rates, measures, and parameter statistics under `resultfolder`.
```

## Model Interpretation

The fitted model provides:
- Gene state transition rates
- Dwell time distributions
- Model likelihood and fit statistics

## Advanced Analysis

### ON/OFF Histograms
```julia
# Generate ON/OFF histograms
write_ONOFFhistograms(fits, "results/")
```

### State Residency
```julia
# Calculate state residency probabilities
write_residency_G_folder("results/my-run")
```

### Burst Analysis
```julia
# For burst summaries, run fits with `burst=true` and inspect the burst output
# files written by `fit`.
```

## Next Steps

- Try fitting models with multiple states
- Experiment with different transition configurations
- Compare results across different genes

For more advanced examples, see:
- Coupled States
- Multiple Conditions
- Hierarchical Models
