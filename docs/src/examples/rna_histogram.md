# RNA Histogram Analysis

This example demonstrates how to analyze RNA count distributions (histograms) from single-molecule FISH (smFISH) or single-cell RNA sequencing (scRNA-seq) data.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("histogram_example")
cd("histogram_example")

# Generate example data using test_fit_simrna
fitted_rates, target_rates = test_fit_simrna(
    rtarget=[0.33, 0.19, 2.5, 1.0],  # Target rates for simulation
    transitions=([1, 2], [2, 1]),    # Simple two-state model
    G=2,                             # 2 gene states
    nRNA=100,                        # Number of RNA molecules
    nalleles=2,                      # Number of alleles
    fittedparam=[1, 2, 3],          # Parameters to fit
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

The `counts.csv` file should contain:
- Cell IDs
- RNA counts per cell
- Optional: Additional cell metadata

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

Now we can fit the model to our RNA count data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
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

### Histogram Analysis

```julia
# Generate and plot histograms
plot_histograms(fits, "results/histograms/")

# Calculate histogram statistics
hist_stats = analyze_histograms(fits)
println(hist_stats)
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

### Multiple Conditions

```julia
# Compare different conditions
conditions = ["CONTROL", "TREATMENT"]
fits_list = []

for cond in conditions
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = transitions,
        datatype = "rna",
        datapath = "data/",
        gene = "MYC",
        datacond = cond
    )
    push!(fits_list, fits)
end

# Compare conditions
compare_conditions(fits_list, conditions, "results/comparison/")
```

### Model Comparison

```julia
# Fit different models
models = [
    (G=2, R=0),  # Basic telegraph
    (G=3, R=0),  # Three-state
    (G=2, R=1)   # With pre-RNA
]

model_fits = []
for (G, R) in models
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = ([1,2], [2,1]),
        datatype = "rna",
        datapath = "data/",
        gene = "MYC",
        datacond = "CONTROL"
    )
    push!(model_fits, (fits, stats))
end

# Compare models
compare_models(model_fits, models, "results/model_comparison/")
```

## Best Practices

1. **Data Quality**
   - Check for outliers
   - Verify count normalization
   - Consider batch effects

2. **Model Selection**
   - Start with simple models
   - Use model selection criteria
   - Validate assumptions

3. **Analysis Pipeline**
   - Document analysis steps
   - Save intermediate results
   - Use consistent naming conventions

## Common Issues and Solutions

### Zero-Inflation
```julia
# Handle zero-inflation
fits = fit(
    G = G,
    R = R,
    transitions = transitions,
    datatype = "rna",
    datapath = "data/",
    gene = "MYC",
    datacond = "CONTROL",
    zero_inflation = true
)
```

### Batch Effects
```julia
# Account for batch effects
fits = fit(
    G = G,
    R = R,
    transitions = transitions,
    datatype = "rna",
    datapath = "data/",
    gene = "MYC",
    datacond = "CONTROL",
    batch_correction = true
)
```

## Next Steps

- Try different model configurations
- Experiment with more complex models
- Compare results across different genes

For more advanced examples, see:
- [Multi-State Models](@ref)
- [Pre-RNA Steps](@ref)
- [Coupled Models](@ref) 