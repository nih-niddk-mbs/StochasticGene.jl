# Dual Reporter Analysis

This example demonstrates how to analyze data from dual reporter systems (e.g., MS2 and PP7) in live cell imaging experiments.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
rna_setup("dual_reporter_example")
cd("dual_reporter_example")
```

## Data Preparation

Place your dual reporter data in the `data/` directory. The data should be organized as follows:

```
data/
├── gene_name/
│   ├── condition1/
│   │   ├── cell1/
│   │   │   ├── ms2.csv
│   │   │   └── pp7.csv
│   │   ├── cell2/
│   │   │   ├── ms2.csv
│   │   │   └── pp7.csv
│   │   └── ...
│   └── condition2/
│       └── ...
```

Each CSV file should contain:
- Time points (in seconds)
- Fluorescence intensity values
- Optional: Background values

## Model Definition

We'll fit a model with:
- 2 gene states (G=2)
- 2 pre-RNA steps (R=2)
- MS2 reporter at first step
- PP7 reporter at second step

```julia
# Define model parameters
G = 2  # Number of gene states
R = 2  # Number of pre-RNA steps
transitions = ([1,2], [2,1])  # Gene state transitions

# Define reporter positions
reporter_positions = Dict(
    "MS2" => 1,  # First step
    "PP7" => 2   # Second step
)
```

## Fitting the Model

Now we can fit the model to our dual reporter data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    datatype = "dual_trace",
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
# Analyze MS2 data
ms2_stats = analyze_reporter(fits, "MS2")
plot_reporter(ms2_stats, "results/ms2/")

# Analyze PP7 data
pp7_stats = analyze_reporter(fits, "PP7")
plot_reporter(pp7_stats, "results/pp7/")
```

### Cross-Correlation Analysis

```julia
# Calculate cross-correlation
correlation = analyze_correlation(fits)

# Plot correlation
plot_correlation(correlation, "results/correlation/")
```

## Advanced Analysis

### Time Lag Analysis

```julia
# Calculate time lags
time_lags = analyze_time_lags(fits)

# Plot time lags
plot_time_lags(time_lags, "results/time_lags/")
```

### Burst Analysis

```julia
# Analyze bursts for each reporter
ms2_bursts = analyze_bursts(fits, "MS2")
pp7_bursts = analyze_bursts(fits, "PP7")

# Compare burst statistics
compare_bursts(ms2_bursts, pp7_bursts, "results/bursts/")
```

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
        datatype = "dual_trace",
        datapath = "data/",
        gene = "MYC",
        datacond = cond,
        reporter_positions = reporter_positions
    )
    push!(fits_list, fits)
end

# Compare conditions
compare_conditions(fits_list, conditions, "results/comparison/")
```

## Best Practices

1. **Data Quality**
   - Ensure proper channel alignment
   - Check for cross-talk between channels
   - Verify background subtraction for each channel

2. **Model Selection**
   - Consider reporter maturation times
   - Account for different reporter properties
   - Validate model assumptions

3. **Analysis Pipeline**
   - Process each reporter consistently
   - Document analysis steps
   - Save intermediate results

## Next Steps

- Try different reporter configurations
- Experiment with more complex models
- Compare results across different genes

For more advanced examples, see:
- [Time Series Analysis](@ref)
- [Coupled Models](@ref)
- [Hierarchical Models](@ref) 