# Coupled Model Analysis

This example demonstrates how to analyze coupled models of gene expression, where multiple genes interact with each other.

## Setup

First, let's set up our project directory and load the package:

```julia
using StochasticGene

# Create project directory
mkdir("coupled_example")
cd("coupled_example")

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

Place your data in the `data/` directory. The data should be organized as follows:

```
data/
├── gene1/
│   ├── condition1/
│   │   ├── data.csv
│   │   └── metadata.csv
│   └── condition2/
│       ├── data.csv
│       └── metadata.csv
└── gene2/
    ├── condition1/
    │   ├── data.csv
    │   └── metadata.csv
    └── condition2/
        ├── data.csv
        └── metadata.csv
```

## Model Definition

We'll fit a coupled model with:
- Two genes, each with 2 states (G=2)
- No pre-RNA steps (R=0)
- Simple transitions between states
- Coupling between genes

```julia
# Define model parameters
G = (2, 2)  # Number of gene states for each gene
R = (0, 0)  # Number of pre-RNA steps for each gene

# Define state transitions
# Format: ((from_states_gene1, to_states_gene1), (from_states_gene2, to_states_gene2))
transitions = (
    ([1, 2], [2, 1]),  # Gene 1 transitions
    ([1, 2], [2, 1])   # Gene 2 transitions
)

# Define transcription rates for each state
transcription_rates = (
    [0.0, 1.0],  # Gene 1 rates
    [0.0, 1.0]   # Gene 2 rates
)

# Define coupling structure
coupling = (
    (1, 2),           # Coupled genes
    (Int[], [1]),     # Coupling directions
    [2, 0],           # Coupling strengths
    [0, 2],           # Coupling effects
    1                 # Coupling type
)
```

## Fitting the Model

Now we can fit the coupled model to our data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
    transcription_rates = transcription_rates,
    datatype = "coupled",
    datapath = "data/",
    genes = ("MYC", "FOS"),
    datacond = "CONTROL",
    coupling = coupling
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

### Gene-Specific Analysis

```julia
# Analyze first gene
gene1_analysis = analyze_gene(fits, 1)
plot_gene(gene1_analysis, "results/gene1/")

# Analyze second gene
gene2_analysis = analyze_gene(fits, 2)
plot_gene(gene2_analysis, "results/gene2/")
```

### Coupling Analysis

```julia
# Analyze coupling strength
coupling_strength = analyze_coupling_strength(fits)
plot_coupling_strength(coupling_strength, "results/coupling/")

# Calculate coupling effects
coupling_effects = calculate_coupling_effects(fits)
plot_coupling_effects(coupling_effects, "results/coupling_effects/")
```

## Advanced Analysis

### Time Series Analysis

```julia
# Analyze time series
time_series = analyze_time_series(fits)
plot_time_series(time_series, "results/time_series/")

# Calculate cross-correlation
cross_corr = calculate_cross_correlation(fits)
plot_cross_correlation(cross_corr, "results/cross_correlation/")
```

### Model Comparison

```julia
# Compare different coupling configurations
configurations = [
    (coupling=((1, 2), (Int[], [1]), [2, 0], [0, 2], 1), name="Strong coupling"),
    (coupling=((1, 2), (Int[], [1]), [1, 0], [0, 1], 1), name="Weak coupling"),
    (coupling=((1, 2), (Int[], Int[]), [0, 0], [0, 0], 1), name="No coupling")
]

config_fits = []
for (coupling, name) in configurations
    fits, stats, measures, data, model, options = fit(
        G = G,
        R = R,
        transitions = transitions,
        transcription_rates = transcription_rates,
        datatype = "coupled",
        datapath = "data/",
        genes = ("MYC", "FOS"),
        datacond = "CONTROL",
        coupling = coupling
    )
    push!(config_fits, (fits, stats, name))
end

# Compare configurations
compare_configurations(config_fits, "results/config_comparison/")
```

## Best Practices

1. **Model Selection**
   - Start with simple coupling
   - Use model selection criteria
   - Validate coupling assumptions

2. **Parameter Estimation**
   - Check parameter identifiability
   - Verify convergence
   - Consider parameter correlations

3. **Interpretation**
   - Relate coupling to biological mechanisms
   - Consider experimental validation
   - Document assumptions

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

- Try different coupling configurations
- Experiment with parameter priors
- Compare results across different gene pairs

For more advanced examples, see:
- [Hierarchical Models](@ref)
- [Joint Analysis](@ref)
- [Time Series Analysis](@ref) 