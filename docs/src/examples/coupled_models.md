# Coupled Model Analysis

!!! warning "Legacy conceptual example"
    This page is retained as a conceptual sketch. For current coupled CSV,
    combined-rate, and key-based coupled workflows, see
    [Cluster and batch workflows](../cluster_batch_workflows.md).

This example demonstrates how to analyze coupled models of gene expression, where multiple genes interact with each other.

!!! note "Batch jobs and combined rate files"
    For **generating scheduler commands**, **staging run specs**, and **stacking single-unit `rates_*.txt` files** into combined starts for coupled fits, see the dedicated guide [Cluster and batch workflows](../cluster_batch_workflows.md) (`makeswarm`, `stage_write_run_specs`, `create_combined_file`, etc.).

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

# Unit 1 in state 2 affects transition 1 of unit 2.
coupling = ((1, 2), [(1, 2, 2, 1)])

# Optional sign metadata constrains gamma during fitting.
inhibitory_coupling = ((1, 2), [(1, 2, 2, 1)], :inhibit)
```

Connections use `(source_unit, source_state, target_unit, target_transition)`.
Use `make_coupling` for one-way coupling, `make_coupling_reciprocal` for two-way
coupling, or `make_coupling_hidden_latent` for the supported hidden three-unit
layout.

Reporter coupling has two forms. `Rany` contributes once when any reporter
position is occupied. `Rsum` contributes once for every occupied position;
`make_coupling("R5", G, R)` is the historical `Rsum` shorthand. See
[Units and models](../concepts/units_and_models.md) for the full representation,
rate ordering, and simulator conventions.

## Fitting the Model

Now we can fit the coupled model to our data:

```julia
# Fit the model
fits, stats, measures, data, model, options = fit(
    G = G,
    R = R,
    transitions = transitions,
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

# Fit writes rates, measures, and parameter statistics under `resultfolder`.
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
# For fitted key-based coupled trace folders, use:
write_traces_key("results/my-coupled-trace-run")
write_correlation_functions_key("results/my-coupled-trace-run")
```

To compare theory with empirical curves centered independently within each
finite trace, use the actual trace lengths in frames:

```julia
write_correlation_functions_key(
    "results/my-coupled-trace-run";
    lags=collect(0:5/3:120),
    trace_center=true,
    window_lengths=216,
    window_interval=5/3,
)
```

The same operation is available without rebuilding the HMM through
`write_correlation_functions_centered(path; ...)`, as
long as the saved theoretical lag grid covers the complete trace window.
To process every theoretical correlation file in a result folder concurrently,
launch Julia with multiple threads (for example, `julia -t 10`) and call
`write_correlation_functions_centered_folder("results/my-coupled-trace-run";
window_lengths=216, window_interval=5/3)`.

To validate a coupled model's theoretical ON-ON cross-correlation against
simulation, use `simulate_trials(...; nexperiments=N)` as described in
[Validate theory against simulation](../api/analysis.md#Validate-theory-against-simulation).

### Model Comparison

```julia
# Compare different coupling configurations
configurations = [
    (coupling=((1, 2), [(1, 2, 2, 2)]), name="Strong coupling"),
    (coupling=((1, 2), [(1, 1, 2, 1)]), name="Weak coupling"),
    (coupling=((1, 2), []), name="No coupling")
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
- Hierarchical Models
- Joint Analysis
- Time Series Analysis
