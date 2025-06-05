# Examples

## Basic Model Definition

```julia
using StochasticGene

# Define a simple model
@define_model SimpleModel num_states::Int, rate::Float64

# Register the model
@register_model SimpleModel

# Create an instance
model = create_model(:SimpleModel, 3, 0.1)
```

## Hidden Markov Model

```julia
using StochasticGene

# Create an HMM model
model = create_model(:HMMModel, 5)

# Create some data
data = TimeSeriesData()

# Run the model
run_model(model, data)
```

## Kalman Filter Model

```julia
using StochasticGene

# Create a Kalman filter model
model = create_model(:KalmanModel, 4, 2)  # 4 states, 2 observations

# Create some data
data = TimeSeriesData()

# Run the model
run_model(model, data)
```

## Coupled System Model

```julia
using StochasticGene

# Create a coupled system model
model = create_model(:CoupledSystemModel, 0.5)  # coupling strength = 0.5

# Create some data
data = TimeSeriesData()

# Run the model
run_model(model, data)
```

## Custom Model with Type Dispatch

```julia
using StochasticGene

# Define a custom model
@define_model CustomModel param1::Int, param2::Float64

# Register the model
@register_model CustomModel

# Define dispatch rules
@dispatch_by_type CustomModel TimeSeriesData

# Create model and data
model = create_model(:CustomModel, 10, 0.5)
data = TimeSeriesData()

# Run the model
run_model(model, data)
```

## Model Registry Usage

```julia
using StochasticGene

# Check available models
println("Available models: ", keys(MODEL_REGISTRY))

# Create a model dynamically
model_type = :HMMModel
if haskey(MODEL_REGISTRY, model_type)
    model = create_model(model_type, 5)
    println("Created model: ", model)
end
```

## Basic Two-State Model

Here's a simple example of fitting a two-state telegraph model to RNA histogram data:

```julia
using StochasticGene

# Set up directory structure
rna_setup("scRNA")
cd("scRNA")

# Fit the model with 4 MCMC chains
fits, stats, measures, data, model, options = fit(nchains=4)
```

## Fitting to RNA-seq Data

This example demonstrates fitting to RNA-seq data with a more complex model:

```julia
using StochasticGene

# Set up directory structure
rna_setup("scRNA")
cd("scRNA")

# Fit a G=3, R=2, S=2 model
fits, stats, measures, data, model, options = fit(
    nchains=4,
    datatype="rna",
    datapath="data/HCT116_testdata/",
    cell="HCT116",
    gene="MYC",
    datacond="MOCK",
    resultfolder="HCT_scRNA",
    G=3,
    R=2,
    S=2,
    insertstep=1,
    transitions=([1,2], [2,1], [2,3], [3,1])
)
```

## Live Cell Imaging Analysis

This example shows how to analyze live cell imaging data:

```julia
using StochasticGene

# Create simulated trace data
simulate_trace_data("data/testtraces/")

# Fit a model to the traces
fits, stats, measures, data, model, options = fit(
    nchains=4,
    datatype="trace",
    datapath="data/testtraces",
    cell="TEST",
    gene="test",
    datacond="testtrace",
    resultfolder="trace-test",
    infolder="trace-test",
    traceinfo=(1.0, 1, -1, 1.0),
    transitions=([1, 2], [2, 1], [2, 3], [3, 1]),
    G=3,
    R=2,
    S=2,
    insertstep=1,
    decayrate=1.0,
    noisepriors=[40.0, 20.0, 200.0, 10.0]
)

# Generate predicted traces
write_traces("results/trace-test/", "data/testtraces", "", 1.0)

# Generate ON/OFF histograms
write_ONOFFhistograms("results/trace-test/")

# Generate G state residence time probabilities
write_residency_G_folder("results/trace-test/")
```

## Parallel Processing

This example demonstrates parallel processing with multiple chains:

```julia
# Start Julia with 4 processors
$ julia -p 4

# On Biowulf, you may need to specify threads
$ julia -p 4 -t 1

# Use StochasticGene on all workers
julia> @everywhere using StochasticGene

# Fit with multiple chains
fits, stats, measures, data, model, options = fit(nchains=4)
```

## Simulation Examples

### Basic Simulation

Here's an example of simulating a model:

```julia
# Define rates and transitions
r = [0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231]
transitions = ([1, 2], [2, 1], [2, 3], [3, 1])

# Simulate the model
h = simulator(
    r,
    transitions,
    3,  # G states
    2,  # R steps
    2,  # S splice sites
    1,  # insertstep
    nhist=150,
    bins=[collect(5/3:5/3:200), collect(.1:.1:20)],
    onstates=[Int[], [2,3]],
    nalleles=2
)
```

### Simulating Trace Data

To generate simulated trace data:

```julia
simulate_trace_data(
    "data/testtraces/",
    ntrials=10,
    r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231, 30, 20, 200, 100, 0.8],
    transitions=([1, 2], [2, 1], [2, 3], [3, 1]),
    G=3,
    R=2,
    S=2,
    insertstep=1,
    onstates=Int[],
    interval=1.0,
    totaltime=1000.0
)
```

## Batch Processing on Biowulf

This example shows how to set up batch processing on Biowulf:

```julia
using StochasticGene

# Create swarm files for specific genes
makeswarm(
    ["CENPL", "MYC"],
    cell="HCT116",
    maxtime=600.0,
    nchains=8,
    datatype="rna",
    G=2,
    transitions=([1,2], [2,1]),
    datacond="MOCK",
    resultfolder="HCT_scRNA",
    datapath="HCT116_testdata/",
    root="."
)

# Or for all genes in the data folder
makeswarm_genes(
    cell="HCT116",
    maxtime=600.0,
    nchains=8,
    datatype="rna",
    G=2,
    transitions=([1,2], [2,1]),
    datacond="MOCK",
    resultfolder="HCT_scRNA",
    datapath="HCT116_testdata/",
    nsets=1,
    root="."
)
```

Then run the swarm file:
```bash
[username@biowulf ~]$ swarm -f fit_HCT116-scRNA-ss_MOCK_2.swarm --time 1:00:00 -t 8 -g 16 --merge-output --module julialang
```

## Analyzing Results

After running fits, you can analyze the results:

```julia
# Create data frames of results
write_dataframes_only("results/HCT_scRNAtest", "data/HCT116_testdata")

# Generate augmented summary
write_augmented("results/HCT_scRNAtest/Summary_HCT116-scRNA-ss_MOCK_2.csv", "results/HCT_scRNAtest")
```

For more detailed information about the functions and their parameters, see the [API Reference](@ref). 