# Getting Started

## Directory Structure Setup

StochasticGene requires a specific directory structure where data are stored and results are saved. At the top is the `root` folder (e.g., "scRNA" or "RNAfits") with subfolders `data` and `results`. Inside `data` are two more folders: `alleles` and `halflives`, containing allele numbers and half lives, respectively.

To create this structure:

```julia
using StochasticGene

# Create folder structure in "scRNA" directory
rna_setup("scRNA")

# Or use current directory as root
rna_setup()

# Move into the folder
cd("scRNA")

# Verify current directory
pwd()
```

## Basic Usage

### Fitting a Simple Model

To fit a model, you need to load data and choose a model. Data currently accepted are:
- Stationary histograms (e.g., smFISH or scRNA)
- Intensity time traces (e.g., trk files)
- Dwell time distributions
- Combinations of the above

Here's a basic example fitting a two-state telegraph model:

```julia
using StochasticGene

# Fit with 4 MCMC chains
fits, stats, measures, data, model, options = fit(nchains=4)
```

The `fit` function returns six variables:
- `fits`: MCMC fit results (posterior samples, log-likelihoods, etc.)
- `stats`: Summary statistics for parameters
- `measures`: Diagnostic measures (including WAIC and its standard error)
- `data`, `model`, `options`: The data, model, and options structures used

### Model Types

Models are distinguished by:
- Number of G states
- Transition graph between G states
- Number of R steps
- Reporter insertion step
- Splicing configuration

For example, a G=3, R=2, S=2 model with reporter at step 1:

```julia
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
```

## Advanced Usage

### Parallel Processing

To use multiple processors:

```julia
# Start Julia with 4 processors
$ julia -p 4

# On Biowulf, you may need to specify threads
$ julia -p 4 -t 1

# Verify processor count
julia> nworkers()
4

# Use StochasticGene on all workers
julia> @everywhere using StochasticGene
```

### Batch Processing on Biowulf

For batch processing on Biowulf:

1. Create swarm files:
```julia
using StochasticGene

# For specific genes
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

# For all genes in data folder
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

2. Run the swarm file:
```bash
[username@biowulf ~]$ swarm -f fit_HCT116-scRNA-ss_MOCK_2.swarm --time 1:00:00 -t 8 -g 16 --merge-output --module julialang
```

### Analyzing Results

To analyze the results:

```julia
# Create data frames of results
write_dataframes_only("results/HCT_scRNAtest", "data/HCT116_testdata")

# Generate augmented summary
write_augmented("results/HCT_scRNAtest/Summary_HCT116-scRNA-ss_MOCK_2.csv", "results/HCT_scRNAtest")
```

## Units

- All rates have units of inverse minutes
- Half lives in the `halflives` file are in hours
- For stationary mRNA distributions, rate units are relative (scaling all rates by a constant won't affect results)
- For traces and dwell time histograms, units are in minutes

For more detailed examples and advanced usage, see the [Examples](@ref) section. 