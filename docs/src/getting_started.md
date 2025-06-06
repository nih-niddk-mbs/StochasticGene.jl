# Getting Started

## Folder Structure

StochasticGene expects a specific directory structure for data and results:

```bash
project_root/
├── data/
│   ├── alleles/      # Contains allele count information
│   └── halflives/    # Contains mRNA half-life information
└── results/          # Output directory for analysis results
```

## Setting Up

Use the built-in helper to set up a new project folder with the required data structure and reference files:

```julia
using StochasticGene
rna_setup()           # Sets up folders in the current directory
# or, to set up in a subfolder:
rna_setup("scRNA")    # Sets up folders in ./scRNA
cd("scRNA")           # Move into the project directory (if you want to work there)
```

This will create the necessary folders and download reference CSV files for alleles and halflives. It does not generate arbitrary mock data, but provides the standard structure and example reference files needed to get started.

## Data Preparation

### RNA Count Data
RNA count data should be provided as a text file with the following format:
```text
# gene_condition.txt
count1
count2
count3
...
```

The data will be loaded into an `RNAData` structure with the following fields:
- `histRNA`: Vector of RNA counts
- `gene`: Gene name
- `condition`: Experimental condition
- `nalleles`: Number of alleles

### Live Cell Imaging Data
For live cell imaging data, provide a text file with time series data:
```text
# gene_condition.txt
time1 intensity1
time2 intensity2
time3 intensity3
...
```

## Fitting a Model

### Basic Two-State Model
Fit a simple two-state model to RNA count data:

```julia
fits = fit(
    G = 2,                    # Two gene states
    R = 0,                    # No pre-RNA steps
    transitions = ([1,2], [2,1]),  # Transitions between states
    datatype = "rna",         # RNA count data
    datapath = "data/HCT116_testdata/",
    gene = "MYC",            # Gene name
    datacond = "MOCK"        # Experimental condition
)
```

### Model with Pre-RNA Steps
Fit a model with pre-RNA processing steps:

```julia
fits = fit(
    G = 2,                    # Two gene states
    R = 3,                    # Three pre-RNA steps
    S = 2,                    # Two splice sites
    insertstep = 1,           # Reporter insertion at step 1
    transitions = ([1,2], [2,1]),  # Gene state transitions
    datatype = "trace",       # Live cell imaging data
    datapath = "data/testtraces",
    gene = "MS2",            # Gene name
    datacond = "testtrace",   # Experimental condition
    traceinfo = (1.0, 1., -1, 1.),  # Frame interval, start, end, active fraction
    noisepriors = [40., 20., 200., 10.],  # Noise parameters
    nchains = 4              # Number of MCMC chains
)
```

## Model Components

### Gene States (G)
- Arbitrary number of gene states
- User-specified transitions between states
- One active state for transcription initiation
- Example transitions:
  - Two-state: `([1,2], [2,1])` (telegraph model)
  - Three-state: `([1,2], [2,1], [2,3], [3,1])` (cyclic model)

### Pre-RNA Steps (R)
- Irreversible forward transitions
- mRNA ejection from final R step
- Optional reporter insertion step
- Example: R=3 with reporter at step 1:
  ```julia
  R = 3
  insertstep = 1
  ```

### Splicing (S)
- Up to R splice sites
- PreRNA with or without spliced intron
- Multiple configurations per R step
- Example: S=2 with R=3:
  ```julia
  R = 3
  S = 2
  ```

## Model Types

StochasticGene supports several model types:

### G Models (Telegraph Models)
- Basic two-state model
- Multiple gene states
- Example:
  ```julia
  G = 2
  R = 0
  transitions = ([1,2], [2,1])
  ```

### GR Models (with Pre-RNA Steps)
- Gene states plus pre-RNA processing
- Example:
  ```julia
  G = 2
  R = 3
  transitions = ([1,2], [2,1])
  ```

### GRS Models (with Splicing)
- Gene states, pre-RNA steps, and splicing
- Example:
  ```julia
  G = 2
  R = 3
  S = 2
  insertstep = 1
  ```

### Coupled Models (Multiple Alleles)
- Multiple alleles with coupling
- Example:
  ```julia
  G = (2,2)  # Two alleles, each with 2 states
  R = (3,3)  # Three pre-RNA steps for each
  coupling = (
      (1,1),  # Model indices
      ([2], [1]),  # Source units
      (["G2"], ["G1"]),  # Source states
      ([2], [1])  # Target transitions
  )
  ```

## Data Types

The package can handle:

### RNA Count Data
- Single molecule FISH (smFISH)
- Single cell RNA sequencing (scRNA-seq)
- Example:
  ```julia
  datatype = "rna"
  datapath = "data/rna_counts/"
  ```

### Live Cell Imaging
- MS2 reporter data
- PP7 reporter data
- Example:
  ```julia
  datatype = "trace"
  datapath = "data/ms2_traces/"
  traceinfo = (1.0, 1., -1, 1.)  # Frame interval, start, end, active fraction
  ```

### Dwell Time Analysis
- ON/OFF state durations
- Example:
  ```julia
  datatype = "dwelltime"
  datapath = "data/dwell_times/"
  ```

### Combined Data
- RNA counts with dwell times
- Example:
  ```julia
  datatype = "rnadwelltime"
  datapath = ["data/rna_counts/", "data/dwell_times/"]
  ```

## Next Steps

- Check the [Examples](@ref) section for more complex usage scenarios
- Read the [API Reference](@ref) for detailed function documentation
- Join the [GitHub discussions](https://github.com/nih-niddk-mbs/StochasticGene.jl/discussions) for community support
