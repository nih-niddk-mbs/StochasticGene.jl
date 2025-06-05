# Getting Started

## Folder Structure

StochasticGene expects a specific directory structure for data and results:

```bash
project_root/
├── data/
│   ├── alleles/
│   └── halflives/
└── results/
```

## Setting Up

Use the built-in helper to set up a new project folder with mock data:

```julia
using StochasticGene
rna_setup("scRNA")  # Creates the structure in the "scRNA" folder
cd("scRNA")         # Move into the project directory
```

## Fitting a Model

Fit a simple two-state model to the mock data:

```julia
fits = fit(G=2, R=0, transitions=([1,2],[2,1]), datatype="rna", datapath="data/HCT116_testdata/")
```

This will use the default mock data and model settings. For more advanced usage, see the API and Examples sections.

## Model Components

### Gene States (G)
- Arbitrary number of gene states
- User-specified transitions between states
- One active state for transcription initiation

### Pre-RNA Steps (R)
- Irreversible forward transitions
- mRNA ejection from final R step
- Optional reporter insertion step

### Splicing (S)
- Up to R splice sites
- PreRNA with or without spliced intron
- Multiple configurations per R step

## Model Types

StochasticGene supports several model types:

- G models (telegraph models)
- GR models (with pre-RNA steps)
- GRS models (with splicing)
- Coupled models (multiple alleles)

## Data Types

The package can handle:

- mRNA count distributions (smFISH, scRNA)
- Image intensity traces (live cell imaging)
- Dwell time distributions
- Combined data types

## Basic Usage Examples

### RNA Count Analysis
```julia
# Fit a two-state model to RNA count data
fits = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/",
    gene = "MYC"
)
```

### Live Cell Imaging
```julia
# Analyze MS2 reporter data
fits = fit(
    G = 2,
    R = 1,
    transitions = ([1,2], [2,1]),
    datatype = "trace",
    datapath = "data/",
    gene = "MS2"
)
```

### Dwell Time Analysis
```julia
# Analyze dwell time distributions
fits = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "dwell",
    datapath = "data/",
    gene = "PP7"
)
```

## Next Steps

- Check the [Examples](@ref) section for more complex usage scenarios
- Read the [API Reference](@ref) for detailed function documentation
- Join the [GitHub discussions](https://github.com/nih-niddk-mbs/StochasticGene.jl/discussions) for community support
