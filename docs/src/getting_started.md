# Getting Started

## Setting Up Your Project

StochasticGene requires a specific directory structure:

```bash
root/
├── data/
│   ├── alleles/
│   └── halflives/
└── results/
```

To set up this structure:

```julia
using StochasticGene
rna_setup("scRNA")  # Creates the structure in the "scRNA" folder
```

The `rna_setup()` command creates:
- A `data` folder with `alleles` and `halflives` subfolders
- A `results` folder for storing analysis results
- Mock data in the `data` folder

## Basic Model Fitting

StochasticGene can fit various data types:
- mRNA count distributions (smFISH, scRNA)
- Image intensity traces (live cell imaging)
- Dwell time distributions
- Combined data types

Models can include:
- Arbitrary number of G (gene) states
- Pre-RNA steps (R)
- Splice sites (S)
- Reporter insertion steps
- Multiple alleles
- Coupled genes/alleles

To fit a basic model:

```julia
fits = fit(
    G = 2,  # Number of gene states
    R = 0,  # Number of pre-RNA steps
    transitions = ([1,2], [2,1]),  # Gene state transitions
    datatype = "rna",  # Data type
    datapath = "data/HCT116_testdata/",  # Path to data
    root = pwd()  # Root directory
)

The `fit` function returns six variables:
- Model parameters
- Posterior distributions
- Predicted data
- Goodness of fit metrics
- Model diagnostics
- Analysis results

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
