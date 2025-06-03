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
