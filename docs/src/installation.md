# Installation

## Prerequisites

- Julia version 1.9.3 or higher
- Git (for development installation)

## Basic Installation

To install StochasticGene.jl, open the Julia REPL and run:

```julia
using Pkg
Pkg.add("StochasticGene")
```

## Development Installation

To install the development version:

```julia
using Pkg
Pkg.develop(url="https://github.com/nih-niddk-mbs/StochasticGene.jl")
```

## Required Packages

The following packages will be automatically installed as dependencies:

- Distributions
- StatsBase
- DataFrames
- CSV
- Plots
- ProgressMeter
- MCMCChains
- Turing
- DifferentialEquations
- Optim
- ForwardDiff
- SpecialFunctions
- Random
- Statistics
- LinearAlgebra

## Verification

To verify the installation, run:

```julia
using StochasticGene
```

If no errors occur, the package is installed correctly.

## Troubleshooting

If you encounter any issues during installation:

1. Ensure you have the correct Julia version (1.9.3 or higher)
2. Try updating your package registry:
   ```julia
   using Pkg
   Pkg.update()
   ```
3. If problems persist, check the [GitHub Issues page](https://github.com/nih-niddk-mbs/StochasticGene.jl/issues) for known issues
4. For additional help, open a new issue on GitHub

## Local Installation

To install StochasticGene on your computer:

```julia
julia> ] add StochasticGene
```

To test the installation:

```julia
julia> ] test StochasticGene
```

To update to the latest version:

```julia
julia> ] update StochasticGene
```

## Biowulf Installation (NIH HPC)

1. Start an interactive session:

```bash
[username@biowulf ~]$ sinteractive --constraint=x2695 --mem=64G
```

2. Load Julia:

```bash
[username@biowulf ~]$ module load julialang
[username@biowulf ~]$ julia -t 1
```

3. Install StochasticGene:

```julia
julia> ] add StochasticGene
```

If you encounter Julia crashes after an update, remove your Julia depot and reinstall:

```bash
[username@biowulf ~]$ rm -r --force .julia
```

Then start Julia and re-add StochasticGene as above.


## Requirements

- Julia 1.9.3 or higher
- Required packages will be automatically installed during installation
