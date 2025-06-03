# Installation

## Local Installation

To install StochasticGene on a local computer:

```julia
julia> ] add StochasticGene
```

After installation, you can test it with:

```julia
julia> ] test StochasticGene
```

To update to a new version:

```julia
julia> ] update StochasticGene
```

## Biowulf Installation

To install StochasticGene on NIH Biowulf:

```bash
[username@biowulf ~]$ sinteractive --constraint=x2695 --mem=64G
[username@biowulf ~]$ module load julialang
[username@biowulf ~]$ julia -t 1
```

Then install StochasticGene:

```julia
julia> ] add StochasticGene
```

If you encounter Julia crashes due to updates, try:

```bash
[username@biowulf ~]$ rm -r --force .julia
```

Then restart Julia and reinstall StochasticGene.

## Requirements

- Julia 1.9.3 or higher
- Required packages will be automatically installed during installation
