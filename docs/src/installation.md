# Installation

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
