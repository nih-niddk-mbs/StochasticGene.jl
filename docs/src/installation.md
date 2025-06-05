# Installation

## Requirements

- Julia version 1.9.3 or higher
- Required packages will be automatically installed

## Local Installation

1. Open a terminal and start Julia:
```bash
$ julia
```

2. In the Julia REPL, enter package mode by typing `]` and install StochasticGene:
```julia
julia> ] add StochasticGene
```
(To exit package mode, press Ctrl+C)

3. (Optional) Test the installation:
```julia
julia> ] test StochasticGene
```
Note: This may take 20+ minutes and may show some warnings or errors regarding dependencies. These are from external packages that StochasticGene references but does not rely on, so they can be safely ignored.

4. To update to a new version:
```julia
julia> ] update StochasticGene
```

## Installation on NIH Biowulf

1. Log into Biowulf and start an interactive session:
```bash
[username@biowulf ~]$ sinteractive --constraint=x2695 --mem=64G
```

2. Load Julia:
```bash
[username@biowulf ~]$ module load julialang
```

3. Start Julia (with a single thread):
```bash
[username@biowulf ~]$ julia -t 1
```

4. Install StochasticGene as described in the local installation section.

### Troubleshooting on Biowulf

If Julia crashes after updating StochasticGene, try:
1. Go to your home directory
2. Delete the Julia directory:
```bash
[username@biowulf ~]$ rm -r --force .julia
```
3. Start Julia and reinstall StochasticGene

Note: Julia is periodically updated on Biowulf, and StochasticGene may not automatically work with the latest version. You can either:
- Call previous versions when allocating CPUs
- Reinstall StochasticGene in the new version

## Verifying Installation

To verify your installation, run:
```julia
julia> using StochasticGene
```

If no errors appear, the installation was successful. 