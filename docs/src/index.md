# StochasticGene.jl

A Julia package for simulating and fitting stochastic models of gene transcription to experimental data.

## Overview

StochasticGene.jl is designed to analyze various types of experimental data including:
- Distributions of mRNA counts per cell (e.g., single molecule FISH (smFISH) or single cell RNA sequencing (scRNA-seq) data)
- Image intensity traces from live cell imaging
- Dwell time distributions of reporters (e.g., MS2 or PP7) imaged in live cells
- Combinations of the above data types

The package implements stochastic (continuous Markov) models with:
- Arbitrary number of G (gene) states
- R (pre-RNA) steps
- S splice sites (up to R)
- Reporter insertion step

Key features:
- Support for multiple alleles of the gene in the same cell
- Configurable reporter visibility at specific R steps
- Coupling between genes/alleles
- Support for both exon and intron reporter constructs
- Classic telegraph models with arbitrary numbers of G states
- Multiple simultaneous reporters (e.g., for introns and transcription factors)
- Parallel processing capabilities
- Scalable from laptop to cluster computing

## Quick Start

```julia
using StochasticGene

# Set up directory structure
rna_setup("scRNA")

# Fit a simple two-state model
fits, stats, measures, data, model, options = fit(nchains=4)
```

## Documentation Structure

- [Installation](@ref): How to install StochasticGene.jl
- [Getting Started](@ref): Basic usage and examples
- [API Reference](@ref): Detailed documentation of all functions and types
- [Examples](@ref): More complex usage examples
- [Contributing](@ref): How to contribute to the project

## System Requirements

- Julia version 1.9.3 or higher
- Required packages will be automatically installed

## Support

- Documentation: This site
- GitHub Issues: [Report bugs or request features](https://github.com/nih-niddk-mbs/StochasticGene.jl/issues)
- Maintainers: Contact the package maintainers for specific questions

## Citation

If you use StochasticGene.jl in your research, please cite:
- Rodriguez et al. Cell (2018)
- Wan et al. Cell (2021)
- Trzaskoma et al. Science Advances (2024)

## License

This package is licensed under the MIT License. 