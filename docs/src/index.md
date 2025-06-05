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

## Key Features

- **Flexible Model Specification**
  - Support for multiple alleles of the gene in the same cell
  - Configurable reporter visibility at specific R steps
  - Coupling between genes/alleles
  - Support for both exon and intron reporter constructs
  - Classic telegraph models with arbitrary numbers of G states
  - Multiple simultaneous reporters (e.g., for introns and transcription factors)

- **Advanced Capabilities**
  - Parallel processing capabilities
  - Scalable from laptop to cluster computing
  - Bayesian parameter estimation
  - MCMC fitting
  - Comprehensive result analysis

## Quick Start

```julia
using StochasticGene

# Set up directory structure
rna_setup("scRNA")

# Fit a simple two-state model
fits, stats, measures, data, model, options = fit(nchains=4)
```

## Documentation Structure

- [Installation](installation.md): How to install StochasticGene.jl
- [Getting Started](getting_started.md): Basic usage and examples
- [API Reference](api/index.md): Detailed documentation of all functions and types
- [Examples](examples/index.md): More complex usage examples
- [Contributing](contributing.md): How to contribute to the project

## System Requirements

- Julia version 1.9.3 or higher
- Required packages will be automatically installed
- For large datasets: Multiprocessor system recommended

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