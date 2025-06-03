# StochasticGene.jl

**Package version:** {{pkg_version}}  
**Julia version:** {{julia_version}}

StochasticGene.jl is a Julia package for simulating and fitting stochastic models of gene transcription to experimental data. It supports:

- Arbitrary numbers of gene states (G), pre-RNA steps (R), and splice sites (S)
- mRNA count distributions (smFISH, scRNA), live cell imaging traces, and dwell time distributions
- Bayesian parameter estimation and MCMC fitting
- Parallel and hierarchical model fitting

## Features

- Flexible model specification (G, R, S, alleles, splicing)
- Support for multiple experimental data types
- Tools for simulation, fitting, and analysis
- Designed for both small and large-scale (cluster) use

## Quick Start

```julia
using StochasticGene

# Set up directory structure and mock data
rna_setup("scRNA")
cd("scRNA")

# Fit a simple two-state telegraph model to mock data
fits = fit(
    G = 2,  # Number of gene states
    R = 0,  # Number of pre-RNA steps
    transitions = ([1,2], [2,1]),  # Gene state transitions
    datatype = "rna",  # Data type
    datapath = "data/HCT116_testdata/",  # Path to data
    gene = "MYC",  # Gene name
    datacond = "MOCK"  # Data condition
)
```

## Documentation Structure

- [Installation](installation.md): How to install StochasticGene.jl
- [Getting Started](getting_started.md): Basic usage and examples
- [API Reference](api/index.md): Complete API documentation
- [Examples](examples/index.md): Detailed usage examples
- [Contributing](contributing.md): How to contribute to the project

## System Requirements

- Julia 1.9.3 or higher
- Required packages will be automatically installed
- For large datasets: Multiprocessor system recommended

## Support

For questions or issues, please:

1. Check the [documentation](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/)
2. Open an issue on the [GitHub repository](https://github.com/nih-niddk-mbs/StochasticGene.jl/issues)
3. Contact the maintainers directly

## Citation

If you use StochasticGene.jl in your research, please cite:

- Rodriguez, J., et al. (2018). Cell
- Wan, L., et al. (2021). Cell
- Trzaskoma, M., et al. (2024). Science Advances

## License

StochasticGene.jl is open source software licensed under the MIT License.