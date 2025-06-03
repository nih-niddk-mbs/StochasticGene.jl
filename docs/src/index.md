# StochasticGene.jl

StochasticGene.jl is a Julia package for simulating and fitting stochastic models of gene transcription to experimental data. It supports a variety of data types and model configurations, making it versatile for different biological applications.

## Features

- **Flexible Model Architecture**
  - Arbitrary number of gene states (G)
  - Pre-RNA steps (R)
  - Splice sites (S)
  - Reporter insertion steps
  - Multiple alleles support

- **Data Types Supported**
  - mRNA count distributions (smFISH, scRNA)
  - Image intensity traces (live cell imaging)
  - Dwell time distributions
  - Combined data types

- **Advanced Fitting**
  - Bayesian parameter estimation
  - MCMC sampling
  - Hierarchical modeling
  - Parallel processing support

- **Analysis Tools**
  - Trace simulation
  - ON/OFF dwell time analysis
  - G state residency probabilities
  - Burst size analysis

## Quick Start

To get started with StochasticGene.jl:

```julia
using StochasticGene

# Fit a basic model
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