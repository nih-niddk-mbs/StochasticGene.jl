# StochasticGene.jl
# Author: Carson C. Chow

"""
    module StochasticGene

A Julia package for stochastic modeling of gene transcription and Bayesian inference.

# Overview
StochasticGene.jl is a comprehensive Julia package for simulating and fitting stochastic models of gene transcription. It provides tools for:

1. **Model Simulation**
   - Generalized telegraph models (GRSM)
   - Multiple gene states (G)
   - Pre-RNA steps (R)
   - Splice sites (S)
   - Reporter insertion
   - Multiple alleles
   - Coupled gene models

2. **Parameter Inference**
   - Bayesian parameter estimation via MCMC
   - Maximum likelihood estimation
   - Hierarchical modeling
   - Parallel processing support
   - Adaptive proposal distributions

3. **Data Types**
   - mRNA count distributions (smFISH, scRNA-seq)
   - Live cell imaging traces
   - Dwell time distributions
   - Combined data types
   - Time-series data

4. **Analysis Tools**
   - Model fitting and comparison
   - Hidden Markov model analysis
   - Burst size analysis
   - ON/OFF state analysis
   - Model diagnostics
   - Posterior analysis

# Key Features

## Model Flexibility
- **Gene States (G)**: Arbitrary number of states
- **Pre-RNA Steps (R)**: Multiple elongation steps
- **Splice Sites (S)**: Support for splicing dynamics
- **Coupled Models**: Multiple interacting alleles/genes
- **Hierarchical Models**: Share information across conditions

## Inference Methods
- **MCMC Sampling**: Robust parameter estimation
- **Adaptive Proposals**: Efficient exploration of parameter space
- **Parallel Chains**: Multi-core and distributed computing
- **Convergence Diagnostics**: Ensure reliable results

## Performance
- **Optimized Backend**: For fast simulations
- **Automatic Differentiation**: For gradient-based methods
- **Memory Efficiency**: Handles large datasets

# Quick Start

```julia
using StochasticGene

# Set up a simple two-state model
fits = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/example_data/",
    gene = "MYC"
)
```

# Module Structure

- `transition_rate_make.jl`: Transition rate calculations
- `metropolis_hastings.jl`: MCMC parameter estimation
- `io.jl`: Input/output operations
- `chemical_master.jl`: Chemical master equation solutions
- `utilities.jl`: Common utility functions
- `simulator_coupled.jl`: Stochastic simulation algorithms
- `fit.jl`: Model fitting functions
- `analysis.jl`: Post-fit analysis tools
- `biowulf.jl`: NIH Biowulf cluster support
- `hmm.jl`: Hidden Markov model functions

# Dependencies

- `CSV`: Data handling
- `DataFrames`: Data manipulation
- `Distributed`: Parallel processing
- `Distributions`: Statistical distributions
- `LSODA`: ODE solvers
- `MultivariateStats`: Statistical analysis
- `Optim`: Optimization algorithms
- `Plots`: Visualization
- `ProgressMeter`: Progress tracking
- `StatsBase`: Statistical functions

# Documentation

For detailed usage, see the [documentation](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/).
"""
module StochasticGene
# __precompile__(true)

using CSV
using DataFrames
using Dates
using DelimitedFiles
using Distributed
using Distributions
using Downloads
using FFTW
using JSON
using LinearAlgebra
using LSODA
using MultivariateStats
using Optim
using OrdinaryDiffEq
using Printf
using SparseArrays
using Statistics
using StatsBase
# using CUDA

export
    assemble_all,
    assemble_measures_model,
    CSV,
    DataFrame,
    datapdf,
    digit_vector,
    fix,
    fix_filenames,
    find_best_models,
    fit,
    folder_path,
    folder_setup,
    get_param,
    get_rates,
    GMmodel,
    GRSMmodel,
    AbstractGRSMmodel,
    AbstractGeneTransitionModel,
    AbstractGMmodel,
    large_deviance,
    large_rhat,
    load_data,
    load_model,
    loglikelihood,
    lsoda,
    make_array,
    make_dataframes,
    make_mat,
    makeswarm,
    makeswarm_genes,
    makeswarm_models,
    make_traces,
    make_traces_dataframe,
    mean,
    mean_elongationtime,
    metropolis_hastings,
    ModelArgs,
    norm,
    normalize_histogram,
    num_all_parameters,
    num_rates,
    on_states,
    predictedarray,
    predictedfn,
    prepare_rates,
    plot_empirical_vs_theory,
    prob_Gaussian,
    prob_Gaussian_grid,
    CorrelationAlgorithm,
    CorrelationTrait,
    CorrelationAlgorithm,
    StandardCorrelation,
    WindowedCorrelation,
    MultiTauCorrelation,
    IDLCorrelation,
    DEFAULT_CORRELATION_ALGORITHM,
    hastrait,
    readfile,
    RNAData,
    score_models_from_traces,
    RNADwellTimeData,
    RNAOnOffData,
    readrates,
    rna_setup,
    run_mh,
    run_mcmc_parallel,
    set_indices,
    set_elements_TCoupledUnit,
    simulate_trace,
    simulate_trace_data,
    simulate_trace_vector,
    simulate_trials,
    simulator,
    source_states,
    sparse,
    summarize_model_scores,
    test_compare,
    test_compare_coupling,
    test_fit_rna,
    test_fit_rnadwelltime,
    test_fit_rnaonoff,
    test_fit_simrna,
    test_fit_trace,
    test_fit_trace_hierarchical,
    test_fit_tracejoint,
    test_fit_tracejoint_hierarchical,
    TComponents,
    T_dimension,
    TraceData,
    TraceRNAData,
    Tsit5,
    unbiased_crosscov,
    crosscorrelation_function,
    compute_cov_empirical,
    export_per_trace_crosscorrelations,
    write_augmented,
    write_cov,
    write_cov_empirical,
    write_dataframes,
    write_dataframes_only,
    write_histograms,
    write_ONOFFhistograms,
    write_residency_G_folder,
    write_RNAhistogram,
    write_traces,
    write_traces_coupling,
    write_traces_coupling_spawn,
    write_winners,
    zero_median

    

### Source files

# Type system and common functions
include("common.jl")

# Transition rate matrices of stochastic models defining master equations
include("transition_rate_structures.jl")
include("transition_rate_elements.jl")
include("transition_rate_functions.jl")
include("transition_rate_make.jl")

# Likelihood functions for fitting algorithms
include("likelihoods.jl")

# Metropolis Hastings MCMC for computing posterior distributions of model parameters
include("metropolis_hastings.jl")

# Input output functions
include("io.jl")

# Chemical master equation solutions of stochastic models for likelihood functions in fitting algorithms
include("chemical_master.jl")

# commonly used functions
include("utilities.jl")

# Probability distributions by direct simulation of stochastic models using Gillespie and Gibson-Bruck algorithms
# include("simulator.jl")
include("simulator_coupled.jl")

# functions for fitting models to data
include("fit.jl")

# functions for post fit analysis and plots
include("analysis.jl")

# functions for use on NIH cluster Biowulf
include("biowulf.jl")

# functions for hidden markov models
include("hmm.jl")

# test functions
include("test.jl")


"""
A Julia module for simulation and Bayesian inference of parameters of stochastic models of gene transcription.

API Overview:


"""

end #Module
