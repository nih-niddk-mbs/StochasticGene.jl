# StochasticGene.jl
# Author: Carson C. Chow

"""
    module StochasticGene

A Julia package for stochastic modeling of gene transcription and Bayesian inference.

# Overview
StochasticGene.jl is a comprehensive Julia package for simulating and fitting stochastic models of gene transcription. It provides tools for:

1. **Model Simulation**
   - Generalized telegraph models (GRSM)
   - Multiple gene states
   - Pre-RNA steps
   - Splice sites
   - Reporter insertion
   - Multiple alleles

2. **Parameter Inference**
   - Bayesian parameter estimation
   - MCMC sampling
   - Hierarchical modeling
   - Parallel processing

3. **Data Types**
   - mRNA count distributions
   - Live cell imaging traces
   - Dwell time distributions
   - Combined data types

4. **Analysis Tools**
   - Model fitting
   - Post-fit analysis
   - Hidden Markov models
   - Burst size analysis
   - ON/OFF state analysis

# Key Features

1. **Flexible Model Architecture**
   - Arbitrary number of gene states
   - Pre-RNA steps
   - Splice sites
   - Reporter insertion
   - Multiple alleles

2. **Advanced Fitting**
   - Bayesian parameter estimation
   - MCMC sampling
   - Hierarchical modeling
   - Parallel processing

3. **Data Compatibility**
   - mRNA count distributions
   - Live cell imaging traces
   - Dwell time distributions
   - Combined data types

4. **Analysis Tools**
   - Model fitting
   - Post-fit analysis
   - Hidden Markov models
   - Burst size analysis
   - ON/OFF state analysis

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
- `Turing`: Bayesian inference

# Usage Example

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

# References

- Rodriguez, J., et al. (2018). Cell
- Wan, L., et al. (2021). Cell
- Trzaskoma, M., et al. (2024). Science Advances

"""
module StochasticGene
__precompile__(true)

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
    new_FISHfolder,
    norm,
    normalize_histogram,
    num_all_parameters,
    num_rates,
    on_states,
    plot,
    plot!,
    plot_histogram,
    predictedarray,
    predictedfn,
    prepare_rates,
    prob_Gaussian,
    prob_GaussianMixture,
    prob_GaussianMixture_6,
    prob_Gaussian_grid,
    readfile,
    RNAData,
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
    write_augmented,
    write_cov,
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
