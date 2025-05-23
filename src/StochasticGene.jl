# StochasticGene.jl
# Author: Carson C. Chow

"""
    module StochasticGene

A Julia module for simulation and Bayesian inference of parameters of stochastic models of gene transcription.

# Overview
This module provides a comprehensive set of tools for simulating stochastic models of gene transcription and performing Bayesian inference on model parameters. It includes functions for:

- Transition rate calculations
- Metropolis Hastings MCMC for computing posterior distributions
- Input/output operations
- Chemical master equation solutions for likelihood functions
- Commonly used utility functions
- Probability distributions via direct simulation using Gillespie and Gibson-Bruck algorithms
- Model fitting to data
- Post-fit analysis and plotting
- Functions for use on the NIH cluster Biowulf
- Hidden Markov models
# - GPU-accelerated computations using CUDA (optional)

# Included Files
- `transition_rate_make.jl`: Functions for transition rate calculations.
- `metropolis_hastings.jl`: Metropolis Hastings MCMC for computing posterior distributions of model parameters.
- `io.jl`: Input/output functions.
- `chemical_master.jl`: Chemical master equation solutions of stochastic models for likelihood functions in fitting algorithms.
- `utilities.jl`: Commonly used utility functions.
- `simulator_coupled.jl`: Probability distributions by direct simulation of stochastic models using Gillespie and Gibson-Bruck algorithms.
- `fit.jl`: Functions for fitting models to data.
- `analysis.jl`: Functions for post-fit analysis and plots.
- `biowulf.jl`: Functions for use on NIH cluster Biowulf.
- `hmm.jl`: Functions for hidden Markov models.

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
