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

using Distributions
using StatsBase
using Statistics
using DelimitedFiles
using Dates
using Distributed
using LinearAlgebra
# using Plots
using SparseArrays
using OrdinaryDiffEq
using LSODA
using DataFrames
using FFTW
using Downloads
using CSV
using MultivariateStats
using Optim

export
    assemble_all,
    assemble_measures_model,
    CSV,
    DataFrame,
    datapdf,
    fix,
    fix_filenames,
    fit,
    folder_path,
    folder_setup,
    get_rates,
    GMmodel,
    GRSMcoupledmodel,
    GRSMhierarchicalmodel,
    GRSMmodel,
    large_deviance,
    large_rhat,
    likelihoodarray,
    likelihoodfn,
    load_data,
    load_model,
    loglikelihood,
    make_array,
    make_dataframes,
    make_ONOFFhistograms,
    makeswarm,
    make_traces,
    make_traces_dataframe,
    mean,
    mean_elongationtime,
    metropolis_hastings,
    ModelArgs,
    new_FISHfolder,
    norm,
    normalize_histogram,
    num_rates,
    on_states,
    plot,
    plot!,
    plot_histogram,
    prob_Gaussian,
    prob_GaussianMixture,
    prob_GaussianMixture_6,
    prob_Gaussian_grid,
    RNAData,
    RNADwellTimeData,
    RNAOnOffData,
    rna_setup,
    run_mh,
    simulate_trace,
    simulate_trace_data,
    simulate_trace_vector,
    simulator,
    TraceData,
    TraceRNAData,
    write_augmented,
    write_dataframes,
    write_dataframes_only,
    write_histograms,
    write_traces,
    write_traces,
    write_winners
    

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
