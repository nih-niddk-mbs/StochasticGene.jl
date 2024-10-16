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
using Plots
using SparseArrays
using OrdinaryDiffEq
using DataFrames
using FFTW
using Downloads
using CSV
using MultivariateStats
using Optim

export
    CSV,
    DataFrame,
    rna_setup,
    folder_setup,
    makeswarm,
    RNAData,
    RNAOnOffData,
    RNADwellTimeData,
    TraceData,
    TraceRNAData,
    GMmodel,
    GRSMmodel,
    GRSMhierarchicalmodel,
    GRSMcoupledmodel,
    fit,
    mean_elongationtime,
    run_mh,
    metropolis_hastings,
    prob_Gaussian,
    prob_GaussianMixture,
    prob_GaussianMixture_6,
    prob_Gaussian_grid,
    make_dataframes,
    write_dataframes,
    write_dataframes_only,
    write_winners,
    write_augmented,
    write_traces,
    fix_filenames,
    fix,
    large_deviance,
    large_rhat,
    assemble_all,
    assemble_measures_model,
    write_histograms,
    make_ONOFFhistograms,
    plot_histogram,
    make_traces,
    make_traces_dataframe,
    write_traces,
    plot,
    plot!,
    simulator,
    simulate_trace,
    simulate_trace_vector,
    simulate_trace_data,
    on_states,
    num_rates,
    load_data,
    load_model,
    make_components_MTD,
    make_components_MTAI,
    make_components_MT,
    make_components_M,
    make_components_T,
    make_components_TAI,
    make_components_TRG,
    make_components_Tcoupled,
    likelihoodarray,
    likelihoodfn,
    loglikelihood,
    datapdf,
    normalize_histogram,
    norm,
    make_array,
    folder_path,
    get_rates,
    new_FISHfolder,
    ModelArgs
    

### Source files

# Type system and common functions
include("common.jl")

# Transition rate matrices of stochastic models defining master equations
include("transition_rate_structures.jl")
include("transition_rate_elements.jl")
include("transition_rate_functions.jl")
include("transition_rate_make.jl")

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


"""
A Julia module for simulation and Bayesian inference of parameters of stochastic models of gene transcription.

API Overview:


"""

end #Module
