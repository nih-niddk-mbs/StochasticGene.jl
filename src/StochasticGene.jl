
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
using DifferentialEquations
using LSODA
using DataFrames
using FFTW
using Downloads
using CSV
using MultivariateStats

export
    rna_setup,
    makeswarm,
    fit_rna,
    make_dataframes,
    write_dataframes,
    write_dataframes_only,
    write_winners,
    write_augmented,
    fix_filenames,
    fix,
    large_deviance,
    assemble_all,
    write_histograms,
    plot_histogram,
    simulatorGM

### Source files

# Type system
include("common.jl")

# commonly used functions
include("utilities.jl")

# Metropolis Hastings MCMC for computing posterior distributions of model parameters
include("metropolis_hastings.jl")

# Input output functions
include("io.jl")

# Chemical master equation solutions of stochastic models for likelihood functions in fitting algorithms
include("chemical_master.jl")

# Transition rate matrices of stochastic models defining master equations
include("transition_rate_matrices.jl")

# Probability distributions by direct simulation of stochastic models using Gillespie and Gibson-Bruck algorithms
include("simulator.jl")
include("telegraphsplice.jl")

# functions specific for Gene Trap experiments of Wan et al.
include("genetrap.jl")

# functions for scRNA and FISH experiments
include("rna.jl")

# functions for post fit analysis and plotting
include("analysis.jl")

# functions for use on NIH cluster Biowulf
include("biowulf.jl")


# include("/Users/carsonc/Dropbox/Larson/GeneTrap_analysis/code/GillespieSimulatorExInt.jl")
# include("scRNA.jl")
# include("tff1.jl")



"""
A Julia module for simulation and Bayesian inference of parameters of stochastic models of gene transcription.

API Overview:


"""

end #Module
