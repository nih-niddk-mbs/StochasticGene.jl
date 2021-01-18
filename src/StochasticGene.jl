"""
StochasticGene
"""
module StochasticGene

using Distributions
using StatsBase
using Statistics
using DelimitedFiles
using Dates
using Distributed
using LinearAlgebra
using PyPlot
using SparseArrays
using DifferentialEquations
using LSODA


export
    Model,
    GModel,
    GRModel,
    GRSModel,
    RNALiveCellSplice,
    ExperimentalData,
    HistogramData,
    RNALiveCell,
    runmodel,
    logprior,
    priordistribution,
    likelihoodtuple,
    crossentropy,
    likelihoodfn,
    readparams,
    readFISH,
    readLCPDF,
    runparallel,
    metropolis_hastings,
    sample,
    waicupdate!,
    mhStep!,
    proposal,
    initialproposal,
    sigmalogNormal,
    loglikelihood,
    MHOPtions,
    offtimeCDF,
    telegraphPDF


### Source files

# Type system
include("common.jl")

# commonly used functions
include("utilities.jl")

# Metropolis Hastings MCMC for computing posterior distributions of model parameters
include("metropolis_hastings.jl")

# Chemical master equation solutions of stochastic models for likelihood functions in fitting algorithms
include("chemical_master.jl")

# Transition rate matrices of stochastic models defining master equations
include("transition_rate_matrices.jl")

# Probability distributions by exact and direct simulation of stochastic models using Gillespie algorithms
# include("gillespie.jl")
include("ClassicTelegraph.jl")

# functions specific for Gene Trap experiments of Wan et al.
include("genetrap.jl")

# functions for scRNA and FISH experiments
include("rna.jl")


# include("/Users/carsonc/Dropbox/Larson/GeneTrap_analysis/code/GillespieSimulatorExInt.jl")
# include("scRNA.jl")
# include("tff1.jl")



"""
A Julia module for simulation and Bayesian inference of parameters of stochastic models of gene transcription.

API Overview:


"""

end #Module
