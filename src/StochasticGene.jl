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
include("geneTrap.jl")
include("scrna.jl")


# include("/Users/carsonc/Dropbox/Larson/GeneTrap_analysis/code/GillespieSimulatorExInt.jl")
# include("scRNA.jl")
# include("tff1.jl")



"""
A Julia module for simulation and Bayesian inference of parameters of stochastic models of gene transcription.

API Overview:

"""
# StochasticGene
# """
# #Probability distributions for Bayesian inference:
# y = data, return  = model parameters
#
# prior = p(r)
# likelihood = sum_y p(y|r)
# lprior = log prior = log p(r)
# ll = log likelihood = log p(y|r) =  \sum llC[i]
# llC::Array = log likelihood of components, e.g. LC, FISH, ...
#
# lpost = log posterior = log p(y|theta)p(theta) = ll + lprior
#
# #For WAIC calculations
# lppd = log pointwise predictive density = log E(p(y|theta)) E over MCMC draws
# pWAIC = var(log p(y|theta)) var over MCMC draws
#
# WAIC = -2(llpd - pWAIC)
# AIC = -2(llpd - k), k = number of parameters
# """
# """
# # Output types
# #all draws
# r,llC,lprior,lpost
#
# #draw stats
# (max likelihood r,rlast, E(r), var(r), r quantiles)
# (E(llC), var(llC), llC quantiles)
# accepted, draws
# (E(p), var(log p)   #For WAIC calculation:
# """

end #Module
