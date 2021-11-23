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
    ExperimentalData,
    SampleData,
    HistogramData,
    AbstractRNAData,
    RNAData,
    TransientRNAData,
    RNAMixedData,
    LiveCellData,
    RNALiveCellData,
    MultiRNALiveCellData,
    Model,
    StochasticGRmodel,
    AbstractGMmodel,
    AbstractGRMmodel,
    AbtractGMlossmodel,
    GModel,
    GMrescaledmodel,
    GMmultimodel,
    GMtransientmodel,
    GMdelaymodel,
    GMlossmodel,
    GMfixedeffectslossmodel,
    GMmixedmodel,
    GRMmodel,
    GRSMmodel,
    Options,
    Results,
    datahistogram,
    datapdf,
    likelihoodfn,
    likelihoodarray,
    likelihoodtuple,
    get_rates,
    get_r,
    get_n,
    get_param,
    setr,
    logprior,
    MHOptions,
    Fit,
    Stats,
    run_mh,
    run_chains,
    metropolis_hastings,
    anneal,
    warmup,
    sample,
    mhstep,
    compute_waic,
    aic,
    initial_proposal,
    proposal,
    mhfactor,
    update_waic,
    extract_chain,
    collate_fit,
    collate_waic,
    pooled_waic,
    merge_param,
    merge_ll,
    merge_fit,
    compute_stats,
    find_ml,
    loglikelihood,
    crossentropy,
    hist_entropy,
    deviance

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
include("ClassicTelegraph.jl")
include("simulator.jl")

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
