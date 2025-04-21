# This file is part of StochasticGene.jl   

###  common.jl

# Data types
"""
    AbstractExperimentalData

Abstract type for experimental data
"""
abstract type AbstractExperimentalData end
"""
    AbstractSampleData

Abstract type for data in the form of samples
"""
abstract type AbstractSampleData <: AbstractExperimentalData end
"""
    AbstractHistogramData

Abstract type for data in the form of a histogram (probability distribution)
"""
abstract type AbstractHistogramData <: AbstractExperimentalData end
"""
    AbstractTraceData
    
Abstract type for intensity time series data
"""
abstract type AbstractTraceData <: AbstractExperimentalData end

"""
  AbstractRNAData{hType}

Abstract type for steady state RNA histogram data
"""
abstract type AbstractRNAData{hType} <: AbstractHistogramData end


"""
    AbstractTraceHistogramData
    
Abstract type for intensity time series data with RNA histogram data
"""

abstract type AbstractTraceHistogramData <: AbstractTraceData end

# Data structures 
#
# Do not use underscore "_" in label

"""
    RNAData{nType,hType}

Structure for storing RNA histogram data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::nType`: Length of the histogram.
- `histRNA::hType`: RNA histograms.
"""
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    label::String
    gene::String
    nRNA::nType
    histRNA::hType
end

struct RNACountData <: AbstractRNAData{Vector{Int}}
    label::String
    gene::String
    nRNA::Int
    countsRNA::Vector{Int}
    yieldfactor::Vector{Float64}
end

struct DwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    bins::Vector
    DwellTimes::Vector
    DTtypes::Vector
end

"""
    RNAOnOffData

Structure for storing RNA ON/OFF time data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::Int`: Length of the histogram.
- `histRNA::Vector`: RNA histograms.
- `bins::Vector`: Number of live cell recording time bins.
- `ON::Vector`: ON time probability density.
- `OFF::Vector`: OFF time probability density.
"""
struct RNAOnOffData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Vector
    bins::Vector
    ON::Vector
    OFF::Vector
end
"""
    RNADwellTimeData

Structure for storing RNA dwell time data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::Int`: Length of the histogram.
- `histRNA::Array`: RNA histograms.
- `bins::Vector{Vector}`: Number of live cell recording time bins.
- `DwellTimes::Vector{Vector}`: Dwell times.
- `DTtypes::Vector`: Types of dwell times.
"""
struct RNADwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Array
    bins::Vector{Vector}
    DwellTimes::Vector{Vector}
    DTtypes::Vector
end
"""
    TraceData{labelType,geneType,traceType}

Structure for storing trace data.

# Fields
- `label::labelType`: Label for the data set.
- `gene::geneType`: Gene name (case sensitive).
- `interval::Float64`: Time interval between trace points.
- `trace::traceType`: Trace data.
"""
struct TraceData{labelType,geneType,traceType} <: AbstractTraceData
    label::labelType
    gene::geneType
    interval::Float64
    trace::traceType
end
"""
    TraceRNAData{traceType,hType}

Structure for storing trace RNA histogram data.

# Fields
- `label::String`: Label for the data set.
"""
struct TraceRNAData{traceType,hType} <: AbstractTraceHistogramData
    label::String
    gene::String
    interval::Float64
    trace::traceType
    nRNA::Int
    histRNA::hType
end

# Model structures



"""
    HMMReporter

Structure for reporters.

# Fields
- `n::Int`: Number of noise parameters.
- `per_state::Vector`: Number of reporters per state.
- `probfn::Function`: Noise distribution function, e.g., `prob_GaussianMixture`.
- `weightind::Int`: Index for mixture model bias parameter (restricted to range [0,1]).
- `offstates::Vector{Int}`: Vector of off states.
"""
struct HMMReporter
    n::Int
    per_state::Vector
    probfn::Function
    weightind::Int
    offstates::Vector{Int}
    noiseparams::Vector{Int}
end



"""
    Abstract model types
"""
abstract type AbstractModel end
abstract type AbstractGeneTransitionModel <: AbstractModel end
abstract type AbstractGMmodel <: AbstractGeneTransitionModel end
abstract type AbstractGRSMmodel{TraitType} <: AbstractGeneTransitionModel end



"""
    Model structures

fields:
- `rates`: transition rates
- `Gtransitions`:  tuple of vectors of G state transitions
- `G`: number of G steps
- `R`: number of R steps
- `S`: indicator for splicing, 0 no splicing, > 1 splicing
- `insertstep`: R step where reporter is inserted (first step where reporter is visible)
- `nalleles`: number of alleles producing RNA
- `splicetype`: choices are "", "offeject", "offdecay"
- `rateprior`: prior distribution for rates
- `proposal`: MCMC proposal distribution
- `fittedparam`: indices of rates to be fitted
- `fixedeffects`: indices of rates that are fixed to each other, in the form of a 2 tuple of vectors
    with index 1 the tied index vector and 2 the corresponding fitted index vector
- `fixedeffects`: tuple of vectors of rates that are locked together
- `method`: method option, for nonhierarchical models 1 indicates solving Master equation directly, otherwise by eigendecomposition, 
            for hierarchical models, 2-tuple, where 1st component is same as above and 2nd is Bool where true means rates are fixed for all individuals
-` reporter`: vector of reporters or sojorn states (onstates) or vectors of vectors depending on model and data

"""

"""
    GMmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GM models.

# Fields
- `rates::RateType`: Transition rates.
- `Gtransitions::Tuple`: Tuple of vectors of G state transitions.
- `G::Int`: Number of G steps.
- `nalleles::Int`: Number of alleles producing RNA.
- `rateprior::PriorType`: Prior distribution for rates.
- `proposal::ProposalType`: MCMC proposal distribution.
- `fittedparam::ParamType`: Indices of rates to be fitted.
"""
struct GMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGMmodel
    rates::RateType
    Gtransitions::Tuple
    G::Int
    nalleles::Int
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end



"""
Trait system
each trait keeps track of its indices
model keeps track of number of model rates
reporter keeps track of noise parameter indices


nallparams it total number of parameters per individual of all types including transition rates, noise parameters, coupling parameters, and grid parameters
nrates is number of transition rates per individual
nparams is number of fitted parameters per individual

"""



"""
    HierarchicalTrait

Structure for hierarchical model traits.

# Fields
- `nhypersets::Int`: Number of hyperparameter sets
- `nrates::Int`: Number of transition rates per individual
- `nparams::Int`: Number of fitted parameters per individual
- `nindividuals::Int`: Number of individuals (traces)
- `individualstart::Int`: Starting index for individual parameters
- `paramstart::Int`: Starting index for hyperparameters
- `hyperindices::Vector{Vector}`: Indices for hyperparameters
- `fittedshared::Vector{Int}`: Indices of shared parameters that are fitted
"""
struct HierarchicalTrait
    nhypersets::Int
    nrates::Int
    nparams::Int
    nindividuals::Int
    individualstart::Int
    paramstart::Int
    hyperindices::Vector{Vector}
    fittedshared::Vector{Int}
end

"""
    GridTrait

Structure for grid-based parameter space exploration.

# Fields
- `ngrid::Int`: Number of grid points
- `gridindices::Vector`: Indices for grid parameters
"""
struct GridTrait
    ngrid::Int
    gridindices::Vector
end

"""
    CouplingTrait

Structure for coupled model traits.

# Fields
- `ncoupling::Int`: Number of coupling parameters
- `couplingindices::Vector`: Indices for coupling parameters
"""
struct CouplingTrait
    ncoupling::Int
    couplingindices::Vector
end

"""
    hastrait(model, trait)

Check if a trait is present in a model.

# Arguments
- `model::AbstractModel`: Model to check.
- `trait::Symbol`: Trait to check.

"""
function hastrait(model, trait)
    !isnothing(model.trait) && haskey(model.trait, trait)
end

"""
    GRSMmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM models.

# Fields
- `trait::TraitType`: Trait type.
- `rates::RateType`: Transition rates.
- `nrates::nratesType`: Number of transition rates.
- `Gtransitions::Tuple`: Tuple of vectors of G state transitions.
- `G::Int`: Number of G steps.
- `R::Int`: Number of R steps.
- `S::Int`: Indicator for splicing, 0 means no splicing, > 1 means splicing.
- `insertstep::Int`: R step where reporter is inserted (first step where reporter is visible).
- `nalleles::Int`: Number of alleles producing RNA.
- `splicetype::String`: Choices are "", "offeject", "offdecay".
- `rateprior::PriorType`: Prior distribution for rates.
- `proposal::ProposalType`: MCMC proposal distribution.
- `fittedparam::ParamType`: Indices of rates to be fitted.
- `fixedeffects::Tuple`: Indices of rates that are fixed to each other, in the form of a 2-tuple of vectors with index 1 being the tied index vector and 2 being the corresponding fitted index vector.
- `method::MethodType`: Method option, for non-hierarchical models 1 indicates solving Master equation directly, otherwise by eigendecomposition.
- `components::ComponentType`: Components of the model.
- `reporter::ReporterType`: Vector of reporters or sojourn states (onstates) or vectors of vectors depending on model and data.
"""

struct GRSMmodel{TraitType,RateType,nratesType,GType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{TraitType}
    trait::TraitType
    rates::RateType
    ratetransforms::Vector{Function}
    nrates::nratesType
    Gtransitions::Tuple
    G::GType
    R::GType
    S::GType
    insertstep::GType
    nalleles::Int
    splicetype::String
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end


"""
    print_model(model::AbstractModel)

Print all fields of model

# Arguments
- `model::AbstractModel`: Model to print.

# Returns
- `nothing`
"""
function print_model(model::AbstractModel)
    for fname in fieldnames(model)
        println("$fname =", getfield(model, fname))
    end
end

"""
Abstract Option types for fitting methods
"""
abstract type Options end
abstract type Results end



# === Behavioral trait functions ===
#
# These functions provide soft behavioral tags for data types
# without relying on type inheritance. They help determine which
# methods are valid for a given data type without enforcing hard
# subtyping relationships.
#
# Usage example:
#     if is_deviance_supported(data)
#         println("Deviance: ", deviance(fits, data, model))
#     end

"""
    is_histogram_compatible(data) -> Bool

Return `true` if the data can be treated as a histogram (e.g. for
plotting, normalization, or PDF comparisons). 

Used to prevent histogram methods from being accidentally applied to
non-histogram types like `RNACountData`.
"""
is_histogram_compatible(::AbstractExperimentalData) = false
is_histogram_compatible(::AbstractHistogramData) = true
is_histogram_compatible(::RNACountData) = false




