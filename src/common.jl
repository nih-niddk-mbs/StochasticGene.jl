"""
Abstract Experimental Data types
"""
abstract type ExperimentalData end
"""
Base for data in the form of samples
"""
abstract type SampleData <: ExperimentalData end
"""
Base for data in the form of a distribution
"""
abstract type HistogramData <: ExperimentalData end

"""
Data structures
"""
struct RNAData <: HistogramData
    name::String
    gene::String
    nRNA::Int
    histRNA::Array
end
struct TransientRNAData{nType,tType,hType} <: HistogramData
    name::String
    gene::String
    nRNA::nType
    time::tType
    histRNA::hType
end
struct RNAMixedData <: HistogramData
    name::String
    gene::String
    nRNA::Int
    type::Array
    time::Array
    histRNA::Array
end
struct LiveCellData <: HistogramData
    name::String
    gene::String
    bins::Array
    OFF::Array
    ON::Array
end
struct RNALiveCellData <: HistogramData
    name::String
    gene::String
    bins::Array
    OFF::Array
    ON::Array
    nRNA::Int
    histRNA::Array
end
struct MultiRNALiveCellData <: HistogramData
    name::String
    gene::String
    nRNA::Array
    histRNA::Tuple
    bins::Tuple
    ON::Tuple
    OFF::Tuple
end

"""
Abstract model types
"""
abstract type Model end
abstract type StochasticGRmodel <: Model end
abstract type AbstractGRMmodel <: StochasticGRmodel end
abstract type AbstractGMmodel <: StochasticGRmodel end

"""
Model structures
"""
struct GMmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMLossmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMMixedmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GRMmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGRMmodel
    G::Int
    R::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

# struct GRSMmodel <: StochasticGRModel
#     G::Int
#     R::Int
#     nalleles::Int
#     type::String
#     rates::Vector
#     rateprior::Vector
#     proposal::Vector
#     fittedparam::Vector
#     method::Int
# end


struct GRSMmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGRMmodel
    G::Int
    R::Int
    nalleles::Int
    type::String
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

function write_model(model)


end

"""
Abstract Option types for fitting methods
"""
abstract type Options end

abstract type Results end
