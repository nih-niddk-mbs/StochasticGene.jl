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

abstract type AbstractRNAData{hType} <: HistogramData end

"""
Data structures
"""
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    name::String
    gene::String
    nRNA::nType
    histRNA::hType
end
struct TransientRNAData{nType,tType,hType} <: AbstractRNAData{hType}
    name::String
    gene::String
    nRNA::nType
    time::tType
    histRNA::hType
end
struct RNAMixedData{hType} <: AbstractRNAData{hType}
    name::String
    gene::String
    nRNA::Array
    fish::Array{Bool,1}
    histRNA::hType
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

struct GMrescaledmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMmultimodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMdelaymodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMlossmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMmixedmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
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


# Functions for saving and loading data and models

"""
write_log(file,datafile,data,model)
write all information necessary for rerunning
"""
function save_data(file::String,data::TransientRNAData)
    f = open(file,"w")
    writedlm(f, [typeof(data)])
    writedlm(f,[data.name])
    writedlm(f,[data.gene])
    writedlm(f,[data.nRNA])
    writedlm(f,[data.time])
    writedlm(f,[data.histRNA])
    close(f)
end


function load_data(file::String,model::AbstractGMmodel)


end

function save_model(file::String,model::AbstractGMmodel)
    f = open(file,"w")
    write(f, model.G)
    writedlm(f,model.nalleles)
    writedlm(f,model.ratepriors)
    writedlm(f,model.proposal)
    writedlm(f,model.fittedparam)
    writedlm(f,model.method)
    close(f)

end

function load_model(file::String,model::AbstractGRMmodel)


end
