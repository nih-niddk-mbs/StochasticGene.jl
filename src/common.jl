# This file is part of StochasticGene.jl   

###  common.jl


"""
Abstract Option types for fitting methods
"""
abstract type Options end
abstract type Results end



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
abstract type AbstractTraceHistogramData <: AbstractExperimentalData end

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


# struct Hierarchy
#     nhypersets::Int
#     nrates::Int
#     nparams::Int
#     nindividuals::Int
#     ratestart::Int
#     paramstart::Int
#     hyperindices::Vector{Vector}
# end

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
end

#=
Model Traits
=#

"""
    Abstract type for model traits that define different capabilities.
"""
abstract type ModelTrait end

"""
    Base GRSM trait that all models must have.
"""
struct GRSMTrait <: ModelTrait end



"""
    struct HierarchicalTrait <: ModelTrait

A structure representing hierarchical traits in a model.

# Fields
- `nhypersets::Int`: Number of hyperparameter sets.
- `nrates::Int`: Number of rates (all parameters) for each individual.
- `nparams::Int`: Number of fitted hyperparameters per set.
- `nindividuals::Int`: Number of individuals (traces).
- `ratestart::Int`: Starting index for rates in the parameter vector.
- `paramstart::Int`: Starting index for hyperparameters in the parameter vector.
- `hyperindices::Vector{Vector}`: Indices of hyperparameters for each set.

# Description
The `HierarchicalTrait` struct is used to represent hierarchical traits in a model. It includes information about the number of hyperparameter sets, the number of rates and parameters for each individual, and the indices of hyperparameters for each set. This struct is useful for models that involve hierarchical or multi-level parameters, allowing for more complex and flexible modeling of traits across different individuals or groups.
"""
struct HierarchicalTrait <: GRSMTrait
    nhypersets::Int
    nrates::Int
    nparams::Int
    nindividuals::Int
    ratestart::Int
    paramstart::Int
    hyperindices::Vector{Vector}
end
"""
    GridTrait

Trait for models that use grid-based parameter space exploration.

# Fields
- `rate_range::UnitRange`: Range of rate parameters
- `noise_range::UnitRange`: Range of noise parameters
- `grid_range::UnitRange`: Range of grid points
- `n_grid::Int`: Number of grid points
"""
struct GridTrait <: GRSMTrait
    raterange::UnitRange
    noiserange::UnitRange
    gridrange::UnitRange
    n_grid::Int
end

"""
    Abstract model types
"""
abstract type AbstractModel end
abstract type AbstractGmodel <: AbstractModel end
abstract type AbstractGMmodel <: AbstractGmodel end
abstract type AbstractGRSMmodel{TraitType} <: AbstractGmodel end

# abstract type AbstractGRSMmodel{RateType,ReporterType} <: AbstractGmodel end
# abstract type AbstractGRSMhierarchicalmodel{RateType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType} end



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
    struct GRSMmodel{TraitType, GType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType} <: AbstractGRSMmodel{TraitType}

A structure for GRSM trait models.

# Fields
- `traits::TraitType`: Traits of the model.
- `rates::Vector`: Transition rates.
- `Gtransitions::Tuple`: Tuple of vectors of G state transitions.
- `G::GType`: Number of G steps.
- `R::GType`: Number of R steps.
- `S::GType`: Indicator for splicing, 0 means no splicing, > 1 means splicing.
- `insertstep::GType`: R step where reporter is inserted (first step where reporter is visible).
- `nalleles::Int`: Number of alleles producing RNA.
- `splicetype::String`: Choices are "", "offeject", "offdecay".
- `rateprior::PriorType`: Prior distribution for rates.
- `proposal::ProposalType`: MCMC proposal distribution.
- `fittedparam::ParamType`: Indices of rates to be fitted.
- `fixedeffects::Tuple`: Indices of rates that are fixed to each other, in the form of a 2-tuple of vectors with index 1 being the tied index vector and 2 being the corresponding fitted index vector.
- `method::MethodType`: Method option, for non-hierarchical models 1 indicates solving Master equation directly, otherwise by eigendecomposition.
- `components::ComponentType`: Components of the model.
- `reporter::ReporterType`: Vector of reporters or sojourn states (onstates) or vectors of vectors depending on model and data.

# Description
The `GRSMmodel` struct is used to represent GRSM trait models. It includes information about the traits, transition rates, state transitions, splicing, alleles, priors, proposals, and various components of the model. This struct is useful for modeling gene transcription with stochastic processes and performing Bayesian inference on model parameters.
"""
struct GRSMmodel{TraitType, GType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType} <: AbstractGRSMmodel{TraitType}
    traits::TraitType
    rates::Vector
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

#### obsolete structures


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
# struct GRSMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
#     rates::RateType
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end



# """
#     GRSMhierarchicalmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

# Structure for GRSM hierarchical models.

# # Fields
# - `hierarchy::Hierarchy`: Hierarchy of rates for hierarchical modeling.
# - as in GRSMmodel
# """
# struct GRSMhierarchicalmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMhierarchicalmodel{RateType,ReporterType}
#     rates::RateType
#     hierarchy::Hierarchy
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end

# """
#     GRSMcoupledmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

# Structure for GRSM coupled models.

# # Fields
# - `coupling::CouplingType`: Coupling information for the model.
# """
# struct GRSMcoupledmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
#     rates::RateType
#     coupling::CouplingType
#     Gtransitions::Tuple
#     G::Tuple
#     R::Tuple
#     S::Tuple
#     insertstep::Tuple
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end

# struct GRSMgridmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
#     rates::RateType
#     raterange::UnitRange
#     noiserange::UnitRange
#     gridrange::UnitRange
#     Ngrid::Int
#     coupling::CouplingType
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end

# struct GRSMgridhierarchicalmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMhierarchicalmodel{RateType,ReporterType}
#     rates::RateType
#     hierarchy::Hierarchy
#     raterange::UnitRange
#     noiserange::UnitRange
#     gridrange::UnitRange
#     Ngrid::Int
#     coupling::CouplingType
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end
# """
#     GRSMcoupledgridmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

# Structure for GRSM coupled grid models.

# # Fields
# - `coupling::CouplingType`: Coupling type.
# """
# struct GRSMcoupledgridmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
#     raterange::UnitRange
#     noiserange::UnitRange
#     gridrange::UnitRange
#     Ngrid::Int
#     coupling::CouplingType
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end
# """
#     GRSMcoupledhierarchicalmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

# Structure for GRSM coupled hierarchical models.

# # Fields
# - `hierarchy::Hierarchy`: Hierarchy of models.
# - `coupling::CouplingType`: Coupling type.
# """
# struct GRSMcoupledhierarchicalmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
#     hierarchy::Hierarchy
#     coupling::CouplingType
#     raterange::UnitRange
#     noiserange::UnitRange
#     gridrange::UnitRange
#     Ngrid::Int
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end

# """
#     GRSMcoupledgridhierarchicalmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

# Structure for GRSM coupled grid hierarchical models.

# # Fields
# - `hierarchy::Hierarchy`: Hierarchy of models.
# - `coupling::CouplingType`: Coupling type.
# """
# struct GRSMcoupledgridhierarchicalmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
#     hierarchy::Hierarchy
#     coupling::CouplingType
#     raterange::UnitRange
#     noiserange::UnitRange
#     gridrange::UnitRange
#     Ngrid::Int
#     Gtransitions::Tuple
#     G::Int
#     R::Int
#     S::Int
#     insertstep::Int
#     nalleles::Int
#     splicetype::String
#     rateprior::PriorType
#     proposal::ProposalType
#     fittedparam::ParamType
#     fixedeffects::Tuple
#     method::MethodType
#     components::ComponentType
#     reporter::ReporterType
# end



