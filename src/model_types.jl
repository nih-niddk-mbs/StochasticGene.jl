
# This file is part of StochasticGene.jl   

# model_types.jl

"""
    Pool

Structure for hierarchical model

    # Fields
- `nhyper::Int`: number of hyper parameter sets
- `nparams::Int`: number of fitted hyper params per set = length(fittedparam)
- `nrates`::Int`: number of rates (all params) for each individual
- `nindividualparams::Int`: number of fitted params per individual
- `nindividuals::Int`: number of individuals (traces)
"""
struct Pool
    nhyper::Int
    nrates::Int
    nparams::Int
    nindividuals::Int
    ratestart::Int
    paramstart::Int
    hyperindices::Vector{Vector}
end

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



"""
    Abstract model types
"""
abstract type AbstractModel end
abstract type AbstractGmodel <: AbstractModel end
abstract type AbstractGMmodel <: AbstractGmodel end
abstract type AbstractGRSMmodel{RateType,ReporterType} <: AbstractGmodel end
abstract type AbstractHierarchicalModel{RateType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType} end

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

### To be replaced with new code


"""
    GRSMmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM models.

# Fields
- `rates::RateType`: Transition rates.
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
struct GRSMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    rates::RateType
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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
    GRSMhierarchicalmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM hierarchical models.

# Fields
- `pool::Pool`: Pool of rates for hierarchical modeling.
- as in GRSMmodel
"""
struct GRSMhierarchicalmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    rates::RateType
    pool::Pool
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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
    GRSMcoupledmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM coupled models.

# Fields
- `coupling::CouplingType`: Coupling information for the model.
"""
struct GRSMcoupledmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    rates::RateType
    coupling::CouplingType
    Gtransitions::Tuple
    G::Tuple
    R::Tuple
    S::Tuple
    insertstep::Tuple
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

struct GRSMgridmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    rates::RateType
    raterange::UnitRange
    noiserange::UnitRange
    gridrange::UnitRange
    Ngrid::Int
    coupling::CouplingType
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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

struct GRSMgridhierarchicalmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    rates::RateType
    pool::Pool
    raterange::UnitRange
    noiserange::UnitRange
    gridrange::UnitRange
    Ngrid::Int
    coupling::CouplingType
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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
    GRSMcoupledgridmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM coupled grid models.

# Fields
- `coupling::CouplingType`: Coupling type.
"""
struct GRSMcoupledgridmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    raterange::UnitRange
    noiserange::UnitRange
    gridrange::UnitRange
    Ngrid::Int
    coupling::CouplingType
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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
    GRSMcoupledhierarchicalmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM coupled hierarchical models.

# Fields
- `pool::Pool`: Pool of models.
- `coupling::CouplingType`: Coupling type.
"""
struct GRSMcoupledhierarchicalmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    pool::Pool
    coupling::CouplingType
    raterange::UnitRange
    noiserange::UnitRange
    gridrange::UnitRange
    Ngrid::Int
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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
    GRSMcoupledgridhierarchicalmodel{RateType, CouplingType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM coupled grid hierarchical models.

# Fields
- `pool::Pool`: Pool of models.
- `coupling::CouplingType`: Coupling type.
"""
struct GRSMcoupledgridhierarchicalmodel{RateType,CouplingType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    pool::Pool
    coupling::CouplingType
    raterange::UnitRange
    noiserange::UnitRange
    gridrange::UnitRange
    Ngrid::Int
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
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

