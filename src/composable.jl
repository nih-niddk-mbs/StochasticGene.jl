# This file is part of StochasticGene.jl

"""
    Module for composable type system that allows mixing and matching of different model capabilities.
"""

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
    CouplingTrait

Trait for models that have coupling between units.

# Fields
- `source_models::Tuple`: Model indices for each unit (e.g., (1,1,2) means units 1 and 2 use model 1, unit 3 uses model 2)
- `source_units::Tuple`: Source units for each unit (e.g., ([2,3], [1], Int[]) means unit 1 influenced by 2,3; unit 2 by 1; unit 3 by none)
- `source_states::Tuple`: Source states (e.g., (3,0) means model 1 influences in G state 3, model 2 doesn't influence)
- `target_transitions::Tuple`: Target transitions (e.g., (0,4) means model 1 not influenced, model 2 influenced at G transition 4)
- `n_coupling_params::Int`: Number of coupling parameters
"""
struct CouplingTrait <: ModelTrait 
    source_models::Tuple
    source_units::Tuple
    source_states::Tuple
    target_transitions::Tuple
    n_coupling_params::Int
end

"""
    HierarchicalTrait

Trait for hierarchical models that share parameters across different levels.

# Fields
- `hierarchy::Hierarchy`: Hierarchy of models for hierarchical modeling
"""
struct HierarchicalTrait <: ModelTrait
    hierarchy::Hierarchy
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
struct GridTrait <: ModelTrait
    rate_range::UnitRange
    noise_range::UnitRange
    grid_range::UnitRange
    n_grid::Int
end

#=
Data Traits
=#

"""
    Abstract type for data traits that define different data types.
"""
abstract type DataTrait end

"""
    RNATrait

Trait for RNA histogram data.

# Fields
- `nRNA::Int`: Length of histogram
- `histRNA::Vector`: RNA histogram data
"""
struct RNATrait <: DataTrait
    nRNA::Int
    histRNA::Vector
end

"""
    OnOffTrait

Trait for ON/OFF time data.

# Fields
- `bins::Vector`: Time bins
- `ON::Vector`: ON time distribution
- `OFF::Vector`: OFF time distribution
"""
struct OnOffTrait <: DataTrait
    bins::Vector
    ON::Vector
    OFF::Vector
end

"""
    DwellTimeTrait

Trait for dwell time data.

# Fields
- `bins::Vector{Vector}`: Time bins for each dwell time type
- `dwell_times::Vector{Vector}`: Dwell time distributions
- `dt_types::Vector`: Types of dwell times
"""
struct DwellTimeTrait <: DataTrait
    bins::Vector{Vector}
    dwell_times::Vector{Vector}
    dt_types::Vector
end

"""
    TraceTrait

Trait for time trace data.

# Fields
- `interval::Float64`: Time interval between measurements
- `trace::Any`: Trace data (type varies based on data format)
"""
struct TraceTrait <: DataTrait
    interval::Float64
    trace::Any
end

#=
Option Traits
=#

"""
    Abstract type for option traits that define different fitting options.
"""
abstract type OptionTrait end

"""
    MCMCTrait

Trait for MCMC fitting options.

# Fields
- `samplesteps::Int`: Number of MCMC sampling steps
- `warmupsteps::Int`: Number of warmup steps
- `annealsteps::Int`: Number of annealing steps
- `temp::Float64`: MCMC temperature
- `tempanneal::Float64`: Annealing temperature
"""
struct MCMCTrait <: OptionTrait
    samplesteps::Int
    warmupsteps::Int
    annealsteps::Int
    temp::Float64
    tempanneal::Float64
end

"""
    TimingTrait

Trait for timing options.

# Fields
- `maxtime::Float64`: Maximum wall time for run
"""
struct TimingTrait <: OptionTrait
    maxtime::Float64
end

#=
Composable Types
=#

"""
    ComposableModel{Traits, RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

A flexible model type that can combine different capabilities through traits.

# Fields
- `traits::Traits`: Tuple of model traits defining capabilities
- `rates::RateType`: Transition rates
- `Gtransitions::Tuple`: G state transitions
- `G::Union{Int,Tuple}`: Number of G steps (Int for single model, Tuple for coupled)
- `R::Union{Int,Tuple}`: Number of R steps
- `S::Union{Int,Tuple}`: Number of S steps
- `insertstep::Union{Int,Tuple}`: Reporter insertion step
- `nalleles::Int`: Number of alleles
- `splicetype::String`: Splice action type
- `rateprior::PriorType`: Prior distribution for rates
- `proposal::ProposalType`: MCMC proposal distribution
- `fittedparam::ParamType`: Parameters to be fitted
- `fixedeffects::Tuple`: Fixed effects
- `method::MethodType`: Method type
- `components::ComponentType`: Model components
- `reporter::ReporterType`: Reporter type
"""
struct ComposableModel{Traits, RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    traits::Traits
    rates::RateType
    Gtransitions::Tuple
    G::Union{Int,Tuple}
    R::Union{Int,Tuple}
    S::Union{Int,Tuple} 
    insertstep::Union{Int,Tuple}
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
    ComposableData{Traits}

A flexible data type that can combine different data types through traits.

# Fields
- `traits::Traits`: Tuple of data traits defining capabilities
- `label::String`: Label for the data set
- `gene::String`: Gene name (case sensitive)
"""
struct ComposableData{Traits} <: AbstractExperimentalData
    traits::Traits
    label::String
    gene::String
end

"""
    ComposableOptions{Traits}

A flexible options type that can combine different option types through traits.

# Fields
- `traits::Traits`: Tuple of option traits
"""
struct ComposableOptions{Traits} <: Options
    traits::Traits
end

#=
Helper Functions
=#

# Trait checking
has_trait(model::ComposableModel, T::Type) = any(t isa T for t in model.traits)
get_trait(model::ComposableModel, T::Type) = first(t for t in model.traits if t isa T)

has_data_trait(data::ComposableData, T::Type) = any(t isa T for t in data.traits)
get_data_trait(data::ComposableData, T::Type) = first(t for t in data.traits if t isa T)

has_option_trait(options::ComposableOptions, T::Type) = any(t isa T for t in options.traits)
get_option_trait(options::ComposableOptions, T::Type) = first(t for t in options.traits if t isa T)

# Constructors
"""
    create_model(base_model::GRSMTrait, traits::ModelTrait...; kwargs...)

Create a new composable model with the specified traits and parameters.
"""
function create_model(base_model::GRSMTrait, traits::ModelTrait...; kwargs...)
    all_traits = (base_model, traits...)
    
    # Convert scalar parameters to tuples if coupling is present
    if any(t isa CouplingTrait for t in traits)
        kwargs = convert_to_tuples(kwargs)
    end
    
    ComposableModel(all_traits, kwargs...)
end

"""
    create_data(traits::DataTrait...; label::String, gene::String)

Create a new composable data structure with the specified traits.
"""
function create_data(traits::DataTrait...; label::String, gene::String)
    ComposableData(traits, label, gene)
end

"""
    create_options(traits::OptionTrait...)

Create a new composable options structure with the specified traits.
"""
function create_options(traits::OptionTrait...)
    ComposableOptions(traits)
end

# Trait constructors
coupling(; source_models, source_units, source_states, target_transitions, n_coupling_params) =
    CouplingTrait(source_models, source_units, source_states, target_transitions, n_coupling_params)

hierarchical(; hierarchy) = HierarchicalTrait(hierarchy)

grid(; rate_range, noise_range, grid_range, n_grid) = 
    GridTrait(rate_range, noise_range, grid_range, n_grid)

rna_data(; nRNA::Int, histRNA::Vector) = RNATrait(nRNA, histRNA)

onoff_data(; bins::Vector, ON::Vector, OFF::Vector) = OnOffTrait(bins, ON, OFF)

dwelltime_data(; bins::Vector{Vector}, dwell_times::Vector{Vector}, dt_types::Vector) = 
    DwellTimeTrait(bins, dwell_times, dt_types)

trace_data(; interval::Float64, trace::Any) = TraceTrait(interval, trace)

mcmc_options(; samplesteps::Int, warmupsteps::Int=0, annealsteps::Int=0, 
            temp::Float64=1.0, tempanneal::Float64=100.0) = 
    MCMCTrait(samplesteps, warmupsteps, annealsteps, temp, tempanneal)

timing_options(; maxtime::Float64) = TimingTrait(maxtime)

# Utility functions
function convert_to_tuples(kwargs)
    result = Dict{Symbol,Any}()
    for (k, v) in kwargs
        if k in (:G, :R, :S, :insertstep) && !(v isa Tuple)
            result[k] = (v,)
        else
            result[k] = v
        end
    end
    result
end 