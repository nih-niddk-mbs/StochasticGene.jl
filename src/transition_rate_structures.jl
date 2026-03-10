# This file is part of StochasticGene.jl
#
# transition_rate_structures.jl
#
# This file defines the core data structures used throughout the StochasticGene.jl package
# for representing transition rate matrices and their components. These structures form
# the foundation for building and manipulating the various types of transition matrices
# needed for stochastic gene expression modeling.
#
# Key structures include:
# - Element: Basic building block for matrix elements with row, column, rate index, and sign
# - Indices: Index ranges for different types of transition rates (gamma, nu, eta, decay)
# - AbstractComponents: Abstract types for different component categories
# - TComponents: Basic transition matrix components
# - TAIComponents: Transition matrix components with on/off state separation
# - TDComponents: Transition matrix components for dwell time analysis
# - Coupled components: Various structures for coupled gene systems
# - MComponents: Matrix components for mRNA distribution calculations
# - MTComponents: Combined M and T matrix components
#
# These structures are used by the other transition rate files to construct the actual
# transition matrices and perform various analyses in the StochasticGene.jl framework.

"""
	struct Element

structure for T transition matrix elements
fields: 
- `a`: row, 
- `b`: column
- `index`: rate vector index
- `pm`: sign of elements

"""
struct Element
    a::Int
    b::Int
    index::Int
    pm::Int
end

"""
	struct Indices

index ranges for rates
gamma: G transitions
nu: R transitions
eta: splice transitions
decay: mRNA decay rate

"""
struct Indices
    gamma::Vector{Int}
    nu::Vector{Int}
    eta::Vector{Int}
    decay::Int
end

mutable struct CoupledDTCache
    TC::Dict
    pss::Dict
end

abstract type AbstractComponents end

"""
	AbstractTComponents

abstract type for components of transition matrices 
"""
abstract type AbstractTComponents <: AbstractComponents end

abstract type AbstractTDComponents <: AbstractTComponents end

abstract type AbstractMComponents <: AbstractComponents end

abstract type AbstractMTComponents <: AbstractComponents end

"""
	struct TComponents

fields:
    nT::Int
    elementsT::Vector

"""
struct TComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector
end

# """
# 	struct TRGComponents

# fields: 
#     nT::Int
#     nG::Int
#     nR::Int
#     elementsG::Vector
#     elementsRGbar::Vector
#     elementsRG::Vector

# """
# struct TRGComponents <: AbstractTComponents
#     nT::Int
#     nG::Int
#     nR::Int
#     elementsG::Vector
#     elementsRGbar::Vector
#     elementsRG::Vector
# end


"""
struct TAIComponents{elementType} <: AbstractTComponents

fields:
    nT::Int
    elementsT::Vector{Element}
    elementsTA::Vector{elementType}
    elementsTI::Vector{elementType}

"""
struct TAIComponents{elementType} <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTA::Vector{elementType}
    elementsTI::Vector{elementType}
end
"""
struct TDComponents <: AbstractTComponents

fields:
nT, elementsT, elementsTD:: Vector{Vector{Element}}

"""
struct TDComponents <: AbstractTDComponents
    nT::Int
    nG::Int
    elementsT::Vector{Element}
    elementsG::Vector{Element}
    elementsTD::Vector{Vector{Element}}
    TDdims::Vector{Int}
end

"""
 	TRGCoupledComponents

fields:
    nT::Int
    nG::Int
    nR::Int
    sourceState::Union{Int, Vector{Int}}
    targetTransition::Union{Int, Vector{Int}}
    elementsG::Vector{Element}
    elementsGt::Vector{Element}
    elementsGs::Vector{Element}
    elementsRGbar::Vector{Element}
    elementsRG::Vector{Element}

"""
struct TRGCoupledUnitComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    sourceState::Union{Int, Vector{Int}}
    targetTransition::Union{Int, Vector{Int}}
    elementsG::Vector{Element}
    elementsGt::Vector{Element}
    elementsGs::Vector{Element}
    elementsRGbar::Vector{Element}
    elementsRG::Vector{Element}
end
struct TCoupledUnitComponents{sourceType, targetType} <: AbstractTComponents
    nT::Int
    sourceState::sourceType
    targetTransition::targetType
    elementsT::Vector{Element}
    elementsSource::Vector{Element}
    elementsTarget::Vector{Element}
end

struct TDCoupledUnitComponents <: AbstractTDComponents
    nT::Int
    nG::Int
    elementsSource::Vector{Element}
    elementsTarget::Vector{Element}
    elementsT::Vector{Element}
    elementsG::Vector{Element}
    elementsTD::Vector{Vector{Element}}
    TDdims::Vector{Int}

end

"""
    ConnectionSpec

Canonical representation of a single coupling connection at the matrix-construction level.

Each connection is identified by a 4-tuple `(β, s, α, t)` where:
- `β::Int`: source unit index
- `s::Int`: source state index within unit β
- `α::Int`: target unit index
- `t::Int`: target transition index within unit α

This flat representation is convenient for fast iteration and indexing when
assembling coupled transition matrices; higher-level code is responsible for
mapping model-level coupling specifications into these connection specifications.
The connection list may be empty (uncoupled T is still built); source and target
info per unit is derived via `source_states_for_unit` and `target_transition_for_unit`.
"""
const ConnectionSpec = NTuple{4,Int}

"""
    source_states_for_unit(connections, unit::Int)

Return the source state(s) for which `unit` acts as a source, derived from the
connection list. Empty if the unit is never a source.

Single source of truth: no separate `source_state` aggregate is stored elsewhere.
"""
function source_states_for_unit(connections, unit::Int)
    [Int(s) for (β, s, α, t) in connections if β == unit]
end

"""
    target_transition_for_unit(connections, unit::Int)

Return the target transition index for `unit` when it is a target (first connection
to that unit), or 0 if the unit is never a target.
"""
function target_transition_for_unit(connections, unit::Int)
    for (β, s, α, t) in connections
        α == unit && return t
    end
    return 0
end

"""
    ConnectionRecord

Precomputed data for one coupling interaction, suitable for direct use in
matrix assembly routines.

Fields:
- `β::Int`: source unit index
- `α::Int`: target unit index
- `U::SparseMatrixCSC`: source operator acting on unit β (built from G/R structure)
- `elementsTarget::Vector{Element}`: element specification for the target operator
  in unit α (used with `make_mat` and the rate vector for unit α)
- `nTα::Int`: size of the target transition matrix for unit α

These records are intended to be built once from model-level information and
then reused whenever coupled transition matrices are constructed.
"""
struct ConnectionRecord
    β::Int
    α::Int
    U::SparseMatrixCSC
    elementsTarget::Vector{Element}
    nTα::Int
end

"""
 	TCoupledComponents

fields:
    N, model, sources, modelcomponents: as before.
    connections: optional precomputed per-connection data (e.g. `Vector{ConnectionRecord}`),
        built from model-level coupling `(unit_model, connections::Vector{ConnectionSpec})`.
    Per-unit source/target are derived from the connection list (e.g. `source_states_for_unit`,
    `target_transition_for_unit`) and stored in modelcomponents[α].sourceState and .targetTransition.

Note: Coupling is supplied as `(unit_model, connections)`; empty `connections` is valid (Tc = uncoupled T).
Higher-level code builds per-connection data from these; matrix constructors in `transition_rate_make.jl`
depend only on the pre-built components.
"""
struct TCoupledComponents{ModelType} <: AbstractComponents
    N::Int
    model::Tuple
    sources::Tuple
    modelcomponents::ModelType
    connections::Union{Nothing,Vector{ConnectionRecord}}
end

"""
    TDCoupledSpecComponent

Per-dwell-spec component for the spec-based coupled dwell path. Holds the sojourn
state set in the space that matches the TC used at prediction time: full coupled
T-space for R-space specs (OFF/ON), coupled G-space for G-space specs (OFFG/ONG).
`is_G_space` selects which TC is built (Gm vs T[unit]) so sojourn and matrix match.
"""
struct TDCoupledSpecComponent
    sojourn::Vector{Int}
    is_G_space::Bool
end

"""
    TDCoupledSpecComponents

Spec-based coupled dwell components: one TDCoupledSpecComponent per DwellSpec,
plus uncoupled full-space elements and connection_data for coupling. Each component's
sojourn set matches the TC space (full T or G). TC is built as make_mat(elements_*, rates_full, N)
plus coupling terms (bilinear) added via connection_data.
"""
struct TDCoupledSpecComponents
    specs::Vector{DwellSpec}
    components::Vector{TDCoupledSpecComponent}
    underlying::TCoupledComponents{Vector{TDCoupledUnitComponents}}
    connection_data::Vector{ConnectionRecord}
    elements_T::Vector{Element}
    nT_T::Int
    elements_G::Vector{Element}
    nT_G::Int
end

struct TForcedComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    targets::Tuple
end




"""
	struct MComponents

structure for matrix components for matrix M

fields:
    elementsT::Vector{Element}
    elementsB::Vector{Element}
    nT::Int
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
"""
struct MComponents <: AbstractMComponents
    elementsT::Vector{Element}
    elementsB::Vector{Element}
    nT::Int
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end



struct MTComponents{MType,Ttype} <: AbstractMTComponents
    mcomponents::MType
    tcomponents::Ttype
end

"""
 	MT0Components

 structure for M and T components

fields:
mcomponents::MComponents
tcomponents::TComponents
"""
struct MT0Components <: AbstractMTComponents
    mcomponents::MComponents
    tcomponents::TComponents
end

"""
 	MTAIComponents

 structure for M and TAI components

fields:
mcomponents::MComponents
tcomponents::TAIComponents
"""
struct MTAIComponents{elementType} <: AbstractMTComponents
    mcomponents::MComponents
    tcomponents::TAIComponents{elementType}
end
"""
 	MTDComponents

 structure for M and TD components

fields:
mcomponents::MComponents
tcomponents::TDComponents
"""
struct MTDComponents <: AbstractMTComponents
    mcomponents::MComponents
    tcomponents::TDComponents
end


