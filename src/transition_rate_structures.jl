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
    Element

Single entry of a transition rate matrix, parameterized by position and sign.

# Fields
- `a::Int`: row index.
- `b::Int`: column index.
- `index::Int`: index into the flat rate vector `r`.
- `pm::Int`: sign (`+1` or `-1`); off-diagonal entries are `+1`, diagonal accumulation uses `-1`.
"""
struct Element
    a::Int
    b::Int
    index::Int
    pm::Int
end

"""
    Indices

Index ranges into the flat rate vector `r` for each class of transition.

# Fields
- `gamma::Vector{Int}`: indices for G-state (promoter) transitions.
- `nu::Vector{Int}`: indices for R-state (RNA/elongation) transitions.
- `eta::Vector{Int}`: indices for splice transitions.
- `decay::Int`: index for the mRNA decay rate.
"""
struct Indices
    gamma::Vector{Int}
    nu::Vector{Int}
    eta::Vector{Int}
    decay::Int
end

"""
    CoupledDTCache

Mutable cache for reusing per-unit coupled transition matrices and steady-state
distributions across dwell-time computations within a single likelihood evaluation.

# Fields
- `TC::Dict`: cached coupled transition matrices keyed by unit index.
- `pss::Dict`: cached steady-state distributions keyed by unit index.
"""
mutable struct CoupledDTCache
    TC::Dict
    pss::Dict
end

"""
    AbstractComponents

Root abstract type for all model component structures that store pre-built
matrix elements and metadata for constructing transition matrices.
"""
abstract type AbstractComponents end

"""
    AbstractTComponents

Abstract type for components that produce a single transition matrix T.
"""
abstract type AbstractTComponents <: AbstractComponents end

"""
    AbstractTDComponents

Abstract type for components that produce both a transition matrix T and one or
more dwell-time matrices TD (T restricted to sojourn-state transitions).
"""
abstract type AbstractTDComponents <: AbstractTComponents end

"""
    AbstractMComponents

Abstract type for components that produce an mRNA birth-death matrix M.
"""
abstract type AbstractMComponents <: AbstractComponents end

"""
    AbstractMTComponents

Abstract type for components that produce both an M matrix and a T matrix.
"""
abstract type AbstractMTComponents <: AbstractComponents end

"""
    TComponents

Components for a single GRS transition matrix.

# Fields
- `nT::Int`: dimension of the transition matrix.
- `elementsT::Vector`: element list for building T via `make_mat`.
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
    TAIComponents{elementType}

Transition matrix components with separate element sets for active (A) and
inactive (I) reporter states, used for on/off dwell-time and trace analysis.

# Fields
- `nT::Int`: dimension of the full transition matrix.
- `elementsT::Vector{Element}`: elements for the full T matrix.
- `elementsTA::Vector{elementType}`: elements restricted to active (on) states.
- `elementsTI::Vector{elementType}`: elements restricted to inactive (off) states.
"""
struct TAIComponents{elementType} <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTA::Vector{elementType}
    elementsTI::Vector{elementType}
end
"""
    TDComponents

Transition matrix components for dwell-time analysis. Stores both the full T
elements and pre-filtered TD element sets (one per sojourn/dtype pair).

# Fields
- `nT::Int`: dimension of the T matrix.
- `nG::Int`: dimension of the G (promoter) sub-matrix.
- `elementsT::Vector{Element}`: elements for the full T matrix.
- `elementsG::Vector{Element}`: elements for the G-only sub-matrix (used for G-state dwell times).
- `elementsTD::Vector{Vector{Element}}`: per-dtype filtered element sets; each produces a TD matrix of size `TDdims[i]`.
- `TDdims::Vector{Int}`: matrix dimension for each TD element set.
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
    TRGCoupledUnitComponents

Per-unit components for the RG-stack coupled transition matrix. Stores separate
element sets for the gene, source, target, and RNA-gene sub-matrices.

# Fields
- `nT::Int`: dimension of the full unit transition matrix.
- `nG::Int`: dimension of the G (promoter) sub-matrix.
- `nR::Int`: number of RNA/elongation steps.
- `sourceState::Union{Int, Vector{Int}}`: source state index(es) for coupling.
- `targetTransition::Union{Int, Vector{Int}}`: target transition index(es) for coupling.
- `elementsG::Vector{Element}`: elements for the bare G matrix.
- `elementsGt::Vector{Element}`: elements for the target-modified G matrix.
- `elementsGs::Vector{Element}`: elements for the source indicator matrix.
- `elementsRGbar::Vector{Element}`: elements for the RNA-gene matrix without coupling.
- `elementsRG::Vector{Element}`: elements for the full RNA-gene matrix with coupling.
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
"""
    TCoupledUnitComponents{sourceType, targetType}

Per-unit components for a general coupled transition matrix. Parameterized on the
source/target types to accommodate scalar, vector, and empty coupling specs.

# Fields
- `nT::Int`: dimension of the unit transition matrix.
- `sourceState::sourceType`: source state(s) for coupling (can be `Int`, `Vector{Int}`, etc.).
- `targetTransition::targetType`: target transition(s) for coupling.
- `elementsT::Vector{Element}`: elements for the uncoupled T matrix.
- `elementsSource::Vector{Element}`: elements for the source indicator.
- `elementsTarget::Vector{Element}`: elements for the target transition modifier.
"""
struct TCoupledUnitComponents{sourceType, targetType} <: AbstractTComponents
    nT::Int
    sourceState::sourceType
    targetTransition::targetType
    elementsT::Vector{Element}
    elementsSource::Vector{Element}
    elementsTarget::Vector{Element}
end

"""
    TDCoupledUnitComponents

Per-unit dwell-time components for the RG-stack coupled system. Extends
`TRGCoupledUnitComponents` with filtered TD element sets.

# Fields
- `nT::Int`: dimension of the full unit transition matrix.
- `nG::Int`: dimension of the G sub-matrix.
- `elementsSource::Vector{Element}`: source indicator elements.
- `elementsTarget::Vector{Element}`: target transition elements.
- `elementsT::Vector{Element}`: elements for the full T matrix.
- `elementsG::Vector{Element}`: elements for the G sub-matrix.
- `elementsTD::Vector{Vector{Element}}`: per-dtype filtered element sets.
- `TDdims::Vector{Int}`: matrix dimension for each TD element set.
"""
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
    TCoupledComponents{ModelType}

Top-level coupled transition matrix components for the RG-stack (Kronecker) path.
Holds per-unit component structs and coupling metadata.

# Fields
- `N::Int`: full coupled state-space dimension.
- `model::Tuple`: unit model ordering (indices into the model parameter arrays).
- `sources::Tuple`: per-unit source state lists derived from the connection list.
- `modelcomponents::ModelType`: per-unit component structs (e.g. `Vector{TCoupledUnitComponents}`).
- `connections::Union{Nothing, Vector{ConnectionRecord}}`: optional precomputed per-connection
  data; `nothing` when using the legacy per-unit path.
"""
struct TCoupledComponents{ModelType} <: AbstractComponents
    N::Int
    model::Tuple
    sources::Tuple
    modelcomponents::ModelType
    connections::Union{Nothing,Vector{ConnectionRecord}}
end


"""
    IndexCoupledFull

Typed rate index for the full-space coupled stack. Identifies a rate by its
model and local position within that model's rate vector.

# Fields
- `model::Int`: model index (indexes into the per-unit `rates` vector-of-vectors).
- `localindex::Int`: position within `rates[model]`.
"""
struct IndexCoupledFull
    model::Int
    localindex::Int
end

"""
    ElementCoupledFull

Single entry of a full-space coupled transition matrix, with typed rate indexing.
Used exclusively by `TCoupledFullComponents` and `TDCoupledFullComponents`.

# Fields
- `a::Int`: row index.
- `b::Int`: column index.
- `idx::IndexCoupledFull`: identifies which rate vector and position to use.
- `pm::Int8`: sign (`+1` or `-1`).
"""
struct ElementCoupledFull
    a::Int
    b::Int
    idx::IndexCoupledFull
    pm::Int8
end


"""
    TCoupledFullComponents

Coupled transition matrix represented by full N×N matrix elements (no Kronecker
at compute time), with uncoupled and coupling elements stored separately.

- `N::Int`: Full state dimension.
- `elements_base::Vector{ElementCoupledFull}`: Uncoupled contributions.
- `elements_coupling::Vector{ElementCoupledFull}`: Coupling contributions.
"""
struct TCoupledFullComponents <: AbstractComponents
    N::Int
    elements_base::Vector{ElementCoupledFull}
    elements_coupling::Vector{ElementCoupledFull}
    targets::Vector{Tuple{Int,Int}}
end

"""
    TDCoupledFullComponents

Coupled dwell-time components in full N×N state space. Extends `TCoupledFullComponents`
with dwell-time element subsets: each TD matrix is T with transitions to non-sojourn
states removed, stored as filtered element vectors at construction time.

- `N::Int`: Full state dimension.
- `elements_base::Vector{ElementCoupledFull}`: Uncoupled T elements (model, localindex).
- `elements_coupling::Vector{ElementCoupledFull}`: Coupling T elements (localindex = k).
- `elementsTD_base::Vector{Vector{Vector{ElementCoupledFull}}}`: [unit][dtype] filtered base elements.
- `elementsTD_coupling::Vector{Vector{Vector{ElementCoupledFull}}}`: [unit][dtype] filtered coupling elements.
- `targets::Vector{Tuple{Int,Int}}`: For coupling k: `(model, localindex)` of target base rate;
  used to compute `coupling_rates[k] = γ_k * rates[model][localindex]`.
"""
struct TDCoupledFullComponents <: AbstractComponents
    N::Int
    elements_base::Vector{ElementCoupledFull}
    elements_coupling::Vector{ElementCoupledFull}
    elementsTD_base::Vector{Vector{Vector{ElementCoupledFull}}}
    elementsTD_coupling::Vector{Vector{Vector{ElementCoupledFull}}}
    targets::Vector{Tuple{Int,Int}}
end

"""
    TForcedComponents

Transition matrix components for a model with externally forced transitions.

# Fields
- `nT::Int`: dimension of the transition matrix.
- `elementsT::Vector{Element}`: elements for the base T matrix.
- `targets::Tuple`: target state/transition specifications for the forced transitions.
"""
struct TForcedComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    targets::Tuple
end




"""
    MComponents

Components for the mRNA birth-death matrix M used in steady-state RNA
distribution calculations.

# Fields
- `elementsT::Vector{Element}`: elements for the promoter/gene transition sub-matrix.
- `elementsB::Vector{Element}`: elements for the burst/birth sub-matrix.
- `nT::Int`: dimension of the T sub-matrix.
- `U::SparseMatrixCSC`: full mRNA coupling matrix.
- `Uminus::SparseMatrixCSC`: degradation part of U.
- `Uplus::SparseMatrixCSC`: production part of U.
"""
struct MComponents <: AbstractMComponents
    elementsT::Vector{Element}
    elementsB::Vector{Element}
    nT::Int
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end



"""
    MTComponents{MType, Ttype}

Generic combined M and T components. Parameterized to support any compatible
pair of M and T component types.

# Fields
- `mcomponents::MType`: components for the M (mRNA) matrix.
- `tcomponents::Ttype`: components for the T (transition) matrix.
"""
struct MTComponents{MType,Ttype} <: AbstractMTComponents
    mcomponents::MType
    tcomponents::Ttype
end

"""
    MT0Components

Combined M and T components for the base GRS model (no on/off state separation).

# Fields
- `mcomponents::MComponents`: components for the M matrix.
- `tcomponents::TComponents`: components for the T matrix.
"""
struct MT0Components <: AbstractMTComponents
    mcomponents::MComponents
    tcomponents::TComponents
end

"""
    MTAIComponents{elementType}

Combined M and TAI components for models with separate active/inactive reporter states.

# Fields
- `mcomponents::MComponents`: components for the M matrix.
- `tcomponents::TAIComponents{elementType}`: components for the T matrix with A/I separation.
"""
struct MTAIComponents{elementType} <: AbstractMTComponents
    mcomponents::MComponents
    tcomponents::TAIComponents{elementType}
end

"""
    MTDComponents

Combined M and TD components for dwell-time analysis models.

# Fields
- `mcomponents::MComponents`: components for the M matrix.
- `tcomponents::TDComponents`: components for the T and TD matrices.
"""
struct MTDComponents <: AbstractMTComponents
    mcomponents::MComponents
    tcomponents::TDComponents
end


