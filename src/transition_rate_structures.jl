# This file is part of StochasticGene.jl
#
# transition_rate_structures.jl
#


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
	struct MComponents

structure for matrix components for matrix M
"""
struct MComponents
    elementsT::Vector{Element}
    elementsB::Vector{Element}
    nT::Int
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end

struct M2Components
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRG::Vector
    elementsR::Vector
    elementsB::Vector{Element}
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end

"""
	AbstractTComponents

abstract type for components of transition matrices 
"""
abstract type AbstractTComponents end

"""
	struct TComponents

fields:
nT, elementsT

"""
struct TComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector
end
struct RGComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    elementsG::Vector
    elementsR::Vector
    elementsRG::Vector
end


"""
struct TAIComponents{elementType} <: AbstractTComponents

fields:
nT, elementsT, elementsTA, elementsTI

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
struct TDComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTG::Vector{Element}
    elementsTD::Vector{Vector{Element}}
end
"""
 	ModelCoupledComponents

fields:
    nT::Int
    nG::Int
    nR::Int
    sourceState::Int
    targetTransition::Int
    elementsG::Vector
    elementsGt::Vector
    elementsGs::Vector
    elementsRG::Vector
    elementsRGbar::Vector
"""
struct ModelCoupledComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    sourceState::Int
    targetTransition::Int
    elementsG::Vector
    elementsGt::Vector
    elementsGs::Vector
    elementsRG::Vector
    elementsRGbar::Vector
end
"""
 	TCoupledComponents

fields:
    N::Int: total number of states
    model::Tuple: model index for each trace
    sources::Tuple: source model for each
    modelcomponents::Vector{ModelCoupledComponents}
"""
struct TCoupledComponents
    N::Int
    model::Tuple
    sources::Tuple
    modelcomponents::Vector{ModelCoupledComponents}
end


"""
 	MTComponents

 structure for M and T components

fields:
mcomponents::MComponents
tcomponents::TComponents
"""
struct MTComponents
    mcomponents::MComponents
    tcomponents::TComponents
end

struct MRGComponents
    mcomponents::M2Components
    tcomponents::RGComponents
end


"""
 	MTAIComponents

 structure for M and TAI components

fields:
mcomponents::MComponents
tcomponents::TAIComponents
"""
struct MTAIComponents{elementType}
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
struct MTDComponents
    mcomponents::MComponents
    tcomponents::TDComponents
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


struct Indices_coupled
    gamma::Vector{Int}
    nu::Vector{Int}
    eta::Vector{Int}
    decay::Int
end
