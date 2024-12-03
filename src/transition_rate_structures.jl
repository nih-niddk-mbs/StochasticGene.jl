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

"""
	AbstractTComponents

abstract type for components of transition matrices 
"""
abstract type AbstractTComponents end

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

"""
	struct TRGComponents

fields: 
    nT::Int
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRGbar::Vector
    elementsRG::Vector

"""
struct TRGComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRGbar::Vector
    elementsRG::Vector
end


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
struct TDComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTG::Vector{Element}
    elementsTD::Vector{Vector{Element}}
end

struct TDRGComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRGbar::Vector
    elementsRG::Vector
    elementsGD::Vector
    elementsRGD::Vector
    elementsRGbarD::Vector
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

struct TDRGCoupledUnitComponents <: AbstractTComponents
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
    elementsGD::Vector{Vector{Element}}
    elementsRD::Vector{Vector{Element}}
end
"""
 	TCoupledComponents

fields:
    N::Int: total number of states
    model::Tuple: model index for each trace
    sources::Tuple: source model for each
    modelcomponents::Vector{TRGCoupledComponents}
"""
struct TCoupledComponents{ModelType}
    N::Int
    model::Tuple
    sources::Tuple
    modelcomponents::ModelType
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
struct MComponents
    elementsT::Vector{Element}
    elementsB::Vector{Element}
    nT::Int
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end

"""
    struct MRGComponents

fields:
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRGbar::Vector
    elementsRG::Vector
    elementsB::Vector{Element}
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end
"""
struct MRGComponents
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRGbar::Vector
    elementsRG::Vector
    elementsB::Vector{Element}
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
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

struct MTRGComponents
    mcomponents::MRGComponents
    tcomponents::TRGComponents
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


