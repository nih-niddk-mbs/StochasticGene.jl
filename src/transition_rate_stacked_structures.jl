# This file is part of StochasticGene.jl
# transition_rate_stacked_structures.jl

abstract type AbstractStackedComponents <: AbstractComponents end

"""
    struct StackedComponents

Fields:
    nT::Int
    elementsT::Vector
    buffer_size::Int
"""
struct StackedComponents <: AbstractStackedComponents
    nT::Int
    elementsT::Vector
    buffer_size::Int
end

"""
    StackedComponents(transitions::Tuple, G, R, S, insertstep, buffer_size, splicetype)

Constructor for StackedComponents
"""
function StackedComponents(transitions::Tuple, G, R, S, insertstep, buffer_size, splicetype)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_stacked_TGRS(transitions, G, R, S, insertstep, indices, splicetype, buffer_size)
    StackedComponents(nT, elementsT, buffer_size)
end

"""
    set_elements_stacked_TGRS(transitions, G, R, S, insertstep, indices, splicetype, buffer_size)

Set elements for stacked TGRS matrix
"""
function set_elements_stacked_TGRS(transitions, G, R, S, insertstep, indices, splicetype, buffer_size)
    elementsT = Vector{StackedElement}(undef, 0)
    
    # Set G transitions
    set_elements_G!(elementsT, transitions, G, indices.gamma)
    
    # Set stacked R transitions
    set_elements_stacked_R!(elementsT, R, buffer_size, indices.nu)
    
    # Set splice transitions if needed
    if S > 0
        set_elements_splice!(elementsT, R, S, insertstep, indices.eta, splicetype)
    end
    
    nT = G * (2^R + buffer_size)
    return elementsT, nT
end
