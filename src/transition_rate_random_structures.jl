# This file is part of StochasticGene.jl
# transition_rate_random_structures.jl

abstract type AbstractRandomComponents <: AbstractComponents end

"""
    struct RandomComponents

Fields:
    nT::Int
    elementsT::Vector
    elementsB::Vector
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
    mean_eject::Float64
    var_eject::Float64
"""
struct RandomComponents <: AbstractRandomComponents
    nT::Int
    elementsT::Vector
    elementsB::Vector
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
    mean_eject::Float64
    var_eject::Float64
end

"""
    RandomComponents(transitions::Tuple, G, R, S, insertstep, splicetype, mean_eject, var_eject)

Constructor for RandomComponents
"""
function RandomComponents(transitions::Tuple, G, R, S, insertstep, splicetype, mean_eject, var_eject)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, elementsB, nT = set_elements_random_TGRS(transitions, G, R, S, insertstep, indices, splicetype, mean_eject, var_eject)
    total = 200  # Maximum mRNA number
    U, Uminus, Uplus = make_mat_U_random(total, 1.0, mean_eject, var_eject)
    RandomComponents(nT, elementsT, elementsB, U, Uminus, Uplus, mean_eject, var_eject)
end

"""
    set_elements_random_TGRS(transitions, G, R, S, insertstep, indices, splicetype, mean_eject, var_eject)

Set elements for random TGRS matrix
"""
function set_elements_random_TGRS(transitions, G, R, S, insertstep, indices, splicetype, mean_eject, var_eject)
    elementsT = Vector{Element}(undef, 0)
    elementsB = Vector{Element}(undef, 0)
    
    # Set G transitions
    set_elements_G!(elementsT, transitions, G, indices.gamma)
    
    # Set R transitions
    set_elements_R!(elementsT, R, S, insertstep, indices.nu, indices.eta, splicetype)
    
    # Set boundary elements with random ejection
    set_elements_B_random!(elementsB, R, indices.nu[R+1], mean_eject, var_eject)
    
    nT = G * 2^R
    return elementsT, elementsB, nT
end

"""
    set_elements_B_random!(elements, R, nu_eject, mean_eject, var_eject)

Set boundary elements with random ejection pathways
"""
function set_elements_B_random!(elements, R, nu_eject, mean_eject, var_eject)
    max_eject = ceil(Int, mean_eject + 3*sqrt(var_eject))
    eject_probs = zeros(max_eject + 1)
    
    # Calculate Poisson probabilities
    for k = 0:max_eject
        eject_probs[k+1] = exp(-mean_eject) * (mean_eject^k) / factorial(k)
    end
    
    # Normalize probabilities
    eject_probs ./= sum(eject_probs)
    
    # Set boundary elements
    for k = 0:max_eject
        push!(elements, Element(R, k, nu_eject, 1))
    end
end
