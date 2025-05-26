# This file is part of StochasticGene.jl
# transition_rate_stacked_make.jl

"""
    make_mat_stacked_R(elements, rates, R, buffer_size)

Create stacked R transition matrix
"""
function make_mat_stacked_R(elements::Vector{StackedElement}, rates::Vector, R::Int, buffer_size::Int)
    n = R + buffer_size
    T = sparse(n, n, zeros(n*n))
    
    for e in elements
        T[e.a, e.b] += e.pm * rates[e.index]
    end
    
    return T
end

"""
    make_mat_stacked_T(components, rates)

Create stacked T matrix
"""
function make_mat_stacked_T(components::StackedComponents, rates::Vector)
    T = sparse(components.nT, components.nT, zeros(components.nT*components.nT))
    
    for e in components.elementsT
        T[e.a, e.b] += e.pm * rates[e.index]
    end
    
    return T
end

"""
    make_mat_stacked_TRG(components, rates)

Create stacked TRG matrix
"""
function make_mat_stacked_TRG(components::StackedComponents, rates::Vector)
    T = make_mat_stacked_T(components, rates)
    
    # Extract G, R, and buffer components
    nG = components.nT ÷ (2^components.buffer_size)
    nR = 2^components.buffer_size
    
    # Create submatrices
    GR = sparse(nG, nR, zeros(nG*nR))
    RGbar = sparse(nR, nR, zeros(nR*nR))
    RG = sparse(nR, nR, zeros(nR*nR))
    
    # Fill submatrices
    for e in components.elementsT
        if e.buffer_state > 0  # Buffer state
            GR[e.a, e.b] += e.pm * rates[e.index]
        else  # Regular R state
            RGbar[e.a, e.b] += e.pm * rates[e.index]
            RG[e.a, e.b] += e.pm * rates[e.index]
        end
    end
    
    return T, GR, RGbar, RG
end
