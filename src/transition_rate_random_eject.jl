# This file is part of StochasticGene.jl
# transition_rate_random_eject.jl

"""
    make_mat_U_random(total::Int, decay::Float64, mean_eject::Float64, var_eject::Float64)

Create U matrix with random ejection pathways

# Arguments
- `total::Int`: Total number of mRNA states
- `decay::Float64`: mRNA decay rate
- `mean_eject::Float64`: Mean number of mRNAs ejected
- `var_eject::Float64`: Variance of ejected mRNA numbers

# Returns
- `U::SparseMatrixCSC`: Transition matrix
- `Uminus::SparseMatrixCSC`: Ejection matrix
- `Uplus::SparseMatrixCSC`: Birth matrix
"""
function make_mat_U_random(total::Int, decay::Float64, mean_eject::Float64, var_eject::Float64)
    # Create probability distribution for ejection numbers
    # Using a Poisson distribution as a starting point
    max_eject = ceil(Int, mean_eject + 3*sqrt(var_eject))  # 3 standard deviations
    eject_probs = zeros(max_eject + 1)
    
    # Calculate Poisson probabilities
    for k = 0:max_eject
        eject_probs[k+1] = exp(-mean_eject) * (mean_eject^k) / factorial(k)
    end
    
    # Normalize probabilities
    eject_probs ./= sum(eject_probs)
    
    # Create matrices
    U = sparse(total, total, zeros(total*total))
    Uminus = sparse(total, total, zeros(total*total))
    Uplus = sparse(total, total, zeros(total*total))
    
    # Set decay terms
    for m = 2:total-1
        U[m, m] = -decay * (m - 1)
        Uplus[m, m+1] = decay * m
    end
    U[total, total] = -decay * (total - 1)
    Uplus[total, total+1] = decay * total
    
    # Set ejection terms
    for m = 1:total
        # Calculate total ejection rate
        total_rate = 1.0  # This will be preserved
        
        # Set ejection probabilities
        for k = 0:min(m, max_eject)
            if m - k >= 1
                Uminus[m, m-k] = eject_probs[k+1] * total_rate
            end
        end
    end
    
    return U, Uminus, Uplus
end

"""
    make_mat_B_random(nT::Int, total::Int, mean_eject::Float64, var_eject::Float64)

Create B matrix with random ejection pathways

# Arguments
- `nT::Int`: Number of T states
- `total::Int`: Total number of mRNA states
- `mean_eject::Float64`: Mean number of mRNAs ejected
- `var_eject::Float64`: Variance of ejected mRNA numbers

# Returns
- `B::SparseMatrixCSC`: Boundary matrix
"""
function make_mat_B_random(nT::Int, total::Int, mean_eject::Float64, var_eject::Float64)
    max_eject = ceil(Int, mean_eject + 3*sqrt(var_eject))
    eject_probs = zeros(max_eject + 1)
    
    # Calculate Poisson probabilities
    for k = 0:max_eject
        eject_probs[k+1] = exp(-mean_eject) * (mean_eject^k) / factorial(k)
    end
    
    # Normalize probabilities
    eject_probs ./= sum(eject_probs)
    
    B = sparse(nT, total*nT, zeros(nT*total*nT))
    
    # Set boundary conditions
    for i = 1:nT
        # Calculate total ejection rate
        total_rate = 1.0  # This will be preserved
        
        # Set ejection probabilities
        for k = 0:max_eject
            if k < total
                B[i, k*nT + i] = eject_probs[k+1] * total_rate
            end
        end
    end
    
    return B
end
