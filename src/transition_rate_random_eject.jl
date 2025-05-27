# This file is part of StochasticGene.jl
# transition_rate_random_eject.jl

"""
    make_mat_U(total::Int, decay::Float64, mean_eject::Float64)

Create U matrix with Poisson-distributed ejection

# Arguments
- `total::Int`: Total number of mRNA states
- `decay::Float64`: mRNA decay rate
- `mean_eject::Float64`: Mean number of mRNAs ejected

# Returns
- `U::SparseMatrixCSC`: Transition matrix
- `Uminus::SparseMatrixCSC`: Ejection matrix
- `Uplus::SparseMatrixCSC`: Birth matrix
"""
function make_mat_U(total::Int, decay::Float64, mean_eject::Float64)
    max_eject = ceil(Int, mean_eject + 3*sqrt(mean_eject))  # 3 standard deviations
    
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
    
    # Set ejection terms using Poisson distribution
    for m = 1:total
        # Calculate Poisson probabilities
        for k = 0:min(m, max_eject)
            if m - k >= 1
                prob = exp(-mean_eject) * (mean_eject^k) / factorial(k)
                Uminus[m, m-k] = prob
            end
        end
    end
    
    return U, Uminus, Uplus
end
