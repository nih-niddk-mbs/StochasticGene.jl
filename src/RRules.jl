using ChainRules
using ChainRulesCore
using LinearAlgebra
using SparseArrays
using Zygote # For testing and its pullbacks
using FiniteDifferences # For verifying the rrule

# Custom adjoint for qr on dense matrices
# function ChainRulesCore.rrule(::typeof(qr), M::Matrix)
#     F = qr(M)
    
#     function qr_pullback(Δ)
#         # Extract the gradient for R
#         R̄ = Δ.factors
        
#         # Compute gradient for input matrix
#         M̄ = zeros(size(M))
        
#         # Backpropagate through QR decomposition
#         # For each column j
#         for j in 1:size(M, 2)
#             # For each row i >= j (upper triangular part)
#             for i in j:size(M, 1)
#                 # Gradient contribution from R
#                 M̄[i, j] = R̄[i, j]
#             end
#         end
        
#         return NoTangent(), M̄
#     end
    
#     return F, qr_pullback
# end

# function normalized_nullspace_original(M::SparseMatrixCSC)
#     m = size(M, 1)
#     p = zeros(m)
#     F = qr(M)   #QR decomposition
#     R = F.R
#     # Back substitution to solve R*p = 0
#     p[end] = 1.0
#     for i in 1:m-1
#         p[m-i] = -R[m-i, m-i+1:end]' * p[m-i+1:end] / R[m-i, m-i]
#     end
#     # Permute elements according to sparse matrix result
#     pp = copy(p)
#     for i in eachindex(p)
#         pp[F.pcol[i]] = p[i]
#     end
#     max.(pp / sum(pp), 0)

# end

# # Original function provided by the user
# function normalized_nullspace(M::SparseMatrixCSC)
#     m = size(M, 1)
#     p = zeros(eltype(M), m) # Use eltype(M) for type stability
#     F = qr(M)
#     R = F.R # R is the upper triangular factor
    
#     # Back substitution to solve R*p = 0
#     # This implies that the last row of R is zero or near zero, allowing a non-trivial solution with p[end]=1.
#     p[end] = 1.0
#     for i in (m-1):-1:1 # Iterate backwards from m-1 down to 1
#         # R[i, i+1:end]' * p[i+1:end] is a dot product sum
#         # R[i,i] is the diagonal element
#         p[i] = -R[i, i+1:end]' * p[i+1:end] / R[i, i]
#     end

#     # Permute elements according to sparse matrix result
#     pp = Vector{eltype(M)}(undef, m) # Preallocate for efficiency
#     for i in eachindex(p)
#         pp[F.pcol[i]] = p[i]
#     end

#     # Normalization and non-negativity constraint
#     sum_pp = sum(pp)
#     output = max.(pp / sum_pp, 0)
#     return output
# end

# # Rrule for normalized_nullspace
# function ChainRulesCore.rrule(::typeof(normalized_nullspace), M::SparseMatrixCSC)
#     # --- Forward Pass ---
#     m = size(M, 1)
#     F = qr(M)
#     R = F.R
#     Q = F.Q

#     # Back substitution to solve R*p = 0
#     p = Vector{eltype(M)}(undef, m)
#     p[end] = 1.0
#     for i in (m-1):-1:1
#         p[i] = -R[i, i+1:end]' * p[i+1:end] / R[i, i]
#     end

#     # Permute elements
#     pp = Vector{eltype(M)}(undef, m)
#     for i in eachindex(p)
#         pp[F.pcol[i]] = p[i]
#     end

#     # Normalization and non-negativity
#     sum_pp = sum(pp)
#     output = max.(pp / sum_pp, 0)

#     # --- Pullback Function ---
#     function normalized_nullspace_pullback(ȳ)
#         # 1. Backpropagate through max and normalization
#         # ∂output/∂pp = (pp / sum_pp > 0) / sum_pp
#         # ∂output/∂sum_pp = -(pp / sum_pp > 0) .* pp / sum_pp^2
#         mask = output .> 0
#         pp̄ = ȳ .* mask / sum_pp
#         sum_pp̄ = -dot(ȳ .* mask, pp) / sum_pp^2
#         pp̄ .+= sum_pp̄

#         # 2. Backpropagate through permutation
#         # ∂pp/∂p = permutation matrix P
#         p̄ = similar(pp)
#         for i in eachindex(p)
#             p̄[i] = pp̄[F.pcol[i]]
#         end

#         # 3. Backpropagate through back substitution
#         # ∂p/∂R[i,i] = p[i] / R[i,i]
#         # ∂p/∂R[i,i+1:end] = -p[i+1:end]' / R[i,i]
#         R̄ = zeros(eltype(M), size(R))
#         for i in 1:(m-1)
#             curr_p_idx = m-i
#             # Diagonal element
#             R̄[curr_p_idx, curr_p_idx] = p̄[curr_p_idx] * p[curr_p_idx] / R[curr_p_idx, curr_p_idx]
#             # Off-diagonal elements
#             R̄[curr_p_idx, curr_p_idx+1:end] = -p̄[curr_p_idx] * p[curr_p_idx+1:end]' / R[curr_p_idx, curr_p_idx]
#         end

#         # 4. Backpropagate through QR
#         # ∂M/∂R = Q
#         M̄ = Q * R̄
#         # Apply sparsity pattern
#         M̄ = M̄ .* (M .!= 0)

#         return NoTangent(), M̄
#     end

#     return output, normalized_nullspace_pullback
# end

################
# Gemini 2.5 Pro
################

# function make_Q(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], fittedparam=[], noisepriors=[0.0, 0.1, 1.0, 0.1])
#     # tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials)
#     # trace, tracescale = zero_median(tracer, zeromedian)
#     # nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
#     data = StochasticGene.TraceData{String,String,Tuple}("trace", "gene", 1, ())
#     model = StochasticGene.load_model(data, rtarget, rtarget, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, rtarget, Int[], 1., .1, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, true)
#     StochasticGene.make_mat_T(model.components, model.rates)
# end

# import ChainRulesCore: rrule, NoTangent, dot, @ignore_derivatives
# using LinearAlgebra, SparseArrays

# # This helper function performs the forward pass and returns all
# # intermediate values needed for the pullback. It is only called from the rrule.
# function _normalized_nullspace_with_intermediates(M::SparseMatrixCSC)
#     m = size(M, 1)

#     # Perform the sparse QR factorization.
#     # We use @ignore_derivatives to signal that the permutation (pcol)
#     # is non-differentiable and should be treated as a constant by AD.
#     F = @ignore_derivatives qr(M)
#     Q, R, pcol = F.Q, F.R, F.pcol

#     # Solve for the nullspace vector `p` of the permuted system (R*p=0).
#     p = Vector{eltype(M)}(undef, m)
#     p[end] = 1.0

#     R_upper = R[1:m-1, 1:m-1]
#     r_last_col = R[1:m-1, m]

#     # Solve the upper triangular system R_upper * p_upper = -r_last_col
#     # We must convert the view of the sparse R to a dense matrix for `det` to work.
#     if abs(det(Matrix(R_upper))) > 1e-12
#         p[1:m-1] = UpperTriangular(R_upper) \ -r_last_col
#     else
#         p[1:m-1] .= 0.0
#     end

#     # Un-permute the nullspace vector to match the original matrix's column order.
#     pp = Vector{eltype(M)}(undef, m)
#     pp[pcol] = p

#     # Normalize the final vector (L1-style normalization)
#     sum_pp = sum(pp)
#     output = abs(sum_pp) > 1e-12 ? max.(pp / sum_pp, 0.0) : zeros(eltype(M), m)

#     return output, Q, R, p, pp, pcol, sum_pp
# end


# # The user-facing function is a simple wrapper around the internal implementation.
# function normalized_nullspace(M::SparseMatrixCSC)
#     output, _, _, _, _, _, _ = _normalized_nullspace_with_intermediates(M)
#     return output
# end


# # This single, comprehensive rrule for the top-level function should now be
# # correctly dispatched by Zygote.
# function ChainRulesCore.rrule(::typeof(normalized_nullspace), M::SparseMatrixCSC)

#     # --- Forward Pass ---
#     # Run the helper to get the primal output and all intermediate values.
#     (output, Q, R, p_of_R, pp, pcol, sum_pp) = _normalized_nullspace_with_intermediates(M)
#     m = size(M, 1)

#     # --- Pullback Function ---
#     function normalized_nullspace_pullback(ȳ)
#         # 1. Backpropagate through normalization (max.(pp/sum_pp, 0)) and non-negativity.
#         pp̄ = if abs(sum_pp) > 1e-12
#             # Apply mask for the max(..., 0) part.
#             masked_grad = ȳ .* (output .> 0)

#             # *** FINAL CORRECTED GRADIENT FOR NORMALIZATION ***
#             # For y = x/sum(x), the vector-Jacobian product is x̄ = (ȳ - dot(ȳ, y)) / sum(x).
#             # The result of the dot product is a scalar, which Julia correctly broadcasts.
#             (masked_grad .- dot(masked_grad, output)) ./ sum_pp
#         else
#             zeros(eltype(M), m)
#         end

#         # 2. Backpropagate through permutation to get the gradient for `p_of_R`.
#         p_of_R_bar = pp̄[pcol]

#         # 3. Project the gradient to be orthogonal to p_of_R.
#         # This is necessary because the choice of scale for p_of_R (p_of_R[end]=1) was arbitrary.
#         p_dot_p = dot(p_of_R, p_of_R)
#         if p_dot_p > 1e-12
#             p_of_R_bar .-= (dot(p_of_R_bar, p_of_R) / p_dot_p) .* p_of_R
#         end

#         # 4. Backpropagate through back-substitution (R*p=0) to get R̄.
#         p̄_upper = p_of_R_bar[1:m-1]
#         R_upper = R[1:m-1, 1:m-1]

#         # We solve R_upper' * λ = p̄_upper for λ.
#         # Convert R_upper' to a dense LowerTriangular matrix for the solver.
#         λ = LowerTriangular(Matrix(R_upper')) \ p̄_upper
#         p_of_R_upper = p_of_R[1:m-1]

#         # Construct the dense gradient R̄ from λ and p_of_R.
#         R̄ = zeros(eltype(M), m, m)
#         R̄[1:m-1, 1:m-1] = -λ * p_of_R_upper'
#         R̄[1:m-1, m] = -λ * p_of_R[end]

#         # 5. Backpropagate through QR decomposition (M_perm ≈ Q*R).
#         # The result M̄_perm will be a dense matrix.
#         M̄_perm = Q * R̄

#         # 6. Un-permute the gradient columns to get the final dense gradient M̄.
#         M̄ = zeros(eltype(M), m, m)
#         M̄[:, pcol] = M̄_perm

#         return NoTangent(), M̄
#     end

#     return output, normalized_nullspace_pullback
# end




# --- Testing the rrule ---

# A 3x3 example matrix with a null space: M = [1 2 3; 4 5 6; 7 8 9] (rank 2)
# M_square = sparse([1, 1, 1, 2, 2, 2, 3, 3, 3],
#                   [1, 2, 3, 1, 2, 3, 1, 2, 3],
#                   [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], 3, 3)

# println("Original M: \n", Matrix(M_square))

# # Calculate the null vector using the function
# null_vec_calculated = normalized_nullspace(M_square)
# println("Calculated null vector: ", null_vec_calculated)
# println("M * null_vec (should be near zero): ", Matrix(M_square) * null_vec_calculated)

# # Compute Jacobian using Zygote (which will use our rrule)
# println("\n--- Zygote Jacobian ---")
# jac_M_zygote = Zygote.jacobian(normalized_nullspace, M_square)[1]
# println("Jacobian of M (Zygote):\n", Matrix(jac_M_zygote))

# # Verify with FiniteDifferences.jl
# println("\n--- FiniteDifferences Jacobian ---")
# fdm = FiniteDifferences.central_fdm(5, 1) # 5th order central difference method

# # Wrapper function for FiniteDifferences to convert dense to sparse
# function normalized_nullspace_dense_wrapper(M_dense::Matrix)
#     M_sparse = sparse(M_dense)
#     return normalized_nullspace(M_sparse)
# end

# M_dense_for_fd = Matrix(M_square) # Convert to dense for FD
# fd_jac_M = FiniteDifferences.jacobian(fdm, normalized_nullspace_dense_wrapper, M_dense_for_fd)[1]
# println("Jacobian of M (FiniteDifferences):\n", fd_jac_M)

# # Compare the jacobians numerically
# diff_jac = Matrix(jac_M_zygote) - fd_jac_M
# println("\nDifference (Zygote - FiniteDifferences):\n", diff_jac)
# println("Norm of the difference: ", norm(diff_jac))

# # Example with a simpler 2x2 matrix that clearly has a null space
# # M = [1 1; 1 1] -> null space is [-1, 1] (or [1, -1])
# M_has_null = sparse([1, 2, 1, 2], [1, 2, 2, 1], [1.0, 1.0, 1.0, 1.0], 2, 2)
# println("\nOriginal M (rank deficient): \n", Matrix(M_has_null))

# null_vec_has_null = normalized_nullspace(M_has_null)
# println("Calculated null vector: ", null_vec_has_null)
# println("M * null_vec (should be near zero): ", Matrix(M_has_null) * null_vec_has_null)

# println("\n--- Zygote Jacobian for M_has_null ---")
# jac_M_has_null_zygote = Zygote.jacobian(normalized_nullspace, M_has_null)[1]
# println("Jacobian of M_has_null (Zygote):\n", Matrix(jac_M_has_null_zygote))

# println("\n--- FiniteDifferences Jacobian for M_has_null ---")
# M_dense_has_null_for_fd = Matrix(M_has_null)
# fd_jac_M_has_null = FiniteDifferences.jacobian(fdm, normalized_nullspace_dense_wrapper, M_dense_has_null_for_fd)[1]
# println("Jacobian of M_has_null (FiniteDifferences):\n", fd_jac_M_has_null)

# diff_jac_has_null = Matrix(jac_M_has_null_zygote) - fd_jac_M_has_null
# println("\nDifference (Zygote - FiniteDifferences) for M_has_null:\n", diff_jac_has_null)
# println("Norm of the difference for M_has_null: ", norm(diff_jac_has_null))

# Test function for gradient of v^2 where v is nullspace
# function test_nullspace_squared_gradient(M)
#     # Function to compute v^2
#     function nullspace_squared(M)
#         v = normalized_nullspace(M)
#         return sum(v.^2)  # sum of squares
#     end
    
#     # Compute gradient using Zygote
#     grad_zygote = Zygote.gradient(nullspace_squared, M)[1]
    
#     # Compute gradient using FiniteDifferences
#     fdm = FiniteDifferences.central_fdm(5, 1)
#     function nullspace_squared_dense(M_dense)
#         M_sparse = sparse(M_dense)
#         return nullspace_squared(M_sparse)
#     end
#     grad_fd = FiniteDifferences.grad(fdm, nullspace_squared_dense, Matrix(M))[1]
    
#     # Compare results
#     println("Original matrix:\n", Matrix(M))
#     println("\nNullspace vector:\n", normalized_nullspace(M))
#     println("\nGradient (Zygote):\n", Matrix(grad_zygote))
#     println("\nGradient (FiniteDifferences):\n", grad_fd)
#     println("\nDifference norm: ", norm(Matrix(grad_zygote) - grad_fd))
    
#     return grad_zygote, grad_fd
# end



###################

# --- Step 1: Define the Custom Adjoint (rrule) ---
# This code defines a single, comprehensive differentiation rule for `normalized_nullspace`.
# It explicitly handles all steps of the gradient calculation, minimizing reliance on
# Zygote's internal behaviors.

import ChainRulesCore: rrule, NoTangent, dot
using LinearAlgebra, SparseArrays, Zygote

# This internal helper function performs the forward pass and returns all
# intermediate values needed for the pullback. It's good practice to keep the
# main function clean and separate the logic needed only for differentiation.
function _normalized_nullspace_with_intermediates(M::SparseMatrixCSC)
    m = size(M, 1)

    # Perform the sparse QR factorization to get the permutation and factors.
    # The permutation itself is treated as non-differentiable.
    F = qr(M)
    Q, R, pcol = F.Q, F.R, F.pcol

    # Solve for the nullspace vector `p` of the permuted system (R*p=0).
    p = Vector{eltype(M)}(undef, m)
    p[end] = 1.0

    R_upper = R[1:m-1, 1:m-1]
    r_last_col = R[1:m-1, m]

    # We must convert the view of the sparse R to a dense matrix for `det` to work.
    if abs(det(Matrix(R_upper))) > 1e-12
        p[1:m-1] = UpperTriangular(R_upper) \ -r_last_col
    else
        p[1:m-1] .= 0.0
    end

    # Un-permute the nullspace vector to match the original matrix's column order.
    pp = Vector{eltype(M)}(undef, m)
    pp[pcol] = p

    # Normalize the final vector (L1-style normalization).
    sum_pp = sum(pp)
    output = abs(sum_pp) > 1e-12 ? max.(pp / sum_pp, 0.0) : zeros(eltype(M), m)

    return output, Q, R, p, pp, pcol, sum_pp
end

# This is the function you will call in your code. Its rrule is defined below.
function normalized_nullspace(M::SparseMatrixCSC)
    # The forward pass just needs the final output.
    output, _, _, _, _, _, _ = _normalized_nullspace_with_intermediates(M)
    return output
end

# This single rrule for `normalized_nullspace` contains all the hand-coded logic.
function ChainRulesCore.rrule(::typeof(normalized_nullspace), M::SparseMatrixCSC)

    # --- Forward Pass ---
    # Run the helper to get the primal output and all intermediate values for the pullback.
    (output, Q, R, p_of_R, pp, pcol, sum_pp) = _normalized_nullspace_with_intermediates(M)
    m = size(M, 1)

    # --- Pullback Function ---
    # This function defines how to go from the gradient of the output (ȳ)
    # back to the gradient of the input matrix (M̄).
    function normalized_nullspace_pullback(ȳ)
        # 1. Backpropagate through normalization and non-negativity.
        pp̄ = if abs(sum_pp) > 1e-12
            masked_grad = ȳ .* (output .> 0)
            # Correct adjoint for y = x/sum(x) is (ȳ - dot(ȳ, y)) / sum(x).
            (masked_grad .- dot(masked_grad, output)) ./ sum_pp
        else
            zeros(eltype(M), m)
        end

        # 2. Backpropagate through the un-permutation step.
        # Forward: pp[pcol] = p_of_R. Adjoint: p_of_R_bar = pp̄[pcol].
        p_of_R_bar = pp̄[pcol]

        # 3. Project the gradient to be orthogonal to p_of_R.
        # This correctly handles the arbitrary scaling choice (p[end]=1) made earlier.
        p_dot_p = dot(p_of_R, p_of_R)
        if p_dot_p > 1e-12
            p_of_R_bar .-= (dot(p_of_R_bar, p_of_R) / p_dot_p) .* p_of_R
        end

        # 4. Backpropagate through the back-substitution (R*p=0) to get R̄.
        p̄_upper = p_of_R_bar[1:m-1]
        R_upper = R[1:m-1, 1:m-1]
        λ = LowerTriangular(Matrix(R_upper')) \ p̄_upper
        p_of_R_upper = p_of_R[1:m-1]

        R̄ = zeros(eltype(M), m, m)
        R̄[1:m-1, 1:m-1] = -λ * p_of_R_upper'
        R̄[1:m-1, m] = -λ * p_of_R[end]

        # 5. Backpropagate through the QR decomposition (M_perm ≈ Q*R).
        M̄_perm = Q * R̄

        # 6. Backpropagate through the column permutation to get the final dense gradient M̄.
        M̄ = zeros(eltype(M), m, m)
        M̄[:, pcol] = M̄_perm

        return NoTangent(), M̄
    end

    return output, normalized_nullspace_pullback
end


# --- Step 2: Set up the Gradient Calculation using your Real Module ---
# This part of the script assumes your `StochasticGene.jl` module is loaded
# and that your functions (`load_model`, `make_mat_T`, etc.) are available.

using StochasticGene

# IMPORTANT ASSUMPTION:
# For this to work, your `StochasticGene.make_mat_T` function must internally
# call the `normalized_nullspace` function that we defined above.

# Your `make_Q` function, unchanged. It will now call the real StochasticGene functions.
function make_Q(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], nsamples=10000000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1:num_rates(transitions, R, S, insertstep)-1; num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+1], propcv=0.01, cv=100.0, noisepriors=[0.0, 0.1, 1.0, 0.1], zeromedian=true, maxtime=100.0, initprior=0.1)
    data = TraceData{String,String,Tuple}("trace", "gene", 1, ())
    model = load_model(data, rtarget, rtarget, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, rtarget, Int[], 1., 0.1, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing, zeromedian)
    StochasticGene.make_mat_T(model.components, model.rates)
end


# Define a scalar loss function. AD requires a single number to differentiate.
function model_loss(rtarget_params)
    # This now calls your real model code
    T = make_Q(rtarget=rtarget_params)
    # A simple example loss function
    return sum(T)
end

# Define the initial parameters you want to find the gradient for.
initial_rtarget = [0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50.0, 5.0, 50.0, 5.0]

println("Calculating gradient of loss with respect to rtarget...")

# Use Zygote to get the gradient of the entire process!
# The `[1]` is necessary to extract the gradient tuple from Zygote's output.
grad_rtarget = Zygote.gradient(model_loss, initial_rtarget)[1]

println("\nInitial rtarget values:")
println(initial_rtarget)
println("\nGradient with respect to rtarget:")
println(grad_rtarget)
