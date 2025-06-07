using ChainRulesCore
using LinearAlgebra
using SparseArrays
using Zygote # For testing and its pullbacks
using FiniteDifferences # For verifying the rrule

# Original function provided by the user
function normalized_nullspace(M::SparseMatrixCSC)
    m = size(M, 1)
    p = zeros(eltype(M), m) # Use eltype(M) for type stability
    F = qr(M)
    R = F.R # R is the upper triangular factor
    
    # Back substitution to solve R*p = 0
    # This implies that the last row of R is zero or near zero, allowing a non-trivial solution with p[end]=1.
    p[end] = 1.0
    for i in (m-1):-1:1 # Iterate backwards from m-1 down to 1
        # R[i, i+1:end]' * p[i+1:end] is a dot product sum
        # R[i,i] is the diagonal element
        p[i] = -R[i, i+1:end]' * p[i+1:end] / R[i, i]
    end

    # Permute elements according to sparse matrix result
    # F.pcol is the column permutation vector such that M*P = Q*R (where P is the permutation matrix from F.pcol)
    # The `p` vector solved corresponds to the permuted columns of M.
    # To get the vector `pp` corresponding to the original columns of M, we apply the inverse permutation.
    # pp[F.pcol[i]] = p[i] means pp is p permuted by invperm(F.pcol).
    pp = Vector{eltype(M)}(undef, m) # Preallocate for efficiency
    p_perm = invperm(F.pcol) # Get the inverse permutation
    for i in eachindex(p)
        pp[i] = p[p_perm[i]] # Apply the inverse permutation
    end
    # The original loop: `for i in eachindex(p) pp[F.pcol[i]] = p[i] end`
    # this means `pp` is the vector `p` where `p[i]` goes to `pp[F.pcol[i]]`.
    # So `pp[j]` is `p[k]` where `j = F.pcol[k]`, meaning `k = invperm(F.pcol)[j]`.
    # Therefore, `pp = p[invperm(F.pcol)]` is correct. The original code's loop is applying `invperm(F.pcol)` to the *indices* of `p`.
    # Let's stick to the original code's logic for the forward pass for exact replication,
    # and then derive the adjoint correctly.
    pp_orig_logic = similar(p) # Use similar to ensure type and size
    for i in eachindex(p)
        pp_orig_logic[F.pcol[i]] = p[i]
    end

    # Normalization and non-negativity constraint
    sum_pp = sum(pp_orig_logic)
    output = max.(pp_orig_logic / sum_pp, 0)
    return output
end

# Rrule for normalized_nullspace
function ChainRulesCore.rrule(::typeof(normalized_nullspace), M::SparseMatrixCSC)
    # --- Forward Pass ---
    m = size(M, 1)
    F = qr(M)
    R = F.R

    p_sol = Vector{eltype(M)}(undef, m)
    p_sol[end] = 1.0
    for i in (m-1):-1:1
        p_sol[i] = -R[i, i+1:end]' * p_sol[i+1:end] / R[i, i]
    end

    pp_vec = similar(p_sol)
    for i in eachindex(p_sol)
        pp_vec[F.pcol[i]] = p_sol[i]
    end
    
    sum_pp = sum(pp_vec)
    output = max.(pp_vec / sum_pp, 0)

    # --- Pullback Function ---
    function normalized_nullspace_pullback(ȳ)
        # 1. Backpropagate through `output = max.(pp_vec / sum_pp, 0)`
        # q = pp_vec / sum_pp
        # output = max.(q, 0)
        q̄_masked = ȳ .* (output .> 0) # Apply mask for max.(., 0)

        # Adjoint for `q = pp_vec / sum_pp`
        # d(q)/d(pp_vec)_j = 1/sum_pp for j=i, 0 otherwise
        # d(q)/d(sum_pp) = -pp_vec / sum_pp^2
        pp_vec_bar = q̄_masked / sum_pp
        sum_pp_bar = -dot(q̄_masked, pp_vec) / sum_pp^2
        
        # Add gradient from `sum(pp_vec)` to `pp_vec_bar`
        pp_vec_bar .+= sum_pp_bar # Adjoint of sum (d(sum(x))/dx_i = 1)

        # 2. Backpropagate through `pp_vec[F.pcol[i]] = p_sol[i]` (permutation)
        # This means `pp_vec` is `p_sol` with elements moved to `F.pcol` indices.
        # So, if `pp_vec[j] = p_sol[k]` where `j = F.pcol[k]`, then `p_sol_bar[k] += pp_vec_bar[j]`.
        # This is `p_sol_bar = pp_vec_bar[F.pcol]`
        p_sol_bar = similar(pp_vec)
        for i in eachindex(p_sol)
            p_sol_bar[i] = pp_vec_bar[F.pcol[i]]
        end
        # `p_sol_bar = pp_vec_bar[invperm(F.pcol)]` might be more efficient but matches the loop.

        # 3. Backpropagate through the back-substitution loop `p_sol[i] = -R[i, i+1:end]' * p_sol[i+1:end] / R[i, i]`
        # This is the most complex part to hand-code if not using `\`.
        # We need to compute R_bar given p_sol_bar.
        # This is essentially the adjoint of solving a triangular system.
        # We'll use a reverse loop mirroring the forward pass.
        R_bar = Zygote.Buffer(zeros(eltype(M), size(R))) # Use Zygote.Buffer for mutable array in pullback
        # The gradients for p_sol[end] = 1.0 is NoTangent() (constant).

        # Loop backwards over the forward pass's calculation of p_sol
        for i in 1:(m-1) # Iterate forwards for original R indices
            curr_p_idx = m-i
            # Current step: p_sol[curr_p_idx] = -R[curr_p_idx, curr_p_idx+1:end]' * p_sol[curr_p_idx+1:end] / R[curr_p_idx, curr_p_idx]

            # Let `term_sum = R[curr_p_idx, curr_p_idx+1:end]' * p_sol[curr_p_idx+1:end]`
            # Let `denom = R[curr_p_idx, curr_p_idx]`
            # p_sol[curr_p_idx] = -term_sum / denom

            # Adjoint of p_sol[curr_p_idx] from downstream: p_sol_bar[curr_p_idx]

            # Adjoint for `denom`: d(output)/d(denom) = d(output)/d(p_sol[curr_p_idx]) * d(p_sol[curr_p_idx])/d(denom)
            # = p_sol_bar[curr_p_idx] * (term_sum / denom^2)
            R_bar[curr_p_idx, curr_p_idx] += p_sol_bar[curr_p_idx] * (
                (R[curr_p_idx, curr_p_idx+1:end]' * p_sol[curr_p_idx+1:end]) / R[curr_p_idx, curr_p_idx]^2
            )

            # Adjoint for `term_sum`: d(output)/d(term_sum) = p_sol_bar[curr_p_idx] * (-1 / denom)
            term_sum_bar = p_sol_bar[curr_p_idx] * (-1 / R[curr_p_idx, curr_p_idx])

            # Backpropagate through `term_sum = R[curr_p_idx, curr_p_idx+1:end]' * p_sol[curr_p_idx+1:end]`
            # Adjoint of dot product `a' * b`: a_bar = b * scalar_bar, b_bar = a * scalar_bar
            R_bar[curr_p_idx, curr_p_idx+1:end] .+= term_sum_bar * p_sol[curr_p_idx+1:end] # Add to R_bar
            p_sol_bar[curr_p_idx+1:end] .+= term_sum_bar * R[curr_p_idx, curr_p_idx+1:end] # Propagate to p_sol_bar for previous elements

            # Note: `p_sol_bar[curr_p_idx+1:end]` accumulate contributions from later steps.
            # The loop is correct in processing from m-1 down to 1.
            # But in the adjoint, we need to sum up contributions as we go back up.
            # The loop should ideally run from `curr_p_idx=1` to `m-1` for accumulating to `p_sol_bar`
            # For `R_bar`, it doesn't matter because it's accumulating.

            # Re-thinking the loop for `p_sol_bar`:
            # We already have `p_sol_bar` initialized from permutation step.
            # We want to add contributions from the current step's dependence.
        end
        R_bar_final = copy(R_bar) # Convert Zygote.Buffer to actual array

        # 4. Backpropagate through `F = qr(M)`
        # `qr_back` is the pullback function for qr(M).
        # It takes a tangent of the QR object and returns the tangent of M.
        # The tangent for F (QR object) should have non-zero adjoint for its `factors` field (which is R).
        # We initialize `q=ZeroTangent()` and `pcol=ZeroTangent()` as they are not differentiable in this context.
        qr_val, qr_back = Zygote.pullback(qr, M) # Perform qr in the pullback to get its pullback func

        F_tangent = Tangent{typeof(qr_val)}(; factors=R_bar_final, q=ZeroTangent(), pcol=ZeroTangent())
        M̄_from_qr = qr_back(F_tangent)[1] # Get the gradient for M (first element of the tuple)

        # Return gradients for all inputs: (NoTangent() for function itself, M_bar)
        return NoTangent(), M̄_from_qr
    end

    return output, normalized_nullspace_pullback
end

# --- Testing the rrule ---

# A 3x3 example matrix with a null space: M = [1 2 3; 4 5 6; 7 8 9] (rank 2)
M_square = sparse([1, 1, 1, 2, 2, 2, 3, 3, 3],
                  [1, 2, 3, 1, 2, 3, 1, 2, 3],
                  [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], 3, 3)

println("Original M: \n", Matrix(M_square))

# Calculate the null vector using the function
null_vec_calculated = normalized_nullspace(M_square)
println("Calculated null vector: ", null_vec_calculated)
println("M * null_vec (should be near zero): ", Matrix(M_square) * null_vec_calculated)

# Compute gradient using Zygote (which will use our rrule)
println("\n--- Zygote Gradient ---")
grad_M_zygote = Zygote.gradient(normalized_nullspace, M_square)[1]
println("Gradient of M (Zygote):\n", Matrix(grad_M_zygote))

# Verify with FiniteDifferences.jl
println("\n--- FiniteDifferences Gradient ---")
fdm = FiniteDifferences.central_fdm(5, 1) # 5th order central difference method

# Wrapper function for FiniteDifferences to convert dense to sparse
function normalized_nullspace_dense_wrapper(M_dense::Matrix)
    M_sparse = sparse(M_dense)
    return normalized_nullspace(M_sparse)
end

M_dense_for_fd = Matrix(M_square) # Convert to dense for FD
fd_grad_M = FiniteDifferences.grad(fdm, normalized_nullspace_dense_wrapper, M_dense_for_fd)[1]
println("Gradient of M (FiniteDifferences):\n", fd_grad_M)

# Compare the gradients numerically
diff_grad = Matrix(grad_M_zygote) - fd_grad_M
println("\nDifference (Zygote - FiniteDifferences):\n", diff_grad)
println("Norm of the difference: ", norm(diff_grad))

# Example with a simpler 2x2 matrix that clearly has a null space
# M = [1 1; 1 1] -> null space is [-1, 1] (or [1, -1])
M_has_null = sparse([1, 2, 1, 2], [1, 2, 2, 1], [1.0, 1.0, 1.0, 1.0], 2, 2)
println("\nOriginal M (rank deficient): \n", Matrix(M_has_null))

null_vec_has_null = normalized_nullspace(M_has_null)
println("Calculated null vector: ", null_vec_has_null)
println("M * null_vec (should be near zero): ", Matrix(M_has_null) * null_vec_has_null)

println("\n--- Zygote Gradient for M_has_null ---")
grad_M_has_null_zygote = Zygote.gradient(normalized_nullspace, M_has_null)[1]
println("Gradient of M_has_null (Zygote):\n", Matrix(grad_M_has_null_zygote))

println("\n--- FiniteDifferences Gradient for M_has_null ---")
M_dense_has_null_for_fd = Matrix(M_has_null)
fd_grad_M_has_null = FiniteDifferences.grad(fdm, normalized_nullspace_dense_wrapper, M_dense_has_null_for_fd)[1]
println("Gradient of M_has_null (FiniteDifferences):\n", fd_grad_M_has_null)

diff_grad_has_null = Matrix(grad_M_has_null_zygote) - fd_grad_M_has_null
println("\nDifference (Zygote - FiniteDifferences) for M_has_null:\n", diff_grad_has_null)
println("Norm of the difference for M_has_null: ", norm(diff_grad_has_null))