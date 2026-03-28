# This file is part of StochasticGene.jl
#
# chemical_master.jl
#
# Functions to solve the chemical master equation
# operates on the transpose of the Markov process transition rate matrix


# Functions to compute steady state mRNA histograms
"""
    steady_state_vector(M::AbstractMatrix; solver::Symbol=:default)

Return a normalized steady-state vector `p` with `M * p ≈ 0` and entries summing to 1.

- `solver=:default` or `:fast` — [`normalized_nullspace`](@ref) (sparse QR, SVD fallback).
- `solver=:augmented` or `:autodiff` — [`normalized_nullspace_augmented`](@ref) (linear solve + AD rules).

Non-sparse `M` is converted with `sparse(M)`.
"""
function steady_state_vector(M::SparseMatrixCSC; solver::Symbol=:default)
    if solver === :default || solver === :fast
        return normalized_nullspace(M)
    elseif solver === :augmented || solver === :autodiff
        return normalized_nullspace_augmented(M)
    else
        throw(
            ArgumentError(
                "steady_state_vector: unknown solver $(repr(solver)); use :default, :fast, :augmented, or :autodiff",
            ),
        )
    end
end

steady_state_vector(M::AbstractMatrix; kwargs...) = steady_state_vector(SparseMatrixCSC(M); kwargs...)

"""
    steady_state(M, nT, nalleles, nhist; steady_state_solver=:default)
    steady_state(M, nT, nalleles; steady_state_solver=:default)

Return steady-state mRNA histogram for G and GR models via [`steady_state_vector`](@ref) and marginalization.

# Arguments
- `M`: Truncated full transition rate matrix including transcription state transitions (GR) and mRNA birth and death
- `nT`: Dimension of state transitions
- `nalleles`: Number of alleles (for convolution)
- `nhist`: Length of histogram to return (four-argument form only)
- `steady_state_solver`: passed to [`steady_state_vector`](@ref) (`:default` or `:augmented` / `:autodiff` for gradients)

# Returns
- `steady_state(M, nT, nalleles, nhist)`: Vector of length `nhist`
- `steady_state(M, nT, nalleles)`: Full marginal histogram convolved over alleles
"""
steady_state(M, nT, nalleles, nhist; kwargs...) = steady_state(M, nT, nalleles; kwargs...)[1:nhist]

function steady_state(M, nT, nalleles; steady_state_solver::Symbol=:default)
    P = steady_state_vector(M; solver=steady_state_solver)
    mhist = marginalize(P, nT)
    allele_convolve(mhist, nalleles)
end



# Functions to compute dwell time distributions

"""
    dwelltimeCDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=Tsit5())

  return dwell time CDF (cumulative distribution function) to reach barrier states starting from Sinit in live states
  live U barrier = all states

  First passage time calculation from live states to barrier states

- `tin`: vector of time bins for dwell time distribution
- `Td`: dwell time rate matrix (within non barrier states)
- `barrier`: set of barrier states
- `Sinit`: initial state
- `method`: 1 for directly solving ODE otherwise use eigendecomposition
"""
function dwelltimeCDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=Tsit5())
    t = [max(2 * tin[1] - tin[2], 0.); tin]
    S = time_evolve(t, Td, Sinit, method)
    return vec(sum(S[:, barrier], dims=1))
end

"""
    dwelltimePDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=Tsit5())

return dwell time PDF (probability density function)
"""
dwelltimePDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=Tsit5()) = pdf_from_cdf(dwelltimeCDF(tin, Td, barrier, Sinit, method))


ontimePDF(tin::Vector, TA::AbstractMatrix, offstates, SAinit::Vector, method=Tsit5()) = dwelltimePDF(tin, TA, offstates, SAinit, method)

offtimePDF(tin::Vector, TI::AbstractMatrix, onstates, SIinit::Vector, method=Tsit5()) = dwelltimePDF(tin, TI, onstates, SIinit, method)

function offonPDF(t::Vector, r::Vector, T::AbstractMatrix, TA::AbstractMatrix, TI::AbstractMatrix, nT::Int, elementsT::Vector, onstates::Vector; steady_state_solver::Symbol=:default)
    pss = steady_state_vector(T; solver=steady_state_solver)
    nonzeros = nonzero_rows(TI)
    offtimePDF(t, TI[nonzeros, nonzeros], nonzero_states(onstates, nonzeros), init_SI(r, onstates, elementsT, pss, nonzeros)), ontimePDF(t, TA, off_states(nT, onstates), init_SA(r, onstates, elementsT, pss))
end

"""
pdf_from_cdf(S)

return PDF (derivative (using finite difference) of CDF)

- `S`: dwell time CDF
"""
function pdf_from_cdf(S)
    P = diff(S)
    P / sum(P)
end


"""
    init_S(r::Vector,livestates::Vector,elements::Vector,pss)

return initial distribution for dwell time distribution

given by transition probability of entering live states in steady state
(probability of barrier state multiplied by transition rate to live state)

- `r`: transition rates
- `sojourn`: set of sojourn states (e.g. onstates for ON Time distribution)
- `elements`: vector of Elements (transition rate matrix element structures)
- `pss`: steady state distribution

"""
function init_S(r::Vector, sojourn::Vector, elements::Vector, pss)
    Sinit = zeros(length(pss))
    for e in elements
        if e.b != e.a && (e.a ∈ sojourn && e.b ∉ sojourn)
            Sinit[e.a] += pss[e.b] * r[e.index]
        end
    end
    s = sum(Sinit)
    (iszero(s) || !isfinite(s)) ? Sinit : Sinit / s
end

function init_S(r::Vector, sojourn::Vector, elements::Vector, pss, nonzeros)
    Sinit = zeros(length(pss))
    for e in elements
        if e.b != e.a && (e.a ∈ sojourn && e.b ∉ sojourn)
            Sinit[e.a] += pss[e.b] * r[e.index]
        end
    end
    Sinit = Sinit[nonzeros]
    # Sinit / sum(Sinit)
end

"""
    init_S(sojourn::Vector, T::SparseMatrixCSC, pss)

Initial distribution for dwell time: **ON** = entrance to sojourn (flux from barrier into sojourn).
**OFF** = stationary distribution of exit states of the previous sojourn, mapped to entrance of
current sojourn (same net formula: flux from complement into sojourn). T(i,j) = rate from j to i.
"""
function init_S(sojourn::Vector, T::SparseMatrixCSC, pss)
    Sinit = zeros(length(pss))
    rows, cols, vals = findnz(T)
    for i in eachindex(rows)
        if cols[i] != rows[i] && (rows[i] ∈ sojourn && cols[i] ∉ sojourn)
            Sinit[rows[i]] += pss[cols[i]] * vals[i]
        end
    end
    s = sum(Sinit)
    (iszero(s) || !isfinite(s)) ? Sinit : Sinit / s
end

"""
    init_S(sojourn::Vector, T::SparseMatrixCSC, pss, dttype::String)

Dwell-time initial distribution by type.
- ON/ONG: entrance to sojourn = flux from barrier into sojourn. Sinit[to] += pss[from] * T(to, from).
- OFF/OFFG: entrance to sojourn (OFF) = exit-state distribution of previous sojourn (ON) mapped by
  transition; same net formula (flux from complement into sojourn). T(i,j) = rate from j to i.
"""
function init_S(sojourn::Vector, T::SparseMatrixCSC, pss, dttype::String)
    Sinit = zeros(length(pss))
    rows, cols, vals = findnz(T)
    for i in eachindex(rows)
        cols[i] == rows[i] && continue
        if rows[i] ∈ sojourn && cols[i] ∉ sojourn
            # Entrance to sojourn: flux from col (outside) to row (inside). T(row,col) = rate col→row.
            Sinit[rows[i]] += pss[cols[i]] * vals[i]
        end
    end
    s = sum(Sinit)
    (iszero(s) || !isfinite(s)) ? Sinit : Sinit / s
end


function init_S(sojourn::Vector, T::SparseMatrixCSC, pss, nonzeros)
    Sinit = init_S(sojourn, T, pss)[nonzeros]
    s = sum(Sinit)
    (iszero(s) || !isfinite(s)) ? Sinit : Sinit / s
end

"""
init_SA(r::Vector,onstates::Vector,elements::Vector,pss::Vector)

return initial distribution for ON time distribution
"""
init_SA(r::Vector, onstates::Vector, elements::Vector, pss::Vector) = init_S(r, onstates, elements, pss)

"""
    init_SI(r::Vector,onstates::Vector,elements::Vector,pss,nonzeros)

return nonzero states of initial distribution for OFF time distribution

"""
function init_SI(r::Vector, onstates::Vector, elements::Vector, pss, nonzeros)
    Sinit = zeros(length(pss))
    for e in elements
        if e.b != e.a && (e.b ∈ onstates && e.a ∉ onstates)
            Sinit[e.a] += pss[e.b] * r[e.index]
        end
    end
    Sinit = Sinit[nonzeros]
    Sinit / sum(Sinit)
end


"""
marginalize(p::Vector,nT,nhist)
marginalize(p::Vector,nT)
marginalize(P::Matrix)

Marginalize over G states
"""
function marginalize(p::AbstractVector, nT, nhist)
    return [sum(@view p[(m-1)*nT+1:m*nT]) for m in 1:nhist]
end

function marginalize(p::AbstractVector, nT)
    nhist = div(length(p), nT)
    marginalize(p, nT, nhist)
end

marginalize(P::Matrix; dims=1) = sum(P, dims=dims)



"""
    unfold(P::Matrix)

reshape matrix into a 1D array
"""
unfold(P::Matrix) = reshape(P, length(P))


"""
time_evolve(t,M::Matrix,Sinit::Vector)

Eigenvalue solution of Linear ODE with rate matrix T and initial vector Sinit
"""
function time_evolve(t, Q::AbstractMatrix, S0::Vector, method)
    if isnothing(method)
        return time_evolve_eig(t, Q, S0)
    else
        return time_evolve_diff(t, Q, S0, method)
    end
end
"""
time_evolve_diff(t,M::Matrix,P0)

Solve master equation problem using DifferentialEquations.jl
"""
function time_evolve_diff(t, Q::SparseMatrixCSC, P0, method=Rosenbrock23())
    tspan = (t[1], t[end])
    prob = ODEProblem(fevolve!, P0, tspan, Q)
    sol = solve(prob, method, saveat=t)
    return sol'
end

"""
    fevolve_inplace(du, u::Vector, p, t)

in place update of du of ODE system for DifferentialEquations,jl
"""
function fevolve!(du, u::Vector, p, t)
    du .= p * u
end

function fevolve(u::Vector, p, t)
    return p * u
end


"""
    time_evolve_eig(t, M::AbstractMatrix, Sinit::Vector)

Solve master equation problem using eigen decomposition
"""
function time_evolve_eig(t, M::AbstractMatrix, Sinit::Vector)
    vals, vects = eig_decompose(M)
    weights = solve_vector(vects, Sinit)
    time_evolve_eig(t, vals, vects, weights)
end
"""
    time_evolve_eig(t::Float64, vals::Vector, vects::Matrix, weights::Vector)



"""
function time_evolve_eig(t::Float64, vals::Vector, vects::Matrix, weights::Vector)
    n = length(vals)
    S = zeros(n)
    for j = 1:n
        for i = 1:n
            S[j] += real(weights[i] * vects[j, i] * exp.(vals[i] * t))
        end
    end
    return S
end
"""
    time_evolve_eig(t::Vector, vals::Vector, vects::Matrix, weights::Vector)

"""
function time_evolve_eig(t::Vector, vals::Vector, vects::Matrix, weights::Vector)
    ntime = length(t)
    n = length(vals)
    S = Array{Float64,2}(undef, ntime, n)
    for j = 1:n
        Sj = zeros(ntime)
        for i = 1:n
            Sj += real(weights[i] * vects[j, i] * exp.(vals[i] * t))
        end
        S[:, j] = Sj
    end
    return S
end

"""
normalized_nullspace_qr(M::SparseMatrixCSC)

Fast QR-based nullspace used historically in this package.
Assumes an irreducible rate matrix of rank n-1; may be fragile
when the matrix is close to singular beyond rank n-1.
"""
function normalized_nullspace_qr(M::SparseMatrixCSC)
    m = size(M, 1)
    p = zeros(m)
    F = qr(M)   # QR decomposition
    R = Matrix(F.R)  # convert to dense once to avoid repeated sparse slicing in back substitution
    # Back substitution to solve R*p = 0
    p[end] = 1.0
    for i in 1:m-1
        p[m-i] = -dot(@view(R[m-i, m-i+1:end]), @view(p[m-i+1:end])) / R[m-i, m-i]
    end
    # Permute elements according to sparse matrix result
    pp = copy(p)
    for i in eachindex(p)
        pp[F.pcol[i]] = p[i]
    end
    max.(pp / sum(pp), 0)
end

"""
normalized_nullspace_svd(M::SparseMatrixCSC)

Robust nullspace via SVD. Computes the right singular vector
associated with the smallest singular value and normalizes it.
More expensive than QR but numerically safer; useful as a
fallback and for benchmarking.
"""
function normalized_nullspace_svd(M::SparseMatrixCSC)
    Md = Matrix(M)
    Fsvd = svd(Md)
    v = Fsvd.V[:, end]  # smallest singular value (svd sorts descending)
    s = sum(v)
    if s != 0 && isfinite(s)
        v ./= s
    end
    max.(v, 0)
end

"""
normalized_nullspace(M::SparseMatrixCSC)

Default steady-state solver. Tries the fast QR-based method first
and falls back to SVD if the QR path produces a non-finite or
ill-normalized vector.
"""
function normalized_nullspace(M::SparseMatrixCSC)
    v = normalized_nullspace_qr(M)
    if all(isfinite, v) && (s = sum(v); isfinite(s) && s != 0)
        return v
    end
    normalized_nullspace_svd(M)
end

"""
    normalized_nullspace_augmented(M::SparseMatrixCSC)

Compute a steady-state vector `p` with `M * p ≈ 0` and `sum(p) == 1` by solving the
`n×n` linear system formed from the first `n-1` rows of `M` and the constraint
`sum(p) = 1`. This avoids sparse QR pivoting and SVD, so it is suitable for
forward-mode AD (e.g. `ForwardDiff` through the linear solve; can be fragile for sparse `\\`).
Reverse mode: `ChainRulesCore.rrule` is defined for this function, so Zygote-style AD
will use the hand-coded pullback. You can also call [`pullback_normalized_nullspace_augmented`](@ref) directly.

The last row of `M` is not used; for a CTMC generator with row sums zero, that row
is linearly dependent on the others.

Output is clipped and renormalized like [`normalized_nullspace`](@ref): negative
entries from roundoff are projected with `max(·,0)` and the vector is rescaled to
sum to 1. If the augmented system is singular or non-finite, falls back to
[`normalized_nullspace_svd`](@ref) (then [`pullback_normalized_nullspace_augmented`](@ref) errors).

# See also
- [`normalized_nullspace`](@ref) — fast QR default for plain `Float64` use.
"""
function normalized_nullspace_augmented(M::SparseMatrixCSC{T}) where T
    n = size(M, 1)
    n == 0 && return T[]
    n == 1 && return ones(T, 1)
    r = _augmented_nullspace_try(M, n)
    return r === nothing ? normalized_nullspace_svd(M) : r.out
end

"""
    _augmented_nullspace_try(M, n) -> Union{Nothing, NamedTuple}

On success, returns `(; out, A, p, s)` where `out` is the normalized steady state,
`A` and `p` satisfy `A * p = b` with `b = [0…0, 1]`, and `s = sum(max.(p,0))`.
Returns `nothing` if the augmented system fails (singular / non-finite).
"""
function _augmented_nullspace_try(M::SparseMatrixCSC{T}, n::Int) where T
    A = vcat(M[1:n-1, :], sparse(ones(T, 1, n)))
    # Dense RHS required: sparse `b` makes `A \\ b` error on some Julia versions
    b = vcat(zeros(T, n - 1), one(T))
    p = try
        A \ b
    catch
        return nothing
    end
    any(!isfinite, p) && return nothing
    v = max.(p, zero(T))
    s = sum(v)
    (!(isfinite(s)) || s == 0) && return nothing
    out = v ./ s
    return (; out, A, p, s)
end

"""Returns `(output, augmented_path_ok)`; if `!augmented_path_ok`, output is from SVD."""
function _normalized_nullspace_augmented_core(M::SparseMatrixCSC, n::Int)
    r = _augmented_nullspace_try(M, n)
    r === nothing ? (normalized_nullspace_svd(M), false) : (r.out, true)
end

"""
    pullback_normalized_nullspace_augmented(M::SparseMatrixCSC, ȳ::AbstractVector)
    pullback_normalized_nullspace_augmented(M, ȳ, r)

Hand-coded reverse-mode map for the augmented steady-state: returns primal
`p = normalized_nullspace_augmented(M)` when the augmented linear system succeeds,
and sparse `M̄` with the same structure as `M`, representing the gradient of
`dot(ȳ, p)` with respect to the entries of `M` (only rows `1:n-1` receive nonzero
sensitivities; the last row of `M` is unused in the augmented formulation).

If the augmented solve would fall back to SVD, throws `ArgumentError` (call
[`normalized_nullspace`](@ref) for a primal-only robust path).

The three-argument form reuses the named tuple `r` from [`_augmented_nullspace_try`](@ref)
so a `ChainRulesCore.rrule` can avoid repeating the forward solve.

Use this from `ChainRulesCore.rrule`, Zygote `@adjoint`, or a custom optimizer.
"""
function pullback_normalized_nullspace_augmented(M::SparseMatrixCSC{T}, ȳ::AbstractVector) where T
    n = size(M, 1)
    n == 0 && return T[], spzeros(T, 0, 0)
    n == 1 && return ones(T, 1), spzeros(T, 1, 1)
    r = _augmented_nullspace_try(M, n)
    if r === nothing
        throw(
            ArgumentError(
                "normalized_nullspace_augmented fell back to SVD; pullback is undefined. " *
                "Use a well-conditioned augmented system or call normalized_nullspace for primal only.",
            ),
        )
    end
    return pullback_normalized_nullspace_augmented(M, ȳ, r)
end

function pullback_normalized_nullspace_augmented(
    M::SparseMatrixCSC{T},
    ȳ::AbstractVector,
    r::NamedTuple,
) where T
    n = size(M, 1)
    (; out, A, p, s) = r
    output = out
    ȳ = ȳ
    v̄ = (ȳ .- dot(ȳ, output)) ./ s
    p̄ = v̄ .* (p .> 0)
    λ = A' \ p̄
    Ā = -λ * p'
    G = zeros(T, n, n)
    G[1:n-1, :] .= Ā[1:n-1, :]
    Ii, Jj, _ = findnz(M)
    M̄_vals = map((i, j) -> G[i, j], Ii, Jj)
    M̄ = sparse(Ii, Jj, M̄_vals, n, n)
    return output, M̄
end

function ChainRulesCore.rrule(::typeof(normalized_nullspace_augmented), M::SparseMatrixCSC{T}) where T
    n = size(M, 1)
    if n == 0
        return T[], _ -> (ChainRulesCore.NoTangent(), ChainRulesCore.ZeroTangent())
    end
    if n == 1
        return ones(T, 1), _ -> (ChainRulesCore.NoTangent(), ChainRulesCore.ZeroTangent())
    end
    r = _augmented_nullspace_try(M, n)
    if r === nothing
        out = normalized_nullspace_svd(M)
        function normalized_nullspace_augmented_pullback_svd(ȳ)
            ChainRulesCore.unthunk(ȳ)
            throw(
                ArgumentError(
                    "normalized_nullspace_augmented fell back to SVD; reverse-mode AD is undefined for this input.",
                ),
            )
        end
        return out, normalized_nullspace_augmented_pullback_svd
    end
    p = r.out
    function normalized_nullspace_augmented_pullback_rr(ȳ)
        ȳ = ChainRulesCore.unthunk(ȳ)
        _, M̄ = pullback_normalized_nullspace_augmented(M, ȳ, r)
        return ChainRulesCore.NoTangent(), M̄
    end
    return p, normalized_nullspace_augmented_pullback_rr
end

function normalized_nullspace_ad(M::SparseMatrixCSC)
    m = size(M, 1)
    F = qr(M)
    R = F.R
    # Back substitution to solve R*p = 0, with p[end] = 1.0
    function backsub(R)
        p = [1.0]
        for i in m-1:-1:1
            val = -R[i, i+1:end]' * reverse(p) / R[i, i]
            p = [val; p]
        end
        p
    end
    p = backsub(R)
    # AD-friendly permutation
    pp = [p[findfirst(==(i), F.pcol)] for i in 1:length(p)]
    max.(pp / sum(pp), 0)
end
"""
    allele_convolve(mhist,nalleles)

    Convolve to compute distribution for contributions from multiple alleles
"""
function _allele_convolve_step(prev::AbstractVector, mhist::AbstractVector)
    nhist = length(mhist)
    return [
        sum(prev[m-m2+1] * mhist[m2+1] for m2 in 0:min(nhist - 1, m)) for m in 0:nhist-1
    ]
end

function allele_convolve(mhist, nalleles)
    nhist = length(mhist)
    a = float.(collect(mhist))
    for _ in 2:nalleles
        a = _allele_convolve_step(a, mhist)
    end
    return a
end
"""
    allele_deconvolve(mhist,nalleles)

    Deconvolve to compute distribution of one allele from contributions of multiple alleles
"""
allele_deconvolve(mhist, nalleles) = irfft((rfft(mhist)) .^ (1 / nalleles), length(mhist))


"""
nhist_loss(nhist, yieldfactor; threshold=0.999)

Compute length of pre-loss histogram using a principled probabilistic approach.

Given an observed histogram with `nhist` bins (support 0 to nhist-1), this function
estimates the true histogram size needed to account for binomial observation loss.

# Arguments
- `nhist::Int`: Size of observed histogram (number of bins)
- `yieldfactor::Float64`: Detection efficiency (0-1)
- `threshold::Float64`: Cumulative probability threshold (default 0.99)
  - Finds the smallest true count where P(observe ≤ nhist-1 | true = j) ≥ threshold
  - Lower values (e.g., 0.95) are less conservative and reduce computational cost

# Returns
- `Int`: Estimated true histogram size

# Method
Uses cumulative probability: finds the smallest true count `j` such that the
probability of observing at most `nhist-1` given true count `j` is at least `threshold`.
This ensures we capture at least `threshold` fraction of the probability mass.

# Example
```julia
# If we observe up to 20 counts with 5% yield, what's the true range?
nhist_loss(21, 0.05)  # Returns ~300-350 (with threshold=0.99)
nhist_loss(21, 0.05, threshold=0.95)  # Returns ~250-300 (less conservative)
```
"""
function nhist_loss(nhist::Int, yieldfactor::Float64; threshold::Float64=0.99)
    if yieldfactor >= 1.0
        return nhist
    end
    
    max_observed = nhist - 1  # 0-indexed: observed counts are 0 to nhist-1
    
    # Binary search for the smallest true count j where 
    # P(observe ≤ max_observed | true = j) ≥ threshold
    # Start with simple scaling as lower bound
    lower = max(nhist, round(Int, max_observed / yieldfactor))
    # Upper bound: use a conservative multiplier
    upper = max(lower + 1, round(Int, max_observed / yieldfactor * 3))
    
    # Binary search
    while upper - lower > 1
        mid = (lower + upper) ÷ 2
        d = Binomial(mid, clamp(yieldfactor, 0.0, 1.0))
        cumprob = cdf(d, max_observed)
        
        if cumprob >= threshold
            upper = mid
        else
            lower = mid
        end
    end
    
    # Check final bounds
    d_lower = Binomial(lower, clamp(yieldfactor, 0.0, 1.0))
    if cdf(d_lower, max_observed) >= threshold
        return lower
    else
        return upper
    end
end

"""
technical_loss(mhist,yieldfactor)

Reduce counts due to technical loss
"""
function technical_loss!(mhist::Vector{<:Vector}, yieldfactor)
    for i in eachindex(mhist)
        mhist[i] = technical_loss(mhist[i], yieldfactor, length(mhist[i]))
    end
end
"""
technical_loss(mhist,yieldfactor,nhist)

Reduce counts due to technical loss using Binomial sampling
"""
function technical_loss(mhist::Vector, yieldfactor, nhist)
    p = zeros(nhist)
    for m in eachindex(mhist)
        d = Binomial(m - 1, clamp(yieldfactor, 0.0, 1.0))
        for c in 1:m+1
            p[c] += mhist[m] * pdf(d, c - 1)
        end
    end
    normalize_histogram(p)
end

function technical_loss_at_k(k::Int, mhist, yieldfactor, nhist_true)
    pmf = normalize_histogram(mhist)
    p = 0
    # Sum over all true counts m >= k up to nhist_true-1
    # nhist_true is the inflated size (nRNA_true) when yield < 1.0
    for m in k:nhist_true-1
        d = Binomial(m, clamp(yieldfactor, 0.0, 1.0))
        p += pmf[m+1] * pdf(d, k)
    end
    return p
end
"""
technical_loss_poisson(mhist,yieldfactor,nhist)

Reduce counts due to technical loss using Poisson sampling
"""
function technical_loss_poisson(mhist, yieldfactor, nhist)
    p = zeros(nhist)
    for m in eachindex(mhist)
        d = Poisson(yieldfactor * (m - 1))
        for c in 1:nhist
            p[c] += mhist[m] * pdf(d, c - 1)
        end
    end
    normalize_histogram(p)
end

"""
additive_noise(mhist,noise,nhist)

Add Poisson noise to histogram
"""
function additive_noise(mhist, noise, nhist)
    p = zeros(nhist)
    d = Poisson(noise)
    for m in 1:nhist
        for n in 1:m
            p[m] += mhist[n] * pdf(d, m - n)
        end
    end
    normalize_histogram(p)
end
"""
threshold_noise(mhist,noise,yieldfactor,nhist)

Add Poisson noise to histogram then reduce counts due to technical loss

"""
function threshold_noise(mhist, noise, yieldfactor, nhist)
    h = additive_noise(mhist, noise, nhist)
    technical_loss(h, yieldfactor, nhist)
end

"""
    make_loss_matrix(nRNA_observed::Int, nRNA_true::Int, yieldfactor::Float64)

Create a loss matrix L where L[i+1, j+1] = P(observe i | true count j)
using binomial coefficients.

# Arguments
- `nRNA_observed::Int`: Size of observed histogram (number of bins)
- `nRNA_true::Int`: Size of true histogram (number of bins, typically larger when yieldfactor < 1.0)
- `yieldfactor::Float64`: Probability of observing each mRNA (0-1). 
  This is the detection efficiency or yield factor.

# Returns
- `Matrix{Float64}`: Loss matrix of size (nRNA_observed × nRNA_true)
  - Rows (i+1): observed mRNA count i (0 to nRNA_observed-1)
  - Columns (j+1): true mRNA count j (0 to nRNA_true-1)
  - L[i+1, j+1] = P(observe i | true count j) = Binomial(j, yieldfactor).pdf(i)

# Notes
- The loss matrix accounts for observation noise where each mRNA has probability
  `yieldfactor` of being detected.
- Columns are normalized to sum to 1 (accounting for truncation at nRNA_observed).
- If yieldfactor = 1.0, the matrix is identity (no loss).
- If yieldfactor < 1.0, the matrix spreads probability from higher true counts
  to lower observed counts.

# Example
```julia
# 5% yield: observe up to 20, true up to ~400
nRNA_true = nhist_loss(21, 0.05)  # ~300-350
L = make_loss_matrix(21, nRNA_true, 0.05)
# Apply to predicted histogram: p_observed = L * p_true
```
"""
function make_loss_matrix(nRNA_observed::Int, nRNA_true::Int, yieldfactor::Float64)
    L = zeros(nRNA_observed, nRNA_true)
    for j in 0:nRNA_true-1  # true count (column)
        for i in 0:min(j, nRNA_observed-1)  # observed count (row, can't observe more than true)
            d = Binomial(j, clamp(yieldfactor, 0.0, 1.0))
            L[i+1, j+1] = pdf(d, i)
        end
        # Normalize column to account for truncation (probabilities for i > nRNA_observed-1 are lost)
        col_sum = sum(L[:, j+1])
        if col_sum > 0
            L[:, j+1] ./= col_sum
        end
    end
    return L
end

# Backward compatibility: if only nRNA provided, assume square matrix (old behavior)
make_loss_matrix(nRNA::Int, yieldfactor::Float64) = make_loss_matrix(nRNA, nRNA, yieldfactor)

"""
solve_vector(A::Matrix,b::vector)
solve A x = b
If matrix divide has error higher than tol
use SVD and pseudoinverse with threshold
"""
function solve_vector(A::Matrix, b::Vector, th=1e-16, tol=1e-1)
    x = A \ b
    if norm(b - A * x, Inf) > tol
        M = svd(A)
        Sv = M.S
        Sv[abs.(Sv).<th] .= 0.0
        Sv[abs.(Sv).>=th] = 1 ./ Sv[abs.(Sv).>=th]
        x = M.V * diagm(Sv) * M.U' * b
    end
    return x[:, 1] # return as vector
end


