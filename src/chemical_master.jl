# This file is part of StochasticGene.jl
#
# chemical_master.jl
#
# Functions to solve the chemical master equation
# operates on the transpose of the Markov process transition rate matrix

using ForwardDiff

"""Sparse `\\` / sparse `qr` / `svd` use LAPACK and do not support `ForwardDiff.Dual`; abstract `Real` eltypes can also hide dual numbers, so route both cases through dense generic linear algebra."""
_use_dense_ad_eltype(::Type{T}) where {T<:Real} = !isconcretetype(T) || T <: ForwardDiff.Dual

function _concrete_value_eltype(vals)
    isempty(vals) && return Float64
    T = typeof(first(vals))
    for x in Iterators.drop(vals, 1)
        T = promote_type(T, typeof(x))
    end
    return T
end

function _dense_generic_matrix(M::SparseMatrixCSC)
    Tv = _concrete_value_eltype(nonzeros(M))
    D = zeros(Tv, size(M))
    rows, cols, vals = findnz(M)
    @inbounds for k in eachindex(vals)
        D[rows[k], cols[k]] = vals[k]
    end
    return D
end

function _dense_generic_vector(v::AbstractVector)
    Tv = _concrete_value_eltype(v)
    out = Vector{Tv}(undef, length(v))
    @inbounds for i in eachindex(v)
        out[i] = v[i]
    end
    return out
end

function _nullspace_robust_fallback(M::SparseMatrixCSC{T}) where T
    if _use_dense_ad_eltype(T)
        return normalized_nullspace_qr(M)
    end
    try
        v = normalized_nullspace_qr(M)
        if all(isfinite, v)
            s = sum(v)
            if isfinite(s) && s != 0
                return v
            end
        end
    catch
    end
    return normalized_nullspace_svd(M)
end

function _augmented_linear_solve(A::SparseMatrixCSC{T}, b::AbstractVector) where T
    _use_dense_ad_eltype(T) && return _dense_generic_matrix(A) \ _dense_generic_vector(b)
    return A \ b
end

# Functions to compute steady state mRNA histograms
"""
    steady_state_vector(M::AbstractMatrix; solver::Symbol=:default)

Return a normalized steady-state vector `p` with `M * p ≈ 0` and entries summing to 1.

- `solver=:default` or `:fast` — [`normalized_nullspace`](@ref) (sparse QR, SVD fallback).
- `solver=:augmented` or `:autodiff` — [`normalized_nullspace_augmented`](@ref) (linear solve + AD rules).

Non-sparse `M` is converted with `sparse(M)`.
"""
function steady_state_vector(M::SparseMatrixCSC; solver::Symbol=:default, nT::Union{Nothing,Int}=nothing)
    if solver === :default || solver === :fast
        return normalized_nullspace(M)
    elseif solver === :augmented || solver === :autodiff
        return normalized_nullspace_augmented(M)
    elseif solver === :closure_thomas
        nT === nothing && throw(
            ArgumentError(
                "steady_state_vector: solver=:closure_thomas requires nT keyword argument",
            ),
        )
        return normalized_nullspace_closure_thomas(M, nT)
    else
        throw(
            ArgumentError(
                "steady_state_vector: unknown solver $(repr(solver)); use :default, :fast, :augmented, :autodiff, or :closure_thomas",
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
    P = steady_state_vector(M; solver=steady_state_solver, nT=nT)
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

# Zero vector for dwell-time entrance-flux accumulation; must promote element type so
# ForwardDiff/Zygote can propagate (e.g. rates and/or pss carry Dual numbers).
function _vector_zero_sinit_like(pss::AbstractVector, others::AbstractVector...)
    T = eltype(pss)
    for v in others
        T = promote_type(T, eltype(v))
    end
    return zeros(T, length(pss))
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
    Sinit = _vector_zero_sinit_like(pss, r)
    for e in elements
        if e.b != e.a && (e.a ∈ sojourn && e.b ∉ sojourn)
            Sinit[e.a] += pss[e.b] * r[e.index]
        end
    end
    s = sum(Sinit)
    (iszero(s) || !isfinite(s)) ? Sinit : Sinit / s
end

function init_S(r::Vector, sojourn::Vector, elements::Vector, pss, nonzeros)
    Sinit = _vector_zero_sinit_like(pss, r)
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
    rows, cols, vals = findnz(T)
    Sinit = _vector_zero_sinit_like(pss, vals)
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
    rows, cols, vals = findnz(T)
    Sinit = _vector_zero_sinit_like(pss, vals)
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
    Sinit = _vector_zero_sinit_like(pss, r)
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

Implementation avoids in-place `setindex!` so reverse-mode AD (Zygote) can differentiate
when this path is used; for gradients prefer [`normalized_nullspace_augmented`](@ref) with
`steady_state_solver=:augmented`.
"""
function _backsub_upper_tri_one(R::AbstractMatrix, m::Int)
    Tr = eltype(R)
    m == 1 && return [one(Tr)]
    foldr(1:m-1; init=[one(Tr)]) do k, p_tail
        pk = -dot(@view(R[k, (k+1):m]), p_tail) / R[k, k]
        [pk; p_tail]
    end
end

function normalized_nullspace_qr(M::SparseMatrixCSC{T}) where T
    m = size(M, 1)
    Md = _use_dense_ad_eltype(T) ? _dense_generic_matrix(M) : M
    F = qr(Md)   # dense QR for ForwardDiff.Dual / abstract Real paths
    R = Matrix(F.R)  # convert to dense once to avoid repeated sparse slicing in back substitution
    p = _backsub_upper_tri_one(R, m)
    pp = hasproperty(F, :pcol) ? p[invperm(getproperty(F, :pcol))] : p
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
    Md = !isconcretetype(eltype(M)) ? _dense_generic_matrix(M) : Matrix(M)
    Fsvd = try
        svd(Md)
    catch
        return normalized_nullspace_qr(M)
    end
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
    _nullspace_robust_fallback(M)
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
[`normalized_nullspace_svd`](@ref) for `Float32`/`Float64`, or dense [`normalized_nullspace_qr`](@ref)
for `ForwardDiff.Dual` (SVD is unavailable). [`pullback_normalized_nullspace_augmented`](@ref) errors on any fallback.

# See also
- [`normalized_nullspace`](@ref) — fast QR default for plain `Float64` use.
"""
function normalized_nullspace_augmented(M::SparseMatrixCSC{T}) where T
    n = size(M, 1)
    n == 0 && return T[]
    n == 1 && return ones(T, 1)
    r = _augmented_nullspace_try(M, n)
    return r === nothing ? _nullspace_robust_fallback(M) : r.out
end

"""
    _augmented_nullspace_try(M, n) -> Union{Nothing, NamedTuple}

On success, returns `(; out, A, p, s)` where `out` is the normalized steady state,
`A` and `p` satisfy `A * p = b` with `b = [0…0, 1]`, and `s = sum(max.(p,0))`.
Returns `nothing` if the augmented system fails (singular / non-finite).
"""
function _augmented_nullspace_try(M::SparseMatrixCSC{T}, n::Int) where T
    A = vcat(M[1:n-1, :], sparse(fill(one(T), 1, n)))
    # Dense RHS required: sparse `b` makes `A \\ b` error on some Julia versions
    b = vcat(zeros(T, n - 1), one(T))
    p = try
        _augmented_linear_solve(A, b)
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

"""Returns `(output, augmented_path_ok)`; if `!augmented_path_ok`, output is from SVD or dense QR (Dual)."""
function _normalized_nullspace_augmented_core(M::SparseMatrixCSC, n::Int)
    r = _augmented_nullspace_try(M, n)
    r === nothing ? (_nullspace_robust_fallback(M), false) : (r.out, true)
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
        out = _nullspace_robust_fallback(M)
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
    nhist_loss(nhist::Int, yieldfactor::Float64; threshold::Float64=0.99, max_nRNA::Int=500)

Estimate the true histogram size from observed histogram size and yield factor.

## Theory

After yield sampling:
- Input `nhist` = observed histogram size (bin containing 99% cumulative probability)
- `yieldfactor` = observed_mean / true_mean
- Since means and quantiles scale together: `true_size ≈ observed_size / yieldfactor`

## Parameters

- `nhist::Int`: Observed histogram size (99th percentile point)
- `yieldfactor::Float64`: Yield factor = observed_mean / true_mean
- `threshold::Float64`: Unused; kept for API compatibility
- `max_nRNA::Int`: Hard cap on expansion (default 500)

## Returns

Estimated true histogram size: `round(nhist / yieldfactor)`, capped at `max_nRNA`.
"""
function nhist_loss(nhist::Int, yieldfactor::Float64; threshold::Float64=0.99, max_nRNA::Int=500)
    if yieldfactor >= 0.99 || nhist == 0
        return nhist
    end
    
    if yieldfactor <= 0.0
        return nhist
    end
    
    # Simple inversion: true ≈ observed / yield
    nRNA_true = round(Int, nhist / yieldfactor)
    
    # Apply cap
    return min(nRNA_true, max_nRNA)
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
        for c in 1:min(m+1, nhist)
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

**MH / primal path:** mutating implementation (preallocated matrix + column normalization in place).
For reverse-mode AD (Zygote / NUTS with `get_rates_ad`), use [`make_loss_matrix_ad`](@ref).

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
        col_sum = sum(@view L[:, j+1])
        if col_sum > 0
            L[:, j+1] ./= col_sum
        end
    end
    return L
end

"""
    make_loss_matrix_ad(nRNA_observed::Int, nRNA_true::Int, yieldfactor::Float64)

Same mathematics as [`make_loss_matrix`](@ref), but **no in-place column writes** (builds columns and
`hcat`). Use with Zygote / [`get_rates_ad`](@ref) / NUTS. For MH, keep [`make_loss_matrix`](@ref).
"""
function make_loss_matrix_ad(nRNA_observed::Int, nRNA_true::Int, yieldfactor::Float64)
    y = clamp(yieldfactor, 0.0, 1.0)
    cols = map(0:nRNA_true-1) do j
        col = [pdf(Binomial(j, y), i) for i in 0:nRNA_observed-1]
        s = sum(col)
        s > 0 ? col ./ s : col
    end
    return reduce(hcat, cols)::Matrix{Float64}
end

# Backward compatibility: if only nRNA provided, assume square matrix (old behavior)
make_loss_matrix(nRNA::Int, yieldfactor::Float64) = make_loss_matrix(nRNA, nRNA, yieldfactor)
make_loss_matrix_ad(nRNA::Int, yieldfactor::Float64) = make_loss_matrix_ad(nRNA, nRNA, yieldfactor)

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

# ============================================================================
# Closure-Thomas steady-state solver
#
# Block-tridiagonal solver derived from the steady-state generating-function
# ODE `(A + zB) G(z) + mu (1-z) G'(z) = 0`. Matching coefficients of `z^n`
# gives the three-term recurrence
#   mu (m+1) P_{m+1} - (mu m I - A) P_m + B P_{m-1} = 0,    P_{-1} = 0,
# a block-tridiagonal system in `m` with `n_T x n_T` blocks. We exploit this
# structure directly via Thomas elimination, avoiding the `O((n_T (M+1))^3)`
# dense / sparse factorization of the full joint generator.
#
# The closure-Thomas path matches the reflective-boundary correction that
# `make_mat_M` applies (the extra `+B` at the cap block to conserve total
# probability); the boundary block becomes `B R_{Mcap-1} + (A + B) - mu Mcap I`,
# which is the same singular operator whose null space the existing
# augmented LU solver targets. Cost: `O((Mcap+1) n_T^3)` for the forward sweep
# plus an `O(n_T^3)` augmented solve at the cap. Beats the existing sparse LU
# significantly once `n_T >= ~12` and the mRNA cap is moderate.
# ============================================================================

"""
    _extract_AB_mu_from_M(M, nT)

Reverse-engineer the closure-form `(A, B, mu, nhist)` from the assembled joint
generator built by [`make_mat_M`](@ref). Assumes the standard m-major Kronecker
layout: block `(m_julia, m_julia')` of `M` occupies rows
`(m_julia - 1) * nT + 1 : m_julia * nT` and similar for columns.

The Julia-side `m_julia = 1` block carries no decay term, so it equals the
non-productive part `A = T - B`. The Julia-side super-diagonal block at
`(1, 2)` is `decay * I_nT`, and the sub-diagonal block at `(2, 1)` is `B`.
"""
function _extract_AB_mu_from_M(M::SparseMatrixCSC, nT::Int)
    N = size(M, 1)
    nhist = div(N, nT)
    nhist * nT == N || throw(ArgumentError("size(M, 1) = $N is not divisible by nT = $nT"))
    Tv = _concrete_value_eltype(nonzeros(M))
    A = Matrix{Tv}(undef, nT, nT)
    @inbounds for j in 1:nT, i in 1:nT
        A[i, j] = M[i, j]
    end
    mu = nhist >= 2 ? M[1, nT + 1] : zero(Tv)
    B = Matrix{Tv}(undef, nT, nT)
    if nhist >= 2
        @inbounds for j in 1:nT, i in 1:nT
            B[i, j] = M[nT + i, j]
        end
    else
        fill!(B, zero(Tv))
    end
    return A, B, mu, nhist
end

"""
    normalized_nullspace_closure_thomas(M::SparseMatrixCSC, nT::Int)

Steady-state vector of the truncated joint generator `M` via block-Thomas
elimination on the closure recurrence. `nT` is the hidden-state dimension
(the block size). The returned vector is in the same m-major linear layout
as `M`, with nonnegative entries and unit total mass.

The reflective boundary correction that [`make_mat_M`](@ref) installs at the
cap is honored: the cap block uses `(A + B) - mu * Mcap * I`, matching the
existing augmented-LU steady state.

Forward-mode AD through `ForwardDiff.Dual` works out of the box (every step
is a dense `n_T x n_T` LU); a reverse-mode `rrule` is not yet provided, so
Zygote-style gradients should continue to use the `:augmented` solver.
"""
function normalized_nullspace_closure_thomas(M::SparseMatrixCSC, nT::Int)
    A, B, mu, nhist = _extract_AB_mu_from_M(M, nT)
    Tv = eltype(A)
    if nhist == 1
        # No mRNA dimension; just the finite-state stationary distribution.
        return _nullspace_robust_fallback(M)
    end
    Mcap = nhist - 1
    Imat = Matrix{Tv}(I, nT, nT)
    R = Vector{Matrix{Tv}}(undef, Mcap)
    # R[1] = R_0: P_0 = R_0 P_1, with A P_0 + mu P_1 = 0.
    R[1] = A \ (-mu .* Imat)
    @inbounds for m in 1:(Mcap - 1)
        Mm = B * R[m] .+ A .- (mu * m) .* Imat
        R[m + 1] = Mm \ (-(mu * (m + 1)) .* Imat)
    end
    # Cap block with reflective boundary: M_Mcap = B R_{Mcap-1} + (A + B) - mu Mcap I.
    M_cap = B * R[Mcap] .+ A .+ B .- (mu * Mcap) .* Imat
    # M_cap is singular by mass conservation; augment with sum-of-block = 1.
    aug = copy(M_cap)
    @inbounds for j in 1:nT
        aug[nT, j] = one(Tv)
    end
    rhs = zeros(Tv, nT)
    rhs[nT] = one(Tv)
    P_Mcap = aug \ rhs
    # Back-substitute P_{m-1} = R_{m-1} P_m.
    P_blocks = Vector{Vector{Tv}}(undef, nhist)
    P_blocks[nhist] = P_Mcap
    @inbounds for j in (nhist - 1):-1:1
        P_blocks[j] = R[j] * P_blocks[j + 1]
    end
    # Flatten in m-major layout matching M.
    p = Vector{Tv}(undef, size(M, 1))
    @inbounds for m in 1:nhist, a in 1:nT
        p[(m - 1) * nT + a] = P_blocks[m][a]
    end
    # Global normalization; flip sign if the kernel direction came out negative.
    s = sum(p)
    if isfinite(s) && s != 0
        if s < 0
            @inbounds for k in eachindex(p)
                p[k] = -p[k]
            end
            s = -s
        end
        @inbounds for k in eachindex(p)
            p[k] /= s
        end
        return max.(p, zero(Tv))
    end
    # Augmented solve degenerated; fall back to the robust nullspace path.
    return _nullspace_robust_fallback(M)
end


