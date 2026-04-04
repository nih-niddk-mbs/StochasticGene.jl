# This file is part of StochasticGene.jl
#
# chemical_master.jl
#
# Functions to solve the chemical master equation
# operates on the transpose of the Markov process transition rate matrix


# Functions to compute steady state mRNA histograms
"""
    steady_state(M, nT, nalleles, nhist)
    steady_state(M, nT, nalleles)

Return steady-state mRNA histogram for G and GR models by computing the null space of the
truncated full transition matrix.

# Arguments
- `M`: Truncated full transition rate matrix including transcription state transitions (GR) and mRNA birth and death
- `nT`: Dimension of state transitions
- `nalleles`: Number of alleles (for convolution)
- `nhist`: Length of histogram to return (four-argument form only)

# Returns
- `steady_state(M, nT, nalleles, nhist)`: Vector of length `nhist`
- `steady_state(M, nT, nalleles)`: Full marginal histogram convolved over alleles
"""
steady_state(M, nT, nalleles, nhist) = steady_state(M, nT, nalleles)[1:nhist]

function steady_state(M, nT, nalleles)
    P = normalized_nullspace(M)
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

function offonPDF(t::Vector, r::Vector, T::AbstractMatrix, TA::AbstractMatrix, TI::AbstractMatrix, nT::Int, elementsT::Vector, onstates::Vector)
    pss = normalized_nullspace(T)
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
function marginalize(p::Vector, nT, nhist)
    mhist = zeros(nhist)
    for m in 1:nhist
        i = (m - 1) * nT
        mhist[m] = sum(p[i+1:i+nT])
    end
    return mhist
end

function marginalize(p::Vector, nT)
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
function allele_convolve(mhist, nalleles)
    nhist = length(mhist)
    mhists = Array{Array{Float64,1}}(undef, nalleles)
    mhists[1] = float.(mhist)
    for i = 2:nalleles
        mhists[i] = zeros(nhist)
        for m = 0:nhist-1
            for m2 = 0:min(nhist - 1, m)
                mhists[i][m+1] += mhists[i-1][m-m2+1] * mhist[m2+1]
            end
        end
    end
    return mhists[nalleles]
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


