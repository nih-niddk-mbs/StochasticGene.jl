# This file is part of StochasticGene.jl
#
# chemical_master.jl
#
# Functions to solve the chemical master equation
# operates on the transpose of the Markov process transition rate matrix


# Functions to compute steady state mRNA histograms
"""
steady_state(M,nT,nalleles,nhist)
steady_state(M,nT,nalleles)

return steady state mRNA histogram for G and GR models

computes null space of the truncated full transition matrix.

-`M`: Truncated full transition rate matrix including transcription state transitions (GR) and mRNA birth and death
-`nT`: dimension of state transitions

"""
steady_state(M, nT, nalleles, nhist) = steady_state(M, nT, nalleles)[1:nhist]

function steady_state(M, nT, nalleles)
    P = normalized_nullspace(M)
    mhist = marginalize(P, nT)
    allele_convolve(mhist, nalleles)
end



# Functions to compute dwell time distributions

"""
    dwelltimeCDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=lsoda())

  return dwell time CDF (cumulative distribution function) to reach barrier states starting from Sinit in live states
  live U barrier = all states

  First passage time calculation from live states to barrier states

- `tin`: vector of time bins for dwell time distribution
- `Td`: dwell time rate matrix (within non barrier states)
- `barrier`: set of barrier states
- `Sinit`: initial state
- `method`: 1 for directly solving ODE otherwise use eigendecomposition
"""
function dwelltimeCDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=lsoda())
    t = [max(2 * tin[1] - tin[2], 0.); tin]
    S = time_evolve(t, Td, Sinit, method)
    return vec(sum(S[:, barrier], dims=1))
end

"""
    dwelltimePDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=lsoda())

return dwell time PDF (probability density function)
"""
dwelltimePDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method=lsoda()) = pdf_from_cdf(dwelltimeCDF(tin, Td, barrier, Sinit, method))


ontimePDF(tin::Vector, TA::AbstractMatrix, offstates, SAinit::Vector, method=lsoda()) = dwelltimePDF(tin, TA, offstates, SAinit, method)

offtimePDF(tin::Vector, TI::AbstractMatrix, onstates, SIinit::Vector, method=lsoda()) = dwelltimePDF(tin, TI, onstates, SIinit, method)

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
    Sinit / sum(Sinit)
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

function init_S(sojourn::Vector, T::SparseMatrixCSC, pss)
    Sinit = zeros(length(pss))
    rows, cols, vals = findnz(T)
    for i in eachindex(rows)
        if cols[i] != rows[i] && (rows[i] ∈ sojourn && cols[i] ∉ sojourn)
            Sinit[rows[i]] += pss[cols[i]] * vals[i]
        end
    end
    Sinit / sum(Sinit)
end


function init_S(sojourn::Vector, T::SparseMatrixCSC, pss, nonzeros) 
    Sinit = init_S(sojourn, T, pss)[nonzeros]
    Sinit / sum(Sinit)
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
    fevolve!(du,u::Vector, p, t)

in place update of du of ODE system for DifferentialEquations,jl
"""
function fevolve!(du, u::Vector, p, t)
    du .= p * u
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
normalized_nullspace(M::SparseMatrixCSC)
Compute the normalized null space of a nxn matrix
of rank n-1 using QR decomposition with pivoting
"""
function normalized_nullspace(M::SparseMatrixCSC)
    m = size(M, 1)
    p = zeros(m)
    F = qr(M)   #QR decomposition
    R = F.R
    # Back substitution to solve R*p = 0
    p[end] = 1.0
    for i in 1:m-1
        p[m-i] = -R[m-i, m-i+1:end]' * p[m-i+1:end] / R[m-i, m-i]
    end
    # Permute elements according to sparse matrix result
    pp = copy(p)
    for i in eachindex(p)
        pp[F.pcol[i]] = p[i]
    end
    pp / sum(pp)

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
nhist_loss(nhist,yieldfactor)

Compute length of pre-loss histogram
"""
nhist_loss(nhist, yieldfactor) = round(Int, nhist / yieldfactor)

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

function technical_loss_at_k(k::Int, mhist, yieldfactor, nhist)
    pmf = normalize_histogram(mhist)
    p = 0
    for m in k:nhist-1
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


