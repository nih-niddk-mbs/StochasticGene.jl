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
    dwelltimeCDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method::Int=1)

  return dwell time CDF (cumulative distribution function) to reach barrier states starting from Sinit in live states
  live U barrier = all states

  First passage time calculation from live states to barrier states

- `tin`: vector of time bins for dwell time distribution
- `Td`: dwell time rate matrix (within non barrier states)
- `barrier`: set of barrier states
- `Sinit`: initial state
- `method`: 1 for directly solving ODE otherwise use eigendecomposition
"""
function dwelltimeCDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method::Int=1)
    t = [2*tin[1]-tin[2]; tin]
    S = time_evolve(t, Td, Sinit, method)
    return sum(S[:, barrier], dims=2)[:, 1]
end

"""
    dwelltimePDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method::Int=1)

return dwell time PDF (probability density function)
"""
dwelltimePDF(tin::Vector, Td::AbstractMatrix, barrier, Sinit::Vector, method::Int=1) = pdf_from_cdf(dwelltimeCDF(tin, Td, barrier, Sinit, method))


ontimePDF(tin::Vector, TA::AbstractMatrix, offstates, SAinit::Vector, method::Int=1) = dwelltimePDF(tin, TA, offstates, SAinit, method)

offtimePDF(tin::Vector, TI::AbstractMatrix, onstates, SIinit::Vector, method::Int=1) = dwelltimePDF(tin, TI, onstates, SIinit, method)

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

marginalize(P::Matrix) = sum(P, dims=1)



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
    if method == 1
        return time_evolve_diff(t, Q, S0)
    else
        return time_evolve_eig(t, Q, S0)
    end
end
"""
time_evolve_diff(t,M::Matrix,P0)

Solve master equation problem using DifferentialEquations.jl
"""
function time_evolve_diff(t, Q::SparseMatrixCSC, P0)
    tspan = (t[1], t[end])
    prob = ODEProblem(fevolve!, P0, tspan, Q)
    sol = solve(prob, lsoda(), saveat=t)
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
    time_evolve(t, vals, vects, weights)
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
function technical_loss!(mhist::Vector, yieldfactor)
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
        for c in 1:nhist
            p[c] += mhist[m] * pdf(d, c - 1)
        end
    end
    normalize_histogram(p)
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




### Legacy code no longer compatible with new format

# """
#     time_evolve_delay(t, r0, r1, delay, n, P0)

# TBW
# """
# function time_evolve_delay(t, r0, r1, delay, n, P0)
#     tspan = (0.0, t[end])
#     nhist = Int(length(P0) / (n + 1)) - 2
#     p = [r0; r1; delay; n; nhist]
#     P0 = [P0; 1.0]
#     prob = ODEProblem(fdelay!, P0, tspan, p)
#     sol = solve(prob, lsoda(), saveat=t)
#     return sol'
# end

# """
#     fdelay!(du, u, p, t)

# TBW
# """
# function fdelay!(du, u, p, t)
#     n = Int(p[end-1])
#     nhist = Int(p[end])
#     delay = p[end-3]
#     r0 = p[1:2*(n+1)]
#     r1 = p[2*(n+1)+1:4*(n+1)]
#     r = r0 * u[end] + (1 - u[end]) * r1
#     M = mat_GM(r, n, nhist)
#     du[1:end-1] = M * u[1:end-1]
#     du[end] = -delay * u[end]
# end

# """
# transient(ts::Vector,r,n,nhist,nalleles,P0)

# Compute mRNA pmf of GM model at times in vector ts starting
# with initial condition P0
# """

# """
# transient(t,r,yieldfactor,n,nhist,nalleles,P0::Vector,method)
# transient(t,n,nhist,nalleles,P0,Mvals,Mvects)

# Compute mRNA pmf of GM model at time t given initial condition P0
# and eigenvalues and eigenvectors of model transition rate matrix
# method = 1 uses JuliaDifferentialEquations.jl
# method != 1 uses eigenvalue decomposition
# """
# function transient(t, r, yieldfactor, n, nalleles, P0::Vector, method)
#     mhist = transient(t, r, n, nalleles, P0, method)
#     technical_loss!(mhist, yieldfactor)
#     return mhist
# end
# function transient(t::Float64, r, n, nalleles, P0::Vector, method)
#     P = transient(t, r, n, P0, method)
#     mh = marginalize(P, n)
#     allele_convolve(mh, nalleles)
# end
# function transient(t::Vector, r, n, nalleles, P0::Vector, method)
#     P = transient(t, r, n, P0, method)
#     # P += abs.(rand(size(P)))
#     mhist = Array{Array,1}(undef, length(t))
#     for i in 1:size(P, 1)
#         mh = marginalize(P[i, :], n)
#         mhist[i] = allele_convolve(mh, nalleles)
#     end
#     return mhist
# end
# function transient(t, r, n, P0::Vector, method)
#     nhist = Int(length(P0) / (n + 1)) - 2
#     M = mat_GM(r, n, nhist)
#     time_evolve(t, M, P0, method)
# end

# function transient_delay(t::Vector, r0::Vector, r1::Vector, delay::Float64, n::Int, nalleles, P0::Vector)
#     P = time_evolve_delay(t, r0, r1, delay, n, P0)
#     mhist = Array{Array,1}(undef, length(t))
#     for i in 1:size(P, 1)
#         mh = marginalize(P[i, 1:end-1], n)
#         mhist[i] = allele_convolve(mh, nalleles)
#     end
#     return mhist
# end



# """
# steady_state(r,transitions,G,R,nhist,nalleles,splicetype="")

# return steady state mRNA histogram for G,R model
#     calls steady_state(M,nT,nalleles)

# -`r`: transition rates
# -`transitions`: tuple of G state transitions
# -`G`: number of G states
# -`R`: number of R states
# -`nhists`: number of mRNA histogram bins
# -`nalleles`: number of alleles (alleles are assumed to be uncoupled)
# -`splicetype`: parameter to specify type of kinetics if necessary

# """
# function steady_state(r, transitions, G, R, nhist, nalleles, splicetype="")
#     components = make_components(transitions, G, R, r, nhist + 2, splicetype, set_indices(length(transitions), R))
#     M = make_mat_M(mcomponents.mcomponents, r)
#     steady_state(M, components.mcomponents.nT, nalleles, nhist)
# end
# """
# steady_state(r,transitions,G,nhist,nalleles)

# returns steaday state mRNA histogram for G models
#     calls steady_state(M,nT,nalleles)

# """
# function steady_state(r, transitions, G, nhist, nalleles)
#     components = make_components_M(transitions, G, nhist + 2, r[end])
#     M = make_mat_M(components, r)
#     steady_state(M, G, nalleles, nhist)
# end


# """
# ontimeCDF(tin::Vector,G::Int,R::Int,TA::AbstractMatrix,pss::Vector,method)

# returns ON dwell time histogram (PDF) of GRS model
# Found by computing accumulated probability into OFF states
# where transitions out of OFF states are zeroed, starting from first instance of ON state
# weighted by steady state distribution (i.e. solving a first passage time problem)
# Chemical master equation is solved  using DifferentialEquations.jl LSODA method

# Type of model is specified by TA matrix

# -`tin` end time of dwell time histogram
# -`G`: number of G states
# -`R`: number of R states

# """
# function ontimeCDF(tin::Vector, G::Int, R::Int, TA::AbstractMatrix, pss::Vector, method)
#     t = [0; tin]
#     SAinit = init_SA(pss, G - 1, R)
#     SA = time_evolve(t, TA, SAinit, method)
#     accumulate(SA, G - 1, R)  # accumulated prob into OFF states
# end

# """
# ontimeCDF(tin::Vector,G::Int,TA::AbstractMatrix,method)

# returns ON time histogram for G model
#     ON states can be any subset of G states

# """
# function ontimeCDF(tin::Vector, TA::AbstractMatrix, offstates, SAinit::Vector, method::Int=1)
#     t = [0; tin]
#     SA = time_evolve(t, TA, SAinit, method)
#     return sum(SA[:, offstates], dims=2)[:, 1]
# end

# ontimePDF(tin::Vector, TA::AbstractMatrix, offstates, SAinit::Vector, method::Int=1) = pdf_from_cdf(dwelltimeCDF(tin, TA, offstates, SAinit, method))
# """
# offtimeCDF(tin::Vector,r::Vector,G::Int,R::Int,TI::AbstractMatrix,pss::Vector,method)

# returns OFF dwell time distributions of GRS model
# using same algorithm as for ON time

# """
# function offtimeCDF(tin::Vector, r::Vector, G::Int, R::Int, TI::AbstractMatrix, pss::Vector, method::Int=1)
#     t = [0; tin]
#     nonzerosI = nonzero_rows(TI)  # only keep nonzero rows to reduce singularity of matrix
#     TI = TI[nonzerosI, nonzerosI]
#     SIinit = init_SI(pss, r, G - 1, R, nonzerosI)
#     SI = time_evolve(t, TI, SIinit, method)
#     accumulate(SI, G - 1, R, nonzerosI, length(pss)) # accumulated prob into ON states
# end
# """
# offtimeCDF(tin::Vector,r::Vector,G::Int,TI::AbstractMatrix,SIinit::Vector,method)

# returns OFF time histogram
# """
# function offtimeCDF(tin::Vector, TI::AbstractMatrix, onstates, SIinit::Vector, method::Int=1)
#     t = [0; tin]
#     SI = time_evolve(t, TI, SIinit, method)
#     return sum(SI[:, onstates], dims=2)[:, 1]
# end

# offtimePDF(tin::Vector, TI::AbstractMatrix, onstates, SIinit::Vector, method::Int=1) = pdf_from_cdf(dwelltimeCDF(tin, TI, onstates, SIinit))



# """
# offonPDF(T, TA, TI, t::Vector, r::Vector, G::Int, R::Int, method::Int=1)

# returns OFF and ON dwell time histograms for GR(S) model
# Takes difference of ON and OFF time CDF to produce PDF
# method = 1 uses DifferentialEquations.jl
# """
# function offonPDF(T, TA, TI, t::Vector, r::Vector, G::Int, R::Int, method::Int=1)
#     pss = normalized_nullspace(T)
#     SA = ontimeCDF(t, G, R, TA, pss, method)
#     SI = offtimeCDF(t, r, G, R, TI, pss, method)
#     return pdf_from_cdf(SI), pdf_from_cdf(SA)
# end

# """
#     offonPDF(T,TA, TI, t::Vector, r::Vector, G::Int, transitions::Tuple, onstates::Vector, method=1)

# returns OFF and ON dwell time histograms for G model
# """
# function offonPDF(TA, TI, t::Vector, r::Vector, G::Int, transitions::Tuple, onstates::Vector, method=1)
#     offstates = off_states(G, onstates)
#     SA = ontimeCDF(t, TA, offstates, init_SA(G, onstates, transitions), method)
#     if length(off) > 1
#         SI = offtimeCDF(t, TI, onstates, init_SI(G, onstates, transitions, r), method)
#     else
#         SI = offtimeCDF(t, TI, onstates, init_SI(G, onstates, transitions), method)
#     end
#     return pdf_from_cdf(SI), pdf_from_cdf(SA)
# end

# function gt_histograms(r, transitions, G, R, nhist, nalleles, bins, onstates, method=1, splicetype="")
#     ntransitions = length(transitions)
#     if R > 0
#         components = make_components(transitions, G, R, r, nhist + 2, splicetype, Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), collect(ntransitions+R+2:ntransitions+2*R+1), ntransitions + 2 * R + 2))
#     else
#         components = make_components(transitions, G, r, nhist + 2, Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), collect(ntransitions+R+2:ntransitions+2*R+1), ntransitions + 2 * R + 2), onstates)
#     end

#     # if splicetype == "offdecay"
#     #     r[end-1] *= survival_fraction(nu,eta,model.R)
#     # end
#     T = make_mat_T(components.tcomponents, r)
#     TA = make_mat_TA(components.tcomponents, r)
#     TI = make_mat_TI(components.tcomponents, r)
#     if R > 0
#         modelOFF, modelON = offonPDF(T, TA, TI, bins, r, G, R, method)
#     else
#         modelOFF, modelON = offonPDF(TA, TI, bins, G, transitions, onstates)
#     end
#     M = make_mat_M(components.mcomponents, r)
#     histF = steady_state(M, components.mcomponents.nT, nalleles, nhist)
#     return modelOFF, modelON, histF
# end


# # """
# #     offstates(G,onstates)

# #     returns OFF G states (complement set of G ON states)
# # """
# # offstates(G, onstates) = setdiff(collect(1:G), onstates)

# """
# init_SA(pss,n,nr)
# Initial condition for first passage time calculation of active (ON) state
# (Sum over all states with first R step occupied weighted by equilibrium state
# of all states with all R steps not occupied with reporters)
# """
# function init_SA(pss::Vector, n::Int, nr::Int)
#     l = length(pss)
#     base = findbase(l, n, nr)
#     SAinit = zeros(l)
#     if base == 3
#         for z = 1:2^(nr-1)
#             aSAinit = (n + 1) + (n + 1) * decimal(vcat(2, digits(z - 1, base=2, pad=nr - 1)))
#             apss = (n + 1) + (n + 1) * decimal(vcat(0, digits(z - 1, base=2, pad=nr - 1)))
#             SAinit[aSAinit] = pss[apss]
#         end
#     else
#         for z = 1:2^(nr-1)
#             aSAinit = (n + 1) + (n + 1) * decimal(vcat(1, zeros(Int, nr - 1)))
#             apss = n + 1
#             SAinit[aSAinit] = pss[apss]
#         end

#     end
#     SAinit / sum(SAinit)
# end
# function init_SA(G::Int, onstates::Vector, transitions::Tuple)
#     start = Int[]
#     for t in transitions
#         if t[1] ∉ onstates && t[2] ∈ onstates
#             push!(start, t[2])
#         end
#     end
#     Sinit = zeros(G)
#     Sinit[minimum(start)] = 1
#     return Sinit
# end


# """
# init_SI(pss,r,n,nr,nonzeros)
# Initial condition for first passage time calculation of inactive (OFF) state
# (Sum over all states with no R step introns occupied weighted by equilibrium state
# of all states with a single R step intron)
# """
# function init_SI(pss::Vector, r, n::Int, nr, nonzeros)
#     l = length(pss)
#     base = findbase(l, n, nr)
#     SIinit = zeros(l)
#     nu = get_nu(r, n, nr)
#     if base == 3
#         eta = get_eta(r, n, nr)
#         # Start of OFF state by ejection
#         for i in 1:n+1, z in 1:2^(nr-1)
#             exons = digits(z - 1, base=2, pad=nr - 1)
#             ainit = i + (n + 1) * decimal(vcat(exons, 0))
#             apss = i + (n + 1) * decimal(vcat(exons, 2))
#             SIinit[ainit] += pss[apss] * nu[nr+1]
#         end
#         # Start of OFF state by splicing
#         for i in 1:n+1, z in 1:2^(nr)
#             exons = digits(z - 1, base=2, pad=nr)
#             intronindex = findall(exons .== 1)
#             ainit = i + (n + 1) * decimal(exons)
#             for j in intronindex
#                 introns = copy(exons)
#                 introns[j] = 2
#                 apss = i + (n + 1) * decimal(introns)
#                 SIinit[ainit] += pss[apss] * eta[j]
#             end
#         end
#     else
#         # Start of OFF state by ejection
#         for i in 1:n+1
#             ainit = i
#             apss = i + (n + 1) * decimal(vcat(zeros(Int, nr - 1), 1), 2)
#             SIinit[ainit] += pss[apss] * nu[nr+1]
#         end
#     end
#     SIinit = SIinit[nonzeros]
#     SIinit / sum(SIinit)
# end


# function init_SI(G::Int, onstates, transitions)
#     start = Int[]
#     for t in transitions
#         if t[1] ∈ onstates && t[2] ∉ onstates
#             push!(start, t[2])
#         end
#     end
#     Sinit = zeros(G)
#     Sinit[minimum(start)] = 1
#     return Sinit
# end

# function init_SI(G::Int, onstates, transitions, r::Vector)
#     start = Int[]
#     tindex = Int[]
#     for (i, t) in enumerate(transitions)
#         if t[1] ∈ onstates && t[2] ∉ onstates
#             push!(start, t[2])
#             push!(tindex, i)
#         end
#     end
#     Sinit = zeros(G)
#     for (i, s) in enumerate(start)
#         Sinit[s] = r[tindex[i]]
#     end
#     Sinit / sum(Sinit)
# end



# """
# accumulate(SA::Matrix,n,nr)
# Sum over all probability vectors accumulated into OFF states
# """
# function accumulate(SA, n, nr)
#     l, p = size(SA)
#     base = findbase(p, n, nr)
#     SAj = zeros(size(SA, 1))
#     for i = 1:n+1, z = 1:base^nr
#         zdigits = digits(z - 1, base=base, pad=nr)
#         if ~any(zdigits .> base - 2)
#             a = i + (n + 1) * (z - 1)
#             SAj += SA[:, a]
#         end
#     end
#     return SAj
# end
# """
# accumulate(SI::Matrix,n,nr,nonzeros)
# Sum over all probability vectors accumulated into ON states
# """
# function accumulate(SI, n, nr, nonzeros, p)
#     # Sum over all probability vectors accumulated into ON states
#     l = size(SI)[1]
#     base = findbase(p, n, nr)
#     SIj = zeros(size(SI)[1])
#     for i = 1:n+1, z = 1:base^nr
#         zdigits = digits(z - 1, base=base, pad=nr)
#         if any(zdigits .> base - 2)
#             a = i + (n + 1) * (z - 1)
#             if a in nonzeros
#                 SIj += SI[:, findfirst(a .== nonzeros)]
#             end
#         end
#     end
#     return SIj
# end

# """
# initial_pmf(T,ejectrate,n,nr,nhist)

# Tensor product of nullspace of T and Poisson density
# with rate = mean ejection rate of mRNA
# Used for initial condition of iteration algorithm to
# compute steady state pmf of GRM Master Equation
# """
# function initial_pmf(T, ejectrate, n, nr, nhist)
#     t0 = normalized_nullspace(T)
#     ejectrate *= sum(t0[(n+1)*2^(nr-1)+1:end])
#     initial_pmf(t0, ejectrate, nhist)
# end
# function initial_pmf(T, ejectrate, n, nhist)
#     t0 = normalized_nullspace(T)
#     initial_pmf(t0, ejectrate, nhist)
# end
# function initial_pmf(t0, ejectrate, nhist)
#     d = Poisson(ejectrate)
#     nT = length(t0)
#     total = nhist + 2
#     P = Array{Float64,2}(undef, nT, total)
#     for m = 1:total
#         P[:, m] = t0 * pdf(d, m)
#     end
#     P
# end
### end of legacy code
