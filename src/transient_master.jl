# This file is part of StochasticGene.jl
#
# transient_master.jl
#
# Transient RNA master-equation solvers that sit beside the existing steady-state
# CME stack. These routines implement a closure/operator-splitting method for the
# case where mature RNA production depends on the finite configuration but does
# not change that configuration.

"""
    TransientMasterProblem(T, beta, decay)

Finite-modulator transient RNA problem for the split master equation.

`T` is the finite configuration generator with column sums zero,
`beta[a]` is the mature-RNA production rate in configuration `a`, and
`decay` is the mature-RNA decay rate per molecule. The joint distribution is
stored as `P[a, m+1]`.

This problem type represents the affine production/linear-decay case used by
[`transient_master_strang`](@ref). It is intentionally separate from the
existing steady-state CME matrix machinery.
"""
struct TransientMasterProblem{Tmat,Tbeta,Tdecay}
    T::Tmat
    beta::Tbeta
    decay::Tdecay
end

function TransientMasterProblem(T::AbstractMatrix, beta::AbstractVector, decay)
    size(T, 1) == size(T, 2) || throw(DimensionMismatch("T must be square"))
    length(beta) == size(T, 1) || throw(DimensionMismatch("length(beta) must match size(T, 1)"))
    return TransientMasterProblem{typeof(T),typeof(beta),typeof(decay)}(T, beta, decay)
end

transient_master_problem(T::AbstractMatrix, beta::AbstractVector, decay) =
    TransientMasterProblem(T, beta, decay)

"""
    transient_master_problem(components::MComponents, rates)

Build a [`TransientMasterProblem`](@ref) from existing `MComponents`.

This adapter uses the existing transition-matrix machinery to build the finite
configuration generator `T` and productive-event matrix `B`. It currently
supports the diagonal-production case where `B` has no off-diagonal entries,
so production increases mature RNA but leaves the finite configuration
unchanged. If `B` contains off-diagonal productive transitions, the exact
closed-form A-flow used here is not valid and this method throws an
`ArgumentError`.

Important: GR/GRSM models with `R > 0` usually have a last-R-step ejection
transition that both changes the finite R configuration and increments mature
mRNA. That marked transition is represented by off-diagonal entries in `B`.
It is therefore detected here and rejected rather than silently approximated
with the diagonal-production closure. Supporting that case requires a
generalized marked-transition A-flow or a fallback to the full joint transient
CME.
"""
function transient_master_problem(components::MComponents, rates::AbstractVector)
    T = make_mat(components.elementsT, rates, components.nT)
    B = make_mat(components.elementsB, rates, components.nT)
    rows, cols, vals = findnz(B)
    @inbounds for k in eachindex(vals)
        if rows[k] != cols[k] && !iszero(vals[k])
            throw(ArgumentError(
                "transient_master_problem currently supports only diagonal productive events; " *
                "found an off-diagonal B entry. Use the full CME transient solver for this model.",
            ))
        end
    end
    beta = zeros(eltype(vals), components.nT)
    @inbounds for k in eachindex(vals)
        beta[rows[k]] += vals[k]
    end
    return TransientMasterProblem(T, beta, _decay_from_U(components.U))
end

function _decay_from_U(U::AbstractMatrix)
    size(U, 1) >= 2 || return zero(eltype(U))
    return -U[2, 2]
end

function _log_factorials(::Type{T}, M::Int) where T
    out = zeros(T, M + 1)
    @inbounds for k in 1:M
        out[k + 1] = out[k] + log(T(k))
    end
    return out
end

function _binomial_thinning(Pa::AbstractVector, p, log_factorial::AbstractVector)
    M = length(Pa) - 1
    T = promote_type(eltype(Pa), typeof(p), eltype(log_factorial))
    out = zeros(T, M + 1)
    if iszero(p)
        out[1] = sum(Pa)
        return out
    elseif p == one(p)
        return T.(Pa)
    end
    log_p = log(p)
    log_q = log1p(-p)
    @inbounds for m_old in 0:M
        pa = Pa[m_old + 1]
        iszero(pa) && continue
        for k in 0:m_old
            log_w = log_factorial[m_old + 1]
            log_w -= log_factorial[k + 1] + log_factorial[m_old - k + 1]
            log_w += k * log_p + (m_old - k) * log_q
            out[k + 1] += pa * exp(log_w)
        end
    end
    return out
end

function _poisson_convolve(Pa::AbstractVector, lambda)
    M = length(Pa) - 1
    T = promote_type(eltype(Pa), typeof(lambda))
    iszero(lambda) && return T.(Pa)
    poi = zeros(T, M + 1)
    poi[1] = exp(-lambda)
    @inbounds for k in 1:M
        poi[k + 1] = poi[k] * lambda / k
    end
    out = zeros(T, M + 1)
    @inbounds for j in 0:M
        pj = Pa[j + 1]
        iszero(pj) && continue
        for k in 0:(M - j)
            out[j + k + 1] += pj * poi[k + 1]
        end
    end
    return out
end

"""
    transient_A_flow!(P, problem, dt, log_factorial)

Apply the exact mature-RNA birth/death flow for duration `dt`, with the finite
configuration frozen. Existing molecules survive with probability
`exp(-decay * dt)`, and new molecules are added by Poisson convolution with
mean `(beta[a] / decay) * (1 - exp(-decay * dt))` for each configuration `a`.
"""
function transient_A_flow!(
    P::AbstractMatrix,
    problem::TransientMasterProblem,
    dt,
    log_factorial::AbstractVector,
)
    p_survive = exp(-problem.decay * dt)
    @inbounds for a in axes(P, 1)
        Pa = _binomial_thinning(view(P, a, :), p_survive, log_factorial)
        lambda = (problem.beta[a] / problem.decay) * (1 - p_survive)
        Pa = _poisson_convolve(Pa, lambda)
        P[a, :] .= Pa
    end
    return P
end

"""
    transient_B_flow!(P, expT_dt)

Apply the finite-configuration flow `exp(T * dt)` to each RNA-count slice.
"""
function transient_B_flow!(P::AbstractMatrix, expT_dt::AbstractMatrix)
    P .= expT_dt * P
    return P
end

"""
    transient_master_strang(problem, P0, t_final, n_steps; return_joint=true)

Evolve a transient RNA distribution using Strang splitting:
`A(dt/2) B(dt) A(dt/2)`.

`P0` is an `nT × nhist` joint distribution with `P0[a, m+1]`. When
`return_joint=false`, returns the marginal mature-RNA histogram instead.
"""
function transient_master_strang(
    problem::TransientMasterProblem,
    P0::AbstractMatrix,
    t_final::Real,
    n_steps::Integer;
    return_joint::Bool=true,
)
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    size(P0, 1) == length(problem.beta) || throw(DimensionMismatch("P0 rows must match length(beta)"))
    M = size(P0, 2) - 1
    Twork = promote_type(eltype(P0), eltype(problem.T), eltype(problem.beta), typeof(problem.decay), typeof(t_final))
    P = Matrix{Twork}(P0)
    log_factorial = _log_factorials(Twork, M)
    dt = t_final / n_steps
    expT_dt = exp(Matrix{Twork}(problem.T) * dt)
    for _ in 1:n_steps
        transient_A_flow!(P, problem, dt / 2, log_factorial)
        transient_B_flow!(P, expT_dt)
        transient_A_flow!(P, problem, dt / 2, log_factorial)
    end
    return_joint && return P
    return vec(sum(P; dims=1))
end

"""
    transient_master_strang(problem, P0, times; steps_per_interval=1, return_joint=true)

Evolve through increasing `times` and return one snapshot per time. The initial
distribution `P0` is interpreted at time zero.
"""
function transient_master_strang(
    problem::TransientMasterProblem,
    P0::AbstractMatrix,
    times::AbstractVector;
    steps_per_interval::Integer=1,
    return_joint::Bool=true,
)
    steps_per_interval > 0 || throw(ArgumentError("steps_per_interval must be positive"))
    isempty(times) && return return_joint ? Matrix{eltype(P0)}[] : Vector{eltype(P0)}[]
    any(diff(times) .< 0) && throw(ArgumentError("times must be nondecreasing"))

    Twork = promote_type(eltype(P0), eltype(problem.T), eltype(problem.beta), typeof(problem.decay), eltype(times))
    P = Matrix{Twork}(P0)
    snapshots = return_joint ? Matrix{eltype(P)}[] : Vector{eltype(P)}[]
    t_prev = zero(eltype(times))
    for t in times
        dt_total = t - t_prev
        if dt_total != 0
            P = transient_master_strang(problem, P, dt_total, steps_per_interval; return_joint=true)
        end
        push!(snapshots, return_joint ? copy(P) : vec(sum(P; dims=1)))
        t_prev = t
    end
    return snapshots
end

"""
    transient_master_initial(config_distribution, nhist; m0=0)

Create an `nT × nhist` joint initial distribution concentrated at RNA count
`m0` with finite-configuration probabilities `config_distribution`.
"""
function transient_master_initial(config_distribution::AbstractVector, nhist::Integer; m0::Integer=0)
    0 <= m0 < nhist || throw(ArgumentError("m0 must satisfy 0 <= m0 < nhist"))
    P0 = zeros(eltype(config_distribution), length(config_distribution), nhist)
    P0[:, m0 + 1] .= config_distribution
    return P0
end

transient_master_marginal(P::AbstractMatrix) = vec(sum(P; dims=1))
