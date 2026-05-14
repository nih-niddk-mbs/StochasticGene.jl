# This file is part of StochasticGene.jl
#
# transient_master.jl
#
# Transient RNA master-equation solvers that sit beside the existing steady-state
# CME stack. The split master equation handled here is
#
#   dP_m/dt = A P_m + B P_{m-1} + decay * ((m + 1) P_{m+1} - m P_m),
#
# where the finite configuration generator decomposes as `A + B` (non-productive
# and productive transitions) and `decay` is the mature-RNA degradation rate.
# Productive transitions may be diagonal (production does not change the finite
# configuration) or off-diagonal (e.g. GR/GRSM ejection that simultaneously
# changes the R configuration and emits a transcript).
#
# Usage:
#
#     using StochasticGene
#
#     transitions = ([1, 2], [2, 1])
#     G, R, nhist = 2, 0, 60
#     decay = 1.0
#     rates = [0.3, 0.5, 8.0]
#     components = StochasticGene.MComponents(transitions, G, R, nhist, decay, "")
#
#     # Diagonal-B adapter (rejects off-diagonal productive events).
#     problem = transient_master_problem(components, rates)
#     # General adapter (accepts off-diagonal B, e.g. R-step ejection).
#     # problem = transient_master_problem_general(components, rates)
#
#     P0 = transient_master_initial([0.625, 0.375], nhist)
#     P  = transient_master_closure(problem, P0, 10.0; ode_steps = 200)
#     marginal = transient_master_marginal(P)
#
# Solver selection guide (in increasing order of overhead):
#
#   transient_master_strang              -- Delta-t^2, diagonal-B only, cheapest
#   transient_master_strang_richardson   -- Delta-t^4, diagonal-B only
#   transient_master_strang_purebd       -- Delta-t^2, any B, both halves close
#   transient_master_closure             -- exact in-window (RK4 ODE error only)
#   transient_master_closure_taylor      -- order-R Taylor in time (tight tol)
#   transient_master_closure_exp_taylor  -- exp on linear part + order-R Taylor
#                                           on the time-varying coupling (stiff
#                                           hidden-state spectra)
#
# State-dependent decay (mu varies by configuration) is handled by
# `TransientMasterProblemStateMu` and the `_statemu` solver variants.

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

# ---------------------------------------------------------------------------
# Matrix-valued composition-multiplier closure and higher-order splittings.
#
# These methods generalize the diagonal-production Strang split of
# `transient_master_strang` to (a) productive transitions that also change the
# finite configuration (off-diagonal `B`, e.g. GR/GRSM ejection steps) and
# (b) higher-order accuracy than the Strang Delta-t^2 baseline. The vector
# generating function `Z(z, t) = sum_m P_m(t) z^m` for the split master equation
# `dP_m/dt = A P_m + B P_{m-1} + decay * ((m+1) P_{m+1} - m P_m)` satisfies the
# first-order matrix PDE `d_t Z = (A + z B) Z + decay (1 - z) d_z Z`, which
# admits the exact composition/multiplier form `Z(z, t) = K_t(z) Z_0(Phi_t(z))`,
# `Phi_t(z) = 1 - (1 - z) exp(-decay t)`. Expanding `K_t(z)` on a finite output
# window `{0, ..., M}` gives a closed lower-triangular matrix ODE that does not
# truncate the mRNA dynamics during propagation; the output cap only restricts
# the requested window.
# ---------------------------------------------------------------------------

"""
    TransientMasterProblemAB(A, B, decay)

Transient RNA problem with general (possibly off-diagonal) productive transitions.

`A` is the part of the finite-configuration generator that does not create
transcripts, `B` is the part that creates one transcript per transition, and
`A + B` is the full finite generator (column-sum zero, nonnegative off-diagonal).
`decay` is the mature-RNA decay rate. Off-diagonal entries `B[i, j]` with
`i != j` correspond to transitions that both change the finite configuration and
produce one transcript (e.g. the R-step ejection of GR/GRSM models).

Use [`transient_master_closure`](@ref) or [`transient_master_strang_purebd`](@ref)
to evolve this problem. The diagonal-production case can also be expressed as a
[`TransientMasterProblem`](@ref), which is the lower-overhead form used by the
existing splitting routines.
"""
struct TransientMasterProblemAB{TA,TB,Tdecay}
    A::TA
    B::TB
    decay::Tdecay
end

function TransientMasterProblemAB(A::AbstractMatrix, B::AbstractMatrix, decay)
    size(A) == size(B) || throw(DimensionMismatch("A and B must have the same size"))
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("A and B must be square"))
    return TransientMasterProblemAB{typeof(A),typeof(B),typeof(decay)}(A, B, decay)
end

"""
    transient_master_problem_general(components::MComponents, rates)

Build a [`TransientMasterProblemAB`](@ref) from existing `MComponents`, allowing
off-diagonal productive transitions. This is the general analogue of
[`transient_master_problem`](@ref), which only handles diagonal `B`.
"""
function transient_master_problem_general(components::MComponents, rates::AbstractVector)
    T = make_mat(components.elementsT, rates, components.nT)
    B = make_mat(components.elementsB, rates, components.nT)
    A = Matrix(T) .- Matrix(B)
    return TransientMasterProblemAB(A, Matrix(B), _decay_from_U(components.U))
end

_problem_AB(p::TransientMasterProblemAB) = (Matrix{Float64}(p.A), Matrix{Float64}(p.B))

function _problem_AB(p::TransientMasterProblem)
    Bmat = Matrix{Float64}(Diagonal(Float64.(p.beta)))
    Amat = Matrix{Float64}(p.T) .- Bmat
    return (Amat, Bmat)
end

# Build the multiplier RHS into preallocated 3D arrays.
# K, dK have shape (n_T, n_T, M+1) and store K[m] as the (:, :, m+1) slab.
function _multiplier_rhs_batched!(dK::Array{Float64,3},
                                  K::Array{Float64,3},
                                  t::Float64,
                                  t_final::Float64,
                                  A::Matrix{Float64},
                                  B::Matrix{Float64},
                                  decay::Float64,
                                  tmpC::Matrix{Float64},
                                  tmpD::Matrix{Float64})
    n_T = size(K, 1)
    M_plus_1 = size(K, 3)
    p_decay = exp(-decay * (t_final - t))
    q = 1.0 - p_decay
    @inbounds for j in 1:n_T, i in 1:n_T
        tmpC[i, j] = A[i, j] + q * B[i, j]
        tmpD[i, j] = p_decay * B[i, j]
    end
    K_flat = reshape(K, n_T, n_T * M_plus_1)
    dK_flat = reshape(dK, n_T, n_T * M_plus_1)
    mul!(dK_flat, tmpC, K_flat)
    if M_plus_1 >= 2
        K_src = reshape(view(K, :, :, 1:(M_plus_1 - 1)), n_T, n_T * (M_plus_1 - 1))
        dK_tgt = reshape(view(dK, :, :, 2:M_plus_1), n_T, n_T * (M_plus_1 - 1))
        mul!(dK_tgt, tmpD, K_src, 1.0, 1.0)
    end
    return dK
end

function _integrate_multiplier_batched(A::Matrix{Float64},
                                       B::Matrix{Float64},
                                       decay::Float64,
                                       M::Int,
                                       t_final::Float64,
                                       ode_steps::Int)
    n_T = size(A, 1)
    M_plus_1 = M + 1
    K = zeros(n_T, n_T, M_plus_1)
    @inbounds for i in 1:n_T
        K[i, i, 1] = 1.0
    end
    if t_final == 0.0
        return K
    end
    ode_steps > 0 || throw(ArgumentError("ode_steps must be positive"))
    h = t_final / ode_steps
    k1 = similar(K); k2 = similar(K); k3 = similar(K); k4 = similar(K)
    Ktmp = similar(K)
    tmpC = zeros(n_T, n_T); tmpD = zeros(n_T, n_T)
    t = 0.0
    @inbounds for _ in 1:ode_steps
        _multiplier_rhs_batched!(k1, K, t, t_final, A, B, decay, tmpC, tmpD)
        @. Ktmp = K + 0.5 * h * k1
        _multiplier_rhs_batched!(k2, Ktmp, t + 0.5 * h, t_final, A, B, decay, tmpC, tmpD)
        @. Ktmp = K + 0.5 * h * k2
        _multiplier_rhs_batched!(k3, Ktmp, t + 0.5 * h, t_final, A, B, decay, tmpC, tmpD)
        @. Ktmp = K + h * k3
        _multiplier_rhs_batched!(k4, Ktmp, t + h, t_final, A, B, decay, tmpC, tmpD)
        @. K = K + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += h
    end
    return K
end

# Compose initial distribution with the binomial-thinning kernel from mRNA decay
# over the full interval. Returns the n_T x (M_out+1) coefficient matrix.
function _compose_initial_after_death(P0::Matrix{Float64}, p_survive::Float64, M_out::Int)
    n_T, M_in_plus_1 = size(P0)
    M_in = M_in_plus_1 - 1
    C = zeros(n_T, M_out + 1)
    if p_survive <= 0.0
        C[:, 1] .= vec(sum(P0; dims=2))
        return C
    elseif p_survive >= 1.0
        for n in 0:min(M_in, M_out)
            C[:, n + 1] .= P0[:, n + 1]
        end
        return C
    end
    q = 1.0 - p_survive
    @inbounds for m_old in 0:M_in
        col = P0[:, m_old + 1]
        all(iszero, col) && continue
        w = q^m_old
        C[:, 1] .+= w .* col
        if m_old > 0
            ratio = p_survive / q
            for n in 0:(min(m_old, M_out) - 1)
                w *= ((m_old - n) / (n + 1)) * ratio
                C[:, n + 2] .+= w .* col
            end
        end
    end
    return C
end

# P[:, n+1] = sum_{j=0..n} K[:, :, j+1] * C[:, n-j+1]
function _apply_multiplier_batched(K::Array{Float64,3}, C::Matrix{Float64})
    n_T, M_plus_1 = size(C)
    M = M_plus_1 - 1
    P = zeros(n_T, M_plus_1)
    @inbounds for j in 0:M
        Kj = view(K, :, :, j + 1)
        C_block = view(C, :, 1:(M - j + 1))
        P_block = view(P, :, (j + 1):(M + 1))
        mul!(P_block, Kj, C_block, 1.0, 1.0)
    end
    return P
end

"""
    transient_master_closure(problem, P0, t_final; M_out, ode_steps, return_joint=true)

Evolve a transient RNA distribution by the matrix-valued composition/multiplier
closure: exact (up to ODE error) in-window solver with no operator-splitting
error. Works on either [`TransientMasterProblem`](@ref) or
[`TransientMasterProblemAB`](@ref); the latter is required when productive
transitions are off-diagonal.

The closure integrates a lower-triangular matrix ODE for the multipliers
`K[m](t)` on `m = 0:M_out`, composes the initial distribution with binomial
thinning by survival probability `exp(-decay * t_final)`, and applies the
multipliers by Cauchy convolution. The mRNA dynamics are not truncated during
propagation; truncation only enters through the requested output window and
through any tail mass omitted from `P0`.

`ode_steps` controls the RK4 step count for the multiplier ODE; the default
scales with `t_final`. The cost per RK4 step is `O(M_out * n_T^3)` via batched
GEMM, an order of magnitude faster than a naive per-`m` loop at small `n_T`.
"""
function transient_master_closure(
    problem::Union{TransientMasterProblem,TransientMasterProblemAB},
    P0::AbstractMatrix,
    t_final::Real;
    M_out::Int = size(P0, 2) - 1,
    ode_steps::Int = max(200, ceil(Int, 100 * Float64(t_final))),
    return_joint::Bool = true,
)
    A, B = _problem_AB(problem)
    size(P0, 1) == size(A, 1) || throw(DimensionMismatch("P0 row count must match hidden-state dimension"))
    M_out >= 0 || throw(ArgumentError("M_out must be nonnegative"))
    decay = Float64(problem.decay)
    tf = Float64(t_final)
    P0f = Matrix{Float64}(P0)
    K = _integrate_multiplier_batched(A, B, decay, M_out, tf, ode_steps)
    p_survive = exp(-decay * tf)
    C = _compose_initial_after_death(P0f, p_survive, M_out)
    P = _apply_multiplier_batched(K, C)
    return_joint && return P
    return vec(sum(P; dims=1))
end

"""
    transient_master_strang_richardson(problem, P0, t_final, n_steps_base; return_joint=true)

Richardson extrapolation of [`transient_master_strang`](@ref) at step counts
`n_steps_base` and `2 * n_steps_base`. The leading Strang `Delta-t^2` commutator
term cancels, leaving `Delta-t^4` error at roughly three times the cost of a
single base-step Strang run. This is the dissipative-semigroup-safe higher-order
upgrade: Yoshida-Forest 4th order would require a negative substep, which would
break the binomial-thinning kernel in [`transient_A_flow!`](@ref).
"""
function transient_master_strang_richardson(
    problem::TransientMasterProblem,
    P0::AbstractMatrix,
    t_final::Real,
    n_steps_base::Integer;
    return_joint::Bool = true,
)
    n_steps_base > 0 || throw(ArgumentError("n_steps_base must be positive"))
    P_K = transient_master_strang(problem, P0, t_final, n_steps_base; return_joint=true)
    P_2K = transient_master_strang(problem, P0, t_final, 2 * n_steps_base; return_joint=true)
    P = (4.0 .* P_2K .- P_K) ./ 3.0
    return_joint && return P
    return vec(sum(P; dims=1))
end

# L_A half of the pure-BD split: per-state binomial thinning by mRNA decay only,
# no immigration. Block-diagonal in the finite configuration; closed form.
function _pure_death_flow!(P::AbstractMatrix, decay, dt, log_factorial::AbstractVector)
    p_survive = exp(-decay * dt)
    @inbounds for a in axes(P, 1)
        Pa = _binomial_thinning(view(P, a, :), p_survive, log_factorial)
        P[a, :] .= Pa
    end
    return P
end

# L_B half of the pure-BD split: dot P_m = A P_m + B P_{m-1}, no death term.
# Lower-triangular in m: production from level M to level M+1 is the only
# leakage and is invisible to the in-window levels 0..M. Step is RK4 on the
# whole n_T(M+1)-dim system.
function _apply_LB_rk4!(P::AbstractMatrix, A::Matrix{Float64}, B::Matrix{Float64}, dt::Float64)
    rhs! = function(dP, P_in)
        mul!(dP, A, P_in)
        @views mul!(dP[:, 2:end], B, P_in[:, 1:end-1], 1.0, 1.0)
    end
    k1 = similar(P); rhs!(k1, P)
    Ptmp = P .+ (0.5 * dt) .* k1
    k2 = similar(P); rhs!(k2, Ptmp)
    Ptmp .= P .+ (0.5 * dt) .* k2
    k3 = similar(P); rhs!(k3, Ptmp)
    Ptmp .= P .+ dt .* k3
    k4 = similar(P); rhs!(k4, Ptmp)
    P .+= (dt / 6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
    return P
end

"""
    transient_master_strang_purebd(problem, P0, t_final, n_steps; return_joint=true)

Strang split with `L_A` = pure mRNA degradation (per-state binomial thinning,
closed form) and `L_B` = finite-configuration transitions plus production
(lower-triangular in mRNA index, RK4 on the truncated system). Unlike the
diagonal-`B` Strang of [`transient_master_strang`](@ref), production may couple
configurations arbitrarily, so this splitting handles the GR/GRSM ejection case.
Both halves close: `L_A` by tensor-product binomial thinning, `L_B` by
level-by-level RK4 with no upper-cap feedback because the death coupling is
already in `L_A`.
"""
function transient_master_strang_purebd(
    problem::Union{TransientMasterProblem,TransientMasterProblemAB},
    P0::AbstractMatrix,
    t_final::Real,
    n_steps::Integer;
    return_joint::Bool = true,
)
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    A, B = _problem_AB(problem)
    size(P0, 1) == size(A, 1) || throw(DimensionMismatch("P0 row count must match hidden-state dimension"))
    decay = Float64(problem.decay)
    P = Matrix{Float64}(P0)
    M = size(P, 2) - 1
    log_factorial = _log_factorials(Float64, M)
    dt = Float64(t_final) / n_steps
    for _ in 1:n_steps
        _pure_death_flow!(P, decay, dt / 2, log_factorial)
        _apply_LB_rk4!(P, A, B, dt)
        _pure_death_flow!(P, decay, dt / 2, log_factorial)
    end
    return_joint && return P
    return vec(sum(P; dims=1))
end

# ---------------------------------------------------------------------------
# Order-R Taylor integration of the multiplier ODE.
#
# The multiplier ODE is `dot K[m] = (A + q(t) B) K[m] + p(t) B K[m-1]` with
# `p(t) = exp(-decay (t_final - t))` and `q(t) = 1 - p(t)`. Expanding `K(t + s)`
# in `s` and applying Leibniz on the time-dependent operator gives
#
#   K[m]^(k+1) = (A + B) K[m]^(k)
#                + p(t) B sum_{j=0..k} C(k, j) decay^j (K[m-1]^(k-j) - K[m]^(k-j))
#
# (the `q`-dependent piece collapses because q + p = 1). For m = 0 the
# `K[m-1]` terms vanish. The step update is `K[m] <- sum_{k=0..R} h^k/k! K[m]^(k)`.
# Per-step cost: O(R^2 * M * n_T^2) for the Leibniz accumulator plus
# O(R * M * n_T^3) for the matrix products. Beats RK4 when target accuracy is
# tight enough to offset the per-step overhead with a far larger step size.
# ---------------------------------------------------------------------------

function _taylor_step!(K::Vector{Matrix{Float64}},
                       derivs::Vector{Vector{Matrix{Float64}}},
                       ApB::Matrix{Float64},
                       B::Matrix{Float64},
                       decay::Float64,
                       p_t::Float64,
                       h::Float64,
                       accum::Matrix{Float64},
                       R::Int)
    M_plus_1 = length(K)
    @inbounds for m in 1:M_plus_1
        copyto!(derivs[1][m], K[m])
    end

    @inbounds for k in 0:(R - 1)
        for m in 1:M_plus_1
            out = derivs[k + 2][m]
            mul!(out, ApB, derivs[k + 1][m])

            fill!(accum, 0.0)
            pow_mu = 1.0
            for j in 0:k
                c = binomial(k, j) * pow_mu
                Km_kj = derivs[k - j + 1][m]
                if m >= 2
                    Kmm1_kj = derivs[k - j + 1][m - 1]
                    @. accum += c * (Kmm1_kj - Km_kj)
                else
                    @. accum -= c * Km_kj
                end
                pow_mu *= decay
            end
            mul!(out, B, accum, p_t, 1.0)
        end
    end

    h_pow = 1.0
    fact_k = 1.0
    @inbounds for k in 1:R
        h_pow *= h
        fact_k *= k
        c = h_pow / fact_k
        for m in 1:M_plus_1
            @. K[m] += c * derivs[k + 1][m]
        end
    end
    return K
end

function _integrate_multiplier_taylor(A::Matrix{Float64},
                                      B::Matrix{Float64},
                                      decay::Float64,
                                      M::Int,
                                      t_final::Float64,
                                      n_steps::Int,
                                      R::Int)
    n_T = size(A, 1)
    K = [zeros(n_T, n_T) for _ in 0:M]
    @inbounds for i in 1:n_T
        K[1][i, i] = 1.0
    end
    if t_final == 0.0
        return K
    end
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))

    ApB = A .+ B
    derivs = [[zeros(n_T, n_T) for _ in 0:M] for _ in 0:R]
    accum = zeros(n_T, n_T)

    h = t_final / n_steps
    t = 0.0
    @inbounds for _ in 1:n_steps
        p_t = exp(-decay * (t_final - t))
        _taylor_step!(K, derivs, ApB, B, decay, p_t, h, accum, R)
        t += h
    end
    return K
end

function _apply_multiplier_vec(K::Vector{Matrix{Float64}}, C::Matrix{Float64})
    n_T, M_plus_1 = size(C)
    M = M_plus_1 - 1
    P = zeros(n_T, M_plus_1)
    @inbounds for n in 0:M
        for j in 0:n
            P[:, n + 1] .+= K[j + 1] * C[:, n - j + 1]
        end
    end
    return P
end

"""
    transient_master_closure_taylor(problem, P0, t_final;
                                    M_out, n_steps, R=8, return_joint=true)

Order-`R` Taylor integration of the multiplier ODE in place of the RK4 path
used by [`transient_master_closure`](@ref). Typically beats RK4 for target
accuracy below ~1e-9, because the order-`R` truncation lets `n_steps` shrink
while the per-step cost grows only as `R^2` in the Leibniz accumulator.
"""
function transient_master_closure_taylor(
    problem::Union{TransientMasterProblem,TransientMasterProblemAB},
    P0::AbstractMatrix,
    t_final::Real;
    M_out::Int = size(P0, 2) - 1,
    n_steps::Int = max(40, ceil(Int, 20 * Float64(t_final))),
    R::Int = 8,
    return_joint::Bool = true,
)
    A, B = _problem_AB(problem)
    size(P0, 1) == size(A, 1) || throw(DimensionMismatch("P0 row count must match hidden-state dimension"))
    M_out >= 0 || throw(ArgumentError("M_out must be nonnegative"))
    decay = Float64(problem.decay)
    tf = Float64(t_final)
    P0f = Matrix{Float64}(P0)
    K = _integrate_multiplier_taylor(A, B, decay, M_out, tf, n_steps, R)
    p_survive = exp(-decay * tf)
    C = _compose_initial_after_death(P0f, p_survive, M_out)
    P = _apply_multiplier_vec(K, C)
    return_joint && return P
    return vec(sum(P; dims=1))
end

# ---------------------------------------------------------------------------
# Exponential-Taylor integration of the multiplier ODE.
#
# Split as `dot K[m] = L K[m] + p(t) B (K[m-1] - K[m])` with `L := A + B`
# constant. L is column-stochastic, so `exp(L h)` is a contraction --
# matrix-exponential treatment removes the linear stiffness limit that bounds
# RK4 and pure Taylor. Variation of parameters gives
#
#   K[m](t+h) = phi_0(Lh) K[m](t)
#             + sum_{r=0..R-1} h^{r+1} phi_{r+1}(Lh) N[m]^(r)(t),
#
# with `N[m] = p(t) B (K[m-1] - K[m])`. phi-functions use the downward-stable
# recursion `phi_j(Z) = Z * phi_{j+1}(Z) + I/j!`, seeded by a truncated series.
# ---------------------------------------------------------------------------

function _compute_phi_matrices(L::Matrix{Float64}, h::Float64, R_max::Int)
    n_T = size(L, 1)
    Z = L .* h
    Id = Matrix{Float64}(I, n_T, n_T)
    p_top = R_max + 8

    fact_top = 1.0
    @inbounds for k in 1:p_top
        fact_top *= k
    end
    phi_p = Id ./ fact_top
    term = copy(Id)
    fact_inc = fact_top
    @inbounds for l in 1:80
        term = term * Z
        fact_inc *= (p_top + l)
        contrib = term ./ fact_inc
        phi_p .+= contrib
        maximum(abs, contrib) < 1e-18 && break
    end

    phi = Vector{Matrix{Float64}}(undef, R_max + 1)
    current = phi_p
    fact_j = fact_top
    j = p_top
    @inbounds while j > 0
        fact_j /= j
        j -= 1
        current = Z * current + Id ./ fact_j
        if j <= R_max
            phi[j + 1] = copy(current)
        end
    end
    return phi
end

function _exp_taylor_step!(K::Vector{Matrix{Float64}},
                           K_derivs::Vector{Vector{Matrix{Float64}}},
                           N_derivs::Vector{Vector{Matrix{Float64}}},
                           L::Matrix{Float64},
                           B::Matrix{Float64},
                           decay::Float64,
                           p_t::Float64,
                           h::Float64,
                           phi_funcs::Vector{Matrix{Float64}},
                           accum::Matrix{Float64},
                           new_K::Matrix{Float64},
                           R::Int)
    M_plus_1 = length(K)
    @inbounds for m in 1:M_plus_1
        copyto!(K_derivs[1][m], K[m])
    end

    @inbounds for k in 0:(R - 1)
        for m in 1:M_plus_1
            fill!(accum, 0.0)
            pow_mu = 1.0
            for j in 0:k
                c = binomial(k, j) * pow_mu
                Km_kj = K_derivs[k - j + 1][m]
                if m >= 2
                    Kmm1_kj = K_derivs[k - j + 1][m - 1]
                    @. accum += c * (Kmm1_kj - Km_kj)
                else
                    @. accum -= c * Km_kj
                end
                pow_mu *= decay
            end
            Nmk = N_derivs[k + 1][m]
            mul!(Nmk, B, accum, p_t, 0.0)
            if k + 1 < R
                Knext = K_derivs[k + 2][m]
                mul!(Knext, L, K_derivs[k + 1][m])
                @. Knext += Nmk
            end
        end
    end

    @inbounds for m in 1:M_plus_1
        mul!(new_K, phi_funcs[1], K[m])
        h_pow = h
        for r in 0:(R - 1)
            mul!(new_K, phi_funcs[r + 2], N_derivs[r + 1][m], h_pow, 1.0)
            h_pow *= h
        end
        copyto!(K[m], new_K)
    end
    return K
end

function _integrate_multiplier_exp_taylor(A::Matrix{Float64},
                                          B::Matrix{Float64},
                                          decay::Float64,
                                          M::Int,
                                          t_final::Float64,
                                          n_steps::Int,
                                          R::Int)
    n_T = size(A, 1)
    K = [zeros(n_T, n_T) for _ in 0:M]
    @inbounds for i in 1:n_T
        K[1][i, i] = 1.0
    end
    if t_final == 0.0
        return K
    end
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))

    L = A .+ B
    h = t_final / n_steps
    phi_funcs = _compute_phi_matrices(L, h, R)

    K_derivs = [[zeros(n_T, n_T) for _ in 0:M] for _ in 0:(R - 1)]
    N_derivs = [[zeros(n_T, n_T) for _ in 0:M] for _ in 0:(R - 1)]
    accum = zeros(n_T, n_T)
    new_K = zeros(n_T, n_T)

    t = 0.0
    @inbounds for _ in 1:n_steps
        p_t = exp(-decay * (t_final - t))
        _exp_taylor_step!(K, K_derivs, N_derivs, L, B, decay, p_t, h,
                          phi_funcs, accum, new_K, R)
        t += h
    end
    return K
end

"""
    transient_master_closure_exp_taylor(problem, P0, t_final;
                                        M_out, n_steps, R=8, return_joint=true)

Exponential-Taylor integration of the multiplier ODE: matrix exponential on
the constant linear part `L = A + B` and order-`R` Taylor on the time-varying
coupling `p(t) B (K[m-1] - K[m])`. Removes the stability cliff that pure Taylor
and RK4 face in stiff hidden-state regimes (large rate-matrix spectrum), at the
cost of one phi-function bundle per step size (amortized across all `n_steps`
on a uniform grid).
"""
function transient_master_closure_exp_taylor(
    problem::Union{TransientMasterProblem,TransientMasterProblemAB},
    P0::AbstractMatrix,
    t_final::Real;
    M_out::Int = size(P0, 2) - 1,
    n_steps::Int = max(20, ceil(Int, 10 * Float64(t_final))),
    R::Int = 8,
    return_joint::Bool = true,
)
    A, B = _problem_AB(problem)
    size(P0, 1) == size(A, 1) || throw(DimensionMismatch("P0 row count must match hidden-state dimension"))
    M_out >= 0 || throw(ArgumentError("M_out must be nonnegative"))
    decay = Float64(problem.decay)
    tf = Float64(t_final)
    P0f = Matrix{Float64}(P0)
    K = _integrate_multiplier_exp_taylor(A, B, decay, M_out, tf, n_steps, R)
    p_survive = exp(-decay * tf)
    C = _compose_initial_after_death(P0f, p_survive, M_out)
    P = _apply_multiplier_vec(K, C)
    return_joint && return P
    return vec(sum(P; dims=1))
end

# ---------------------------------------------------------------------------
# State-dependent decay (mu varies by configuration).
#
# When the mature-RNA decay rate depends on the finite configuration, the
# scalar-drift closure no longer applies directly. Two routes are exposed:
#  - `transient_master_strang_statewise`: per-state mRNA flow + global finite-
#    state matrix exponential. Cheap and exact within each half, with the usual
#    Strang Delta-t^2 commutator error from non-diagonal B.
#  - `transient_master_strang_shifted`: closure on a uniform `mu_0` plus a
#    dense matrix-exponential correction for the residual (mu - mu_0). The
#    correction has cost O((n_T (M+1))^3) per step, so it is only viable when
#    the joint dimension is moderate.
# ---------------------------------------------------------------------------

"""
    TransientMasterProblemStateMu(A, B, decay)

Transient RNA problem with state-dependent mature-RNA decay. `decay[a]` is the
decay rate in finite configuration `a`. The non-productive and productive
finite-state generators are `A` and `B` as in [`TransientMasterProblemAB`](@ref).
"""
struct TransientMasterProblemStateMu{TA,TB,Tdecay}
    A::TA
    B::TB
    decay::Tdecay
end

function TransientMasterProblemStateMu(A::AbstractMatrix, B::AbstractMatrix, decay::AbstractVector)
    size(A) == size(B) || throw(DimensionMismatch("A and B must have the same size"))
    size(A, 1) == size(A, 2) || throw(DimensionMismatch("A and B must be square"))
    length(decay) == size(A, 1) || throw(DimensionMismatch("length(decay) must match size(A, 1)"))
    all(>=(0), decay) || throw(ArgumentError("decay rates must be nonnegative"))
    return TransientMasterProblemStateMu{typeof(A),typeof(B),typeof(decay)}(A, B, decay)
end

function _statewise_mrna_flow!(P::AbstractMatrix,
                               beta::AbstractVector,
                               decay::AbstractVector,
                               s::Float64,
                               log_factorial::AbstractVector)
    n_T = size(P, 1)
    @inbounds for a in 1:n_T
        mu_a = Float64(decay[a])
        Pa = view(P, a, :)
        p_survive = exp(-mu_a * s)
        Pa_new = _binomial_thinning(Pa, p_survive, log_factorial)
        lam = mu_a > 0 ? (Float64(beta[a]) / mu_a) * (1.0 - p_survive) : Float64(beta[a]) * s
        Pa_new = _poisson_convolve(Pa_new, lam)
        P[a, :] .= Pa_new
    end
    return P
end

"""
    transient_master_strang_statewise(problem, P0, t_final, n_steps; return_joint=true)

Strang split for [`TransientMasterProblemStateMu`](@ref): per-state mRNA
birth/death (closed form, with state-specific decay) sandwiched around a
full-step `exp((A + B) * dt)` on the finite-state generator. Valid when `B` is
diagonal (production does not change configuration); the productive transitions
contribute via `beta[a] = B[a, a]`.
"""
function transient_master_strang_statewise(
    problem::TransientMasterProblemStateMu,
    P0::AbstractMatrix,
    t_final::Real,
    n_steps::Integer;
    return_joint::Bool = true,
)
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    A = Matrix{Float64}(problem.A)
    B = Matrix{Float64}(problem.B)
    # Statewise Strang absorbs production into the per-state mRNA flow, so
    # off-diagonal productive transitions would be silently dropped. Reject.
    @inbounds for j in axes(B, 2), i in axes(B, 1)
        if i != j && !iszero(B[i, j])
            throw(ArgumentError(
                "transient_master_strang_statewise requires diagonal B; use the " *
                "shifted/closure variants for off-diagonal productive transitions",
            ))
        end
    end
    beta = [B[a, a] for a in axes(B, 1)]
    decay = Float64.(problem.decay)
    P = Matrix{Float64}(P0)
    M = size(P, 2) - 1
    log_factorial = _log_factorials(Float64, M)
    dt = Float64(t_final) / n_steps
    T_hidden = A .+ B
    expT = exp(T_hidden .* dt)
    for _ in 1:n_steps
        _statewise_mrna_flow!(P, beta, decay, dt / 2, log_factorial)
        P .= expT * P
        _statewise_mrna_flow!(P, beta, decay, dt / 2, log_factorial)
    end
    return_joint && return P
    return vec(sum(P; dims=1))
end

# Dense joint generator for the state-mu problem on m = 0:M. Used by the
# shifted Strang variant for the residual half.
function _build_statemu_residual_generator(decay::AbstractVector, decay0::Float64, M::Int)
    n_T = length(decay)
    N = n_T * (M + 1)
    L = zeros(N, N)
    idx(a, m) = (a - 1) * (M + 1) + (m + 1)
    @inbounds for m in 0:M
        for a in 1:n_T
            d = Float64(decay[a]) - decay0
            L[idx(a, m), idx(a, m)] -= d * m
            if m < M
                L[idx(a, m), idx(a, m + 1)] += d * (m + 1)
            end
        end
    end
    return L
end

"""
    transient_master_strang_shifted(problem, P0, t_final, n_steps; mu0, return_joint=true)

Strang split for [`TransientMasterProblemStateMu`](@ref) using a uniform
reference decay `mu0`: a closure half-step at `mu0` is sandwiched around a
dense matrix-exponential step on the residual `(decay[a] - mu0) m P_a(m)`
correction. The residual operator is block-diagonal in finite configuration
(only acts on the mRNA index), so the matrix exponential costs `O((n_T (M+1))^3)`
once per step (amortizable across uniform grids). This is only viable for
moderate joint dimensions; for large joint state spaces, prefer the statewise
splitting.
"""
function transient_master_strang_shifted(
    problem::TransientMasterProblemStateMu,
    P0::AbstractMatrix,
    t_final::Real,
    n_steps::Integer;
    mu0::Real = minimum(Float64.(problem.decay)),
    return_joint::Bool = true,
)
    n_steps > 0 || throw(ArgumentError("n_steps must be positive"))
    mu0 > 0 || throw(ArgumentError("mu0 must be positive"))
    A = Matrix{Float64}(problem.A)
    B = Matrix{Float64}(problem.B)
    closure_problem = TransientMasterProblemAB(A, B, Float64(mu0))
    P = Matrix{Float64}(P0)
    n_T, M_plus_1 = size(P)
    M = M_plus_1 - 1
    dt = Float64(t_final) / n_steps
    L_residual = _build_statemu_residual_generator(problem.decay, Float64(mu0), M)
    E_residual = exp(L_residual .* dt)
    ode_steps_per_half = max(50, ceil(Int, 100 * dt))
    pvec = zeros(n_T * M_plus_1)
    idx(a, m) = (a - 1) * M_plus_1 + (m + 1)
    for _ in 1:n_steps
        P = transient_master_closure(closure_problem, P, dt / 2;
                                     M_out=M, ode_steps=ode_steps_per_half)
        @inbounds for a in 1:n_T, m in 0:M
            pvec[idx(a, m)] = P[a, m + 1]
        end
        pvec = E_residual * pvec
        @inbounds for a in 1:n_T, m in 0:M
            P[a, m + 1] = pvec[idx(a, m)]
        end
        P = transient_master_closure(closure_problem, P, dt / 2;
                                     M_out=M, ode_steps=ode_steps_per_half)
    end
    return_joint && return P
    return vec(sum(P; dims=1))
end