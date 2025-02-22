# This file is part of StochasticGene.jl   

### hmm.jl
### Fit discrete HMMs and continuous hidden Markov process models directly to observations (e.g. intensity traces)
###
### Notation in discrete HMM algorithms follows Rabier, 1989
###
### Functions for forward, backward, and Viterbi HMM algorihms
### For continuous processes, numerically solve forward Kolmogorov equation to obtain transition probability matrix
###

"""
    fkf!(du,u::Matrix, p, t)

in place update of du of Kolmogorov forward equation for DifferentialEquations.jl

"""
function fkf!(du, u::Matrix, p, t)
    du .= u * p
end

"""
    fkb!(du,u::Matrix, p, t)

in place update of du of Kolmogorov backward equation for DifferentialEquations.jl
"""
function fkb!(du, u::Matrix, p, t)
    du .= -u * p
end
"""
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_forward(Q, interval, method=Tsit5(), save=false)
    tspan = (0.0, interval)
    prob = ODEProblem(fkf!, Matrix(I, size(Q)), tspan, Q)
    solve(prob, method, save_everystep=save)[:, 2]
end
"""
kolmogorov_backward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_backward(Q, interval, method=Tsit5(), save=false)
    tspan = (0.0, interval)
    prob = ODEProblem(fkb!, Matrix(I, size(Q)), tspan, Q)
    solve(prob, method, save_everystep=save)[:, 2]
end


function prob_Gaussian(par, reporters)
    Normal(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2))
end
"""
    prob_Gaussian(par, reporters_per_state, N)

return Gaussian Distribution 
mean = background + number of reporters_per_state * reporter mean
variance = sum of variances of background and reporters_per_state

- `par`: 4 dimemsional vector of mean and std parameters
- `reporters_per_state`: number of reporters per HMM state
-`N`: number of HMM states
"""
function prob_Gaussian(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_Gaussian(par, reporters_per_state[i])
    end
    d
end

"""
    prob_GaussianMixture(par,reporters_per_state,N)

return Gaussian Mixture distribution with 4 Gaussian parameters and 1 weight parameter

"""
function prob_GaussianMixture(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture(par, reporters_per_state[i])
    end
    d
end
function prob_GaussianMixture(par, reporters)
    MixtureModel(Normal, [(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2)), (par[1], par[2])], [par[5], 1 - par[5]])
end

"""
    prob_GaussianMixture_6(par, reporters_per_state, N)

Gaussian Mixture distribution with 6 Gaussian parameters and 1 weight parameter
"""
function prob_GaussianMixture_6(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture_6(par, reporters_per_state[i])
    end
    d
end
function prob_GaussianMixture_6(par, reporters)
    MixtureModel(Normal, [(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2)), (par[5], par[6])], [par[7], 1 - par[7]])
end


"""
    prob_Gaussian_ind(par, reporters_per_state, N)

TBW
"""
function prob_Gaussian_ind(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_Gaussian_ind(par, reporters_per_state[i])
    end
    d
end
function prob_Gaussian_ind(par, reporters)
    if reporters > 0
        return Normal(reporters * par[3], sqrt(reporters) * par[4])
    else
        return Normal(par[1], par[2])
    end
end


"""
    prob_Gaussian_grid(par, reporters_per_state, Nstate, Ngrid, f::Function=kronecker_delta)

Generates a 3D array of Normal distributions based on the given parameters and reporters per state.

# Arguments
- `par`: Parameters for the Gaussian distribution.
- `reporters_per_state`: Number of reporters per state.
- `Nstate`: Number of states.
- `Ngrid`: Number of positions.
- `f::Function`: Function to use for Kronecker delta (default is `kronecker_delta`).

# Returns
- `Array{Distribution{Univariate,Continuous}}`: A 3D array of Normal distributions.
"""
function prob_Gaussian_grid(par, reporters_per_state, Nstate, Ngrid, f::Function=kronecker_delta)
    d = Array{Distribution{Univariate,Continuous}}(undef, Nstate, Ngrid, Ngrid)
    for j in 1:Nstate
        for k in 1:Ngrid
            for l in 1:Ngrid
                σ = sqrt(par[2]^2 + reporters_per_state[j] * par[4]^2 * f(k, l))
                d[j, k, l] = Normal(par[1] + reporters_per_state[j] * par[3] * f(k, l), σ)
            end
        end
    end
    return d
end



"""
    make_ap(r, interval, components::TComponents)

Return computed discrete HMM transition probability matrix a and equilibrium state probability p0
a is computed by numerically integrating Kolmogorov Forward equation for the underlying stochastic continuous time Markov process behind the GRSM model
p0 is left nullspace of transition rate matrix Q (right nullspace of Q')

Arguments:
- `r`: transition rates
- `interval`: time interval between intensity observations (frame interval)
- `components`: T matrix components

Qtr is the transpose of the Markov process transition rate matrix Q

"""
function make_ap(r, interval, components::TComponents, method=Tsit5())
    Qtr = make_mat_T(components, r) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

"""
    make_ap_coupled(r, couplingStrength, interval, components)


"""
function make_ap_coupled(r, couplingStrength, interval, components, method=Tsit5())
    Qtr = make_mat_TC(components, r, couplingStrength)
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

"""
    make_ap(r, interval, elementsT::Vector, N)

Arguments:
- `r`: transition rates
- `interval`: time interval between intensity observations (frame interval)
- `elementsT`: vector of T matrix elements
- `N`: number of HMM states

"""
function make_ap(r, interval, elementsT::Vector, N, method=Tsit5())
    Qtr = make_mat(elementsT, r, N) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end
"""
    make_logap(r, transitions, interval, G)

return log of a and p0
"""
function make_logap(r, interval, elementsT, N)
    a, p0 = make_ap(r, interval, elementsT, N)
    log.(max.(a, 0)), log.(max.(p0, 0))
end

"""
    make_a_grid(param, Ngrid)

TBW
"""
function make_a_grid(param, Ngrid)
    as = zeros(Ngrid, Ngrid)
    # d = zeros(Ngrid, Ngrid)
    for i in 1:Ngrid
        for j in 1:Ngrid
            as[i, j] = exp(-grid_distance(i, j, round(Int, sqrt(Ngrid)))^2 / (2 * param^2))
            # d[i, j] = grid_distance(i, j, div(Ngrid, 2))
        end
    end
    as ./ sum(as, dims=2)
end

"""
    set_d(noiseparams::Vector, reporters_per_state::Vector, probfn::Vector, N::Int)

TBW
"""
function set_d(noiseparams::Vector, reporters_per_state::Vector, probfn::Vector, N::Int)
    d = Vector[]
    for i in eachindex(noiseparams)
        push!(d, probfn[i](noiseparams[i], reporters_per_state[i], N))
    end
    return d
end

function set_d(noiseparams, reporter::Vector, N)
    ps = [r.per_state for r in reporter]
    pf = [r.probfn for r in reporter]
    set_d(noiseparams, ps, pf, N)
end
"""
    set_b(trace, params, reporters_per_state, probfn::Function=prob_Gaussian)

returns matrix b = P(Observation_i | State_j) for Gaussian distribution

-`trace`: Tx2 matrix of intensities.  Col 1 = time, Col 2 = intensity
-`N`: number of hidden states
-`T`: number of observations
"""
function set_b(trace, params, reporters_per_state, probfn::Function, N)
    # N = length(reporters_per_state)
    d = probfn(params, reporters_per_state, N)
    b = Matrix{Float64}(undef, N, length(trace))
    t = 1
    for obs in trace
        for j in 1:N
            b[j, t] = pdf(d[j], obs)
        end
        t += 1
    end
    return b
end

function set_b(trace, d, N)
    b = Matrix{Float64}(undef, N, length(trace))
    for (t, obs) in enumerate(trace)
        for j in 1:N
            b[j, t] = pdf(d[j], obs)
        end
    end
    return b
end

function set_b(trace::Matrix, d::Vector{Vector}, N)
    b = ones(N, size(trace, 1))
    t = 1
    for obs in eachrow(trace)
        for j in 1:N
            for i in eachindex(d)
                b[j, t] *= pdf(d[i][j], obs[i])
            end
        end
        t += 1
    end
    return b
end



"""
    set_b_coupled(trace, d, N)
    set_b_coupled(trace, params, reporter::Vector{HMMReporter}, N)
    set_b_coupled(trace, params, rep_per_state::Vector, probfn::Vector , N)


returns matrix b for coupled system
"""
function set_b_coupled(trace, d::Vector{Vector}, N)
    b = ones(N, size(trace, 1))
    t = 1
    for obs in eachrow(trace)
        for j in 1:N
            for i in eachindex(d)
                b[j, t] *= pdf(d[i][j], obs[i])
            end
        end
        t += 1
    end
    return b
end


function set_b_coupled(trace, params, reporters_per_state::Vector, probfn::Vector, N)
    d = set_d(params, reporters_per_state, probfn, N)
    set_b_coupled(trace, d, N)
end

function set_b_coupled(trace, params, reporter::Vector{HMMReporter}, N)
    d = Vector[]
    for i in eachindex(params)
        rep = reporter[i]
        push!(d, rep.probfn(params[i], rep.per_state, N))
    end
    set_b_coupled(trace, d, N)
end



"""
    set_b_grid(trace, params, reporters_per_state, probfn::Function, Nstate, Ngrid)

Calculates the probability matrix `b` for given trace data using a specified probability function.

# Arguments
- `trace`: The trace data.
- `params`: Parameters for the probability function.
- `reporters_per_state`: Number of reporters per state.
- `probfn::Function`: The probability function to use.
- `Nstate`: Number of states.
- `Ngrid`: Number of positions.

# Returns
- `Matrix`: The probability matrix `b`.
"""
function set_b_grid(trace, params, reporters_per_state, probfn::Function, Nstate, Ngrid)
    d = probfn(params, reporters_per_state, Nstate, Ngrid)
    set_b_grid(trace, d, Nstate, Ngrid)
end
function set_b_grid(trace, d, Nstate, Ngrid)
    b = ones(Nstate, Ngrid, size(trace, 2))
    t = 1
    for obs in eachcol(trace)
        for j in 1:Nstate
            for k in 1:Ngrid
                for l in 1:Ngrid
                    b[j, k, t] *= pdf(d[j, k, l], obs[l])
                end
            end
        end
        t += 1
    end
    return b
end

"""
    set_logb(trace, params, reporters_per_state, probfn=prob_Gaussian)

    returns log of matrix b
"""
function set_logb(trace, params, reporters_per_state, probfn::Function, N)
    # N = length(reporters_per_state)
    d = probfn(params, reporters_per_state, N)
    logb = Matrix{Float64}(undef, N, length(trace))
    t = 1
    for obs in trace
        for j in 1:N
            logb[j, t] = logpdf(d[j], obs)
        end
        t += 1
    end
    return logb
end

"""
    set_logb_coupled(trace, params, reporter, N)


"""
function set_logb_coupled(trace, params, reporter, N)
    d = Vector[]
    for i in eachindex(params)
        rep = reporter[i]
        push!(d, rep.probfn(params[i], rep.per_state, N))
    end
    logb = zeros(N, size(trace, 1))
    t = 1
    for obs in eachrow(trace)
        for j in 1:N
            for i in eachindex(d)
                logb[j, t] += logpdf(d[i][j], obs[i])
            end
        end
        t += 1
    end
    return logb
end

function set_b_background(obs::Float64, d::Vector{Distribution{Univariate,Continuous}})
    b = Array{Float64}(undef, size(d))
    for j in CartesianIndices(d)
        b[j] = pdf(d[j], obs)
    end
    return b
end

function set_b_background(obs::Float64, d::Vector{Vector}, k::Int, N)
    b = ones(N)
    for j in 1:N
        b[j] *= pdf(d[k][j], obs)
    end
    return b
end

"""
    forward_inner_operation!(α, a, b::Vector, i, j, t)

TBW
"""
function forward_inner_operation!(α, a, b::Vector, i, j, t)
    α[j, t] += α[i, t-1] * a[i, j] * b[j]
end

function forward_inner_operation!(α, a, b::Matrix, i, j, t)
    α[j, t] += α[i, t-1] * a[i, j] * b[j, t]
end

"""
forward(a, b, p0, N, T)

returns forward variable α, and scaling parameter array C using scaled forward algorithm
α[i,t] = P(O1,...,OT,qT=Si,λ)
Ct = Prod_t 1/∑_i α[i,t]

"""
function forward(a::Matrix, b, p0, N, T)
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    α[:, 1] = p0 .* b[:, 1]
    C[1] = 1 / sum(α[:, 1])
    α[:, 1] *= C[1]
    for t in 2:T
        for j in 1:N
            for i in 1:N
                forward_inner_operation!(α, a, b, i, j, t)
            end
        end
        C[t] = 1 / sum(α[:, t])
        α[:, t] *= C[t]
    end
    return α, C
end

function forward(a::Matrix, b, p0)
    N, T = size(b)
    forward(a, b, p0, N, T)
end

function forward(atuple::Tuple, b::Array, p0)
    a, a_grid = atuple
    Nstate, Ngrid, T = size(b)
    α = zeros(Nstate, Ngrid, T)
    C = Vector{Float64}(undef, T)
    α[:, :, 1] = p0 .* b[:, :, 1]
    C[1] = 1 / sum(α[:, :, 1])
    α[:, :, 1] *= C[1]
    for t in 2:T
        for l in 1:Ngrid, k in 1:Nstate
            for j in 1:Ngrid, i in 1:Nstate
                α[i, j, t] += α[k, l, t-1] * a[k, i] * a_grid[l, j] * b[i, j, t]
            end
        end
        C[t] = 1 / sum(α[:, :, t])
        α[:, :, t] *= C[t]
    end
    return α, C
end

"""
    forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)

TBW
"""
function forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    α = zeros(Nstate, Ngrid, T)
    C = Vector{Float64}(undef, T)
    α[:, :, 1] = p0 .* b[:, :, 1]
    C[1] = 1 / sum(α[:, :, 1])
    α[:, :, 1] *= C[1]
    for t in 2:T
        for l in 1:Ngrid, k in 1:Nstate
            for j in 1:Ngrid, i in 1:Nstate
                α[i, j, t] += α[k, l, t-1] * a[k, i] * a_grid[l, j] * b[i, j, t]
            end
        end
        C[t] = 1 / sum(α[:, :, t])
        α[:, :, t] *= C[t]
    end
    return α, C
end

"""
    forward_matrixmult(a, b, p0, N, T)

TBW
"""
function forward_matrixmult(a, b, p0, N, T)
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    α[:, 1] = p0 .* b[:, 1]
    C[1] = 1 / sum(α[:, 1])
    α[:, 1] *= C[1]
    for t in 2:T
        α[:, t] = α[:, t-1] .* a' * b[:, t]
        C[t] = 1 / sum(α[:, t])
        α[:, t] *= C[t]
    end
    return α, C
end

"""
forward_log(a, b, p0, N, T)
forward_log!(ϕ, ψ, loga, logb, logp0, N, T)

returns log α

(computations are numerically stable)

"""
function forward_log(loga, logb, logp0, N, T)
    ψ = zeros(N)
    ϕ = Matrix{Float64}(undef, N, T)
    forward_log!(ϕ, ψ, loga, logb, logp0, N, T)
    return ϕ
end
function forward_log!(ϕ, ψ, loga, logb, logp0, N, T)
    ϕ[:, 1] = logp0 + logb[:, 1]
    for t in 2:T
        for k in 1:N
            for j in 1:N
                ψ[j] = ϕ[j, t-1] + loga[j, k] + logb[k, t]
            end
            ϕ[k, t] = logsumexp(ψ)
        end
    end
end


"""
backward_scaled(a,b)

return backward variable β using scaled backward algorithm

β[i,T] = P(O[t+1]...O[t] | qT = Si,λ)

"""
function backward(a, b, C, N, T)
    β = ones(N, T)
    β[:, T] /= C[T]
    for t in T-1:-1:1
        for i in 1:N
            for j in 1:N
                β[i, t] += a[i, j] * b[j, t+1] * β[j, t+1]
            end
        end
        β[:, t] /= C[t]
    end
    return β
end

"""
backward_log(a, b, N, T)

return log β

"""
function backward_log(a, b, N, T)
    loga = log.(a)
    ψ = zeros(N)
    ϕ = Matrix{Float64}(undef, N, T)
    ϕ[:, T] = [0.0, 0.0]
    for t in T-1:-1:1
        for i in 1:N
            for j in 1:N
                ψ[j] = ϕ[j, t+1] + loga[i, j] + log.(b[j, t+1])
            end
            ϕ[i, t] = logsumexp(ψ)
        end
    end
    return ϕ
end

"""
    ll_background(a::Matrix, b::Vector, p0, N, T)

TBW
"""
function ll_background(obs::Float64, d::Vector{Distribution{Univariate,Continuous}}, a::Matrix, p0, nstates, nframes, weight)
    _, C = forward(a, set_b_background(obs, d), p0, nstates, nframes)
    weight * sum(log.(C))
end

function ll_background(obs::Vector, d::Vector{Distribution{Univariate,Continuous}}, a::Matrix, p0, nstates, nframes, weight)
    _, C = forward(a, set_b(obs, d, nstates), p0, nstates, nframes)
    weight * sum(log.(C))
end

function ll_background(obs::Vector, d::Vector{Vector}, a::Matrix, p0, nstates, nframes, weight)
    l = 0
    for i in eachindex(obs)
        b = set_b_background(obs[i], d, i, nstates)
        _, C = forward(a, b, p0, nstates, nframes)
        l += weight[i] * sum(log.(C))
    end
    l
end


p_off(a, p0, offstates, nframes) = sum(p0[offstates]' * a[offstates, offstates]^nframes)

"""
    ll_off(a, p0, offstates, poff, nframes)

L ∝ - log P(O | r) - p_inactive/p_active log (P(off | r))
"""
function ll_off(a, p0, offstates, poff, nframes)
    p = sum(p0[offstates]' * a[offstates, offstates]^nframes)
    -(1 - poff) * log(1 - p) - poff * log(p)
end
"""
    ll_off_coupled(a, p0, offstates, weight::Vector, nframes)

TBW
"""
function ll_off_coupled(a, p0, offstates, weight::Vector, nframes)
    l = 0
    for i in eachindex(weight)
        l += ll_off(a, p0, offstates[i], weight[i], nframes)
    end
    l
end

"""
    ll_hmm(a::Matrix, p0::Vector, d, traces, nT)

Calculates the log-likelihood of observed traces given the transition matrix `a`, 
the equilibrium probabilities `p0`, and the observation distributions `d`.

# Arguments
- `a::Matrix`: The transition probability matrix.
- `p0::Vector`: The initial state probabilities.
- `d`: The observation distributions for each state.
- `traces`: The observed traces for which the log-likelihood is calculated.
- `nT`: The number of time points.

# Returns
- `Float64`: The total log-likelihood of the observed traces.
- `Array{Float64}`: An array of log-likelihoods for each trace.
"""
function ll_hmm(a::Matrix, p0::Vector, d, traces, nT)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        _, C = forward(a, set_b(traces[i], d, nT), p0, nT, size(traces[i], 1))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm(r::Vector, a::Matrix, p0::Vector, n_noiseparams, reporters_per_state, probfn, traces, nT)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        b = set_b(traces[i], r[end-n_noiseparams+1:end, i], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, size(traces[i], 1))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm(noiseparams::Vector, a::Matrix, p0::Vector, reporter, traces, nT)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter, nT)
        b = set_b(traces[i], d, nT)
        _, C = forward(a, b, p0, nT, size(traces[i], 1))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm(r::Vector, interval::Float64, components::AbstractComponents, n_noiseparams::Int, reporters_per_state, probfn, traces, nT, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[:, i], interval, components, method)
        b = set_b(traces[i], r[end-n_noiseparams+1:end, i], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, size(traces[i], 1))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm(rin, reporter, interval, traces, nT)
    r, couplingStrength, noiseparams = rin
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap_coupled(r[:, i], couplingStrength[:, i], interval, components)
        d = set_d(noiseparams[i], reporter, nT)
        b = set_b(traces[i], d, nT)
        _, C = forward(a, b, p0, nT, size(traces[i], 1))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end


### Trait models
function ll_hmm_trait(r::Vector, nstates::Int, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)
    a, p0 = make_ap(r, interval, components, method)
    d = probfn(r[reporters.noiseparams], reporters.per_state, nstates)
    lb = trace[3] > 0.0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(a, p0, d, trace[1], nstates)
    ll + lb, logpredictions
end

function ll_hmm_trait(r::@NamedTuple{rshared, rindividual}, nstates::Int, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)
    rshared, rindividual = r
    a, p0 = make_ap(rshared[:, 1], interval, components, method[1])
    d = probfn(rshared[reporters.noiseparams, 1], reporters.per_state, nstates)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    if method[2]
        ll, logpredictions = ll_hmm(rindividual, a, p0, reporters.n, reporters.per_state, reporters.probfn, trace[1], nstates)
    else
        ll, logpredictions = ll_hmm(rindividual, interval, components, reporters.n, reporters.per_state, reporters.probfn, trace[1], nstates, method[1])
    end
    ll + lb, logpredictions
end

function ll_hmm_trait(r::@NamedTuple{r, couplingStrength, noiseparams}, nstates::Int, components::TCoupledComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method)
    rates, couplingStrength, noiseparams = r
    a, p0 = make_ap_coupled(rates, couplingStrength, interval, components, method[1])
    d = set_d(noiseparams, reporter, nT)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background([n[1] for n in noiseparams], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(a, p0, d, trace[1], nstates)
    ll + lb, logpredictions
end


function ll_hmm_trait(r::@NamedTuple{rshared, rindividual, couplingStrength, noiseparams}, nstates, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method)
    rshared, rindividual, couplingStrength, noiseparams = r
    a, p0 = make_ap_coupled(rshared, couplingStrength, interval, components, method)
    d = set_d(noiseparams.shared, reporter, nT)
    lb = any(trace[3] .> 0.0) ? length(trace[1]) * ll_background([n[1] for n in noiseparams], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    if method[2]
        ll, logpredictions = ll_hmm(noiseparams.individual, a, p0, reporters.n, reporters.per_state, reporters.probfn, trace[1], nstates)
    else
        ll, logpredictions = ll_hmm(rindividual, interval, components, reporters.n, reporters.per_state, reporters.probfn, trace[1], nstates, method[1])
    end
    ll + lb, logpredictions
end

function ll_hmm_trait(r::@NamedTuple{r::Vector{Float64}, noiseparams::Vector{Float64}, prid::Float64}, N::Tuple, components::TComponents, reporters, interval, trace)
    rates, noiseparams, pgrid = r
    Nstate, Ngrid = N
    a_grid = make_a_grid(pgrid, Ngrid)
    a, p0 = make_ap(rates, interval, components)
    d = probfn(noiseparams, reportersper_state, Nstate, Ngrid)
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = size(t, 2)
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm_trait(r::@NamedTuple{rshared::Vector{Float64}, rindividual::Vector{Float64}, noiseparams::Vector{Float64}, pgrid::Float64}, N::Tuple, components::TComponents, reporters, interval, trace)
    rshared, rindividual, noiseparams, pgrid = r
    Nstate, Ngrid = N
    a_grid = make_a_grid(pgrid, Ngrid)
    a, p0 = make_ap(rshared, interval, components)
    logpredictions = Array{Float64}(undef, length(trace[1]))
    for t in trace[1]
        d = probfn(rindividual[noiseparams, i], reporters_per_state, Nstate, Ngrid)
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, size(t, 2))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end


"""
    ll_hmm_trait(r::@NamedTuple{rates::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Float64}, 
                 nstates::Tuple, components::TCoupledComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method)

Compute log-likelihood for coupled + grid trait models.
"""
function ll_hmm_trait(r::@NamedTuple{rates::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Float64},
    nstates::Tuple, components::TCoupledComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method)
    Nstate, Ngrid = nstates
    a_grid = make_a_grid(r.pgrid, Ngrid)
    a, p0 = make_ap_coupled(r.rates, r.couplingStrength, interval, components)

    d = set_d(r.noiseparams, reporter, Nstate)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background([n[1] for n in r.noiseparams], d, a, p0, Nstate, trace[4], trace[3]) : 0.0

    logpredictions = Array{Float64}(undef, length(trace[1]))
    for (i, t) in enumerate(trace[1])
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, size(t, 2))
        @inbounds logpredictions[i] = sum(log.(C))
    end

    sum(logpredictions) + lb, logpredictions
end

"""
    ll_hmm_trait(r::@NamedTuple{rshared::Vector, rindividual::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Float64},
                 nstates::Tuple, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method)

Compute log-likelihood for coupled + hierarchical + grid trait models.
"""
function ll_hmm_trait(r::@NamedTuple{rshared::Vector, rindividual::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Float64},
    nstates::Tuple, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method)
    Nstate, Ngrid = nstates
    a_grid = make_a_grid(r.pgrid, Ngrid)
    a, p0 = make_ap_coupled(r.rshared, r.couplingStrength, interval, components)

    d = set_d(r.noiseparams.shared, reporter, Nstate)
    lb = any(trace[3] .> 0.0) ? length(trace[1]) * ll_background([n[1] for n in r.noiseparams], d, a, p0, Nstate, trace[4], trace[3]) : 0.0

    logpredictions = Array{Float64}(undef, length(trace[1]))
    for (i, t) in enumerate(trace[1])
        d = set_d(r.noiseparams.individual[i], reporter[i], Nstate)
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, size(t, 2))
        @inbounds logpredictions[i] = sum(log.(C))
    end

    sum(logpredictions) + lb, logpredictions
end

#####
"""
    ll_hmm(r, nstates, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)

return total loglikelihood of traces with reporter noise and loglikelihood of each trace
"""
function ll_hmm(r::Vector, nstates::Int, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)
    a, p0 = make_ap(r, interval, components)
    d = probfn(r[end-n_noiseparams+1:end], reporters_per_state, nstates)
    lb = trace[3] > 0.0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(a, p0, d, trace[1], nstates)
    ll + lb, logpredictions
end

"""
    ll_hmm_hierarchical(rshared, rindividual::Matrix, nT, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)

TBW
"""
function ll_hmm_hierarchical(rshared, rindividual::Matrix, nT, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)
    a, p0 = make_ap(rshared[:, 1], interval, components)
    d = probfn(rshared[end-n_noiseparams+1:end, 1], reporters_per_state, nT)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nT, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(rindividual, interval::Float64, components, n_noiseparams, reporters_per_state, probfn, trace[1], nT)
    return ll + lb, vcat(logpredictions, lhp)
end

"""
    ll_hmm_hierarchical_rateshared(rshared, r::Matrix, nT, components::TComponents, n_noiseparams, reporters_per_state, probfn, offstates, interval, trace)

TBW
"""
function ll_hmm_hierarchical_rateshared(rshared, rindividual::Matrix, nT, components::TComponents, n_noiseparams, reporters_per_state, probfn, interval, trace)
    a, p0 = make_ap(rshared[:, 1], interval, components)
    d = probfn(rshared[end-n_noiseparams+1:end, 1], reporters_per_state, nT)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nT, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(rindividual, a, p0, n_noiseparams, reporters_per_state, probfn, trace[1], nT)
    ll + lb, logpredictions
end

"""
    ll_hmm_coupled(r, couplingStrength, noiseparams, components, reporter, interval, traces)

TBW
"""
function ll_hmm_coupled(r, couplingStrength, noiseparams::Vector, components, reporter::Vector{HMMReporter}, interval, trace)
    nT = components.N
    a, p0 = make_ap_coupled(r, couplingStrength, interval, components)
    d = set_d(noiseparams, reporter, nT)
    lb = any(trace[3] .> 0.0) ? length(trace[1]) * ll_background([n[1] for n in noiseparams], d, a, p0, nT, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(a, p0, d, trace[1], nT)
    ll + lb, logpredictions
end

"""
    ll_hmm_grid(r, p, Nstate, Ngrid, components::StochasticGene.TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)

TBW
"""
function ll_hmm_grid(r, noiseparams, pgrid, Nstate, Ngrid, components::TComponents, reporters_per_state, probfn, interval, trace)
    a_grid = make_a_grid(pgrid, Ngrid)
    a, p0 = make_ap(r, interval, components)
    d = probfn(noiseparams, reporters_per_state, Nstate, Ngrid)
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = size(t, 2)
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm_grid_hierarchical(rshared, rindividual, pgrid, Nstate, Ngrid, components::TComponents, reporters_per_state, probfn, interval, trace)
    a_grid = make_a_grid(pgrid, Ngrid)
    a, p0 = make_ap(rshared, interval, components)
    logpredictions = Array{Float64}(undef, length(trace[1]))
    for t in trace[1]
        d = probfn(rindividual[end-n_noiseparams+1:end, i], reporters_per_state, Nstate, Ngrid)
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, size(t, 2))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end




"""
expected_transitions(α, a, b, β, N, T)

returns ξ and γ 
ξ[i,j,t] = P(q[t] = S[i], q[t+1] = S[j] | O, λ)
γ[i,t] = ∑_j ξ[i,j,t]
"""
function expected_transitions(α, a, b, β, N, T)
    ξ = Array{Float64}(undef, N, N, T - 1)
    γ = Array{Float64}(undef, N, T - 1)
    for t in 1:T-1
        for j = 1:N
            for i = 1:N
                ξ[i, j, t] = α[i, t] * a[i, j] * b[j, t+1] * β[j, t+1]
            end
        end
        S = sum(ξ[:, :, t])
        ξ[:, :, t] = S == 0.0 ? zeros(N, N) : ξ[:, :, t] / S
        γ[:, t] = sum(ξ[:, :, t], dims=2)
    end
    return ξ, γ
end

"""
    expected_transitions_log(logα, a, b, logβ, N, T)

TBW
"""
function expected_transitions_log(logα, a, b, logβ, N, T)
    ξ = Array{Float64}(undef, N, N, T - 1)
    γ = Array{Float64}(undef, N, T - 1)
    for t in 1:T-1
        for j = 1:N
            for i = 1:N
                ξ[i, j, t] = logα[i, t] + log(a[i, j]) + log(b[j, t+1]) + logβ[j, t+1]
            end
        end
        S = logsumexp(ξ[:, :, t])
        ξ[:, :, t] .-= S
        for i in 1:N
            γ[i, t] = logsumexp(ξ[i, :, t])
        end
    end
    return ξ, γ
end
"""
expected_a(a, b, p0, N, T)
expected_a(ξ, γ, N)

returns the expected probability matrix a
"""
function expected_a(a, b, p0, N, T)
    α, C = forward(a, b, p0, N, T)
    β = backward(a, b, C, N, T)
    ξ, γ = expected_transitions(α, a, b, β, N, T)
    expected_a(ξ, γ, N)
end
function expected_a(ξ, γ, N::Int)
    a = zeros(N, N)
    ξS = sum(ξ, dims=3)
    γS = sum(γ, dims=2)
    for i in 1:N, j in 1:N
        a[i, j] = ξS[i, j] / γS[i]
    end
    return a
end
function expected_a_log(a, b, p0, N, T)
    α = forward_log(a, b, p0, N, T)
    β = backward_log(a, b, N, T)
    ξ, γ = expected_transitions_log(α, a, b, β, N, T)
    expected_a_log(ξ, γ, N)
end

function expected_a_log(ξ, γ, N::Int)
    a = zeros(N, N)
    ξS = zeros(N, N)
    γS = zeros(N)
    for i in 1:N
        for j in 1:N
            ξS[i, j] = logsumexp(ξ[i, j, :])
        end
        γS[i] = logsumexp(γ[i, :])
    end
    for i in 1:N, j in 1:N
        a[i, j] = ξS[i, j] - γS[i]
    end
    return a
end

function expected_a_loop(a, b, p0, N, T)
    α = forward_loop(a, b, p0, N, T)
    β = backward_loop(a, b, N, T)
    ξ, γ = expected_transitions(α, a, b, β, N, T)
    expected_rate(ξ, γ, N)
end

"""
    viterbi_log(loga, logb, logp0, N, T)

returns maximum likelihood state path using Viterbi algorithm
"""
function viterbi_log(loga, logb, logp0, N, T)
    ϕ = similar(logb)
    ψ = similar(ϕ)
    q = Vector{Int}(undef, T)
    ϕ[:, 1] = logp0 + logb[:, 1]
    ψ[:, 1] .= 0
    for t in 2:T
        for j in 1:N
            m, ψ[j, t] = findmax(ϕ[:, t-1] + loga[:, j])
            ϕ[j, t] = m + logb[j, t]
        end
    end
    q[T] = argmax(ϕ[:, T])
    for t in T-1:-1:1
        q[t] = ψ[q[t+1], t+1]
    end
    return q
end
"""
    viterbi(a, b, p0, N, T)

returns maximum likelihood state path using Viterbi algorithm
"""
function viterbi(a, b, p0, N, T)
    loga = log.(max.(a, 0.0))
    logb = log.(max.(b, 0.0))
    logp0 = log.(max.(p0, 0.0))
    viterbi_log(loga, logb, logp0, N, T)
end

function viterbi_grid_log(loga, loga_grid, logb, logp0, Nstate, Ngrid, T)
    ϕ = similar(logb)
    ψ = similar(ϕ)
    q = Vector{Int}(undef, T)
    ϕ[:, :, 1] = logp0 .+ logb[:, :, 1]
    ψ[:, :, 1] .= 0
    for t in 2:T
        for j in 1:Ngrid, i in 1:Nstate
            m, ψ[i, j, t] = findmax(ϕ[:, :, t-1] .+ loga[:, i] .+ loga_grid[:, j])
            ϕ[i, j, t] = m + logb[i, j, t]
        end
        q[T] = argmax(ϕ[:, :, T])
        for t in T-1:-1:1
            q[t] = ψ[q[t+1], t+1]
        end
    end
    return q
end
function viterbi_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    loga = log.(max.(a, 0.0))
    loga_grid = log.(max.(a_grid, 0.0))
    logb = log.(max.(b, 0.0))
    logp0 = log.(max.(p0, 0.0))
    viterbi_grid_log(loga, loga_grid, logb, logp0, Nstate, Ngrid, T)
end


"""
    covariance_functions(rin, transitions, G::Tuple, R, S, insertstep, interval, probfn, coupling, lags::Vector)

TBW
"""
function covariance_functions(rin, transitions, G::Tuple, R, S, insertstep, interval, probfn, coupling, lags::Vector)
    components2 = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    components = TRGCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    sourceStates = [c.sourceState for c in components.modelcomponents]
    r, couplingStrength, noiseparams = prepare_rates(rin, sourceStates, transitions, G, R, S, insertstep, [4, 4])
    num_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
    mean_intensity = Vector[]
    max_intensity = Float64[]
    for i in eachindex(noiseparams)
        mi = mean.(probfn(noiseparams[i], num_per_state[i], components.N))
        mmax = max(maximum(mi), 1.0)
        push!(max_intensity, mmax)
        push!(mean_intensity, mi / mmax)
    end
    a, p0 = make_ap_coupled(r, couplingStrength, interval, components)
    m1 = mean_hmm(p0, mean_intensity[1])
    m2 = mean_hmm(p0, mean_intensity[2])
    cc12 = (crosscov_hmm(a, p0, mean_intensity[1], mean_intensity[2], lags) .- m1 .* m2) * max_intensity[1] * max_intensity[2]
    cc21 = (crosscov_hmm(a, p0, mean_intensity[2], mean_intensity[1], lags) .- m1 .* m2) * max_intensity[1] * max_intensity[2]
    ac1 = (crosscov_hmm(a, p0, mean_intensity[1], mean_intensity[1], lags) .- m1 .^ 2) * max_intensity[1]^2
    ac2 = (crosscov_hmm(a, p0, mean_intensity[2], mean_intensity[2], lags) .- m2 .^ 2) * max_intensity[2]^2
    m1, m2, cc12, cc21, ac1, ac2, vcat(reverse(cc21), cc12[2:end]), vcat(-reverse(lags), lags[2:end])
end

function autocov_hmm(r, transitions, G, R, S, insertstep, interval, probfn, lags::Vector)
    components = TComponents(transitions, G, R, S, insertstep, "")
    mean_intensity = mean.(probfn(r[end-3:end], num_reporters_per_state(G, R, S, insertstep), components.nT))
    a, p0 = make_ap(r, interval, components)
    crosscov_hmm(a, p0, mean_intensity, mean_intensity, lags) .- mean_hmm(p0, mean_intensity) .^ 2
end

function crosscov_hmm(a, p0, meanintensity1, meanintensity2, lags)
    cc = zeros(length(lags))
    m1 = meanintensity1 .* p0
    al = a^lags[1]
    as = a^(lags[2] - lags[1])
    for l in eachindex(lags)
        for i in eachindex(meanintensity1)
            for j in eachindex(meanintensity2)
                cc[l] += m1[i] * al[i, j] * meanintensity2[j]
            end
        end
        al *= as
    end
    cc
end

function mean_hmm(p0, meanintensity)
    sum(p0 .* meanintensity)
end


"""
    predicted_statepath(r::Vector, N::Int, elementsT, noiseparams, reporters_per_state, probfn, T::Int, interval)
    predicted_statepath(r, tcomponents, reporter, T, interval)

return predicted state path using Viterbi algorithm
"""
# function predicted_statepath(trace, interval, r::Vector, N::Int, elementsT, noiseparams, reporters_per_state, probfn)
#     loga, logp0 = make_logap(r, interval, elementsT, N)
#     logb = set_logb(trace, r[end-noiseparams+1:end], reporters_per_state, probfn)
#     viterbi_log(loga, logb, logp0, N, length(trace))
# end

# function predicted_statepath(trace, interval, r, tcomponents, reporter)
#     predicted_statepath(trace, interval, r, tcomponents.nT, tcomponents.elementsT, reporter.n, reporter.per_state, reporter.probfn)
# end

# function predicted_statepath(trace, interval, model::AbstractGmodel)
#     tcomponents = tcomponent(model)
#     predicted_statepath(trace, interval, model.rates, tcomponents.nT, tcomponents.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn)
# end

function predicted_state(r, N, components, reporter, interval, trace)
    a, p0 = make_ap(r, interval, components)
    b = set_logb(trace, r[end-reporter.n+1:end], reporter.per_state, reporter.probfn, N)
    viterbi(a, b, p0, N, length(trace))
end

function predicted_states(r::Vector, nT, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, traces)
    states = Vector{Int}[]
    observation_dist = Vector[]
    a, p0 = make_ap(r, interval, components)
    d = probfn(r[end-n_noiseparams+1:end], reporters_per_state, nT)
    for t in traces
        T = length(t)
        b = set_b(t, r[end-n_noiseparams+1:end], reporters_per_state, probfn, nT)
        spath = viterbi(a, b, p0, nT, T)
        push!(states, spath)
        # push!(observation_dist, [mean(d[s]) for s in spath])
        push!(observation_dist, [d[s] for s in spath])
    end
    states, observation_dist
end


function predicted_states(r::Matrix, nT, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, traces)
    states = Vector{Int}[]
    observation_dist = Vector[]
    a, p0 = make_ap(r[:, 1], interval, components)
    for (i, t) in enumerate(traces)
        T = length(t)
        b = set_b(t, r[end-n_noiseparams+1:end, i], reporters_per_state, probfn, nT)
        spath = viterbi(a, b, p0, nT, T)
        push!(states, spath)
        d = probfn(r[end-n_noiseparams+1:end, i], reporters_per_state, nT)
        push!(observation_dist, [d[s] for s in spath])
    end
    states, observation_dist
end


function predicted_states(rates, coupling, transitions, G::Tuple, R, S, insertstep, components, n_noise, reporters_per_state, probfn, interval, traces)
    sourceStates = coupling[3]
    r, couplingStrength, noiseparams = prepare_rates(rates, sourceStates, transitions, G, R, S, insertstep, n_noise)
    nT = components.N
    a, p0 = make_ap_coupled(r, couplingStrength, interval, components)
    states = Array[]
    d = []
    for i in eachindex(noiseparams)
        push!(d, probfn[i](noiseparams[i], reporters_per_state[i], nT))
    end
    for t in traces
        T = size(t, 1)
        b = set_b_coupled(t, noiseparams, reporters_per_state, probfn, nT)
        push!(states, viterbi(a, b, p0, nT, T))
    end
    units = Vector[]
    observation_dist = Vector[]
    for s in states
        push!(units, [unit_state(i, G, R, S, coupling[1]) for i in s])
        push!(observation_dist, [[d[i] for d in d] for i in s])
    end
    units, observation_dist
end

function predicted_states_grid(r::Vector, Nstates, Ngrid, components::TComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, traces)
    states = Vector{Int}[]
    observation_dist = Vector[]
    a, p0 = make_ap(r, interval, components)
    a_grid = make_a_grid(r[end], Ngrid)
    d = probfn(r[end-n_noiseparams:end-1], reporters_per_state, Nstates, Ngrid)
    for t in traces
        T = length(t)
        b = set_b_grid(t, d, Nstates, Ngrid)
        spath = viterbi_grid(a, a_grid, b, p0, Nstates, Ngrid, T)
        push!(states, spath)
        # push!(observation_dist, [mean(d[s]) for s in spath])
        push!(observation_dist, [d[s] for s in spath])
    end
    states, observation_dist
end

### NEW TRAIT-BASED IMPLEMENTATION (EXPERIMENTAL) ###
#=
The following code implements a new trait-based approach for handling 
different model combinations. This is currently separate from the main
implementation to allow for testing and validation.
=#

"""
    ll_hmm_trait(r::Vector, nstates::Int, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)

Compute log-likelihood for basic model without traits (trait-based implementation).
"""
function ll_hmm_trait(r::Vector, nstates::Int, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)
    a, p0 = make_ap(r, interval, components, method)
    d = probfn(r[reporters.noiseparams], reporters.per_state, nstates)
    lb = trace[3] > 0.0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(a, p0, d, trace[1], nstates)
    ll + lb, logpredictions
end

"""
    ll_hmm_trait(r::@NamedTuple{rates::Vector, couplingStrength::Vector, noiseparams::Vector}, nstates::Int, components::TCoupledComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method)

Compute log-likelihood for coupled trait models (trait-based implementation).
"""
function ll_hmm_trait(r::@NamedTuple{rates::Vector, couplingStrength::Vector, noiseparams::Vector}, nstates::Int, components::TCoupledComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method)
    a, p0 = make_ap_coupled(r.rates, r.couplingStrength, interval, components, method[1])
    d = set_d(r.noiseparams, reporter, nstates)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background([n[1] for n in r.noiseparams], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    ll, logpredictions = ll_hmm(a, p0, d, trace[1], nstates)
    ll + lb, logpredictions
end

"""
    ll_hmm_trait(r::@NamedTuple{rshared::Vector, rindividual::Vector}, nstates::Int, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)

Compute log-likelihood for hierarchical trait models (trait-based implementation).
"""
function ll_hmm_trait(r::@NamedTuple{rshared::Vector, rindividual::Vector}, nstates::Int, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)
    a, p0 = make_ap(r.rshared[:, 1], interval, components, method[1])
    d = probfn(r.rshared[reporters.noiseparams, 1], reporters.per_state, nstates)
    lb = trace[3] > 0 ? length(trace[1]) * ll_background(trace[2], d, a, p0, nstates, trace[4], trace[3]) : 0.0
    if method[2]
        ll, logpredictions = ll_hmm(r.rindividual, a, p0, reporters.n, reporters.per_state, reporters.probfn, trace[1], nstates)
    else
        ll, logpredictions = ll_hmm(r.rindividual, interval, components, reporters.n, reporters.per_state, reporters.probfn, trace[1], nstates, method[1])
    end
    ll + lb, logpredictions
end

"""
    ll_hmm_trait(r::@NamedTuple{rates::Vector, pgrid::Float64}, nstates::Tuple, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)

Compute log-likelihood for grid trait models (trait-based implementation).
"""
function ll_hmm_trait(r::@NamedTuple{rates::Vector, pgrid::Float64}, nstates::Tuple, components::TComponents, reporters::HMMReporter, interval::Float64, trace::Tuple, method)
    Nstate, Ngrid = nstates
    a_grid = make_a_grid(r.pgrid, Ngrid)
    a, p0 = make_ap(r.rates, interval, components)
    d = probfn(r.rates[end-reporters.noiseparams+1:end], reporters.per_state, Nstate, Ngrid)
    logpredictions = Array{Float64}(undef, length(trace[1]))
    for (i, t) in enumerate(trace[1])
        b = set_b_grid(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, size(t, 2))
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

### END NEW TRAIT-BASED IMPLEMENTATION ###
