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
    ll_hmm(r, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, offstates, interval, trace)

return total loglikelihood of traces with reporter noise and loglikelihood of each trace
"""
function ll_hmm(r, nT, components::TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, offstates, interval, trace)
    a, p0 = make_ap(r, interval, components)
    lb = trace[3] > 0.0 ? length(trace[1]) * ll_background(a, p0, offstates, trace[3], trace[4]) : 0.0
    ll, lp = ll_hmm(r, nT, n_noiseparams, reporters_per_state, probfn, trace[1], a, p0)
    return ll + lb, lp
end

function ll_hmm(r, nT, n_noiseparams::Int, reporters_per_state, probfn, traces, a, p0)
    logpredictions = Array{Float64}(undef, 0)
    for t in traces
        T = length(t)
        b = set_b(t, r[end-n_noiseparams+1:end], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions), logpredictions
end

function ll_hmm_2(r, nT, components::TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, offstates, interval, trace)
    a, p0 = make_ap(r, interval, components)
    d = probfn(r[end-n_noiseparams+1:end], reporters_per_state, nT)
    lb = trace[3] > 0.0 ? length(trace[1]) * ll_background(a, p0, offstates, trace[3], trace[4]) : 0.0
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = length(t)
        b = set_b(t, d, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions) + lb, logpredictions
end

"""
    ll_hmm_coupled(r, couplingStrength, noiseparams, components, reporter, interval, traces)

TBW
"""
function ll_hmm_coupled_2(r, couplingStrength, noiseparams::Vector, components, reporter::Vector{HMMReporter}, interval, trace)
    nT = components.N
    a, p0 = make_ap_coupled(r, couplingStrength, interval, components)
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = size(t, 1)
        b = set_b_coupled(t, noiseparams, reporter, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    offstates = [r.offstates for r in reporter]
    lb = prod(trace[3]) > 0.0 ? length(trace[1]) * ll_background_coupled(a, p0, offstates, trace[3], trace[4]) : 0.0
    sum(logpredictions) + lb, logpredictions
end

function ll_hmm_coupled(r, couplingStrength, noiseparams::Vector, components, reporter::Vector{HMMReporter}, interval, trace)
    nT = components.N
    a, p0 = make_ap_coupled(r, couplingStrength, interval, components)
    ps = [r.per_state for r in reporter]
    pf = [r.probfn for r in reporter]
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = size(t, 1)
        b = set_b_coupled(t, noiseparams, ps, pf, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    offstates = [r.offstates for r in reporter]
    lb = prod(trace[3]) > 0.0 ? length(trace[1]) * ll_background_coupled(a, p0, offstates, trace[3], trace[4]) : 0.0
    sum(logpredictions) + lb, logpredictions
end

"""
    ll_background(a, p0, offstates, weight, n)

L ∝ - log P(O | r) - p_inactive/p_active log (P(off | r))
"""
function ll_background(a, p0, offstates, poff, nframes)
    p = sum(p0[offstates]' * a[offstates, offstates]^nframes)
    l = -(1 - poff) * log(1 - p) - poff * log(p)
    l
end

p_off(a, p0, offstates, nframes) = sum(p0[offstates]' * a[offstates, offstates]^nframes)

"""
    ll_background_coupled(a, p0, offstates, weight, n)

TBW
"""
function ll_background_coupled(a, p0, offstates, weight::Vector, nframes)
    l = 0
    for i in eachindex(weight)
        l += ll_background(a, p0, offstates[i], weight[i], nframes)
    end
    l
end


### Obsolete

# function ll_hmm(r, nT, elementsT::Vector, noiseparams, reporters_per_state, probfn, offstates, interval, trace)
#     a, p0 = make_ap(r, interval, elementsT, nT)
#     lb = trace[3] > 0.0 ? ll_background(a, p0, offstates, trace[3], trace[4]) : 0.0
#     ll, lp = ll_hmm(r, nT, noiseparams, reporters_per_state, probfn, trace[1], log.(max.(a, 0)), log.(max.(p0, 0)))
#     return ll + lb, lp
# end

"""
    ll_hmm_log(r, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, offstates, interval, trace)

"""
function ll_hmm_log(r, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, offstates, interval, trace)
    a, p0 = make_ap(r, interval, components)
    lb = trace[3] > 0.0 ? ll_background(a, p0, offstates, trace[3], trace[4]) : 0.0
    ll, lp = ll_hmm_log(r, nT, noiseparams, reporters_per_state, probfn, trace[1], log.(max.(a, 0)), log.(max.(p0, 0)))
    return ll + lb, lp
end
"""
    ll_hmm_log(r, nT, noiseparams::Int, reporters_per_state, probfn, traces, loga, logp0)

"""
function ll_hmm_log(r, nT, noiseparams::Int, reporters_per_state, probfn, traces, loga, logp0)
    logpredictions = Array{Float64}(undef, 0)
    for t in traces
        T = length(t)
        logb = set_logb(t, r[end-noiseparams+1:end], reporters_per_state, probfn, nT)
        l = forward_log(loga, logb, logp0, nT, T)
        push!(logpredictions, logsumexp(l[:, T]))
    end
    -sum(logpredictions), -logpredictions
end

"""
    ll_hmm_hierarchical(r::Matrix, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, interval, trace)

TBW
"""
function ll_hmm_hierarchical(r::Matrix, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, interval, trace)
    logpredictions = Array{Float64}(undef, 0)
    for (i, t) in enumerate(trace[1])
        T = length(t)
        a, p0 = make_ap(r[:, i], interval, components)
        b = set_b(t, r[end-noiseparams+1:end, i], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm_hierarchical_rateshared(r::Matrix, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, interval, trace)

TBW
"""
function ll_hmm_hierarchical_rateshared(r::Matrix, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, interval, trace)
    logpredictions = Array{Float64}(undef, 0)
    a, p0 = make_ap(r[:, 1], interval, components)
    for (i, t) in enumerate(trace[1])
        T = length(t)
        b = set_b(t, r[end-noiseparams+1:end, i], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm_hierarchical_rateshared_background(r::Matrix, nT, elementsT::Vector, noiseparams, reporters_per_state, probfn, offstates, interval, trace)

TBW
"""
function ll_hmm_hierarchical_rateshared_background(r::Matrix, nT, elementsT::Vector, noiseparams, reporters_per_state, probfn, offstates, interval, trace)
    logpredictions = Array{Float64}(undef, 0)
    a, p0 = make_ap(r[:, 1], interval, elementsT, nT)
    lb = trace[3] > 0 ? ll_background(a, p0, offstates, trace[3], trace[4]) : 0.0
    for (i, t) in enumerate(trace[1])
        T = length(t)
        b = set_b(t, r[end-noiseparams+1:end, i], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions) + lb, logpredictions
end

"""
    ll_hmm_hierarchical_rateshared_background(r::Matrix, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, offstates, interval, trace)

TBW
"""
function ll_hmm_hierarchical_rateshared_background(r::Matrix, nT, components::TRGComponents, noiseparams, reporters_per_state, probfn, offstates, interval, trace)
    logpredictions = Array{Float64}(undef, 0)
    a, p0 = make_ap(r[:, 1], interval, components)
    lb = trace[3] > 0 ? ll_background(a, p0, offstates, trace[3], trace[4]) : 0.0
    for (i, t) in enumerate(trace[1])
        T = length(t)
        b = set_b(t, r[end-noiseparams+1:end, i], reporters_per_state, probfn, nT)
        _, C = forward(a, b, p0, nT, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions) + lb, logpredictions
end


"""
    ll_hmm_grid(r, p, Nstate, Ngrid, components::StochasticGene.TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)

TBW
"""
function ll_hmm_grid(r, p, Nstate, Ngrid, components::StochasticGene.TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)
    a_grid = make_a_grid(p, Ngrid)
    a, p0 = StochasticGene.make_ap(r, interval, components)
    d = probfn(r[end-n_noiseparams+1:end], reporters_per_state, Nstate, Ngrid)
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = length(t)
        b = set_b_grid_v(t, d, Nstate, Ngrid)
        _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
        push!(logpredictions, sum(log.(C)))
    end
    sum(logpredictions), logpredictions
end

"""
    make_ap(r, interval, components::TRGComponents)

Return computed discrete HMM transition probability matrix a and equilibrium state probability p0
a is computed by numerically integrating Kolmogorov Forward equation for the underlying stochastic continuous time Markov process behind the GRSM model
p0 is left nullspace of transition rate matrix Q (right nullspace of Q')

Arguments:
- `r`: transition rates
- `interval`: time interval between intensity observations (frame interval)
- `components`: TRG matrix components

Qtr is the transpose of the Markov process transition rate matrix Q

"""
function make_ap(r, interval, components::TRGComponents)
    Qtr = make_mat_TRG(components, r) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(Qtr', interval), normalized_nullspace(Qtr)
end

"""
    make_ap_coupled(r, couplingStrength, interval, components)


"""
function make_ap_coupled(r, couplingStrength, interval, components)
    Qtr = make_mat_TC(components, r, couplingStrength)
    kolmogorov_forward(Qtr', interval), normalized_nullspace(Qtr)
end

"""
    make_ap(r, interval, elementsT::Vector, N)

Arguments:
- `r`: transition rates
- `interval`: time interval between intensity observations (frame interval)
- `elementsT`: vector of T matrix elements
- `N`: number of HMM states

"""
function make_ap(r, interval, elementsT::Vector, N)
    Qtr = make_mat(elementsT, r, N) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(Qtr', interval), normalized_nullspace(Qtr)
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
    d = zeros(Ngrid, Ngrid)
    for i in 1:Ngrid
        for j in 1:Ngrid
            as[i, j] = exp(-grid_distance(i, j, div(Ngrid, 2))^2 / (2 * param^2))
            d[i, j] = grid_distance(i, j, div(Ngrid, 2))
        end
    end
    as ./ sum(as, dims=2)
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
    set_b_coupled(trace, params, reporter::Vector{HMMReporter}, N)
    set_b_coupled(trace, params, rep_per_state::Vector, probfn::Vector , N)
    set_b_coupled(trace, d, N)

returns matrix b for coupled system
"""

function set_b_coupled(trace, params, rep_per_state::Vector, probfn::Vector, N)
    d = Vector[]
    for i in eachindex(params)
        push!(d, probfn[i](params[i], rep_per_state[i], N))
    end
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

function set_b_coupled(trace, d, N)
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
function set_b_grid_v(trace, params, reporters_per_state, probfn::Function, Nstate, Ngrid)
    d = probfn(params, reporters_per_state, Nstate, Ngrid)
    set_b_grid_v(trace, d, Nstate, Ngrid)
end
function set_b_grid_v(trace, d, Nstate, Ngrid)
    b = ones(Nstate, Ngrid, length(trace))
    t = 1
    for obs in eachcol(trace)
        for j in 1:Nstate
            for k in 1:Ngrid
                for l in 1:Ngrid
                    b[j, k, t] *= StochasticGene.pdf(d[j, k, l], obs[l])
                end
            end
        end
        t += 1
    end
    return b
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
function prob_Gaussian(par, reporters)
    Normal(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2))
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
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_forward(Q, interval, save=false, method=Tsit5())
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
function kolmogorov_backward(Q, interval, save=false, method=Tsit5())
    tspan = (0.0, interval)
    prob = ODEProblem(fkb!, Matrix(I, size(Q)), tspan, Q)
    solve(prob, method, save_everystep=save)[:, 2]
end
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
forward(a, b, p0, N, T)

returns forward variable α, and scaling parameter array C using scaled forward algorithm
α[i,t] = P(O1,...,OT,qT=Si,λ)
Ct = Prod_t 1/∑_i α[i,t]

"""
function forward(a, b, p0, N, T)
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    α[:, 1] = p0 .* b[:, 1]
    C[1] = 1 / sum(α[:, 1])
    α[:, 1] *= C[1]
    for t in 2:T
        for j in 1:N
            for i in 1:N
                α[j, t] += α[i, t-1] * a[i, j] * b[j, t]
            end
        end
        C[t] = 1 / sum(α[:, t])
        α[:, t] *= C[t]
    end
    return α, C
end

function forward_mm(a, b, p0, N, T)
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
forward_loop(a, b, p0, N, T)

return α using unscaled forward algorithm
(numerically unstable for large T)

"""
function forward_loop(a, b, p0, N, T)
    α = zeros(N, T)
    α[:, 1] = p0 .* b[:, 1]
    for t in 2:T
        for j in 1:N
            for i in 1:N
                α[j, t] += α[i, t-1] * a[i, j] * b[j, t]
            end
        end
    end
    return α
end

function forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    α = zeros(Nstate, Ngrid, T)
    C = Vector{Float64}(undef, T)
    α[:, :, 1] = p0 .* b[:, :, 1]
    C[1] = 1 / sum(α[:, :, 1])
    α[:, :, 1] *= C[1]
    for t in 2:T
        for l in 1:Ngrid
            for k in 1:Nstate
                for j in 1:Ngrid
                    for i in 1:Nstate
                        α[i, j, t] += α[k, l, t-1] * a[k, i] * a_grid[l, j] * b[i, j, t]
                    end
                end
            end
        end
        C[t] = 1 / sum(α[:, :, t])
        α[:, :, t] *= C[t]
    end
    return α, C
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
backward_loop(a, b, N, T)

returns β using unscaled backward algorithm
(numerically unstable for large T)
"""
function backward_loop(a, b, N, T)
    β = zeros(N, T)
    β[:, T] = [1.0, 1.0]
    for t in T-1:-1:1
        for i in 1:N
            for j in 1:N
                β[i, t] += a[i, j] * b[j, t+1] * β[j, t+1]
            end
        end
    end
    return β
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
    viterbi(loga, logb, logp0, N, T)

returns maximum likelihood state path using Viterbi algorithm
"""
function viterbi(loga, logb, logp0, N, T)
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
    viterbi_exp(a, b, p0, N, T)

returns maximum likelihood state path using Viterbi algorithm
"""
function viterbi_exp(a, b, p0, N, T)
    loga = log.(max.(a, 0.0))
    logb = log.(max.(b, 0.0))
    logp0 = log.(max.(p0, 0.0))
    viterbi(loga, logb, logp0, N, T)
end

function covariance_functions(rin, transitions, G::Tuple, R, S, insertstep, interval, probfn, coupling, lags::Vector)
    components = make_components_Tcoupled(coupling, transitions, G, R, S, insertstep, "")
    r, couplingStrength, noiseparams = prepare_rates(rin, coupling[2], transitions, G, R, S, insertstep, [4,4])
    mean_intensity = Vector[]
    for i in eachindex(params)
        push!(mean_intensity, mean.(probfn(noiseparams[i], num_reporters_per_state(G, R, S, insertstep, coupling[1]), components.N)))
    end
    a, p0 = make_ap_coupled(r, couplingStrength, interval, components)
    m1 = mean_hmm(p0, mean_intensity[1])
    m2 = mean_hmm(p0, mean_intensity[2])
    cc12 = crosscov_hmm(a, p0, mean_intensity[1], mean_intensity[2], lags) - m1 .* m2
    cc21 = crosscov_hmm(a, p0, mean_intensity[2], mean_intensity[1], lags) - m1 .* m2
    ac1 = crosscov_hmm(a, p0, mean_intensity[1], mean_intensity[1], lags) - m1.^2
    ac2 = crosscov_hmm(a, p0, mean_intensity[2], mean_intensity[2], lags) - m2.^2
end

function autocov_hmm(r, transitions, G, R, S, insertstep, interval, probfn, lags::Vector)
    components = make_components_TRG(transitions, G, R, S, insertstep, "")
    mean_intensity = mean.(probfn(r[end-3:end], num_reporters_per_state(G, R, S, insertstep), components.nT))
    a, p0 = make_ap(r, interval, components)
    crosscov_hmm(a, p0, mean_intensity, mean_intensity, lags) .- mean_hmm(p0, mean_intensity).^2
end

function crosscov_hmm(a, p0, meanintensity1, meanintensity2, lags)
    cc = zeros(length(lags))
    for lag in lags
        al = a^lag
        for i in eachindex(meanintensity1)
            for j in eachindex(meanintensity2)
                cc[lag] += meanintensity1[i] * p0[i] * al[i, j] * meanintensity2[j]
            end
        end
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
#     viterbi(loga, logb, logp0, N, length(trace))
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
    viterbi_exp(a, b, p0, N, length(trace))
end

function predicted_states(r::Vector, nT, components::TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, traces)
    states = Vector{Int}[]
    observation_dist = Vector[]
    a, p0 = make_ap(r, interval, components)
    d = probfn(r[end-n_noiseparams+1:end], reporters_per_state, nT)
    for t in traces
        T = length(t)
        b = set_b(t, r[end-n_noiseparams+1:end], reporters_per_state, probfn, nT)
        spath = viterbi_exp(a, b, p0, nT, T)
        push!(states, spath)
        # push!(observation_dist, [mean(d[s]) for s in spath])
        push!(observation_dist, [d[s] for s in spath])
    end
    states, observation_dist
end


function predicted_states(r::Matrix, nT, components::TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, traces)
    states = Vector{Int}[]
    observation_dist = Vector[]
    a, p0 = make_ap(r[:, 1], interval, components)
    for (i, t) in enumerate(traces)
        T = length(t)
        b = set_b(t, r[end-n_noiseparams+1:end, i], reporters_per_state, probfn, nT)
        spath = viterbi_exp(a, b, p0, nT, T)
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
        push!(states, viterbi_exp(a, b, p0, nT, T))
    end
    units = Vector[]
    observation_dist = Vector[]
    for s in states
        push!(units, [unit_state(i, G, R, S, coupling[1]) for i in s])
        push!(observation_dist, [[d[i] for d in d] for i in s])
    end
    units, observation_dist
end


# """
#     predicted_trace(statepath, noise_dist)
#     predicted_trace(statepath, r, reporter, nstates)

# return predicted trace from state path
# """
# function predicted_trace(statepath, noise_dist)
#     [mean(noise_dist[state]) for state in statepath]
# end

# function predicted_trace(statepath, r, reporter, nstates)
#     d = model.reporter.probfn(r[end-model.reporter.n+1:end], reporter.per_state, nstates)
#     predicted_trace(statepath, d)
# end



# """
#     predicted_trace_state(trace, interval, r::Vector, tcomponents, reporter, noise_dist)

# return predicted trace and state path
# """
# function predicted_trace_state(trace, interval, r::Vector, tcomponents, reporter, noise_dist)
#     t = predicted_state(r, tcomponents.nT, tcomponents, reporter, interval, trace)
#     tp = predicted_trace(t, noise_dist)
#     return tp, t
# end

# function predicted_trace_state(trace, interval, model)
#     d = model.reporter.probfn(model.rates[end-model.reporter.n+1:end], model.reporter.per_state, tcomponent(model).nT)
#     predicted_trace_state(trace, interval, model.rates, model.tcomponents, model.reporter, d)
# end

# """
#     predicted_states(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model::AbstractGmodel)

# return vector of predicted state vectors
# """
# function predicted_states(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model::AbstractGmodel)
#     ts = Vector{Int}[]
#     for t in data.trace[1]
#         push!(ts, predicted_statepath(t, data.interval, model))
#     end
#     ts
# end

# """
#     predicted_traces(ts::Vector{Vector}, model)
#     predicted_traces(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model)

# return vector of predicted traces and vector of state paths
# """
# function predicted_traces(ts::Vector, model)
#     tp = Vector{Float64}[]
#     d = model.reporter.probfn(model.rates[end-model.reporter.n+1:end], model.reporter.per_state, tcomponent(model).nT)
#     for t in ts
#         push!(tp, [mean(d[state]) for state in t])
#     end
#     tp
# end
# function predicted_traces(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model)
#     predicted_traces(predicted_states(data, model), model)
# end

# function predicted_traces(ts::Vector, model)
#     tp = Vector{Float64}[]
#     d = model.reporter.probfn(model.rates[end-model.reporter.n+1:end], model.reporter.per_state, tcomponent(model).nT)
#     for t in ts
#         push!(tp, [mean(d[state]) for state in t])
#     end
#     tp, ts
# end