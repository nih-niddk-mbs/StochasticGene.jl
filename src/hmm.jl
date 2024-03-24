# This file is part of StochasticGene.jl   

### hmm.jl
### Fit discrete HMMs and continuous hidden Markov process models directly to intensity traces
###
### Notation in discrete HMM algorithms follows Rabier, 1989
###
### Functions for forward, backward, and Viterbi HMM algorihms
### For continuous processes, numerically solve forward Kolmogorov equation to obtain transition probability matrix
###

"""
    ll_hmm(r, nT, reporters, elementsT, interval, trace)

return total loglikelihood of traces with reporter noise and loglikelihood of each trace
"""
function ll_hmm_hierarchical(r::Matrix, nT, elementsT::Vector, noiseparams, reporters_per_state, probfn, interval, trace)
    logpredictions = Array{Float64}(undef, 0)
    for i in eachindex(trace)
        T = length(trace[i])
        loga, logp0 = make_logap(r[:, i], interval, elementsT, nT)
        logb = set_logb(trace[i], nT, r[end-noiseparams+1:end, i], reporters_per_state, probfn)
        l = forward_log(loga, logb, logp0, nT, T)
        push!(logpredictions, logsumexp(l[:, T]))
    end
    -sum(logpredictions), -logpredictions
end
function ll_hmm(r, nT, elementsT::Vector, noiseparams, reporters_per_state, probfn, interval, trace)
    loga, logp0 = make_logap(r, interval, elementsT, nT)
    ll_hmm(r, nT, noiseparams, reporters_per_state, probfn, trace, loga, logp0)
end

function ll_hmm(r, nT, noiseparams::Int, reporters_per_state, probfn, trace, loga, logp0)
    logpredictions = Array{Float64}(undef, 0)
    for t in trace[1]
        T = length(t)
        logb = set_logb(t, nT, r[end-noiseparams+1:end], reporters_per_state, probfn)
        l = forward_log(loga, logb, logp0, nT, T)
        push!(logpredictions, logsumexp(l[:, T]))
    end
    # for t in trace[2]
    #     T = length(t)
    #     logb = set_logb(t, nT, r[end-noiseparams+1:end], reporters_per_state, probfn)
    #     l = forward_log(loga, logb, logp0, nT, T)
    #     push!(logpredictions, trace[3] / length(trace[2]) * logsumexp(l[:, T]))
    # end
    -sum(logpredictions), -logpredictions
end

function ll_hmm(r, nT, elementsT::Vector, noiseparams, reporters_per_state, probfn, interval, trace, nascent)
    a, p0 = make_ap(r, interval, elementsT, nT)
    lln = ll_nascent(p0, reporters_per_state, nascent)
    ll, logpredictions = ll_hmm(r, nT, noiseparams, reporters_per_state, probfn, trace, log.(max.(a, 0)), log.(p0))
    push!(logpredictions, lln)
    ll += lln
    ll, logpredictions
end

function ll_nascent(p0, reporters_per_state, nascent)
    pn = sum(p0[reporters_per_state.>0])
    d = Binomial(nascent[2], pn)
    -logpdf(d, nascent[1])
end

"""
make_ap(r, interval, elementsT, nT )

Return computed discrete HMM transition probability matrix a and equilibrium state probability p0
a is computed by numerically integrating Kolmogorov Forward equation for the underlying stochastic continuous time Markov process behind the GRSM model
p0 is left nullspace of transition rate matrix Q (right nullspace of Q')

Arguments:
- `r`: transition rates
- `interval`: time interval between intensity observations
- `elementsT`: structure of T matrix elements
- `N`: number of HMM states

Qtr is the transpose of the Markov process transition rate matrix Q

"""
function make_ap(r, interval, elementsT, N)
    Qtr = make_mat(elementsT, r, N) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(sparse(Qtr'), interval), normalized_nullspace(Qtr)
end

make_p0(r, elementsT, N) = normalized_nullspace(make_mat(elementsT, r, N))
"""
    make_logap(r, transitions, interval, G)

return log of a and p0
"""
function make_logap(r, interval, elementsT, N)
    a, p0 = make_ap(r, interval, elementsT, N)
    log.(max.(a, 0)), log.(max.(p0, 0))
end

"""
set_logb(trace, N, params, reporters)

returns matrix logb = P(Observation_i | State_j) for Gaussian distribution

-`trace`: Tx2 matrix of intensities.  Col 1 = time, Col 2 = intensity
-`N`: number of hidden states
-`T`: number of observations

"""
function set_logb(trace, N, params, reporters, probfn=prob_Gaussian)
    d = probfn(params, reporters, N)
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

function set_logb_discrete(trace, N, onstates)
    logb = Matrix{Float64}(undef, N, length(trace))
    t = 1
    for obs in trace
        if (obs > 0.5 && state ∈ onstates) || (obs < 0.5 && state ∉ onstates)
            logb[j, t] = 0.0
        else
            logb[j, t] = -Inf
        end
    end
end

"""
    prob_Gaussian(par, reporters, N)

return Gaussian Distribution 
mean = background + number of reporters * reporter mean
variance = sum of variances of background and reporters

- `par`: 4 dimemsional vector of mean and std parameters
- `reporters`: number of reporters per HMM state
-`N`: number of HMM states
"""
function prob_Gaussian(par, reporters, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = Normal(par[1] + reporters[i] * par[3], sqrt(par[2]^2 + reporters[i] * par[4]^2))
    end
    d
end
"""
    prob_GaussianMixture(par,reporters,N)

return Gaussian Mixture distribution with 4 Gaussian parameters and 1 weight parameter

"""
function prob_GaussianMixture(par, reporters, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = MixtureModel(Normal, [(par[1] + reporters[i] * par[3], sqrt(par[2]^2 + reporters[i] * par[4]^2)), (par[1], par[2])], [par[5], 1 - par[5]])
    end
    d
end


"""
    prob_GaussianMixture_6(par, reporters, N)

Gaussian Mixture distribution with 6 Gaussian parameters and 1 weight parameter
"""
function prob_GaussianMixture_6(par, reporters, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = MixtureModel(Normal, [(par[1] + reporters[i] * par[3], sqrt(par[2]^2 + reporters[i] * par[4]^2)), (par[5], par[6])], [par[7], 1 - par[7]])
    end
    d
end


"""
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_forward(Q, interval, method=Tsit5())
    tspan = (0.0, interval)
    prob = ODEProblem(fkf!, Matrix(I, size(Q)), tspan, Q)
    solve(prob, method, save_everystep=false)[:, 2]
end
"""
    fkf!(du,u::Matrix, p, t)

in place update of du of ODE system for DifferentialEquations,jl
"""
function fkf!(du, u::Matrix, p, t)
    du .= u * p
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
    loga = log.(a)
    logb = log.(b)
    logp0 = log.(p0)
    viterbi(loga, logb, logp0, N, T)
end


"""
    predicted_statepath(r::Vector, N::Int, elementsT, noiseparams, reporters_per_state, probfn, T::Int, interval)
    predicted_statepath(model::AbstractGmodel, T, interval)
    predicted_statepath(r, tcomponents, reporter, T, interval)

return predicted state path using Viterbi algorithm
"""
function predicted_statepath(trace, interval, r::Vector, N::Int, elementsT, noiseparams, reporters_per_state, probfn)
    loga, logp0 = make_logap(r, interval, elementsT, N)
    logb = set_logb(trace, N, r[end-noiseparams+1:end], reporters_per_state, probfn)
    viterbi(loga, logb, logp0, N, length(trace))
end

function predicted_statepath(trace, interval, model::AbstractGmodel)
    tcomponents = tcomponent(model)
    predicted_statepath(trace, interval, model.rates, tcomponents.nT, tcomponents.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn)
end

function predicted_statepath(trace, interval, r, tcomponents, reporter)
    predicted_statepath(trace, interval, r, tcomponents.nT, tcomponents.elementsT, reporter.n, reporter.per_state, reporter.probfn)
end


"""
    predicted_states(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model::AbstractGmodel)

return vector of predicted state vectors
"""
function predicted_states(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model::AbstractGmodel)
    ts = Vector{Int}[]
    for t in data.trace[1]
        push!(ts, predicted_statepath(t, data.interval, model))
    end
    ts
end



"""
    predicted_trace(statepath, noise_dist)
    predicted_trace(statepath, r, reporter, nstates)

return predicted trace from state path
"""
function predicted_trace(statepath, noise_dist)
    [mean(noise_dist[state]) for state in statepath]
end

function predicted_trace(statepath, r, reporter, nstates)
    d = model.reporter.probfn(r[end-model.reporter.n+1:end], reporter.per_state, nstates)
    predicted_trace(statepath, d)
end

"""
    predicted_trace_state(trace, interval, r::Vector, tcomponents, reporter, noise_dist)

return predicted trace and state path
"""
function predicted_trace_state(trace, interval, r::Vector, tcomponents, reporter, noise_dist)
    t = predicted_statepath(trace, interval, r, tcomponents, reporter)
    tp = predicted_trace(t, noise_dist)
    return tp, t
end

function predicted_trace_state(trace, interval, model)
    d = model.reporter.probfn(model.rates[end-model.reporter.n+1:end], model.reporter.per_state, tcomponent(model).nT)
    predicted_trace_state(trace, interval, model.rates, model.tcomponents, model.reporter, d)
end

# function predicted_trace_state(trace, interval, model)
#     t = predicted_statepath(trace, interval, model)
#     d = model.reporter.probfn(model.rates[end-model.reporter.n+1:end], model.reporter.per_state, tcomponent(model).nT)
#     tp = [mean(d[state]) for state in t]
#     return tp, t
# end



"""
    predicted_traces(ts::Vector{Vector}, model)
    predicted_traces(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model)

return vector of predicted traces and vector of state paths
"""
function predicted_traces(ts::Vector, model)
    tp = Vector{Float64}[]
    d = model.reporter.probfn(model.rates[end-model.reporter.n+1:end], model.reporter.per_state, tcomponent(model).nT)
    for t in ts
        push!(tp, [mean(d[state]) for state in t])
    end
    tp, ts
end
function predicted_traces(data::Union{AbstractTraceData,AbstractTraceHistogramData}, model)
    predicted_traces(predicted_states(data, model), model)
end

"""
    tcomponent(model)

return tcomponent of model
"""
tcomponent(model) = typeof(model.components) == TComponents ? model.components : model.components.tcomponents