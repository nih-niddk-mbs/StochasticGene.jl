### hmm.jl
### Fit Markov models directly to intensity traces
###
### Notation in discrete HMM algorithms follows Rabier, 1989

function loglikelihood(param,data::AbstractTraceData,model)
	r = get_rates(param,model)
	loglikelihood(r,model.Gtransitions,data.interval,data.trace)
end

function loglikelihood(r, G, onstates, transitions, interval, trace)
    logpredictions = Array{Float64}(undef,0)
    for t in trace
        T = length(t)
        a, p0 = make_ap(r, transitions, interval, G)
        b = set_b(t)
        l = forward_log(a, b, p0, G, T)
        push!(logpredictions,logsumexp(l[:, T]))
    end
    -logsumexp(logpredictions), -logpredictions
end


"""
make_ap(r, transitions, interval, G, R=0, S=0)

Return computed discrete HMM transition probability matrix a and equilibrium state probability p0
a is computed by numerically integrating Kolmogorov Forward equation for the underlying stochastic continuous time Markov process behind the GM model
p0 is nullspace of transition rate matrix Q

Arguments:
- `r`: transition rates
- `transitions`: Tuple of G state transitions
- `interval`: time interval between intensity observations
- `G`: number of G states
- `R`: number of R states
- `S`: number of S states

Q is the transpose of the Markov process transition rate matrix

"""
function make_logap(r, transitions, interval, G, R=0, S=0)
    Q = make_mat(set_elements_T(transitions, collect(1:length(transitions))), r, G)
    log.(kolmogorov_forward(sparse(Q'), interval)[2]), log.(normalized_nullspace(Q))
end

"""
set_logb(trace, N, T)

returns matrix logb = P(Observation_i | State_j)

-`trace`: Tx2 matrix of intensities.  Col 1 = time, Col 2 = intensity
-`N`: number of hidden states
-`T`: number of observations

"""

function set_logb(trace, N, T, prob, params,onstates)
    d = prob(params, onstates, N)
    logb = Matrix{Float64}(undef, N, T)
    t = 1
    for obs in trace
        for j in 1:N
            logb[j, t] = logpdf(d[j],obs)
        end
        t += 1
    end
    return logb
end

function prob_GaussianMixture(par,onstates,N)
    d = Array{Distribution{Univariate, Continuous}}(undef,N)
    for i in 1:N
        if i ∈ onstates
            d[i] = MixtureModel(Normal, [(par[1], par[2]), (par[3], par[4])], [par[5], 1-par[5]])
        else
            d[i] = Normal(par[1],par[2])
        end
    end
    d
end

function prob_nonoise(obs, state::Int, onstates)
    if (obs > 0.5 && state ∈ onstates) || (obs < 0.5 && state ∉ onstates)
        return 1.0
    else
        return 0.0
    end
end

function prob_Poisson(obs, state, lambda)
    d = Poisson(lambda[state])
    pdf(d, obs)
end



function set_b_og(trace)
    T, M = size(trace)
    b = Matrix{Float64}(undef, M, T)
    t = 1
    for obs in eachrow(trace)
        b[:, t] = [mod(obs[2] + 1, 2), obs[2]]
        t += 1
    end
    return b
end

"""
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_forward(Q, interval)
    global Q_kf = copy(Q)
    tspan = (0.0, interval)
    prob = ODEProblem(fkf, Matrix(I, size(Q)), tspan)
    # sol = solve(prob,saveat=t, lsoda(),abstol = 1e-4, reltol = 1e-4)
    sol = solve(prob, lsoda(), save_everystep=false)
    return sol
end
"""
fkf(u::Matrix,p,t)

"""
fkf(u::Matrix, p, t) = u * Q_kf

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
viterbi(a, b, p0, N, T)

returns maximum likelihood state path using Viterbi algorithm

"""
function viterbi(a, b, p0, N, T)
    loga = log.(a)
    logb = log.(b)
    ϕ = similar(logb)
    ψ = similar(ϕ)
    q = Vector{Int}(undef, T)
    ϕ[:, 1] = log.(p0) + logb[:, 1]
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

