### hmm.jl
### Fit models directly to intensity traces


"""
make_a(r, transitions, interval, G, R=0, S=0)

Return discrete hmm transition matrix a 
Computed by numerically integrating Kolmogorov Forward equation for the underlying stochastic continuous time Markov process behind the GM model

"""
function make_ap(r, transitions, interval, G, R=0, S=0)
    Q = make_mat(set_elements_T(transitions, collect(1:length(transitions))), r, G)
    kolmogorov_forward(Q, interval)[2],normalized_nullspace(sparse(Q'))
end

"""


"""
function set_b(trace)
    b = Vector[]
    for t in eachrow(trace)
        push!(b,[mod(t[2]+1,2),t[2]])
    end
    return b
end


"""
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation at time = interval
- `T`: transition rate matrix
"""
function kolmogorov_forward(Q, interval)
    global Q_global = copy(Q)
    tspan = (0.0, interval)
    prob = ODEProblem(fkf, Matrix(I, size(Q)), tspan)
    # sol = solve(prob,saveat=t, lsoda(),abstol = 1e-4, reltol = 1e-4)
    sol = solve(prob, lsoda(), save_everystep=false)
    return sol
end

function expected_transitions(α, a, b, β, T)
    ξ = Array{Matrix{Float64}}(undef, T)
    γ = Array{Vector{Float64}}(undef, T)
    for t in eachindex(ξ)
        ξ[t] = α[t]' * a * (b[t+1] .* β[t+1])
        ξ[t] ./= sum(ξ[t])
        γ[t] = ξ[t] ./ sum(ξ, dims=2)
    end
    return ξ, γ
end

expected_rate(ξ, γ) = sum(ξ) ./ sum(γ)




"""
fkf(u::Matrix,p,t)

"""
fkf(u::Matrix, p, t) = u * Q_global

"""
forward(a,b,p0)

"""
function forward(a, b, p0, T)
    α = [(p0 .* b[1])']
    for t in 1:T-1
        push!(α,α[t] * (a .* b[t+1]))
    end
    return α, sum(α)
end

function forward_loop(a, b, p0, N, T)
    α = [p0 .* b[1]]
    m = zeros(N)
    αt = zeros(N)
    for t in 1:T-1
        for j in 1:N
            m[j] = 0
            for i in 1:N
                m[j] += α[t][i] * a[i, j]
            end
            αt[j] = m[j] * b[t+1][j]
        end
        push!(α,αt)
    end
    return α, sum(α)
end

forward(a, b, p0, N, T) = forward_log(log.(a), log.(b), log.(p0), N, T)

function forward_log(a, b, p0, N, T)
    loga = log.(a)
    ϕt = zeros(N)
    ψ = zeros(N)
    ϕ = [log.(p0) .+ log.(b[1])]
    for t in 2:T
        for k in 1:N
            for j in 1:N
                ψ[j] = ϕ[t-1][j] + loga[j, k] + log.(b[t][k])
            end
            ϕt[k] = logsumexp(ψ)
        end
        push!(ϕ,ϕt)
    end
    ϕ, logsumexp(ϕ[T])
end

"""
forward_scaled(a,b,p0)

"""
function forward_scaled(a, b, p0, T)
    c = zeros(T)
    C = similar(c)
    α = [(p0 .* b[1])']
    α̂ = copy(α)
    c[1] = 1 / sum(α[1])
    C[1] = c[1]
    for t in 2:T
        push!(α, α̂[t-1] * (a .* b[t]))
        c[t] = 1 / sum(α[t])
        C[t] *= c[t]
        push!(α̂, C[t] * α[t])
    end
    return α, α̂, c, C
end

"""
backward(a,b)

"""
function backward(a, b, T)
    β = similar(b)
    β[T] = 1
    for t in T-1:-1:1
        β[t] = a * (b[t+1] .* β[t+1])
    end
    return β
end

function backward_log(loga, logb, logp0, N, T)
    ϕ = similar(logb)
    ψ = zeros(N)
    ϕ[T] = logp0 .+ logb[1]
    for t in T-1:-1:1
        for k in 1:N
            for j in 1:N
                ψ[j] = loga[j, k] + logb[t+1][k] + ϕ[t+1][k]
            end
            ϕ[t][k] = logsumexp(ψ)
        end
    end
    logsumexp(ϕ[1])
end

"""
backward_scaled(a,b)

"""
function backward_scaled(a, b, c, T)
    β = similar(b)
    β[T] = 1
    for t in T-1:-1:2
        β[t] = a * (b[t+1] .* β̂[t+1])
        β̂[t] = c[t] * β[t]
    end
    return β, β̂
end

function viterbi(loga, logb, logp0, N, T)
    ϕ = similar(logb)
    ψ = similar(ϕ)
    ϕ[1] = logp0 .+ logb[1]
    ψ[1] = 0
    for t in 2:T
        for j in 1:N
            m[j], ψ[t][j] = findmax(ϕ[t-1][j] + loga[:, j])
            ϕ[t] = m[j] + logb[t][j]
        end
    end
    maximum(ϕ[T])
end

