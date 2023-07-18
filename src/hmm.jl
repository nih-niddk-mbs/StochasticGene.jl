### hmm.jl
### Fit models directly to intensity traces


"""
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation at time = interval
- `T`: transition rate matrix
"""
function kolmogorov_forward(Q::Matrix, interval)
    global Q_global = copy(Q)
    tspan = (0.0, interval)
    prob = ODEProblem(fkf, Matrix(I, size(T)), tspan)
    # sol = solve(prob,saveat=t, lsoda(),abstol = 1e-4, reltol = 1e-4)
    sol = solve(prob, lsoda(), save_everystep=false)
    return sol'
end

"""
fkf(u::Matrix,p,t)

"""
fkf(u::Matrix, p, t) = u * Q_global

"""
forward(a,b,p0)

"""
function forward(a, b, p0, T)
    α = similar(b)
    α[1] = p0 .* b[1]
    for t in 1:T-1
        α[t+1] = α[t] * a .* b[t+1]
    end
    return α, sum(α)
end

function forward_loop(a, b, p0, N, T)
    α = similar(b)
    α[1] = p0 .* b[1]
    for t in 1:T-1
        for j in 1:N
            m[j] = 0
            for i in 1:N
                m[j] += α[t][i] * a[i, j] 
            end
            α[t+1][j] = m[j] * b[t+1][j]
        end
    end
    return α, sum(α)
end

forward(a, b, p0, N, T) = forward_log(log.(a), log.(b), log.(p0), N, T)

function forward_log(loga, logb, logp0, N, T)
    ϕ = similar(logb)
    ψ = zeros(N)
    ϕ[1] = logp0 .+ logb[1]
    for t in 2:T
        for k in 1:N
        for j in 1:N
            ψ[j] = ϕ[t-1][j] + loga[j, k] + logb[t][k]
        end
        ϕ[t][k] = logsumexp(ψ)
    end
    logsumexp(ϕ[T])
end

"""
forward_scaled(a,b,p0)

"""
function forward_scaled(a, b, p0, T)
    α = similar(b)
    α̂ = similar(α)
    c = similar(a)
    C = similar(c)
    α[1] = p0 .* b[1]
    c[1] = 1 / sum(α)
    C[1] = c[1]
    for t in 2:T
        α̂[t-1] = cc * α[t-1]
        α[t] = α̂[t-1] * a .* b[t]
        c[t] = 1 / sum(α[t])
        C[t] *= c[t]
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

