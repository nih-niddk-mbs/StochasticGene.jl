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
set_b(trace) = set_b(trace,2,size(trace)[1])

function set_b(trace,N,T)
    b = Matrix{Float64}(undef, N, T)
    i = 1
    for t in eachrow(trace)
        b[:,i] = [mod(t[2]+1,2),t[2]]
        # push!(b,[mod(t[2]+1,2),t[2]])
        # push!(b,pdf(d,mod(t[2]+1,2)))
        i += 1
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

function expected_transitions(α, a, b, β, N, T)
    ξ = fill(zeros(N,N),T-1)
    γ = fill(zeros(N),T-1)
    for t in 1:T-1
        for i = 1:N
            for j = 1:N
                ξ[t][i,j] = α[t][i] * a[i,j] * b[t+1][j] * β[t+1][j]
            end
            γ[t][i] = sum(ξ[t][i,:])
        end
        ξ[t]
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
        push!(α,α[t] * (a .* b[:,t+1]))
    end
    return α, sum(α)
end

function forward_loop(a, b, p0, N, T)
    α = Matrix{Float64}(undef, N, T)
    α[:, 1] = p0 .* b[:,1]
    for t in 1:T-1
        m = zeros(N)
        for j in 1:N
            for i in 1:N
                m[j] += α[i, t] * a[i, j]
            end
            α[j, t+1] = m[j] * b[j,t+1]
        end
    end
    return α, sum(α)
end

function forward_log(a, b, p0, N, T)
    loga = log.(a)
    ψ = zeros(N)
    ϕ = Matrix{Float64}(undef, N, T)
    ϕ[:,1] = log.(p0) .+ log.(b[:,1])
    for t in 2:T
        for k in 1:N
            for j in 1:N
                ψ[j] = ϕ[j,t-1] + loga[j, k] + log.(b[k,t])
            end
            ϕ[k,t] = logsumexp(ψ)
        end
    end
    ϕ, logsumexp(ϕ[:,T])
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
    β = [[1.,1.]]
    for t in 1:T-1
        push!(β, a * (b[T-t+1] .* β[t]))
    end
    return reverse(β)
end

function backward_log(a, b, N, T)
    loga = log.(a)
    ψ = zeros(N)
    ϕt = zeros(N)
    ϕ = [[0.,0.]]
    for t in 1:T-1
        for k in 1:N
            for j in 1:N
                ψ[j] = loga[j, k] + log.(b[T-t+1][k]) + ϕ[t][k]
            end
            ϕt[k] = logsumexp(ψ)
        end
        push!(ϕ,ϕt)
    end
    ϕ, logsumexp(ϕ[T])
end

"""
backward_scaled(a,b)

"""
function backward_scaled(a, b, c, T)
    β = [c[T] .* [1,1]]
    β̂ = [c[T] .* [1,1]]
    for t in 1:T-1
        push!(β, a * (b[T-t+1].* β̂[1]))
        β̂[t+1] = c[T-t] * β[t+1]
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

