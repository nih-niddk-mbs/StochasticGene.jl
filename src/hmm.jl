### hmm.jl
### Fit models directly to intensity traces


"""
kolmogorov_forward(Q::Matrix,interval)

return the solution of the Kolmogorov forward equation at time = interval
- `T`: transition rate matrix
"""
function kolmogorov_forward(Q::Matrix,interval)
    global Q_global = copy(Q)
    tspan = (0.,interval)
    prob = ODEProblem(fkf,Matrix(I,size(T)),tspan)
    # sol = solve(prob,saveat=t, lsoda(),abstol = 1e-4, reltol = 1e-4)
    sol = solve(prob, lsoda(),save_everystep = false)
    return sol'
end

"""
fkf(u::Matrix,p,t)

"""
fkf(u::Matrix,p,t) = u*Q_global

"""
forward(a,b,p0)

"""
function forward(a,b,p0,T)
    α = similar(b)
    α[1] = p0 .* b[1]
    for t in 1:T-1
        α[t+1] = α[t] * a .* b[t+1]
    end
    return α, sum(α)
end

"""
forward_scaled(a,b,p0)

"""
function forward_scaled(a,b,p0,T)
    α = similar(b)
    α̂ = similar(α)  
    c = similar(a)
    C = similar(c)
    α[1] = p0 .* b[1]
    c[1] = 1/sum(α)
    C[1] = c[1]
    for t in 2:T
        α̂[t-1] = cc*α[t-1]
        α[t] = α̂[t-1] * a .* b[t]
        c[t] = 1/sum(α[t])
        C[t] *= c[t]
    end
    return α, α̂, c, C
end

"""
backward(a,b)

"""
function backward(a,b,T)
    β = similar(b)
    β[T] = 1
    for t in T-1:-1:2
        β[t] = a * (b[t+1].*β[t+1])
    end
    return β
end


"""
backward_scaled(a,b)

"""
function backward_scaled(a,b,c,T)
    β = similar(b)
    β[T] = 1
    for t in T-1:-1:2
        β[t] = a * (b[t+1].*β̂[t+1])
        β̂[t] = c[t] * β[t]
    end
    return β, β̂
end

function viterbi(loga,logb,p0,T)
    ϕ = similar(logb)
    ϕ[1] = log.(p0) + log.(logb[1])
    for t in 2:T
        ϕ[t] = maximum(ϕ[t-1] + loga) + logb[t]
    end
    maximum(ϕ[T])
end

