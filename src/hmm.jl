### hmm.jl
### Fit models directly to intensity traces

function forward()
    α[1] = p0 .* b[1]
    for t = 2:T
        α[t+1] = α[t] * a .* b[t+1]
    end

    sum.(α)

end

backward()

function kolmogorov_forward(T::Matrix,interval)
    global T_global = copy(T)
    tspan = (0.,interval)
    prob = ODEProblem(fkf,Matrix(I,size(T)),tspan)
    # sol = solve(prob,saveat=t, lsoda(),abstol = 1e-4, reltol = 1e-4)
    sol = solve(prob, lsoda(),save_everystep = false)
    return sol'
end

fkf(u::Matrix,p,t) = u*T_global
