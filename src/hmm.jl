# This file is part of StochasticGene.jl   

### hmm.jl
### Fit discrete HMMs and continuous hidden Markov process models directly to observations (e.g. intensity traces)
###
### Notation in discrete HMM algorithms follows Rabier, 1989
###
### Functions for forward, backward, and Viterbi HMM algorihms
### For continuous processes, numerically solve forward Kolmogorov equation to obtain transition probability matrix
###

using CUDA
using CUDA: @allowscalar

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


function prob_Gaussian(par, reporters::Int)
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
function prob_Gaussian(par, reporters_per_state::T, N) where {T<:Vector}
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_Gaussian(par, reporters_per_state[i])
    end
    d
end

function prob_Gaussian(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
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
function prob_GaussianMixture(par, reporters::Int)
    MixtureModel(Normal, [(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2)), (par[1], par[2])], [par[5], 1 - par[5]])
end

function prob_GaussianMixture(par, reporters_per_state::T, N) where {T<:Vector}
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture(par, reporters_per_state[i])
    end
    d
end

function prob_GaussianMixture(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture(par, reporters_per_state[i])
    end
    d
end


"""
    prob_GaussianMixture_6(par, reporters_per_state, N)

Gaussian Mixture distribution with 6 Gaussian parameters and 1 weight parameter
"""
function prob_GaussianMixture_6(par, reporters)
    MixtureModel(Normal, [(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2)), (par[5], par[6])], [par[7], 1 - par[7]])
end

function prob_GaussianMixture_6(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture_6(par, reporters_per_state[i])
    end
    d
end

function prob_GaussianMixture_6(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture_6(par, reporters_per_state[i])
    end
    d
end


"""
    prob_Gaussian_ind(par, reporters_per_state, N)

TBW
"""
function prob_Gaussian_ind(par, reporters::Int)
    if reporters > 0
        return Normal(reporters * par[3], sqrt(reporters) * par[4])
    else
        return Normal(par[1], par[2])
    end
end

function prob_Gaussian_ind(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_Gaussian_ind(par, reporters_per_state[i])
    end
    d
end

function prob_Gaussian_ind(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_Gaussian_ind(par, reporters_per_state[i])
    end
    d
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
function prob_Gaussian_grid(par, reporters_per_state, Nstate::Int, Ngrid::Int, f::Function=kronecker_delta)
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

function prob_Gaussian_grid(par, reporters_per_state, N::Tuple, f::Function=kronecker_delta)
    Nstate, Ngrid = N
    prob_Gaussian_grid(par, reporters_per_state, Nstate, Ngrid, f)
end

function prob_Gaussian_grid(par, reporters_per_state, Ngrid::Int, f::Function=kronecker_delta)
    Nstate = length(reporters_per_state)
    prob_Gaussian_grid(par, reporters_per_state, Nstate, Ngrid, f)
end

"""
    make_ap(r, interval, components, method)

Return computed discrete HMM transition probability matrix a and equilibrium state probability p0
a is computed by numerically integrating Kolmogorov Forward equation for the underlying stochastic continuous time Markov process behind the GRSM model
p0 is left nullspace of transition rate matrix Q (right nullspace of Q')

Arguments:
- `r`: transition rates
- `interval`: time interval between intensity observations (frame interval)
- `components`: T matrix components

Qtr is the transpose of the Markov process transition rate matrix Q

"""
function make_ap(rates, interval, components::TComponents, method=Tsit5())
    Qtr = make_mat_T(components, rates) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

function make_ap(rates, couplingStrength, interval, components::TCoupledComponents, method=Tsit5())
    Qtr = make_mat_TC(components, rates, couplingStrength)
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

function make_ap(r::Tuple, interval, components::TComponents, method=Tsit5())
    r, couplingStrength = r
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

function make_a_grid(parN)
    par, Ngrid = parN
    make_a_grid(par, Ngrid)
end


"""
    set_d(noiseparams, reporters_per_state, probfn, N)

TBW
"""

function set_d(noiseparams, reporters_per_state::Vector{Int}, probfn::T, N::Int) where {T<:Function}
    probfn(noiseparams, reporters_per_state, N)
end

function set_d(noiseparams::Vector{T}, reporters_per_state::Vector{Vector{Int}}, probfn::Vector, N::Int) where {T<:AbstractVector}
    d = Vector{Distribution}[]
    for i in eachindex(noiseparams)
        push!(d, probfn[i](noiseparams[i], reporters_per_state[i], N))
    end
    return d
end

function set_d(noiseparams, reporters_per_state, probfn, Nstate::Int, Ngrid::Int)
    probfn(noiseparams, reporters_per_state, Nstate, Ngrid)
end

function set_d(noiseparams, reporters_per_state, probfn, Ngrid::Int)
    probfn(noiseparams, reporters_per_state, Ngrid)
end


"""
    set_d(noiseparams, reporter, N)


TBW
"""
function set_d(noiseparams, reporter::HMMReporter, N::Int)
    set_d(noiseparams, reporter.per_state, reporter.probfn, N)
end

function set_d(noiseparams, reporter::Vector{HMMReporter}, N::Int)
    ps = [r.per_state for r in reporter]
    pf = [r.probfn for r in reporter]
    set_d(noiseparams, ps, pf, N)
end



"""
    set_d(noiseparams, reporters_per_state, probfn)

TBW
"""

function set_d(noiseparams, reporters_per_state::Vector{Int}, probfn::T) where {T<:Function}
    probfn(noiseparams, reporters_per_state)
end

function set_d(noiseparams::Vector{T}, reporters_per_state::Vector{Vector{Int}}, probfn::Vector) where {T<:AbstractVector}
    d = Vector{Distribution}[]
    for i in eachindex(noiseparams)
        push!(d, probfn[i](noiseparams[i], reporters_per_state[i]))
    end
    return d
end

"""
    set_d(noiseparams, reporter)


TBW
"""
function set_d(noiseparams, reporter::HMMReporter)
    set_d(noiseparams, reporter.per_state, reporter.probfn)
end

function set_d(noiseparams, reporter::Vector{HMMReporter})
    ps = [r.per_state for r in reporter]
    pf = [r.probfn for r in reporter]
    set_d(noiseparams, ps, pf)
end


"""
    set_b(trace, d)

TBW
"""
function set_b(trace::Vector, d::Vector{T}) where {T<:Distribution}
    N = length(d)
    b = Matrix{Float64}(undef, N, length(trace))
    for (t, obs) in enumerate(trace)
        for j in 1:N
            b[j, t] = pdf(d[j], obs)
        end
    end
    return b
end

function set_b(trace, d::Vector{T}) where {T<:Vector}
    N = length(d[1])
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

function set_b(trace, d::Array{T,3}) where {T<:Distribution}
    Nstate, Ngrid, _ = size(d)
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
set_b(trace, d, N)

"""

function set_b(trace::Vector, d::Vector, N)
    b = Matrix{Float64}(undef, N, length(trace))
    for (t, obs) in enumerate(trace)
        for j in 1:N
            b[j, t] = pdf(d[j], obs)
        end
    end
    return b
end




function set_b(trace, d::Vector{T}, N) where {T<:Vector}
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



function set_b(trace, d, Nstate, Ngrid)
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

function set_b(trace, d, N::Tuple)
    Nstate, Ngrid = N
    set_b(trace, d, Nstate, Ngrid)
end

"""
    set_b(trace, params, reporters_per_state, probfn, N)

TBW
"""

function set_b(trace, params, reporters_per_state, probfn, N)
    d = set_d(params, reporters_per_state, probfn, N)
    set_b(trace, d, N)
end

function set_b(trace, params, reporters_per_state::Vector{Int}, probfn::T, N::Int) where {T<:Function}
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

function set_b(trace, params, reporters_per_state, probfn, N::Tuple)
    Nstate, Ngrid = N
    d = probfn(params, reporters_per_state, Nstate, Ngrid)
    set_b(trace, d, Nstate, Ngrid)
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

function set_b_background(obs, d::Vector{Distribution{Univariate,Continuous}})
    b = Array{Float64}(undef, size(d))
    for j in CartesianIndices(d)
        b[j] = pdf(d[j], obs)
    end
    return reshape(b, :, 1)
end

function set_b_background(obs, d::Vector{<:Vector}, k::Int, N)
    b = ones(N)
    for j in 1:N
        b[j] *= pdf(d[k][j], obs)
    end
    return reshape(b, :, 1)
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
    forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)

TBW
"""
function forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    if CUDA.functional() && (Nstate * Nstate * Ngrid * Ngrid * T > 1000)
        return forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)
    else
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
end

"""
    forward_grid(a, a_grid, b, p0)

TBW
"""
function forward_grid(a, a_grid, b, p0)
    Nstate, Ngrid, T = size(b)
    forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
end

"""
forward(a, b, p0, N, T)

returns forward variable α, and scaling parameter array C using scaled forward algorithm
α[i,t] = P(O1,...,OT,qT=Si,λ)
Ct = Prod_t 1/∑_i α[i,t]

# """
function forward(a::Matrix, b, p0, N, T)
    if CUDA.functional() && (N * N * T > 1000)
        return forward_gpu(a, b, p0, N, T)
    else
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
end

"""
    forward(a, b, p0)

TBW
"""
function forward(a::Matrix, b, p0)
    N, T = size(b)
    forward(a, b, p0, N, T)
end

function forward(atuple::Tuple, b::Array, p0)
    a, a_grid = atuple
    forward_grid(a, a_grid, b, p0)
end

"""
    forward_gpu(a::Matrix{Float64}, b::Matrix{Float64}, p0::Vector{Float64}, N::Int64, T::Int64)

GPU-accelerated version of the forward algorithm.
Returns the forward variable α and scaling parameter array C.
"""
function forward_gpu(a::Matrix{Float64}, b::Matrix{Float64}, p0::Vector{Float64}, N::Int64, T::Int64)
    # Move data to GPU
    a_gpu = CuArray(a)
    b_gpu = CuArray(b)
    p0_gpu = CuArray(p0)
    
    # Allocate GPU arrays for α and C
    α_gpu = CUDA.zeros(Float64, N, T)
    C_gpu = CUDA.zeros(Float64, T)
    
    # Initialize first time step
    # Use proper GPU array operations instead of scalar indexing
    α_gpu[:, 1] .= p0_gpu .* b_gpu[:, 1]
    
    # Use CUDA.@allowscalar for operations that require scalar indexing
    CUDA.@allowscalar C_gpu[1] = 1 / sum(α_gpu[:, 1])
    CUDA.@allowscalar α_gpu[:, 1] *= C_gpu[1]
    
    # Compute forward probabilities using matrix multiplication on GPU
    for t in 2:T
        # Use proper GPU array operations
        α_gpu[:, t] .= (a_gpu * α_gpu[:, t-1]) .* b_gpu[:, t]
        
        # Use CUDA.@allowscalar for operations that require scalar indexing
        CUDA.@allowscalar C_gpu[t] = 1 / sum(α_gpu[:, t])
        CUDA.@allowscalar α_gpu[:, t] *= C_gpu[t]
    end
    
    # Return results back to CPU
    return Array(α_gpu), Array(C_gpu)
end

"""
    forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)

GPU-accelerated version of the forward_grid algorithm that computes forward probabilities in parallel.
Returns the forward variable α and scaling parameter array C.

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution
- `Nstate`: Number of states
- `Ngrid`: Number of grid points
- `T`: Number of time steps

# Returns
- Tuple of (α, C) where α is the forward variable and C is the scaling parameter array
"""
function forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)
    # Move data to GPU
    a_gpu = CuArray(a)
    a_grid_gpu = CuArray(a_grid)
    b_gpu = CuArray(b)
    p0_gpu = CuArray(p0)

    # Allocate GPU arrays for results
    α_gpu = CuArray{Float64}(undef, Nstate, Ngrid, T)
    C_gpu = CuArray{Float64}(undef, T)

    # Initialize first time step
    α_gpu[:, :, 1] .= p0_gpu .* b_gpu[:, :, 1]
    C_gpu[1] = sum(α_gpu[:, :, 1])
    α_gpu[:, :, 1] ./= C_gpu[1]

    # Compute forward probabilities using matrix multiplication on GPU
    for t in 2:T
        # Reshape for efficient matrix multiplication
        α_prev = reshape(α_gpu[:, :, t-1], Nstate * Ngrid, 1)
        a_combined = kron(a_grid_gpu, a_gpu)  # Kronecker product for combined transitions
        α_temp = reshape(a_combined * α_prev, Nstate, Ngrid)
        α_gpu[:, :, t] .= α_temp .* b_gpu[:, :, t]

        C_gpu[t] = sum(α_gpu[:, :, t])
        α_gpu[:, :, t] ./= C_gpu[t]
    end

    # Return results back to CPU
    return Array(α_gpu), Array(C_gpu)
end

"""
    forward_grid_gpu(a, a_grid, b, p0)

GPU-accelerated version of forward_grid that automatically determines dimensions.
Returns the forward variable α and scaling parameter array C.

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution

# Returns
- Tuple of (α, C) where α is the forward variable and C is the scaling parameter array
"""
function forward_grid_gpu(a, a_grid, b, p0)
    Nstate, Ngrid, T = size(b)
    forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)
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
    ll_off(obs, d, a, p0, weight)

TBW
"""
# function ll_off(obs::Float64, d::Vector{Distribution{Univariate,Continuous}}, a::Matrix, p0, weight)
#     b = set_b_background(obs, d)
#     _, C = forward(a, b, p0)
#     weight * sum(log.(C))
# end

# function ll_off(obs::Vector, d::Vector{Distribution{Univariate,Continuous}}, a::Matrix, p0, weight)
#     b = set_b(obs, d)
#     _, C = forward(a, b, p0)
#     weight * sum(log.(C))
# end



# function ll_off(obs::Vector, d::Vector{T}, a::Matrix, p0, weight) where {T<:Vector}
#     l = 0
#     for i in eachindex(obs)
#         b = set_b(obs[i], d)
#         _, C = forward(a, b, p0)
#         l += weight[i] * sum(log.(C))
#     end
#     l
# end

# function ll_off(trace, rates, noiseparams, reporter, interval, components, method)
#     a, p0 = make_ap(rates, interval, components.elementsT, components.nT, method)
#     d = set_d(noiseparams, reporter, components.nT)
#     b = set_b(trace[2], d)
#     _, C = forward(a, b, p0)
#     sum(log.(C)) * trace[3] * length(trace[1])
# end

function ll_off(trace::Tuple, d::Vector{Distribution{Univariate,Continuous}}, a::Matrix, p0)
    if trace[3] > 0.0
        b = set_b(trace[2], d)
        _, C = forward(a, b, p0)
        sum(log.(C)) * trace[3] * length(trace[1])
    else
        0.0
    end
end

function ll_off(trace::Tuple, rates, noiseparams, reporter::Vector, interval, components, method)
    l = 0.0
    components = components.modelcomponents
    dims = [components[i].nT for i in eachindex(components)]
    for i in eachindex(rates)
        if typeof(trace[3]) <: AbstractVector && trace[3][i] > 0.0
            a, p0 = make_ap(rates[i], interval, components[i].elementsT, components[i].nT, method)
            rps = reduce_reporters_per_state(reporter[i].per_state, dims, i)
            d = set_d(noiseparams[i], rps, reporter[i].probfn, dims[i])
            b = set_b(trace[2][i], d)
            _, C = forward(a, b, p0)
            l += sum(log.(C)) * trace[3][i]
        end
    end
    l * length(trace[1])
end

    ### Called by trait likelihoods

"""
    ll_hmm(a::Matrix, p0::Vector, d, traces)

"""
function ll_hmm(a::Matrix, p0::Vector, d, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(noiseparams, a::Matrix, p0::Vector, reporter, traces)

"""
function ll_hmm(noiseparams::Vector, a::Matrix, p0::Vector, reporter, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(r::Vector, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function ll_hmm(r::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], interval, components, method)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(r, couplingStrength, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function ll_hmm(r::Vector, couplingStrength::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], couplingStrength[i], interval, components, method)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

### grid trait likelihoods

"""
    ll_hmm(a, a_grid, p0::Vector, d, traces)

"""
function ll_hmm(a::Matrix, a_grid::Matrix, p0::Vector, d, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        _, C = forward_grid(a, a_grid, set_b(traces[i], d), p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(noiseparams, a, a_grid, p0::Vector, reporter, traces)

"""
function ll_hmm(noiseparams::Vector, Ngrid::Int, a::Matrix, a_grid::Matrix, p0::Vector, reporter, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        _, C = forward_grid(a, a_grid, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(r::Vector, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function ll_hmm(r::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], interval, components, method)
        a_grid = make_a_grid(pgrid, Ngrid)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        _, C = forward_grid(a, a_grid, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(r, couplingStrength, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function ll_hmm(r::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], couplingStrength[i], interval, components, method)
        a_grid = make_a_grid(pgrid, Ngrid)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        _, C = forward_grid(a, a_grid, b, p0)
        @inbounds logpredictions[i] = sum(log.(C))
    end
    sum(logpredictions), logpredictions
end


###  likelihoods
"""
    ll_hmm(r, components, reporter, interval, trace, method)

"""
# no traits
function ll_hmm(r::Tuple{T1,T2}, components::TComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2}
    r, noiseparams = r
    a, p0 = make_ap(r, interval, components, method)
    d = set_d(noiseparams, reporter)
    # lb = trace[3] > 0.0 ? length(trace[1]) * ll_off(trace[2], d, a, p0, trace[3]) : 0.0
    lb = ll_off(trace, d, a, p0)
    ll, logpredictions = ll_hmm(a, p0, d, trace[1])
    ll + lb, logpredictions
end

# coupled
function ll_hmm(r::Tuple{T1,T2,T3}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3}
    rates, noiseparams, couplingStrength = r
    a, p0 = make_ap(rates, couplingStrength, interval, components, method)
    d = set_d(noiseparams, reporter)
    # lb = any(trace[3] .> 0.0) ? length(trace[1]) * ll_off(trace[2], d, a, p0, trace[3]) : 0.0
    lb = ll_off(trace, rates, noiseparams, reporter, interval, components, method)
    ll, logpredictions = ll_hmm(a, p0, d, trace[1])
    ll + lb, logpredictions
end

# hierarchical
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6}, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    d = set_d(noiseshared[1], reporter)
    # lb = any(trace[3] .> 0.0) ? length(trace[1]) * ll_off(trace[2], d, a, p0, trace[3]) : 0.0
    lb = ll_off(trace, d, a, p0)
    if method[2]
        ll, logpredictions = ll_hmm(noiseindividual, a, p0, reporter, trace[1])
    else
        ll, logpredictions = ll_hmm(rindividual, noiseindividual, interval, components, reporter, trace[1], method[1])
    end
    lhp = ll_hierarchy(pindividual, rhyper)
    ll + lb + sum(lhp), vcat(logpredictions, lhp)
end

# coupled, hierarchical
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, couplingshared, couplingindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    d = set_d(noiseshared[1], reporter)
    # lb = any(trace[3] .> 0.0) ? length(trace[1]) * ll_off(trace[2], d, a, p0, trace[3]) : 0.0
    lb = ll_off(trace, rshared[1], noiseshared[1], reporter, interval, components, method[1])
    if method[2]
        ll, logpredictions = ll_hmm(noiseindividual, a, p0, reporter, trace[1])
    else
        ll, logpredictions = ll_hmm(rindividual, couplingindividual, noiseindividual, interval, components, reporter, trace[1])
    end
    lhp = ll_hierarchy(pindividual, rhyper)
    ll + lb + sum(lhp), vcat(logpredictions, lhp)
end

### grid trait likelihoods
"""
    ll_hmm(r, pgrid, Ngrid, components, reporter, interval, trace, method)
    
"""
# grid
function ll_hmm(r::Tuple{T1,T2,T3}, Ngrid::Int, components::TComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2,T3}
    r, noiseparams, pgrid = r
    a, p0 = make_ap(r, interval, components, method)
    a_grid = make_a_grid(pgrid[1], Ngrid)
    d = set_d(noiseparams, reporter, Ngrid)
    ll, logpredictions = ll_hmm(a, a_grid, p0, d, trace[1])
    ll, logpredictions
end

# coupled, grid
function ll_hmm(r::Tuple{T1,T2,T3,T4}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3,T4}
    r, noiseparams, couplingStrength, pgrid = r
    a, p0 = make_ap(r, couplingStrength, interval, components, method)
    a_grid = make_a_grid(pgrid[1][1], Ngrid)
    d = set_d(noiseparams, reporter, Ngrid)
    ll, logpredictions = ll_hmm(a, a_grid, p0, d, trace[1])
    ll, logpredictions
end

# hierarchical, grid
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, Ngrid::Int, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, _, noiseindividual, _, _, pgridshared, pgridindividual = r
    if method[2]
        a, p0 = make_ap(rshared[1], interval, components, method[1])
        a_grid = make_a_grid(pgridshared[1][1], Ngrid)
        states, observation_dist = predict_trace(noiseindividual, Ngrid, a, a_grid, p0, reporter, trace[1])
    else
        states, observation_dist = predict_trace(rindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1], method[1])
    end
    lhp = ll_hierarchy(pindividual, rhyper)
    ll + sum(lhp), vcat(observation_dist, lhp)
end

# coupled, hierarchical, grid
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    rshared, rindividual, _, noiseindividual, _, _, couplingshared, couplingindividual, pgridshared, pgridindividual = r
    if method[2]
        a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
        a_grid = make_a_grid(pgridshared[1], Ngrid)
        states, observation_dist = predict_trace(noiseindividual, a, a_grid, p0, reporter, trace[1])
    else
        states, observation_dist = predict_trace(rindividual, couplingindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1], method[1])
    end
    lhp = ll_hierarchy(pindividual, rhyper)
    ll + sum(lhp), vcat(observation_dist, lhp)
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

function viterbi(a, b, p0)
    N, T = size(b)
    viterbi(a, b, p0, N, T)
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

function viterbi_grid(a, a_grid, b, p0)
    Nstate, Ngrid, T = size(b)
    viterbi_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
end

"""
    covariance_functions(rin, transitions, G::Tuple, R, S, insertstep, interval, probfn, coupling, lags::Vector)

TBW
"""
function covariance_functions(rin, transitions, G::Tuple, R, S, insertstep, interval, probfn, coupling, lags::Vector)
    components = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    # components = TRGCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    sourceStates = [c.sourceState for c in components.modelcomponents]
    r, couplingStrength, noiseparams = prepare_rates_coupled(rin, sourceStates, transitions, R, S, insertstep, [4, 4])
    num_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
    mean_intensity = Vector[]
    max_intensity = Float64[]
    for i in eachindex(noiseparams)
        mi = mean.(probfn(noiseparams[i], num_per_state[i], components.N))
        mmax = max(maximum(mi), 1.0)
        push!(max_intensity, mmax)
        push!(mean_intensity, mi / mmax)
    end
    a, p0 = make_ap(r, couplingStrength, interval, components)
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

########################
# Old functions
########################
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

# function predicted_statepath(trace, interval, model::AbstractGeneTransitionModel)
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
    a, p0 = make_ap(r[:, end], interval, components)
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

function prepare_rates_coupled(rates, sourceStates, transitions, R::Tuple, S, insertstep, n_noise)
    r = Vector{Float64}[]
    noiseparams = Vector{Float64}[]
    couplingStrength = Float64[]
    j = 1
    for i in eachindex(R)
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + n_noise[i]
        push!(r, rates[j:j+n-1])
        j += n
    end
    for i in eachindex(R)
        s = sourceStates[i]
        if (s isa Integer && s > 0) || (s isa Vector && !isempty(s))
            push!(couplingStrength, rates[j])
            j += 1
        else
            push!(couplingStrength, 0.0)
        end
    end
    for i in eachindex(r)
        push!(noiseparams, r[i][end-n_noise[i]+1:end])
    end
    return r, couplingStrength, noiseparams
end
function predicted_states(rates::Vector, coupling, transitions, G::Tuple, R, S, insertstep, components, n_noise, reporters_per_state, probfn, interval, traces)
    sourceStates = coupling[3]
    r, couplingStrength, noiseparams = prepare_rates_coupled(rates, sourceStates, transitions, R, S, insertstep, n_noise)
    # r, couplingStrength, noiseparams = prepare_rates_coupled(rates, nrates, reporter, couplingindices)
    nT = components.N
    a, p0 = make_ap(r, couplingStrength, interval, components)
    states = Array[]
    d = []
    for i in eachindex(noiseparams)
        push!(d, probfn[i](noiseparams[i], reporters_per_state[i], nT))
    end
    for t in traces
        T = size(t, 1)
        b = set_b(t, noiseparams, reporters_per_state, probfn, nT)
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

function predicted_states(rates::Tuple, coupling, transitions, G::Tuple, R, S, insertstep, components, n_noise, reporters_per_state, probfn, interval, traces)
    sourceStates = coupling[3]
    nT = components.N
    rshared, noiseparams, couplingStrength = rates
    a, p0 = make_ap(rshared[1], couplingStrength, interval, components)
    states = Array[]
    units = Vector[]
    observation_dist = Vector[]
    for i in eachindex(traces)
        T = size(traces[i], 1)
        b = set_b(traces[i], noiseparams[i], reporters_per_state, probfn, nT)
        spath = viterbi(a, b, p0, nT, T)
        push!(states, spath)
        d = set_d(noiseparams[i], reporters_per_state, probfn, nT)
        push!(units, [unit_state(i, G, R, S, coupling[1]) for i in spath])
        push!(observation_dist, [[d[i] for d in d] for i in spath])
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

########################
# New functions
########################
#### Return states and d (not observation_dist)

"""
    predict_trace(a::Matrix, p0::Vector, d, traces)

"""
function predict_trace(a::Matrix, p0::Vector, d, traces)
    states = Vector{Int}[]
    # observation_dist = Vector[]
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        push!(states, viterbi(a, b, p0))
        # spath = viterbi(a, b, p0)
        # push!(states, [unit_state(i, 3, 3, S, coupling[1]) for i in spath])
        # push!(observation_dist, [[d[i] for d in d] for i in spath])
        # push!(observation_dist, [d[s] for s in spath])
    end
    # states, observation_dist
    states, d
end

"""
    predict_trace(noiseparams, a::Matrix, p0::Vector, reporter, traces)

"""
function predict_trace(noiseparams::Vector, a::Matrix, p0::Vector, reporter, traces)
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        di = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], di)
        push!(states, viterbi(a, b, p0))
        push!(d, di)
    end
    states, d
end

"""
    predict_trace(r::Vector, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function predict_trace(r::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], interval, components, method)
        di = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], di)
        push!(states, viterbi(a, b, p0))
        push!(d, di)
    end
    states, d
end

"""
    predict_trace(r, couplingStrength, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function predict_trace(r::Vector, couplingStrength::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], couplingStrength[i], interval, components, method)
        di = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], di)
        push!(states, viterbi(a, b, p0))
        push!(d, di)
    end
    states, d
end

### grid trait likelihoods

"""
    predict_trace(a, a_grid, p0::Vector, d, traces)

"""
function predict_trace(a::Matrix, a_grid::Matrix, p0::Vector, d, traces)
    states = Vector{Int}[]
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        push!(states, viterbi_grid(a, a_grid, b, p0))
    end
    states, d
end

"""
    predict_trace(noiseparams, a, a_grid, p0::Vector, reporter, traces)

"""
function predict_trace(noiseparams::Vector, Ngrid::Int, a::Matrix, a_grid::Matrix, p0::Vector, reporter, traces)
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        push!(states, viterbi_grid(a, a_grid, b, p0))
        push!(d, d)
    end
    states, d
end

"""
    predict_trace(r::Vector, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function predict_trace(r::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], interval, components, method)
        a_grid = make_a_grid(pgrid, Ngrid)
        di = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], di)
        push!(states, viterbi_grid(a, a_grid, b, p0))
        push!(d, di)
    end
    states, d
end

"""
    predict_trace(r, couplingStrength, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function predict_trace(r::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], couplingStrength[i], interval, components, method)
        a_grid = make_a_grid(pgrid, Ngrid)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        push!(states, viterbi_grid(a, a_grid, b, p0))
    end
    states, d
end



"""
    predict_trace(r, components, reporter, interval, trace, method)

"""
# no traits
function predict_trace(r::Tuple{T1,T2}, components::TComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2}
    r, noiseparams = r
    a, p0 = make_ap(r, interval, components, method)
    d = set_d(noiseparams, reporter)
    predict_trace(a, p0, d, trace[1])
end

# coupled
function predict_trace(r::Tuple{T1,T2,T3}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3}
    rates, noiseparams, couplingStrength = r
    a, p0 = make_ap(rates, couplingStrength, interval, components, method)
    d = set_d(noiseparams, reporter)
    predict_trace(a, p0, d, trace[1])
end

# hierarchical
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6}, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6}
    rshared, rindividual, noiseshared, noiseindividual, _, _ = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    d = set_d(noiseshared[1], reporter)
    if method[2]
        states, observation_dist = predict_trace(noiseindividual, a, p0, reporter, trace[1])
    else
        states, observation_dist = predict_trace(rindividual, noiseindividual, interval, components, reporter, trace[1], method[1])
    end
    states, observation_dist
end

# coupled, hierarchical
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, _, _, couplingshared, couplingindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    d = set_d(noiseshared[1], reporter)
    if method[2]
        states, observation_dist = predict_trace(noiseindividual, a, p0, reporter, trace[1])
    else
        states, observation_dist = predict_trace(rindividual, couplingindividual, noiseindividual, interval, components, reporter, trace[1])
    end
    states, observation_dist
end

### grid trait likelihoods
"""
    predict_trace(r, pgrid, Ngrid, components, reporter, interval, trace, method)
    
"""
# grid
function predict_trace(r::Tuple{T1,T2,T3}, Ngrid::Int, components::TComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2,T3}
    r, noiseparams, pgrid = r
    a, p0 = make_ap(r, interval, components, method)
    a_grid = make_a_grid(pgrid[1], Ngrid)
    d = set_d(noiseparams, reporter, Ngrid)
    predict_trace(a, a_grid, p0, d, trace[1])
end

# coupled, grid
function predict_trace(r::Tuple{T1,T2,T3,T4}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3,T4}
    r, noiseparams, couplingStrength, pgrid = r
    a, p0 = make_ap(r, couplingStrength, interval, components, method)
    a_grid = make_a_grid(pgrid[1][1], Ngrid)
    d = set_d(noiseparams, reporter, Ngrid)
    predict_trace(a, a_grid, p0, d, trace[1])
end

# hierarchical, grid
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, Ngrid::Int, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, _, noiseindividual, _, _, pgridshared, pgridindividual = r
    if method[2]
        a, p0 = make_ap(rshared[1], interval, components, method[1])
        a_grid = make_a_grid(pgridshared[1][1], Ngrid)
        states, observation_dist = predict_trace(noiseindividual, Ngrid, a, a_grid, p0, reporter, trace[1])
    else
        states, observation_dist = predict_trace(rindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1], method[1])
    end
    states, observation_dist
end

# coupled, hierarchical, grid
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    rshared, rindividual, _, noiseindividual, _, _, couplingshared, couplingindividual, pgridshared, pgridindividual = r
    if method[2]
        a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
        a_grid = make_a_grid(pgridshared[1], Ngrid)
        states, observation_dist = predict_trace(noiseindividual, a, a_grid, p0, reporter, trace[1])
    else
        states, observation_dist = predict_trace(rindividual, couplingindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1], method[1])
    end
    states, observation_dist
end


function predict_trace(param, data, model::AbstractGRSMmodel)
    r = prepare_rates(param, model)
    if hastrait(model, :grid)
        predict_trace(r, model.trait.grid.ngrid, model.components, model.reporter, data.interval, data.trace, model.method)
    else
        predict_trace(r, model.components, model.reporter, data.interval, data.trace, model.method)
    end
end
