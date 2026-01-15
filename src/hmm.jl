# This file is part of StochasticGene.jl   

### hmm.jl
### Fit discrete HMMs and continuous hidden Markov process models directly to observations (e.g. intensity traces)
###
### Notation in discrete HMM algorithms follows Rabier, 1989
###
### Functions for forward, backward, and Viterbi HMM algorihms
### For continuous processes, numerically solve forward Kolmogorov equation to obtain transition probability matrix
###

# using CUDA
# using CUDA: @allowscalar

"""
    fkf(u::Matrix, p, t)

update of u of Kolmogorov forward equation for DifferentialEquations.jl
"""
function fkf(u::Matrix, p, t)
    u * p
end

"""
    fkf!(du,u::Matrix, p, t)

in place update of du of Kolmogorov forward equation for DifferentialEquations.jl

"""
function fkf!(du, u::Matrix, p, t)
    du .= u * p
end

"""
    fkb(u::Matrix, p, t)

update of u of Kolmogorov backward equation for DifferentialEquations.jl
"""
function fkb(u::Matrix, p, t)
    -u * p
end

"""
    fkb!(du,u::Matrix, p, t)

in place update of du of Kolmogorov backward equation for DifferentialEquations.jl
"""
function fkb!(du, u::Matrix, p, t)
    du .= -u * p
end

"""
    kolmogorov_forward_ad(Q, interval, method=Tsit5(), save=false)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_forward_ad(Q, interval, method=Tsit5(), save=false)
    tspan = (zero(eltype(Q)), interval)
    u0 = Matrix{eltype(Q)}(I, size(Q)...)
    prob = ODEProblem(fkf, u0, tspan, Q)
    sol = solve(prob, method; save_everystep=save)
    return sol.u[end]
end

"""
    kolmogorov_forward_inplace(Q, interval, method=Tsit5(), save=false)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_forward(Q, interval, method=Tsit5(), save=false)
    tspan = (0.0, interval)
    prob = ODEProblem(fkf!, Matrix(I, size(Q)), tspan, Q)
    sol = solve(prob, method, save_everystep=save)
    sol.u[end]
end


"""
    kolmogorov_backward(Q, interval, method=Tsit5(), save=false)

return the solution of the Kolmogorov backward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_backward_ad(Q, interval, method=Tsit5(), save=false)
    tspan = (0.0, interval)
    u0 = Matrix{eltype(Q)}(I, size(Q)...)
    prob = ODEProblem(fkb, u0, tspan, Q)
    sol = solve(prob, method, save_everystep=save)
    sol.u[end]
end


"""
    kolmogorov_backward_inplace(Q, interval, method=Tsit5(), save=false)

return the solution of the Kolmogorov forward equation 
returns initial condition and solution at time = interval

- `Q`: transition rate matrix
- `interval`: interval between frames (total integration time)
"""
function kolmogorov_backward(Q, interval, method=Tsit5(), save=false)
    tspan = (0.0, interval)
    prob = ODEProblem(fkb!, Matrix(I, size(Q)), tspan, Q)
    sol = solve(prob, method, save_everystep=save)
    sol.u[end]
end


"""
    prob_Gaussian(par, reporters::Int)

Create a Gaussian distribution for a single reporter count.

# Arguments
- `par`: 4-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std]
- `reporters::Int`: Number of reporters

# Returns
- `Normal`: Gaussian distribution with mean = background_mean + reporters * reporter_mean and 
  variance = background_std² + reporters * reporter_std²
"""
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

"""
    prob_Gaussian_ad(par, reporters_per_state::T, N) where {T<:Vector}

Create an array of Gaussian distributions using automatic differentiation-friendly approach.

# Arguments
- `par`: 4-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state
- `N`: Number of HMM states

# Returns
- `Vector{Normal}`: Array of Gaussian distributions, one for each state
"""
function prob_Gaussian_ad(par, reporters_per_state::T, N) where {T<:Vector}
    [prob_Gaussian(par, reporters_per_state[i]) for i in 1:N]
end

"""
    prob_Gaussian(par, reporters_per_state::T) where {T<:Vector}

Create an array of Gaussian distributions for multiple states.

# Arguments
- `par`: 4-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state

# Returns
- `Vector{Normal}`: Array of Gaussian distributions, one for each state
"""
function prob_Gaussian(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    prob_Gaussian(par, reporters_per_state, N)
    # d = Array{Distribution{Univariate,Continuous}}(undef, N)
    # for i in 1:N
    #     d[i] = prob_Gaussian(par, reporters_per_state[i])
    # end
    # d
end

function prob_Gaussian_ad(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    [prob_Gaussian_ad(par, reporters_per_state[i]) for i in 1:N]
end

"""
    prob_GaussianMixture(par, reporters::Int)

Create a Gaussian mixture distribution for a single reporter count.

# Arguments
- `par`: 5-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std, weight]
- `reporters::Int`: Number of reporters

# Returns
- `MixtureModel`: Gaussian mixture with two components:
  1. Background + reporters: mean = background_mean + reporters * reporter_mean, 
     std = sqrt(background_std² + reporters * reporter_std²)
  2. Background only: mean = background_mean, std = background_std
  Weight of first component = par[5], weight of second component = 1 - par[5]
"""
function prob_GaussianMixture(par, reporters::Int)
    MixtureModel(Normal, [(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2)), (par[1], par[2])], [par[5], 1 - par[5]])
end

"""
    prob_GaussianMixture(par, reporters_per_state::T, N) where {T<:Vector}

Create an array of Gaussian mixture distributions for multiple states.

# Arguments
- `par`: 5-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std, weight]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state
- `N`: Number of HMM states

# Returns
- `Vector{MixtureModel}`: Array of Gaussian mixture distributions, one for each state
"""
function prob_GaussianMixture(par, reporters_per_state::T, N) where {T<:Vector}
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture(par, reporters_per_state[i])
    end
    d
end

"""
    prob_GaussianMixture(par, reporters_per_state::T) where {T<:Vector}

Create an array of Gaussian mixture distributions for multiple states.

# Arguments
- `par`: 5-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std, weight]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state

# Returns
- `Vector{MixtureModel}`: Array of Gaussian mixture distributions, one for each state
"""
function prob_GaussianMixture(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture(par, reporters_per_state[i])
    end
    d
end

"""
    prob_GaussianMixture_6(par, reporters)

Create a Gaussian mixture distribution with 6 parameters for a single reporter count.

# Arguments
- `par`: 7-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std, 
         second_mean, second_std, weight]
- `reporters::Int`: Number of reporters

# Returns
- `MixtureModel`: Gaussian mixture with two components:
  1. Background + reporters: mean = background_mean + reporters * reporter_mean, 
     std = sqrt(background_std² + reporters * reporter_std²)
  2. Second component: mean = par[5], std = par[6]
  Weight of first component = par[7], weight of second component = 1 - par[7]
"""
function prob_GaussianMixture_6(par, reporters)
    MixtureModel(Normal, [(par[1] + reporters * par[3], sqrt(par[2]^2 + reporters * par[4]^2)), (par[5], par[6])], [par[7], 1 - par[7]])
end

"""
    prob_GaussianMixture_6(par, reporters_per_state, N)

Create an array of Gaussian mixture distributions with 6 parameters for multiple states.

# Arguments
- `par`: 7-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std, 
         second_mean, second_std, weight]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state
- `N`: Number of HMM states

# Returns
- `Vector{MixtureModel}`: Array of Gaussian mixture distributions, one for each state
"""
function prob_GaussianMixture_6(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture_6(par, reporters_per_state[i])
    end
    d
end

"""
    prob_GaussianMixture_6(par, reporters_per_state::T) where {T<:Vector}

Create an array of Gaussian mixture distributions with 6 parameters for multiple states.

# Arguments
- `par`: 7-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std, 
         second_mean, second_std, weight]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state

# Returns
- `Vector{MixtureModel}`: Array of Gaussian mixture distributions, one for each state
"""
function prob_GaussianMixture_6(par, reporters_per_state::T) where {T<:Vector}
    N = length(reporters_per_state)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_GaussianMixture_6(par, reporters_per_state[i])
    end
    d
end


"""
    prob_Gaussian_ind(par, reporters::Int)

Create an independent Gaussian distribution for a single reporter count.

# Arguments
- `par`: 4-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std]
- `reporters::Int`: Number of reporters

# Returns
- `Normal`: Gaussian distribution
  If reporters > 0: mean = reporters * reporter_mean, std = sqrt(reporters) * reporter_std
  If reporters = 0: mean = background_mean, std = background_std
"""
function prob_Gaussian_ind(par, reporters::Int)
    if reporters > 0
        return Normal(reporters * par[3], sqrt(reporters) * par[4])
    else
        return Normal(par[1], par[2])
    end
end

"""
    prob_Gaussian_ind(par, reporters_per_state, N)

Create an array of independent Gaussian distributions for multiple states.

# Arguments
- `par`: 4-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state
- `N`: Number of HMM states

# Returns
- `Vector{Normal}`: Array of Gaussian distributions, one for each state
  If reporters > 0: mean = reporters * reporter_mean, std = sqrt(reporters) * reporter_std
  If reporters = 0: mean = background_mean, std = background_std
"""
function prob_Gaussian_ind(par, reporters_per_state, N)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = prob_Gaussian_ind(par, reporters_per_state[i])
    end
    d
end

"""
    prob_Gaussian_ind(par, reporters_per_state::T) where {T<:Vector}

Create an array of independent Gaussian distributions for multiple states.

# Arguments
- `par`: 4-dimensional vector of parameters [background_mean, background_std, reporter_mean, reporter_std]
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state

# Returns
- `Vector{Normal}`: Array of Gaussian distributions, one for each state
  If reporters > 0: mean = reporters * reporter_mean, std = sqrt(reporters) * reporter_std
  If reporters = 0: mean = background_mean, std = background_std
"""
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
    for l in 1:Nstate
        for k in 1:Ngrid
            for j in 1:Ngrid
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

# function make_ap(rates, interval, components::TComponents, method=Tsit5())
#     Qtr = make_mat_T(components, rates) ##  transpose of the Markov process transition rate matrix Q
#     kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
# end

function make_ap(rates::Vector, interval, components::TComponents, method=Tsit5())
    Qtr = make_mat_T(components, rates) ##  transpose of the Markov process transition rate matrix Q
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

function make_ap(rates, couplingStrength, interval, components::TCoupledComponents, method=Tsit5())
    Qtr = make_mat_TC(components, rates, couplingStrength)
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

function make_ap(r::Tuple{<:Vector}, interval, components::TComponents, method=Tsit5())
    r, couplingStrength = r
    Qtr = make_mat_TC(components, r, couplingStrength)
    kolmogorov_forward(Qtr', interval, method), normalized_nullspace(Qtr)
end

function make_ap(rates, couplingStrength, interval, components::TForcedComponents, method=Tsit5())
    r = set_rates_forced(rates,couplingStrength,components)
    Q = [make_mat_T(components, r[i]) for i in 1:2]
    a = [kolmogorov_forward(Qtr', interval, method) for Qtr in Q]
    p0 = [normalized_nullspace(Qtr) for Qtr in Q]
    return a, p0
end

function set_rates_forced(rates, couplingStrength, components::TForcedComponents)
    r = [deepcopy(rates) for _ in 1:2]  # Create two independent copies
    r[2][components.targets[2]] *= (1. + couplingStrength[1])
    return r
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

Create a grid transition probability matrix using Gaussian kernel.

# Arguments
- `param`: Standard deviation parameter for the Gaussian kernel
- `Ngrid`: Number of grid points

# Returns
- `Matrix{Float64}`: Normalized transition probability matrix where each row sums to 1
  The transition probability from position i to j is proportional to exp(-distance(i,j)²/(2*param²))
"""
function make_a_grid(param, Ngrid)
    as = zeros(Ngrid, Ngrid)
    # d = zeros(Ngrid, Ngrid)
    for j in 1:Ngrid
        for i in 1:Ngrid
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
    set_d_background(noiseparams, sigma2, reporters_per_state, probfn, N)

Create emission distributions with additional background noise variance.

# Arguments
- `noiseparams`: Vector of noise parameters [background_mean, background_std, reporter_mean, reporter_std, ...]
- `sigma2`: Additional variance to add to background noise
- `reporters_per_state`: Number of reporters per HMM state
- `probfn`: Function to create probability distributions
- `N`: Number of HMM states

# Returns
- `Vector{Distribution}`: Array of emission distributions with modified background variance
  The background standard deviation is modified as sqrt(background_std² + sigma2²)
"""
function set_d_background(noiseparams, sigma2, reporters_per_state, probfn, N)
    n = copy(noiseparams)
    n[2] = sqrt(n[2]^2 + sigma2^2)
    probfn(n, reporters_per_state, N)
end

function set_d_background(noiseparams, sigma2::Matrix, reporters_per_state, probfn, N)
    n = copy(noiseparams)
    n[2] = sqrt(n[2]^2 + sigma2[1,2]^2)
    probfn(n, reporters_per_state, N)
end

"""
    set_d(noiseparams, reporters_per_state, probfn, N)

Create emission distributions for HMM states.

# Arguments
- `noiseparams`: Vector of noise parameters
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state
- `probfn::Function`: Function to create probability distributions
- `N::Int`: Number of HMM states

# Returns
- `Vector{Distribution}`: Array of emission distributions, one for each state
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

function set_d_ad(noiseparams::Vector{T}, reporters_per_state::Vector{Vector{Int}}, probfn::Vector, N::Int) where {T<:AbstractVector}
    [probfn[i](noiseparams[i], reporters_per_state[i], N) for i in eachindex(noiseparams)]
end


"""
    set_d(noiseparams, reporter, N)

Create emission distributions using HMMReporter structure.

# Arguments
- `noiseparams`: Vector of noise parameters
- `reporter::HMMReporter`: Reporter structure containing per_state and probfn
- `N::Int`: Number of HMM states

# Returns
- `Vector{Distribution}`: Array of emission distributions, one for each state
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

Create emission distributions for HMM states (auto-determine N).

# Arguments
- `noiseparams`: Vector of noise parameters
- `reporters_per_state::Vector{Int}`: Number of reporters per HMM state
- `probfn::Function`: Function to create probability distributions

# Returns
- `Vector{Distribution}`: Array of emission distributions, one for each state
"""
function set_d(noiseparams, reporters_per_state::Vector{Int}, probfn::T) where {T<:Function}
    probfn(noiseparams, reporters_per_state)
end

function set_d(noiseparams::Vector{T1}, reporters_per_state::Vector{Vector{Int}}, probfn::Vector{T2}) where {T1<:AbstractVector, T2<:Function}
    d = Vector{Distribution}[]
    for i in eachindex(noiseparams)
        push!(d, probfn[i](noiseparams[i], reporters_per_state[i]))
    end
    return d
end

function set_d_ad(noiseparams::Vector{T}, reporters_per_state::Vector{Vector{Int}}, probfn::Vector) where {T<:AbstractVector}
    [probfn[i](noiseparams[i], reporters_per_state[i]) for i in eachindex(noiseparams)]
end

"""
    set_d(noiseparams, reporter)

Create emission distributions using HMMReporter structure (auto-determine N).

# Arguments
- `noiseparams`: Vector of noise parameters
- `reporter::HMMReporter`: Reporter structure containing per_state and probfn

# Returns
- `Vector{Distribution}`: Array of emission distributions, one for each state
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

Create emission probability matrix from observations and distributions.

# Arguments
- `trace::Vector`: Vector of observations
- `d::Vector{Distribution}`: Array of emission distributions, one for each state

# Returns
- `Matrix{Float64}`: Emission probability matrix b[j,t] = P(observation_t | state_j)
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

function set_b(trace::Matrix, d::Tuple)
    [trace[:,1], set_b(trace[:,2], d[2])]
end

function set_b_ad(trace::Vector, d::Vector{T}) where {T<:Distribution}
    N = length(d)
    Tlen = length(trace)
    [pdf(d[j], trace[t]) for j in 1:N, t in 1:Tlen]
end

function set_b(trace, d::Vector{T}) where {T<:Vector}
    N = length(d[1])
    b = ones(N, size(trace, 1))
    t = 1
    for obs in eachrow(trace)
        for i in eachindex(d)
            for j in 1:N
                b[j, t] *= pdf(d[i][j], obs[i])
            end
        end
        t += 1
    end
    return b
end

function set_b_ad(trace, d::Vector{T}) where {T<:Vector}
    N = length(d[1])
    Tlen = size(trace, 1)
    [prod(pdf(d[i][j], trace[t, i]) for i in eachindex(d)) for j in 1:N, t in 1:Tlen]
end

function set_b(trace, d::Array{T,3}) where {T<:Distribution}
    Nstate, Ngrid, _ = size(d)
    b = ones(Nstate, Ngrid, size(trace, 2))
    t = 1
    for obs in eachcol(trace)
        for l in 1:Ngrid
            for k in 1:Ngrid
                for j in 1:Nstate
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
        for l in 1:Ngrid
            for k in 1:Ngrid
                for j in 1:Nstate
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

function set_logb_coupled_inplace(trace, params, reporter, N)
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
    set_logb_coupled(trace, params, reporter, N)

    returns log of matrix b
"""
function set_logb_coupled(trace, params, reporter, N)
    # Build d as a vector of matrices, one for each param/reporter
    d = [reporter[i].probfn(params[i], reporter[i].per_state, N) for i in eachindex(params)]
    # Compute logb as a matrix using comprehensions
    logb = [sum(logpdf(d[i][j, 1], obs[i]) for i in eachindex(d)) for j in 1:N, obs in eachrow(trace)]
    return logb
end

function set_b_background_inplace(obs, d::Vector{<:Distribution})
    b = Array{Float64}(undef, size(d))
    for j in CartesianIndices(d)
        b[j] = pdf(d[j], obs)
    end
    return reshape(b, :, 1)
end

function set_b_background_inplace(obs, d::Vector{<:Vector}, k::Int, N)
    b = ones(N)
    for j in 1:N
        b[j] *= pdf(d[k][j], obs)
    end
    return reshape(b, :, 1)
end

function set_b_background(obs, d::Vector{<:Distribution})
    b = [pdf(dj, obs) for dj in d]
    return reshape(b, :, 1)
end

function set_b_background(obs, d::Vector{<:Vector}, k::Int, N)
    b = [pdf(d[k][j], obs) for j in 1:N]
    return reshape(b, :, 1)
end

"""
    forward_inner_operation!(α, a, b::Vector, i, j, t)

Update forward variable α[j,t] using the forward algorithm recursion.

# Arguments
- `α`: Forward variable matrix
- `a::Matrix`: Transition probability matrix
- `b::Vector`: Emission probabilities for current observation
- `i, j`: State indices
- `t`: Time index

# Details
Performs the update: α[j,t] += α[i,t-1] * a[i,j] * b[j]
This is the core operation in the forward algorithm recursion.
"""
function forward_inner_operation!(α, a::Matrix, b::Vector, i, j, t)
    α[j, t] += α[i, t-1] * a[i, j] * b[j]
end

function forward_inner_operation_ad(α, a, b::Vector, i, j, t)
    α_new = similar(α)
    α_new[j, t] += α[i, t-1] * a[i, j] * b[j]
    return α_new
end

function forward_inner_operation!(α, a::Matrix, b::Matrix, i, j, t)
    α[j, t] += α[i, t-1] * a[i, j] * b[j, t]
end

function forward_inner_operation_ad(α, a::Matrix, b::Matrix, i, j, t)
    α_new = similar(α)
    α_new[j, t] += α[i, t-1] * a[i, j] * b[j, t]
    return α_new
end

function forward_inner_operation!(α, a::Vector{T1}, b::Vector{T2}, i, j, t) where {T1 <: AbstractArray, T2 <: AbstractArray}
    if ~b[1][j]
        forward_inner_operation!(α, a[1], b[2], i, j, t)
    else
        forward_inner_operation!(α, a[2], b[2], i, j, t)
    end
end


"""
forward(a, b, p0, N, T)

returns forward variable α, and scaling parameter array C using scaled forward algorithm
α[i,t] = P(O1,...,OT,qT=Si,λ)
Ct = Prod_t 1/∑_i α[i,t]

# """
function forward(a, b, p0, N, T)
    # if CUDA.functional() && (N * N * T > 1000)
    #     return forward_gpu(a, b, p0, N, T)
    # else
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    α[:, 1] = p0 .* b[:, 1]
    C[1] = 1 / max(sum(α[:, 1]), eps(Float64))
    α[:, 1] *= C[1]
    for t in 2:T
        for j in 1:N
            for i in 1:N
                forward_inner_operation!(α, a, b, i, j, t)
            end
        end
        C[t] = 1 / max(sum(α[:, t]), eps(Float64))
        α[:, t] *= C[t]
    end
    return α, C
    # end
end

function forward(a::Vector{T1}, b::Vector{T2}, p0, N, T) where {T1 <: AbstractArray, T2 <: AbstractArray}
    # if CUDA.functional() && (N * N * T > 1000)
    #     return forward_gpu(a, b, p0, N, T)
    # else
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    if ~b[1][1]
        α[:, 1] = p0[1] .* b[2][:, 1]
    else
        α[:, 1] = p0[2] .* b[2][:, 1]
    end
    C[1] = 1 / max(sum(α[:, 1]), eps(Float64))
    α[:, 1] *= C[1]
    for t in 2:T
        for j in 1:N
            for i in 1:N
                forward_inner_operation!(α, a, b, i, j, t)
            end
        end
        C[t] = 1 / max(sum(α[:, t]), eps(Float64))
        α[:, t] *= C[t]
    end
    return α, C
    # end
end

function forward(a::Vector{T1}, b::Vector{T2}, p0) where {T1 <: AbstractArray, T2 <: AbstractArray}
    N, T = size(b[2])
    forward(a, b, p0, N, T)
end

function forward_ad(a::Matrix, b, p0, N, T)
    b = max.(b, eps(Float64))
    function step(α_prev, t)
        α_new = [sum(α_prev[i] * a[i, j] * b[j, t] for i in 1:N) for j in 1:N]
        st = sum(α_new)
        if st == 0.0
            α_new = fill(1.0 / N, N)
            st = 1.0
        end
        c = 1 / max(st, eps(Float64))
        α_new = α_new .* c
        return α_new, c
    end
    # Initial step
    α1 = p0 .* b[:, 1]
    s1 = sum(α1)
    if s1 == 0.0
        α1 = fill(1.0 / N, N)
        s1 = 1.0
    end
    c1 = 1 / s1
    α1 = α1 .* c1

    # Recurrence using foldl to avoid mutation
    function recur((αs, cs, α_prev), t)
        α_new, c = step(α_prev, t)
        (vcat(αs, [α_new]), vcat(cs, [c]), α_new)
    end
    αs, cs, _ = foldl(recur, 2:T; init=([α1], [c1], α1))
    α = hcat(αs...)
    return α, cs
end

"""
    forward(a, b, p0)

Compute forward variables using scaled forward algorithm (auto-determine dimensions).

# Arguments
- `a::Matrix`: Transition probability matrix
- `b::Matrix`: Emission probability matrix
- `p0`: Initial state distribution

# Returns
- `Tuple{Matrix{Float64}, Vector{Float64}}`: (α, C) where α is the forward variable matrix and C is the scaling parameter array
"""
function forward(a::Matrix, b::Matrix, p0)
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
# function forward_gpu(a::Matrix{Float64}, b::Matrix{Float64}, p0::Vector{Float64}, N::Int64, T::Int64)
#     # Move data to GPU
#     a_gpu = CuArray(a)
#     b_gpu = CuArray(b)
#     p0_gpu = CuArray(p0)

#     # Allocate GPU arrays for α and C
#     α_gpu = CUDA.zeros(Float64, N, T)
#     C_gpu = CUDA.zeros(Float64, T)

#     # Initialize first time step
#     # Use proper GPU array operations instead of scalar indexing
#     α_gpu[:, 1] .= p0_gpu .* b_gpu[:, 1]

#     # Use CUDA.@allowscalar for operations that require scalar indexing
#     CUDA.@allowscalar C_gpu[1] = 1 / sum(α_gpu[:, 1])
#     CUDA.@allowscalar α_gpu[:, 1] *= C_gpu[1]

#     # Compute forward probabilities using matrix multiplication on GPU
#     for t in 2:T
#         # Use proper GPU array operations
#         α_gpu[:, t] .= (a_gpu * α_gpu[:, t-1]) .* b_gpu[:, t]

#         # Use CUDA.@allowscalar for operations that require scalar indexing
#         CUDA.@allowscalar C_gpu[t] = 1 / sum(α_gpu[:, t])
#         CUDA.@allowscalar α_gpu[:, t] *= C_gpu[t]
#     end

#     # Return results back to CPU
#     return Array(α_gpu), Array(C_gpu)
# end

"""
    forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)

Compute forward variables for grid-based HMM using scaled forward algorithm.

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix of size (Nstate, Ngrid, T)
- `p0`: Initial state distribution
- `Nstate`: Number of states
- `Ngrid`: Number of grid points
- `T`: Number of time steps

# Returns
- `Tuple{Array{Float64,3}, Vector{Float64}}`: (α, C) where α is the forward variable array of size (Nstate, Ngrid, T) and C is the scaling parameter array
"""
function forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    # if CUDA.functional() && (Nstate * Nstate * Ngrid * Ngrid * T > 1000)
    #     return forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)
    # else
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
    # end
end

"""
    forward_grid_ad(a, a_grid, b, p0, Nstate, Ngrid, T)

Compute forward variables for grid-based HMM using automatic differentiation-friendly approach.

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix of size (Nstate, Ngrid, T)
- `p0`: Initial state distribution
- `Nstate`: Number of states
- `Ngrid`: Number of grid points
- `T`: Number of time steps

# Returns
- `Tuple{Array{Float64,3}, Vector{Float64}}`: (α, C) where α is the forward variable array and C is the scaling parameter array
"""
function forward_grid_ad(a, a_grid, b, p0, Nstate, Ngrid, T)
    αs = Vector{Matrix{Float64}}(undef, T)
    C = Vector{Float64}(undef, T)
    αs[1] = p0 .* b[:, :, 1]
    C[1] = 1 / sum(αs[1])
    αs[1] *= C[1]
    for t in 2:T
        α_new = zeros(Nstate, Ngrid)
        for l in 1:Ngrid, k in 1:Nstate
            for j in 1:Ngrid, i in 1:Nstate
                α_new[i, j] += αs[t-1][k, l] * a[k, i] * a_grid[l, j] * b[i, j, t]
            end
        end
        C[t] = 1 / sum(α_new)
        αs[t] = α_new * C[t]
    end
    α = cat(αs..., dims=3)  # Stack along the 3rd dimension
    return α, C
end

"""
    forward_grid(a, a_grid, b, p0)

Compute forward variables for grid-based HMM (auto-determine dimensions).

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix of size (Nstate, Ngrid, T)
- `p0`: Initial state distribution

# Returns
- `Tuple{Array{Float64,3}, Vector{Float64}}`: (α, C) where α is the forward variable array and C is the scaling parameter array
"""
function forward_grid(a, a_grid, b, p0)
    Nstate, Ngrid, T = size(b)
    forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
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
# function forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)
#     # Move data to GPU
#     a_gpu = CuArray(a)
#     a_grid_gpu = CuArray(a_grid)
#     b_gpu = CuArray(b)
#     p0_gpu = CuArray(p0)

#     # Allocate GPU arrays for results
#     α_gpu = CuArray{Float64}(undef, Nstate, Ngrid, T)
#     C_gpu = CuArray{Float64}(undef, T)

#     # Initialize first time step
#     α_gpu[:, :, 1] .= p0_gpu .* b_gpu[:, :, 1]
#     C_gpu[1] = sum(α_gpu[:, :, 1])
#     α_gpu[:, :, 1] ./= C_gpu[1]

#     # Compute forward probabilities using matrix multiplication on GPU
#     for t in 2:T
#         # Reshape for efficient matrix multiplication
#         α_prev = reshape(α_gpu[:, :, t-1], Nstate * Ngrid, 1)
#         a_combined = kron(a_grid_gpu, a_gpu)  # Kronecker product for combined transitions
#         α_temp = reshape(a_combined * α_prev, Nstate, Ngrid)
#         α_gpu[:, :, t] .= α_temp .* b_gpu[:, :, t]

#         C_gpu[t] = sum(α_gpu[:, :, t])
#         α_gpu[:, :, t] ./= C_gpu[t]
#     end

#     # Return results back to CPU
#     return Array(α_gpu), Array(C_gpu)
# end

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
# function forward_grid_gpu(a, a_grid, b, p0)
#     Nstate, Ngrid, T = size(b)
#     forward_grid_gpu(a, a_grid, b, p0, Nstate, Ngrid, T)
# end


"""
    forward_matrixmult_inplace(a, b, p0, N, T)

Compute forward variables using matrix multiplication (in-place version).

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Tuple{Matrix{Float64}, Vector{Float64}}`: (α, C) where α is the forward variable matrix and C is the scaling parameter array
"""
function forward_matrixmult_inplace(a, b, p0, N, T)
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
    forward_matrixmult(a, b, p0, N, T)

Compute forward variables using matrix multiplication (functional version).

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Tuple{Matrix{Float64}, Vector{Float64}}`: (α, C) where α is the forward variable matrix and C is the scaling parameter array
"""
function forward_matrixmult(a, b, p0, N, T)
    αs = Vector{Vector{Float64}}(undef, T)
    C = Vector{Float64}(undef, T)
    αs[1] = p0 .* b[:, 1]
    C[1] = 1 / sum(αs[1])
    αs[1] *= C[1]
    for t in 2:T
        α_new = αs[t-1] .* (a' * b[:, t])
        C[t] = 1 / sum(α_new)
        αs[t] = α_new * C[t]
    end
    α = hcat(αs...)  # Convert vector of vectors to matrix
    return α, C
end

"""
    forward_log(loga, logb, logp0, N, T)

Compute log forward variables using numerically stable log-space algorithm.

# Arguments
- `loga`: Log transition probability matrix
- `logb`: Log emission probability matrix
- `logp0`: Log initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Log forward variable matrix log(α) where α[i,t] = P(O_1,...,O_t, q_t = S_i | λ)
"""
function forward_log(loga, logb, logp0, N, T)
    ψ = zeros(N)
    ϕ = Matrix{Float64}(undef, N, T)
    forward_log!(ϕ, ψ, loga, logb, logp0, N, T)
    return ϕ
end
"""
    forward_log!(ϕ, ψ, loga, logb, logp0, N, T)

Compute log forward variables using numerically stable log-space algorithm (in-place).

# Arguments
- `ϕ`: Output log forward variable matrix (modified in-place)
- `ψ`: Temporary workspace vector
- `loga`: Log transition probability matrix
- `logb`: Log emission probability matrix
- `logp0`: Log initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Details
Computes log(α) where α[i,t] = P(O_1,...,O_t, q_t = S_i | λ) using logsumexp for numerical stability.
"""
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
    forward_log_ad(loga, logb, logp0, N, T)

Compute log forward variables using automatic differentiation-friendly approach.

# Arguments
- `loga`: Log transition probability matrix
- `logb`: Log emission probability matrix
- `logp0`: Log initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Log forward variable matrix log(α) where α[i,t] = P(O_1,...,O_t, q_t = S_i | λ)
"""
function forward_log_ad(loga, logb, logp0, N, T)
    ϕs = Vector{Vector{Float64}}(undef, T)
    ψ = zeros(N)
    ϕs[1] = logp0 + logb[:, 1]
    for t in 2:T
        ϕ_new = zeros(N)
        for k in 1:N
            for j in 1:N
                ψ[j] = ϕs[t-1][j] + loga[j, k] + logb[k, t]
            end
            ϕ_new[k] = logsumexp(ψ)
        end
        ϕs[t] = ϕ_new
    end
    ϕ = hcat(ϕs...)  # Convert vector of vectors to matrix
    return ϕ
end

"""
    backward_scaled(a, b, C, N, T)

Compute backward variable β using scaled backward algorithm.

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `C`: Scaling parameter array from forward algorithm
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Backward variable matrix β where β[i,t] = P(O_{t+1},...,O_T | q_t = S_i, λ)
"""
function backward_scaled(a, b, C, N, T)
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
    backward_scaled_ad(a, b, C, N, T)

Compute backward variable β using scaled backward algorithm (automatic differentiation-friendly).

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `C`: Scaling parameter array from forward algorithm
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Backward variable matrix β where β[i,t] = P(O_{t+1},...,O_T | q_t = S_i, λ)
"""
function backward_scaled_ad(a, b, C, N, T)
    βs = Vector{Vector{Float64}}(undef, T)
    βs[T] = ones(N) / C[T]
    for t in (T-1):-1:1
        β_new = zeros(N)
        for i in 1:N
            for j in 1:N
                β_new[i] += a[i, j] * b[j, t+1] * βs[t+1][j]
            end
        end
        βs[t] = β_new / C[t]
    end
    β = hcat(βs...)  # Convert vector of vectors to matrix
    return β
end

"""
    backward_log(a, b, N, T)

Compute log backward variable using numerically stable log-space algorithm.

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Log backward variable matrix log(β) where β[i,t] = P(O_{t+1},...,O_T | q_t = S_i, λ)
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
    backward_log_ad(a, b, N, T)

Compute log backward variables using automatic differentiation-friendly approach.

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Log backward variable matrix log(β) where β[i,t] = P(O_{t+1},...,O_T | q_t = S_i, λ)
"""
function backward_log_ad(a, b, N, T)
    loga = log.(a)
    ψ = zeros(N)
    ϕs = Vector{Vector{Float64}}(undef, T)
    ϕs[T] = zeros(N)  # or [0.0, 0.0] if N == 2
    for t in (T-1):-1:1
        ϕ_new = zeros(N)
        for i in 1:N
            for j in 1:N
                ψ[j] = ϕs[t+1][j] + loga[i, j] + log(b[j, t+1])
            end
            ϕ_new[i] = logsumexp(ψ)
        end
        ϕs[t] = ϕ_new
    end
    ϕ = hcat(ϕs...)  # Convert vector of vectors to matrix
    return ϕ
end

p_off(a, p0, offstates, nframes) = sum(p0[offstates]' * a[offstates, offstates]^nframes)

"""
    ll_off(a, p0, offstates, poff, nframes)

L ∝ - log P(O | r) - p_inactive/p_active log (P(off | r))
"""
# function ll_off1(a, p0, offstates, poff, nframes)
#     p = sum(p0[offstates]' * a[offstates, offstates]^nframes)
#     -(1 - poff) * log(1 - p) - poff * log(p)
# end

# function ll_off(a, p0, offstates, poff, trace)
#     nframes = length(trace)
#     p = sum(p0[offstates]' * a[offstates, offstates]^nframes)
#    - poff * log(p) * nframes
# end
# """
#     ll_off_coupled(a, p0, offstates, weight::Vector, nframes)

# TBW
# """
# function ll_off_coupled(a, p0, offstates, weight::Vector, nframes)
#     l = 0
#     for i in eachindex(weight)
#         l += ll_off(a, p0, offstates[i], weight[i], nframes)
#     end
#     l
# end

function ll_B(a, p0, offstates, weight::Float64, trace)
    p = 0
    for t in trace[1]
        nframes = length(t)
        p += log(sum(p0[offstates]' * a[offstates, offstates]^nframes)) - nframes * 0.2 * log(2 * π)
    end
    return weight * p  # Convention: return positive (weighted) log-likelihood
end
function ll_B(a, p0, offstates::Vector, weight::Vector, trace)
    l = 0
    for i in eachindex(weight)
        l += ll_B(a, p0, offstates[i], weight[i], trace)
    end
    return l  # Convention: return positive log-likelihood
end

"""
    ll_off(trace, noiseparams, reporter, components, a, p0)

Compute log-likelihood contribution from off-state observations.

# Arguments
- `trace::Tuple`: Trace data containing off-state information
- `noiseparams`: Noise parameters for emission distributions
- `reporter`: Reporter structure
- `components`: Model components
- `a::Matrix`: Transition probability matrix
- `p0`: Initial state distribution

# Returns
- `Float64`: Log-likelihood contribution from off-state observations
"""

# function ll_off(trace::Tuple, d::Vector{Distribution{Univariate,Continuous}}, a::Matrix, p0)
#     if trace[3] > 0.0
#         b = set_b_background(trace[2], d)
#         _, C = forward(a, b, p0)
#         return sum(log.(C)) * trace[4] * trace[3] * length(trace[1])
#     else
#         return 0.0
#     end
# end

function ll_off(trace::Tuple, noiseparams, reporter, components, a::Matrix, p0)
    if trace[3] > 0.0
        d = set_d_background(noiseparams, trace[5], reporter.per_state, reporter.probfn, components.nT)
        b = set_b_background(trace[2], d)
        _, C = forward(a, b, p0)
        return -sum(log.(C)) * trace[3] * trace[4] * length(trace[1])
    else
        return 0.0
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
            d = set_d_background(noiseparams[i], trace[5][i], rps, reporter[i].probfn, dims[i])
            b = set_b_background(trace[2][i], d)
            _, C = forward(a, b, p0)
            l += -sum(log.(C)) * trace[3][i] * trace[4]
        end
    end
    l * length(trace[1])
end

### Called by trait likelihoods

"""
    ll_hmm(a::Matrix, p0::Vector, d, traces)

"""
function _ll_hmm(a, p0, d, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

# function _ll_hmm(a::Vector{T1}, p0::Vector{T2}, d, traces) where {T1 <: Matrix, T2 <: Vector}
#     logpredictions = Array{Float64}(undef, length(traces))
#     for i in eachindex(traces)
#         b = set_b(traces[i], d)
#         _, C = forward(a, b, p0)
#         @inbounds logpredictions[i] = -sum(log.(C))
#     end
#     sum(logpredictions), logpredictions
# end

function _ll_hmm_ad(a::Matrix, p0::Vector, d, traces)
    logpredictions = [
        begin
            b = set_b(traces[i], d)
            _, C = forward(a, b, p0)
            -sum(log.(C))
        end
        for i in eachindex(traces)
    ]
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(noiseparams, a::Matrix, p0::Vector, reporter, traces)

"""
function _ll_hmm(noiseparams::Vector, a, p0::Vector, reporter, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

function _ll_hmm_forced(noiseparams::Vector, a, p0::Vector, reporter, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        d = (1, set_d(noiseparams[i], reporter))
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    _ll_hmm(r::Vector, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _ll_hmm(r::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], interval, components, method)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

function _ll_hmm_ad(r::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = [-sum(log.(last(forward(make_ap(r[i], interval, components, method)..., set_b(traces[i], set_d(noiseparams[i], reporter)), p0)))) for i in eachindex(traces)]
    sum(logpredictions), logpredictions
end

"""
    _ll_hmm(r, couplingStrength, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _ll_hmm(r::Vector, couplingStrength::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], couplingStrength[i], interval, components, method)
        d = set_d(noiseparams[i], reporter)
        b = set_b(traces[i], d)
        _, C = forward(a, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

### grid trait likelihoods

"""
    ll_hmm(a, a_grid, p0::Vector, d, traces)

"""
function _ll_hmm_grid(a::Matrix, a_grid::Matrix, p0::Vector, d, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        _, C = forward_grid(a, a_grid, set_b(traces[i], d), p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(noiseparams, a, a_grid, p0::Vector, reporter, traces)

"""
function _ll_hmm_grid(noiseparams::Vector, Ngrid::Int, a::Matrix, a_grid::Matrix, p0::Vector, reporter, traces)
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        _, C = forward_grid(a, a_grid, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(r::Vector, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _ll_hmm_grid(r::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], interval, components, method)
        a_grid = make_a_grid(pgrid, Ngrid)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        _, C = forward_grid(a, a_grid, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hmm(r, couplingStrength, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _ll_hmm_grid(r::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
    logpredictions = Array{Float64}(undef, length(traces))
    for i in eachindex(traces)
        a, p0 = make_ap(r[i], couplingStrength[i], interval, components, method)
        a_grid = make_a_grid(pgrid, Ngrid)
        d = set_d(noiseparams[i], reporter, Ngrid)
        b = set_b(traces[i], d)
        _, C = forward_grid(a, a_grid, b, p0)
        @inbounds logpredictions[i] = -sum(log.(C))
    end
    sum(logpredictions), logpredictions
end

"""
    ll_hierarchy(pindividual, rhyper)

Loglikelihood for coupled hierarchical model individual parameters.
    lognormal distribution constructed from hyper untransformed noise parameters
"""
function ll_hierarchy(pindividual, rhyper)
    # d = distribution_array(mulognormal(rhyper[1], rhyper[2]), sigmalognormal(rhyper[2]))
    d = distribution_array(rhyper[1], sigmanormal.(rhyper[1], rhyper[2]))
    lhp = Float64[]
    for pc in pindividual
        lhpc = 0
        for i in eachindex(pc)
            lhpc += logpdf(d[i], pc[i])  # Convention: accumulate positive log-likelihoods
        end
        push!(lhp, lhpc)
    end
    lhp
end


###  likelihoods
"""
    ll_hmm(r, components, reporter, interval, trace, method)

"""
# no traits
function ll_hmm(r::Tuple{T1,T2}, components::TComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2}
    rates, noiseparams = r
    a, p0 = make_ap(rates, interval, components, method)
    ll, logpredictions = _ll_hmm(a, p0, set_d(noiseparams, reporter), trace[1])
    lb = ll_off(trace, noiseparams, reporter, components, a, p0)
    ll + lb, logpredictions
end

# coupled
function ll_hmm(r::Tuple{T1,T2,T3}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3}
    rates, noiseparams, couplingStrength = r
    a, p0 = make_ap(rates, couplingStrength, interval, components, method)
    d = set_d(noiseparams, reporter)
    ll, logpredictions = _ll_hmm(a, p0, d, trace[1])
    lb = ll_off(trace, rates, noiseparams, reporter, interval, components, method)
    ll + lb, logpredictions
end

# forced
function ll_hmm(r::Tuple{T1,T2,T3}, components::TForcedComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2,T3}
    rates, noiseparams, couplingStrength = r
    a, p0 = make_ap(rates, couplingStrength, interval, components, method)
    d = (1, set_d(noiseparams, reporter))
    ll, logpredictions = _ll_hmm(a, p0, d, trace[1])
    lb = ll_off(trace, noiseparams, reporter, components, a[1], p0[1])
    ll + lb, logpredictions
end

# hierarchical
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6}, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    if method[2]
        ll, logpredictions = _ll_hmm(noiseindividual, a, p0, reporter, trace[1])
    else
        ll, logpredictions = _ll_hmm(rindividual, noiseindividual, interval, components, reporter, trace[1], method[1])
    end
    lb = ll_off(trace, noiseshared[1], reporter, components, a, p0)
    lhp = ll_hierarchy(pindividual, rhyper)
    ll + lb + sum(lhp), vcat(logpredictions, lhp)
end

# coupled, hierarchical
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, couplingshared, couplingindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    if method[2]
        ll, logpredictions = _ll_hmm(noiseindividual, a, p0, reporter, trace[1])
    else
        ll, logpredictions = _ll_hmm(rindividual, couplingindividual, noiseindividual, interval, components, reporter, trace[1])
    end
    lb = ll_off(trace, rshared[1], noiseshared[1], reporter, interval, components, method[1])
    lhp = ll_hierarchy(pindividual, rhyper)
    ll + lb + sum(lhp), vcat(logpredictions, lhp)
end

# forced, hierarchical
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, components::TForcedComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, couplingshared, couplingindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    if method[2]
        ll, logpredictions = _ll_hmm_forced(noiseindividual, a, p0, reporter, trace[1])
    else
        ll, logpredictions = _ll_hmm(rindividual, couplingindividual, noiseindividual, interval, components, reporter, trace[1])
    end
    lb = ll_off(trace, noiseshared[1], reporter, components, a[1], p0[1])
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
    _ll_hmm_grid(a, a_grid, p0, d, trace[1])
end

# coupled, grid
function ll_hmm(r::Tuple{T1,T2,T3,T4}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3,T4}
    r, noiseparams, couplingStrength, pgrid = r
    a, p0 = make_ap(r, couplingStrength, interval, components, method)
    a_grid = make_a_grid(pgrid[1][1], Ngrid)
    d = set_d(noiseparams, reporter, Ngrid)
    _ll_hmm_grid(a, a_grid, p0, d, trace[1])
end

# hierarchical, grid
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, Ngrid::Int, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, pgridshared, pgridindividual = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    a_grid = make_a_grid(pgridshared[1][1], Ngrid)
    d = set_d(noiseshared[1], reporter)
    if method[2]
        ll, logpredictions = _ll_hmm_grid(noiseindividual, Ngrid, a, a_grid, p0, d, trace[1])
    else
        ll, logpredictions = _ll_hmm_grid(rindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1])
    end
    ll, logpredictions
end

# coupled, hierarchical, grid
function ll_hmm(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    rshared, rindividual, _, noiseindividual, _, _, couplingshared, couplingindividual, pgridshared, pgridindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    a_grid = make_a_grid(pgridshared[1][1], Ngrid)
    d = set_d(noiseshared[1], reporter)
    if method[2]
        ll, logpredictions = _ll_hmm_grid(noiseindividual, Ngrid, a, a_grid, p0, d, trace[1])
    else
        ll, logpredictions = _ll_hmm_grid(rindividual, couplingindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1])
    end
    ll, logpredictions
end


########################
# Predict trace functions
########################
#### Return states and d (not observation_dist)

"""
    _predict_trace(a, p0, d, traces)

Predict most likely state sequences for multiple traces using Viterbi algorithm.

# Arguments
- `a`: Transition probability matrix
- `p0`: Initial state distribution
- `d`: Emission distributions
- `traces`: Vector of observation traces

# Returns
- `Tuple{Vector{Vector{Int}}, Any}`: (states, d) where states contains the most likely state sequence for each trace
"""
function _predict_trace(a, p0, d, traces)
    states = Vector{Int}[]
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        push!(states, viterbi(a, b, p0))
    end
    if d isa Tuple
        return states, d[2]
    else
        return states, d
    end
end

# function _predict_trace_forced(noiseparams::Vector, a, p0::Vector, reporter, traces)
#     states = Vector{Int}[]
#     for i in eachindex(traces)
#         d = (1, set_d(noiseparams[i], reporter))
#         b = set_b(traces[i], d)
#         push!(states, viterbi(a, b, p0))
#     end
#     states, d
# end
"""
    _predict_trace_forced(noiseparams, a, p0, reporter, traces)

Predict most likely state sequences for forced model with multiple traces.

# Arguments
- `noiseparams::Vector`: Vector of noise parameters for each trace
- `a`: Transition probability matrix
- `p0::Vector`: Initial state distribution
- `reporter`: Reporter structure
- `traces`: Vector of observation traces

# Returns
- `Tuple{Vector{Vector{Int}}, Vector}`: (states, d) where states contains the most likely state sequence for each trace and d contains emission distributions
"""
function _predict_trace_forced(noiseparams::Vector, a, p0::Vector, reporter, traces)
    states = Vector{Int}[]
    d = Vector[]
    for i in eachindex(traces)
        di = (1, set_d(noiseparams[i], reporter))
        b = set_b(traces[i], di)
        push!(states, viterbi(a, b, p0))
        push!(d, di[2])
    end
    states, d
end


"""
    _predict_trace(noiseparams, a, p0, reporter, traces)

Predict most likely state sequences using noise parameters and reporter.

# Arguments
- `noiseparams::Vector`: Vector of noise parameters for each trace
- `a::Matrix`: Transition probability matrix
- `p0::Vector`: Initial state distribution
- `reporter`: Reporter structure
- `traces`: Vector of observation traces

# Returns
- `Tuple{Vector{Vector{Int}}, Vector}`: (states, d) where states contains the most likely state sequence for each trace and d contains emission distributions
"""
function _predict_trace(noiseparams::Vector, a::Matrix, p0::Vector, reporter, traces)
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
    _predict_trace(r::Vector, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _predict_trace(r::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
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
    _predict_trace(r, couplingStrength, noiseparams, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _predict_trace(r::Vector, couplingStrength::Vector, noiseparams::Vector, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
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
    _predict_trace(a, a_grid, p0, d, traces)

Predict most likely state sequences for grid-based HMM.

# Arguments
- `a::Matrix`: State transition probability matrix
- `a_grid::Matrix`: Grid transition probability matrix
- `p0::Vector`: Initial state distribution
- `d`: Emission distributions
- `traces`: Vector of observation traces

# Returns
- `Tuple{Vector{Vector{Int}}, Any}`: (states, d) where states contains the most likely state sequence for each trace
"""
function _predict_trace(a::Matrix, a_grid::Matrix, p0::Vector, d, traces)
    states = Vector{Int}[]
    for i in eachindex(traces)
        b = set_b(traces[i], d)
        push!(states, viterbi_grid(a, a_grid, b, p0))
    end
    states, d
end

"""
    _predict_trace(noiseparams, a, a_grid, p0::Vector, reporter, traces)

"""
function _predict_trace(noiseparams::Vector, Ngrid::Int, a::Matrix, a_grid::Matrix, p0::Vector, reporter, traces)
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
    _predict_trace(r::Vector, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _predict_trace(r::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
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
    _predict_trace(r, couplingStrength, noiseparams, pgrid, Ngrid, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())

"""
function _predict_trace(r::Vector, couplingStrength::Vector, noiseparams::Vector, pgrid::Vector, Ngrid::Int, interval::Float64, components::AbstractComponents, reporter, traces, method=Tsit5())
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
    _predict_trace(a, p0, d, trace[1])
end

# coupled
function predict_trace(r::Tuple{T1,T2,T3}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3}
    rates, noiseparams, couplingStrength = r
    a, p0 = make_ap(rates, couplingStrength, interval, components, method)
    d = set_d(noiseparams, reporter)
    _predict_trace(a, p0, d, trace[1])
end

# hierarchical
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6}, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6}
    rshared, rindividual, noiseshared, noiseindividual, _, _ = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    d = set_d(noiseshared[1], reporter)
    if method[2]
        states, observation_dist = _predict_trace(noiseindividual, a, p0, reporter, trace[1])
    else
        states, observation_dist = _predict_trace(rindividual, noiseindividual, interval, components, reporter, trace[1], method[1])
    end
    states, observation_dist
end

# coupled, hierarchical
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, _, _, couplingshared, couplingindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    d = set_d(noiseshared[1], reporter)
    if method[2]
        states, observation_dist = _predict_trace(noiseindividual, a, p0, reporter, trace[1])
    else
        states, observation_dist = _predict_trace(rindividual, couplingindividual, noiseindividual, interval, components, reporter, trace[1])
    end
    states, observation_dist
end

# forced
function predict_trace(r::Tuple{T1,T2,T3}, components::TForcedComponents, reporter::HMMReporter, interval, trace, method=Tsit5()) where {T1,T2,T3}
    rates, noiseparams, couplingStrength = r
    a, p0 = make_ap(rates, couplingStrength, interval, components, method)
    d = (1, set_d(noiseparams, reporter))
    _predict_trace(a, p0, d, trace[1])
end

# forced, hierarchical
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, components::TForcedComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method::Tuple=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, _, _, couplingshared, couplingindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    d = (1, set_d(noiseshared[1], reporter))
    if method[2]
        states, observation_dist = _predict_trace_forced(noiseindividual, a, p0, reporter, trace[1])
    else
        states, observation_dist = _predict_trace(rindividual, couplingindividual, noiseindividual, interval, components, reporter, trace[1])
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
    _predict_trace(a, a_grid, p0, d, trace[1])
end

# coupled, grid
function predict_trace(r::Tuple{T1,T2,T3,T4}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval, trace, method=Tsit5()) where {T1,T2,T3,T4}
    r, noiseparams, couplingStrength, pgrid = r
    a, p0 = make_ap(r, couplingStrength, interval, components, method)
    a_grid = make_a_grid(pgrid[1][1], Ngrid)
    d = set_d(noiseparams, reporter, Ngrid)
    _predict_trace(a, a_grid, p0, d, trace[1])
end

# hierarchical, grid
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8}, Ngrid::Int, components::TComponents, reporter::HMMReporter, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8}
    rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, pgridshared, pgridindividual = r
    a, p0 = make_ap(rshared[1], interval, components, method[1])
    a_grid = make_a_grid(pgridshared[1][1], Ngrid)
    d = set_d(noiseshared[1], reporter)
    if method[2]
        _predict_trace(noiseindividual, a, a_grid, p0, d, trace[1])
    else
        _predict_trace(rindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1])
    end
end

# coupled, hierarchical, grid
function predict_trace(r::Tuple{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}, Ngrid::Int, components::TCoupledComponents, reporter::Vector{HMMReporter}, interval::Float64, trace::Tuple, method=(Tsit5(), true)) where {T1,T2,T3,T4,T5,T6,T7,T8,T9,T10}
    rshared, rindividual, _, noiseindividual, _, _, couplingshared, couplingindividual, pgridshared, pgridindividual = r
    a, p0 = make_ap(rshared[1], couplingshared[1], interval, components, method[1])
    a_grid = make_a_grid(pgridshared[1][1], Ngrid)
    d = set_d(noiseshared[1], reporter)
    if method[2]
        _predict_trace(noiseindividual, a, a_grid, p0, d, trace[1])
    else
        _predict_trace(rindividual, couplingindividual, noiseindividual, pgridindividual, Ngrid, interval, components, reporter, trace[1])
    end
end


function predict_trace(param, data, model::AbstractGRSMmodel)
    r = prepare_rates(param, model)
    if hastrait(model, :grid)
        predict_trace(r, model.trait.grid.ngrid, model.components, model.reporter, data.interval, data.trace, model.method)
    else
        predict_trace(r, model.components, model.reporter, data.interval, data.trace, model.method)
    end
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

Compute expected transition counts using log-space forward and backward variables.

# Arguments
- `logα`: Log forward variable matrix
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `logβ`: Log backward variable matrix
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Tuple{Array{Float64,3}, Array{Float64,2}}`: (ξ, γ) where ξ[i,j,t] = P(q_t = S_i, q_{t+1} = S_j | O, λ) and γ[i,t] = ∑_j ξ[i,j,t]
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

Compute expected transition probability matrix from observations.

# Arguments
- `a`: Current transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Matrix{Float64}`: Expected transition probability matrix computed from forward-backward algorithm
"""
function expected_a(a, b, p0, N, T)
    α, C = forward(a, b, p0, N, T)
    β = backward(a, b, C, N, T)
    ξ, γ = expected_transitions(α, a, b, β, N, T)
    expected_a(ξ, γ, N)
end
"""
    expected_a(ξ, γ, N)

Compute expected transition probability matrix from expected transition counts.

# Arguments
- `ξ`: Expected transition counts array of size (N, N, T-1)
- `γ`: Expected state counts array of size (N, T-1)
- `N::Int`: Number of states

# Returns
- `Matrix{Float64}`: Expected transition probability matrix where a[i,j] = ∑_t ξ[i,j,t] / ∑_t γ[i,t]
"""
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

function viterbi(a::Vector{T1}, b::Vector{T2}, p0, N, T) where {T1<:AbstractArray,T2<:AbstractArray}
    logb = log.(max.(b[2], 0.0))
    ϕ = similar(logb)
    ψ = similar(ϕ)
    q = Vector{Int}(undef, T)
    if ~b[1][1]
        ϕ[:, 1] = log.(max.(p0[1], 0.0)) .+ logb[:, 1]
    else
        ϕ[:, 1] = log.(max.(p0[2], 0.0)) .+ logb[:, 1]
    end
    ψ[:, 1] .= 0
    for t in 2:T
        for j in 1:N
            if !b[1][1]
                m, ψ[j, t] = findmax(ϕ[:, t-1] + log.(max.(a[1][:, j], 0.0)))
            else
                m, ψ[j, t] = findmax(ϕ[:, t-1] + log.(max.(a[2][:, j], 0.0)))
            end
            ϕ[j, t] = m + logb[j, t]
        end
    end
    q[T] = argmax(ϕ[:, T])
    for t in T-1:-1:1
        q[t] = ψ[q[t+1], t+1]
    end
    return q



    # α = zeros(N, T)
    # C = Vector{Float64}(undef, T)
    # if ~b[1][1]
    #     α[:, 1] = p0[1] .* b[2][:, 1]
    # else
    #     α[:, 1] = p0[2] .* b[2][:, 1]
    # end
    # C[1] = 1 / max(sum(α[:, 1]), eps(Float64))
    # α[:, 1] *= C[1]
    # for t in 2:T
    #     for j in 1:N
    #         for i in 1:N
    #             forward_inner_operation!(α, a, b, i, j, t)
    #         end
    #     end
    #     C[t] = 1 / max(sum(α[:, t]), eps(Float64))
    #     α[:, t] *= C[t]
    # end
    # return α, C
    # end
end
"""
    viterbi(a, b, p0, N, T)

Find maximum likelihood state sequence using Viterbi algorithm.

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution
- `N`: Number of states
- `T`: Number of time steps

# Returns
- `Vector{Int}`: Most likely state sequence q* = argmax P(q | O, λ)
"""
function viterbi(a, b, p0, N, T)
    loga = log.(max.(a, 0.0))
    logb = log.(max.(b, 0.0))
    logp0 = log.(max.(p0, 0.0))
    viterbi_log(loga, logb, logp0, N, T)
end

"""
    viterbi(a, b, p0)

Find maximum likelihood state sequence using Viterbi algorithm (auto-determine dimensions).

# Arguments
- `a`: Transition probability matrix
- `b`: Emission probability matrix
- `p0`: Initial state distribution

# Returns
- `Vector{Int}`: Most likely state sequence q* = argmax P(q | O, λ)
"""
function viterbi(a, b, p0)
    N, T = size(b)
    viterbi(a, b, p0, N, T)
end

function viterbi(a::Vector{T1}, b::Vector{T2}, p0) where {T1<:AbstractArray,T2<:AbstractArray}
    N, T = size(b[2])
    viterbi(a, b, p0, N, T)
end

"""
    viterbi_grid_log(loga, loga_grid, logb, logp0, Nstate, Ngrid, T)

Find maximum likelihood state sequence for grid-based HMM using Viterbi algorithm in log-space.

# Arguments
- `loga`: Log state transition probability matrix
- `loga_grid`: Log grid transition probability matrix
- `logb`: Log emission probability matrix
- `logp0`: Log initial state distribution
- `Nstate`: Number of states
- `Ngrid`: Number of grid points
- `T`: Number of time steps

# Returns
- `Vector{Int}`: Most likely state sequence q* = argmax P(q | O, λ)
"""
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
"""
    viterbi_grid(a, a_grid, b, p0, Nstate, Ngrid, T)

Find maximum likelihood state sequence for grid-based HMM using Viterbi algorithm.

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix of size (Nstate, Ngrid, T)
- `p0`: Initial state distribution
- `Nstate`: Number of states
- `Ngrid`: Number of grid points
- `T`: Number of time steps

# Returns
- `Vector{Int}`: Most likely state sequence q* = argmax P(q | O, λ)
"""
function viterbi_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    loga = log.(max.(a, 0.0))
    loga_grid = log.(max.(a_grid, 0.0))
    logb = log.(max.(b, 0.0))
    logp0 = log.(max.(p0, 0.0))
    viterbi_grid_log(loga, loga_grid, logb, logp0, Nstate, Ngrid, T)
end

"""
    viterbi_grid(a, a_grid, b, p0)

Find maximum likelihood state sequence for grid-based HMM using Viterbi algorithm (auto-determine dimensions).

# Arguments
- `a`: State transition probability matrix
- `a_grid`: Grid transition probability matrix
- `b`: Emission probability matrix of size (Nstate, Ngrid, T)
- `p0`: Initial state distribution

# Returns
- `Vector{Int}`: Most likely state sequence q* = argmax P(q | O, λ)
"""
function viterbi_grid(a, a_grid, b, p0)
    Nstate, Ngrid, T = size(b)
    viterbi_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
end

"""
    correlation_functions(rin, transitions, G::Tuple, R, S, insertstep, probfn, coupling, lags::Vector; offset=0.0)

Compute theoretical correlation functions for coupled HMM model.

Computes uncentered cross-correlations and autocorrelations for intensity (with noise), ON states (binary),
and reporter counts (deterministic per state). All correlations are uncentered (E[xy]).

The lag interval is automatically inferred from the `lags` vector as the minimum spacing between consecutive lags.
The transition matrix is constructed to represent transitions over this lag interval, allowing non-integer lags
to be handled via eigendecomposition.

# Arguments
- `rin`: Input rate parameters
- `transitions::Tuple`: Transition definitions for each unit
- `G::Tuple`: Number of gene states for each unit
- `R`: Number of RNA states for each unit (Tuple for coupled models)
- `S`: Number of splice sites for each unit (Tuple for coupled models)
- `insertstep`: Reporter insertion step definitions (Tuple for coupled models)
- `probfn`: Probability function for observation model (e.g., `prob_Gaussian`)
- `coupling::Tuple`: Coupling structure between units
- `lags::Vector{<:Real}`: Time lags for covariance calculation (positive lags only, monotonic increasing from 0)
  - **Lags must be uniform or multiples of a minimal value** (e.g., [0, 5/3, 10/3, 15/3, ...])
  - The lag interval is inferred as the first spacing: `lags[2] - lags[1]`
  - Lags are converted to step indices (number of intervals): `step_indices = [0, 1, 2, ...]`
  - The transition matrix `a` is constructed for one lag interval, then `a^step_indices[i]` is used
- `offset::Float64=0.0`: Offset added to ON states and reporter counts (default: 0.0)
  - Applied as: ON = float(num_per_state .> 0.0) .+ offset
  - Added for compatibility with experimental data processing conventions

# Returns
- `Tuple` containing (in order, total of 22 values):
  1. `tau`: Time lags (SYMMETRIZED: includes negative lags, format: [-max_lag, ..., -1, 0, 1, ..., max_lag])
  2. `cc`: Cross-correlation function for intensity (SYMMETRIZED: includes negative lags, uncentered E[xy])
  3. `ac1`: Autocorrelation function for intensity unit 1 (positive lags only, NOT symmetrized, uncentered E[xx])
  4. `ac2`: Autocorrelation function for intensity unit 2 (positive lags only, NOT symmetrized, uncentered E[xx])
  5. `m1`: Mean intensity for unit 1
  6. `m2`: Mean intensity for unit 2
  7. `v1`: Variance of intensity for unit 1 (E[O²] - E[O]², includes noise variance)
  8. `v2`: Variance of intensity for unit 2 (E[O²] - E[O]², includes noise variance)
  9. `ccON`: Cross-correlation function for ON states (SYMMETRIZED: includes negative lags, uncentered E[xy])
  10. `ac1ON`: Autocorrelation function for ON states unit 1 (positive lags only, NOT symmetrized, uncentered E[xx])
  11. `ac2ON`: Autocorrelation function for ON states unit 2 (positive lags only, NOT symmetrized, uncentered E[xx])
  12. `m1ON`: Mean ON state probability for unit 1 (with offset applied)
  13. `m2ON`: Mean ON state probability for unit 2 (with offset applied)
  14. `v1`: Variance of intensity for unit 1 (duplicate of #7)
  15. `v2`: Variance of intensity for unit 2 (duplicate of #8)
  16. `ccReporters`: Cross-correlation function for reporter counts (SYMMETRIZED: includes negative lags, uncentered E[xy])
  17. `ac1Reporters`: Autocorrelation function for reporter counts unit 1 (positive lags only, NOT symmetrized, uncentered E[xx])
  18. `ac2Reporters`: Autocorrelation function for reporter counts unit 2 (positive lags only, NOT symmetrized, uncentered E[xx])
  19. `m1Reporters`: Mean reporter count for unit 1 (with offset applied)
  20. `m2Reporters`: Mean reporter count for unit 2 (with offset applied)
  21. `v1Reporters`: Variance of reporter counts for unit 1
  22. `v2Reporters`: Variance of reporter counts for unit 2

# Notes
- **Input lags**: Must be positive lags only (monotonic increasing from 0, e.g., `collect(0:1:60)`)
- **Output lags (tau)**: Symmetric lags including negative values: `[-max_lag, ..., -1, 0, 1, ..., max_lag]`
- **Symmetrized vs. Non-symmetrized**:
  - **Symmetrized** (include negative lags): `cc`, `ccON`, `ccReporters`, `tau`
  - **NOT symmetrized** (positive lags only): `ac1`, `ac2`, `ac1ON`, `ac2ON`, `ac1Reporters`, `ac2Reporters`
  - Autocovariances are symmetric functions (ac(τ) = ac(-τ)), so they can be symmetrized by the caller if needed
- **Intensity covariances**: Use `autocorfn_hmm` to account for noise variance at lag 0 (uses second moment)
- **ON state covariances**: Use `crosscorfn_hmm` (binary: ON² = ON, so no separate second moment needed)
- **Reporter covariances**: Use `crosscorfn_hmm` (deterministic per state, variance only from state transitions)
- **Offset**: Applied to ON states and reporter counts to match experimental data processing
  - ON states: `float(num_per_state .> 0.0) .+ offset` (binary with offset)
  - Reporters: `Float64.(num_per_state) .+ offset` (counts with offset)
- All correlations are **uncentered** (E[xy], NOT E[xy] - E[x]E[y])
- To get covariances, subtract the product of means: C_XY(τ) = R_XY(τ) - E[X]E[Y]
- **Lag convention**: Positive τ in input `lags` means first unit leads second unit
- **Cross-covariance symmetry**: `cc` is computed as `cc12` (unit 1 leads) for positive lags and `cc21` (unit 2 leads) for negative lags, then symmetrized
"""
function correlation_functions(rin, transitions, G::Tuple, R, S, insertstep, probfn, coupling, lags::Vector; offset::Float64=0.0)
    components = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    # components = TRGCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    sourceStates = [c.sourceState for c in components.modelcomponents]
    r, couplingStrength, noiseparams = prepare_rates_coupled(rin, sourceStates, transitions, R, S, insertstep, [4, 4])
    num_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
    mean_intensity = Vector[]
    second_moment_intensity = Vector[]  # E[O^2|state] for each state
    ON = Vector[]
    reporters = Vector[]
    for i in eachindex(noiseparams)
        dists = probfn(noiseparams[i], num_per_state[i], components.N)
        mi = mean.(dists)  # E[O|state]
        vi = var.(dists)   # Var(O|state)
        second_moment = vi .+ mi.^2
        push!(mean_intensity, mi)
        push!(second_moment_intensity, second_moment)

        push!(ON, float(num_per_state[i] .> 0.0) .+ offset)
        push!(reporters, Float64.(num_per_state[i]) .+ offset)
    end

    # Infer lag interval from lags vector
    # Lags must be uniform or multiples of a minimal value
    if length(lags) < 2
        lag_interval = lags[1] > 0 ? lags[1] : 1.0
    else
        # Use first spacing as the lag interval (assuming uniform or multiples)
        lag_interval = lags[2] - lags[1]
        if lag_interval <= 0
            error("Lags must be strictly increasing and positive")
        end
        # Verify that all lags are multiples of the interval (within tolerance)
        step_indices = lags ./ lag_interval
        if !all(step -> isapprox(step, round(step), rtol=1e-10), step_indices)
            error("Lags must be uniform or multiples of the lag interval. First spacing: $lag_interval, but found non-multiple lag")
        end
    end

    # Convert lags to step indices (number of intervals)
    # e.g., if lags = [0, 5/3, 10/3, ...] and lag_interval = 5/3, then step_indices = [0, 1, 2, ...]
    step_indices = round.(Int, lags ./ lag_interval)

    # transition matrix a and steady state probabilities p0
    # Construct a to represent transitions over one lag_interval
    a, p0 = make_ap(r, couplingStrength, lag_interval, components)

    # Cross-correlations
    cc12 = crosscorfn_hmm(a, p0, mean_intensity[1], mean_intensity[2], step_indices) 
    cc21 = crosscorfn_hmm(a, p0, mean_intensity[2], mean_intensity[1], step_indices)
    cc = vcat(reverse(cc21), cc12[2:end])

    # Need to account for intensity variance at lag 0 for autocovariance functions
    ac1 = autocorfn_hmm(a, p0, mean_intensity[1], second_moment_intensity[1], step_indices) 
    ac2 = autocorfn_hmm(a, p0, mean_intensity[2], second_moment_intensity[2], step_indices)
    
    # ON state 
    ccON = crosscorfn_hmm(a, p0, ON[2], ON[1], step_indices)
    ccON = vcat(reverse(ccON), ccON[2:end])
    ac1ON = crosscorfn_hmm(a, p0, ON[1], ON[1], step_indices)
    ac2ON = crosscorfn_hmm(a, p0, ON[2], ON[2], step_indices)

    # Reporter
    ccReporters = crosscorfn_hmm(a, p0, reporters[2], reporters[1], step_indices)
    ccReporters = vcat(reverse(ccReporters), ccReporters[2:end])
    ac1Reporters = crosscorfn_hmm(a, p0, reporters[1], reporters[1], step_indices)
    ac2Reporters = crosscorfn_hmm(a, p0, reporters[2], reporters[2], step_indices)

    # Means
    m1 = mean_hmm(p0, mean_intensity[1])
    m2 = mean_hmm(p0, mean_intensity[2])
    m1ON = mean_hmm(p0, ON[1])
    m2ON = mean_hmm(p0, ON[2])
    m1Reporters = mean_hmm(p0, reporters[1])
    m2Reporters = mean_hmm(p0, reporters[2])

    # Variances (at step index 0, which corresponds to lag 0)
    zero_step_idx = findfirst(==(0), step_indices)
    v1 = (ac1[zero_step_idx] - m1^2)
    v2 = (ac2[zero_step_idx] - m2^2)
    v1ON = (ac1ON[zero_step_idx] - m1ON^2)
    v2ON = (ac2ON[zero_step_idx] - m2ON^2)
    v1Reporters = (ac1Reporters[zero_step_idx] - m1Reporters^2)
    v2Reporters = (ac2Reporters[zero_step_idx] - m2Reporters^2)
    
    return vcat(-reverse(lags), lags[2:end]), cc, ac1, ac2, m1, m2, v1, v2, ccON, ac1ON, ac2ON, m1ON, m2ON, v1ON, v2ON, ccReporters, ac1Reporters, ac2Reporters, m1Reporters, m2Reporters, v1Reporters, v2Reporters
end

"""
    crosscorfn_hmm(a, p0, meanintensity1, meanintensity2, step_indices)

Compute uncentered cross-correlation function E[X(t)Y(t+τ)] between two signals.

The transition matrix `a` represents transitions over one lag interval. The `step_indices` vector
contains the number of intervals for each lag (e.g., [0, 1, 2, ...] for uniform spacing).
We compute `a^step_indices[i]` for each step index using integer matrix powers.
"""
function crosscorfn_hmm(a, p0, meanintensity1, meanintensity2, step_indices::Vector{Int})
    cc = zeros(length(step_indices))
    m1 = meanintensity1 .* p0
    al = a^step_indices[1]
    if length(step_indices) > 1
        as = a^(step_indices[2] - step_indices[1])
    else
        as = I  # Identity if only one step
    end
    for l in eachindex(step_indices)
        for i in eachindex(meanintensity1)
            for j in eachindex(meanintensity2)
                cc[l] += m1[i] * al[i, j] * meanintensity2[j]
            end
        end
        if l < length(step_indices)
            # Compute step difference and multiply al by as raised to that power
            delta = step_indices[l+1] - step_indices[l]
            al *= as^delta
        end
    end
    cc
end
"""
    autocorfn_hmm(a, p0, meanintensity, second_moment_intensity, step_indices)

Compute uncentered autocorrelation function E[X(t)X(t+τ)] for HMM model.
Uses second moment at lag 0 to account for noise variance.
"""
function autocorfn_hmm(a, p0, meanintensity, second_moment_intensity, step_indices::Vector{Int})
    ac = crosscorfn_hmm(a, p0, meanintensity, meanintensity, step_indices)
    # Find step index 0 (should be first)
    zero_idx = findfirst(==(0), step_indices)
    if zero_idx !== nothing
        ac[zero_idx] = sum(p0 .* second_moment_intensity)
    else
        ac[1] = sum(p0 .* second_moment_intensity)
    end
    ac
end


"""
    correlation_functions_centered(rin, transitions, G::Tuple, R, S, insertstep, probfn, coupling, lags::Vector; offset=0.0)

Compute theoretical centered correlation functions for coupled HMM model.

Returns centered correlation functions (E[xy] - E[x]E[y]) by calling `correlation_functions` and subtracting means.
See `correlation_functions` for detailed argument descriptions and return value structure.
"""
function correlation_functions_centered(rin, transitions, G::Tuple, R, S, insertstep, probfn, coupling, lags::Vector; offset::Float64=0.0)
    ac1, ac2, cc, ccON, lags, m1, m2, v1, v2, m1ON, m2ON, ac1ON, ac2ON, ccReporters, m1Reporters, m2Reporters, ac1Reporters, ac2Reporters = correlation_functions(rin, transitions, G, R, S, insertstep, probfn, coupling, lags; offset=offset)
    return lags, cc-m1*m2, ac1-m1^2, ac2-m1^2, m1, m2, v1, v2, ccON-m1ON*m2ON, ac1ON-m1ON^2, ac2ON-m2ON^2, m1ON, m2ON, v1ON, v2ON, ccReporters-m1Reporters*m2Reporters,  ac1Reporters-m1Reporters^2, ac2Reporters-m2Reporters^2, m1Reporters, m2Reporters, v1Reporters, v2Reporters
end

"""
    crosscov_hmm(a, p0, meanintensity1, meanintensity2, step_indices, m1, m2)

Compute cross-covariance function E[X(t)Y(t+τ)] - E[X]E[Y] between two signals.
"""
function crosscov_hmm(a, p0, meanintensity1, meanintensity2, step_indices::Vector{Int}, m1, m2)
    crosscorfn_hmm(a, p0, meanintensity1, meanintensity2, step_indices) .- m1 .* m2
end

"""
    crosscov_hmm(a, p0, meanintensity1, meanintensity2, step_indices)

Compute cross-covariance function (auto-computes means from p0).
"""
function crosscov_hmm(a, p0, meanintensity1, meanintensity2, step_indices::Vector{Int})
    crosscov_hmm(a, p0, meanintensity1, meanintensity2, step_indices, mean_hmm(p0, meanintensity1), mean_hmm(p0, meanintensity2))
end

"""
    autocov_hmm(a, p0, meanintensity, second_moment_intensity, lags)

Compute autocovariance function E[X(t)X(t+τ)] - E[X]² for HMM model.
Uses second moment at lag 0 to account for noise variance.
"""
function autocov_hmm(a, p0, meanintensity, second_moment_intensity, step_indices::Vector{Int})
    ac = autocorfn_hmm(a, p0, meanintensity, second_moment_intensity, step_indices)
    ac .- mean_hmm(p0, meanintensity) .^ 2
end

"""
    autocov_hmm(r, transitions, G, R, S, insertstep, probfn, lags::Vector)

Compute autocovariance function from rate parameters and model structure.

The lag interval is automatically inferred from the `lags` vector as the minimum spacing between consecutive lags.
"""
function autocov_hmm(r, transitions, G, R, S, insertstep, probfn, lags::Vector)
    components = TComponents(transitions, G, R, S, insertstep, "")
    dists = probfn(r[end-3:end], num_reporters_per_state(G, R, S, insertstep), components.nT)
    mean_intensity = mean.(dists)
    vi = var.(dists)
    # Second moment: E[O^2|state] = Var(O|state) + E[O|state]^2
    second_moment_intensity = vi .+ mean_intensity.^2
    
    # Infer lag interval from lags vector
    # Lags must be uniform or multiples of a minimal value
    if length(lags) < 2
        lag_interval = lags[1] > 0 ? lags[1] : 1.0
    else
        # Use first spacing as the lag interval (assuming uniform or multiples)
        lag_interval = lags[2] - lags[1]
        if lag_interval <= 0
            error("Lags must be strictly increasing and positive")
        end
        # Verify that all lags are multiples of the interval (within tolerance)
        step_indices = lags ./ lag_interval
        if !all(step -> isapprox(step, round(step), rtol=1e-10), step_indices)
            error("Lags must be uniform or multiples of the lag interval. First spacing: $lag_interval, but found non-multiple lag")
        end
    end
    
    # Convert lags to step indices (number of intervals)
    step_indices = round.(Int, lags ./ lag_interval)
    
    a, p0 = make_ap(r, lag_interval, components)
    autocov_hmm(a, p0, mean_intensity, second_moment_intensity, step_indices)
end

"""
    mean_hmm(p0, meanintensity)

Compute mean of signal: E[X] = Σ p0[i] * meanintensity[i]
"""
function mean_hmm(p0, meanintensity)
    sum(p0 .* meanintensity)
end



# function covariance_functions_scaled(rin, transitions, G::Tuple, R, S, insertstep, interval, probfn, coupling, lags::Vector)
#     components = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
#     # components = TRGCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
#     sourceStates = [c.sourceState for c in components.modelcomponents]
#     r, couplingStrength, noiseparams = prepare_rates_coupled(rin, sourceStates, transitions, R, S, insertstep, [4, 4])
#     num_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
#     mean_intensity = Vector[]
#     second_moment_intensity = Vector[]  # E[O^2|state] for each state
#     max_intensity = Float64[]
#     ON = Vector[]
#     for i in eachindex(noiseparams)
#         dists = probfn(noiseparams[i], num_per_state[i], components.N)
#         mi = mean.(dists)  # E[O|state]
#         vi = var.(dists)   # Var(O|state)
#         # Second moment: E[O^2|state] = Var(O|state) + E[O|state]^2
#         second_moment = vi .+ mi.^2
#         mmax = max(maximum(mi), 1.0)
#         push!(max_intensity, mmax)
#         push!(mean_intensity, mi / mmax)
#         # Normalize second moment by mmax^2 to match normalized mean_intensity
#         push!(second_moment_intensity, second_moment / (mmax^2))
#         push!(ON, float(num_per_state[i] .> 0.0))
#     end
#     a, p0 = make_ap(r, couplingStrength, interval, components)
#     m1 = mean_hmm(p0, mean_intensity[1])
#     m2 = mean_hmm(p0, mean_intensity[2])
#     mON1 = mean_hmm(p0, ON[1])
#     mON2 = mean_hmm(p0, ON[2])

#     cc12 = crosscov_hmm(a, p0, mean_intensity[1], mean_intensity[2], lags, m1, m2) * max_intensity[1] * max_intensity[2]
#     cc21 = crosscov_hmm(a, p0, mean_intensity[2], mean_intensity[1], lags, m1, m2) * max_intensity[1] * max_intensity[2]
#     # For enhancer-leads convention (positive tau means enhancer leads):
#     # Use ON[2], ON[1] to compute <gene(t)enhancer(t+tau)>, which for positive tau means enhancer leads
#     ccON = crosscov_hmm(a, p0, ON[2], ON[1], lags, mON2, mON1)
#     ac1 = crosscov_hmm(a, p0, mean_intensity[1], mean_intensity[1], lags, m1, m1) * max_intensity[1]^2
#     ac2 = crosscov_hmm(a, p0, mean_intensity[2], mean_intensity[2], lags, m2, m2) * max_intensity[2]^2
#     ac1 = autocov_hmm(a, p0, mean_intensity[1], second_moment_intensity[1], lags) * max_intensity[1]^2
#     ac2 = autocov_hmm(a, p0, mean_intensity[2], second_moment_intensity[2], lags) * max_intensity[2]^2
#     # Variance: E[O^2] - E[O]^2, using second moment for E[O^2]
#     v1 = (sum(p0 .* second_moment_intensity[1]) - m1^2) * max_intensity[1]^2
#     v2 = (sum(p0 .* second_moment_intensity[2]) - m2^2) * max_intensity[2]^2
#     # Scale means back to original scale (m1 and m2 are computed on normalized scale)
#     m1_scaled = m1 * max_intensity[1]
#     m2_scaled = m2 * max_intensity[2]
#     # ac1 = vcat(reverse(ac1), ac1[2:end])
#     # ac2 = vcat(reverse(ac2), ac2[2:end])
#     cc = vcat(reverse(cc21), cc12[2:end])
#     ccON = vcat(reverse(ccON), ccON[2:end])
#     ac1, ac2, cc, ccON, vcat(-reverse(lags), lags[2:end]), m1_scaled, m2_scaled, v1, v2
# end

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


