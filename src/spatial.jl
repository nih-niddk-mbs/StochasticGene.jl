# Spatial.jl
#
# Experimental code for including spatial information from images of transcriptional bursting
"""
    kronecker_delta(i, j)

Computes the Kronecker delta of two integers.

# Arguments
- `i`: The first integer.
- `j`: The second integer.

# Returns
- `Int`: Returns 1 if `i` equals `j`, otherwise returns 0.
"""
kronecker_delta(i, j) = i == j ? 1 : 0

"""
    read_tracefiles_spatial(path::String, label::String, start::Int, stop::Int, col=3)

Reads trace files from a specified directory that match a given label and extracts data from them.

# Arguments
- `path::String`: The directory path to search for trace files.
- `label::String`: The label to match in the filenames.
- `start::Int`: The starting index for reading the trace data.
- `stop::Int`: The stopping index for reading the trace data.
- `col::Int`: The column index to read from each trace file (default is 3).

# Returns
- `Vector`: A vector of unique traces read from the files.
"""
function read_tracefiles_spatial(path::String, label::String, start::Int, stop::Int, col=3)
    traces = Vector[]
    if isempty(path)
        return traces
    else
        for (root, dirs, files) in walkdir(path)
            for file in files
                if occursin(label, file) && ~occursin(".DS_Store", file)
                    t = read_tracefile_spatial(joinpath(root, file), start, stop, col)
                    ~isempty(t) && push!(traces, t)
                end
            end
        end
        set = sum.(traces)
        return traces[unique(i -> set[i], eachindex(set))]  # only return unique traces
    end
end

"""
    set_b_spatial(trace, params, reporters_per_state, probfn::Function, Ns, Np)

Calculates the probability matrix `b` for given trace data using a specified probability function.

# Arguments
- `trace`: The trace data.
- `params`: Parameters for the probability function.
- `reporters_per_state`: Number of reporters per state.
- `probfn::Function`: The probability function to use.
- `Ns`: Number of states.
- `Np`: Number of positions.

# Returns
- `Matrix`: The probability matrix `b`.
"""
function set_b_spatial_v(trace, params, reporters_per_state, probfn::Function, Ns, Np)
    d = probfn(params, reporters_per_state, Ns, Np)
    b = Ones(Ns * Np, length(trace))
    t = 1
    for obs in eachcol(trace)
        i = 1
        for j in 1:Ns
            for k in 1:Np
                for l in eachindex(obs)
                    b[i, t] *= pdf(d[j, k, l], obs[l])
                end
            end
            i += 1
        end
        t += 1
    end
    return b
end
function set_b_spatial(trace, params, reporters_per_state, probfn::Function, Ns, Np)
    d = probfn(params, reporters_per_state, Ns, Np)
    b = Ones(Ns, Np, length(trace))
    t = 1
    for obs in eachcol(trace)
        for i in 1:Ns
            for j in 1:Np
                for k in eachindex(obs)
                    b[i, j, t] *= pdf(d[i, j, k], obs[k])
                end
            end
        end
        t += 1
    end
    return b
end
"""
    prob_Gaussian_spatial(par, reporters_per_state, Ns, Np, f::Function=kronecker_delta)

Generates a 3D array of Normal distributions based on the given parameters and reporters per state.

# Arguments
- `par`: Parameters for the Gaussian distribution.
- `reporters_per_state`: Number of reporters per state.
- `Ns`: Number of states.
- `Np`: Number of positions.
- `f::Function`: Function to use for Kronecker delta (default is `kronecker_delta`).

# Returns
- `Array{Distribution{Univariate,Continuous}}`: A 3D array of Normal distributions.
"""
function prob_Gaussian_spatial(par, reporters_per_state, Ns, Np, f::Function=kronecker_delta)
    d = Array{Distribution{Univariate,Continuous}}(undef, Ns, Np, Np)
    for j in 1:Ns
        for k in 1:Np
            for l in 1:Np
                σ = sqrt(par[2]^2 + reporters_per_state[j] * par[4]^2 * f(k, l))
                d[j, k, l] = Normal(par[1] + reporters_per_state[j] * par[3] * f(k, l), σ)
            end
        end
    end
    return d
end

function forward_spatial(as, ap, b, p0, Ns, Np, T)
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    α[:, 1] = p0 .* b[:, 1]
    C[1] = 1 / sum(α[:, 1])
    α[:, 1] *= C[1]
    for t in 2:T
        for i in 1:Ns
            for j in 1:Ns
                for k in 1:Np
                    for l in 1:Np
                        α[j, t] += α[i, t-1] * as[i, j] * ap[k, l] * b[j, t]
                    end
                    α[j, t] += α[i, t-1] * as[i, j] * ap[k, l] * b[j, t]
                end
            end
        end
        C[t] = 1 / sum(α[:, t])
        α[:, t] *= C[t]
    end
    return α, C
end

function make_as(param, Ns)
    as = zeros(Ns, Ns)
    d = zeros(Ns, Ns)
    for i in 1:Ns
        for j in 1:Ns
            as[i, j] = exp(-distance(i, j, div(Ns,2))^2 / (2 * param^2))
            d[i,j] = distance(i, j, div(Ns,2))
        end
    end
    return as ./ sum(as, dims = 2)
end

function distance(i, j, Ns)
    xi = rem(i - 1, Ns)
    yi = div(i - 1, Ns)
    xj = rem(j - 1, Ns)
    yj = div(j - 1, Ns)
    return sqrt((xi-xj)^2 + (yi-yj)^2)
end

function speedtest1(M, n)
    M^n
end

function speedtest2(M, n)
    m = size(M, 1)
    Mout = Diagonal(ones(m))  # Initialize Mout as an identity matrix
    for t in 1:n
        Mtemp = zeros(m, m)  # Temporary matrix to store intermediate results
        for i in 1:m
            for j in 1:m
                for k in 1:m
                    Mtemp[i, j] += Mout[i, k] * M[k, j]
                end
            end
        end
        Mout = Mtemp  # Update Mout with the result of the multiplication
    end
    return Mout
end

function test_fit_trace_spatial(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70])
    trace = StochasticGene.simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("tracespatial", "test", interval, (trace, [], weight, nframes))
    # model = StochasticGene.load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple())
    elongationtime = StochasticGene.mean_elongationtime(rtarget, transitions, R)
    priormean = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, elongationtime)
    model = StochasticGene.load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], 1.0, propcv, "", prob_Gaussian, noisepriors, tuple(), tuple(), 1)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget, fits, stats, measures, model, data
end


# components = StochasticGene.make_components_TRG(transitions, G, R, S, insertstep, splicetype)
# a, p0 = make_ap(r, interval, components)
# reporters_per_state = StochasticGene.num_reporters_per_state(G, R, S, insertstep)