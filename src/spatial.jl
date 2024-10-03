# Spatial.jl
#
# Experimental code for including spatial information from images of transcriptional bursting

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
function set_b_spatial(trace, params, reporters_per_state, probfn::Function, Ns, Np)
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
