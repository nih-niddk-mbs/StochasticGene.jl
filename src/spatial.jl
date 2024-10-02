# Spatial.jl
#
# Experimental code for including spatial information from images of transcriptional bursting


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

function set_b_spatial(trace, params, reporters_per_state, probfn::Function, N)
    # N = length(reporters_per_state)
    d = probfn(params, reporters_per_state, N)
    b = Matrix{Float64}(undef, N, length(trace))
    t = 1
    for obs in trace
        for i in location
            for j in 1:N
                b[i, j, t] = pdf(d[j], obs[i])
            end
        end
        t += 1
    end
    return b
end

function prob_Gaussian_spatial(par, reporters_per_state, position)
    N = length(reporters_per_state)*length(position)
    d = Array{Distribution{Univariate,Continuous}}(undef, N)
    for i in 1:N
        d[i] = Normal(par[1] + reporters_per_state[i] * par[3], sqrt(par[2]^2 + reporters_per_state[i] * par[4]^2) + position[i] * par[5])
    end
    d
end
