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

function set_b_spatial(trace, params, reporters_per_state, probfn::Function, Ns, Np)
    # N = length(reporters_per_state)
    d = probfn(params, reporters_per_state, N)
    b = Ones(Ns * Np, length(trace))
    t = 1
    for obs in trace
        for j in 1:Ns
            for k in 1:Np
                for i in eachindex(obs)
                    b[l, t] *= pdf(d[j, k, i], obs[i])
                end
            end
        end
        t += 1
    end
    return b
end

function prob_Gaussian_spatial(par, reporters_per_state, position)
    Ns = length(reporters_per_state)
    Np = length(position)
    d = Array{Distribution{Univariate,Continuous}}(undef, Ns, Np)
    for j in 1:Ns
        for k in 1:Np
            for i in 1:Np
                σ = sqrt(par[2]^2 + reporters_per_state[i] * par[4]^2 * Int(k==i))
                d[j, k, i] = Normal(par[1] + reporters_per_state[j] * par[3] * Int(k==i), σ)
            end
        end
    end
end

function forward_spatial(a, b, p0, N, T)
    α = zeros(N, T)
    C = Vector{Float64}(undef, T)
    α[:, 1] = p0 .* b[:, 1]
    C[1] = 1 / sum(α[:, 1])
    α[:, 1] *= C[1]
    for t in 2:T
        for j in 1:N
            for i in 1:N
                α[j, t] += α[i, t-1] * a[i, j] * b[j, t]
            end
        end
        C[t] = 1 / sum(α[:, t])
        α[:, t] *= C[t]
    end
    return α, C
end