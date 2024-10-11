# Spatial.jl
#
# Experimental code for including grid information from images of transcriptional bursting
"""
    kronecker_delta(i, j)

Computes the Kronecker delta of two integers.

# Arguments
- `i`: The first integer.
- `j`: The second integer.

# Returns
- `Int`: Returns 1 if `i` equals `j`, otherwise returns 0.
"""
# kronecker_delta(i, j) = i == j ? 1 : 0

# """
#     read_tracefiles_grid(path::String, label::String, start::Int, stop::Int, col=3)

# Reads trace files from a specified directory that match a given label and extracts data from them.

# # Arguments
# - `path::String`: The directory path to search for trace files.
# - `label::String`: The label to match in the filenames.
# - `start::Int`: The starting index for reading the trace data.
# - `stop::Int`: The stopping index for reading the trace data.
# - `col::Int`: The column index to read from each trace file (default is 3).

# # Returns
# - `Vector`: A vector of unique traces read from the files.
# """
# function read_tracefiles_grid(path::String, label::String, start::Int, stop::Int, col=3)
#     traces = Vector[]
#     if isempty(path)
#         return traces
#     else
#         for (root, dirs, files) in walkdir(path)
#             for file in files
#                 if occursin(label, file) && ~occursin(".DS_Store", file)
#                     t = read_tracefile_grid(joinpath(root, file), start, stop, col)
#                     ~isempty(t) && push!(traces, t)
#                 end
#             end
#         end
#         set = sum.(traces)
#         return traces[unique(i -> set[i], eachindex(set))]  # only return unique traces
#     end
# end

# """
#     set_b_grid(trace, params, reporters_per_state, probfn::Function, Nstate, Ngrid)

# Calculates the probability matrix `b` for given trace data using a specified probability function.

# # Arguments
# - `trace`: The trace data.
# - `params`: Parameters for the probability function.
# - `reporters_per_state`: Number of reporters per state.
# - `probfn::Function`: The probability function to use.
# - `Nstate`: Number of states.
# - `Ngrid`: Number of positions.

# # Returns
# - `Matrix`: The probability matrix `b`.
# """
# function set_b_grid_v(trace, params, reporters_per_state, probfn::Function, Nstate, Ngrid)
#     d = probfn(params, reporters_per_state, Nstate, Ngrid)
#     set_b_grid(trace, d, Nstate, Ngrid)
# end
# function set_b_grid_v(trace, d, Nstate, Ngrid)
#     b = ones(Nstate, Ngrid, length(trace))
#     t = 1
#     for obs in eachcol(trace)
#         for j in 1:Nstate
#             for k in 1:Ngrid
#                 for l in eachindex(obs)
#                     b[j, k, t] *= StochasticGene.pdf(d[j, k, l], obs[l])
#                 end
#             end
#         end
#         t += 1
#     end
#     return b
# end


# function forward_grid(as, ap, b, p0, Nstate, Ngrid, T)
#     α = zeros(Nstate, Ngrid, T)
#     C = Vector{Float64}(undef, T)
#     α[:, :, 1] = p0 .* b[:, :, 1]
#     C[1] = 1 / sum(α[:, :, 1])
#     α[:, :, 1] *= C[1]
#     for t in 2:T
#         for l in 1:Ngrid
#             for k in 1:Nstate
#                 for j in 1:Ngrid
#                     for i in 1:Nstate
#                         α[i, j, t] += α[k, l, t-1] * as[k, i] * ap[l, j] * b[i, j, t]
#                     end
#                 end
#             end
#         end
#         C[t] = 1 / sum(α[:, :, t])
#         α[:, :, t] *= C[t]
#     end
#     return α, C
# end

# function distance(i, j, Ngrid)
#     xi = rem(i - 1, Ngrid)
#     yi = div(i - 1, Ngrid)
#     xj = rem(j - 1, Ngrid)
#     yj = div(j - 1, Ngrid)
#     return sqrt((xi - xj)^2 + (yi - yj)^2)
# end

# function make_a_grid(param, Ngrid)
#     as = zeros(Ngrid, Ngrid)
#     d = zeros(Ngrid, Ngrid)
#     for i in 1:Ngrid
#         for j in 1:Ngrid
#             as[i, j] = exp(-grid_distance(i, j, div(Ngrid, 2))^2 / (2 * param^2))
#             d[i, j] = grid_distance(i, j, div(Ngrid, 2))
#         end
#     end
#     as ./ sum(as, dims=2)
# end



# function ll_hmm_grid(r, p, Nstate, Ngrid, components::StochasticGene.TRGComponents, n_noiseparams::Int, reporters_per_state, probfn, interval, trace)
#     a_grid = make_a_grid(p, Ngrid)
#     a, p0 = StochasticGene.make_ap(r, interval, components)
#     d = probfn(r[end-n_noiseparams+1:end], reporters_per_state, Nstate, Ngrid)
#     logpredictions = Array{Float64}(undef, 0)
#     for t in trace[1]
#         T = length(t)
#         b = set_b_grid_v(t, d, Nstate, Ngrid)
#         _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
#         push!(logpredictions, sum(log.(C)))
#     end
#     sum(logpredictions), logpredictions
# end

function test_ll_grid(; r=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], p=0.2, Ngrid=4, transitions=([1, 2], [2, 1]), G=2, R=1, S=1, insertstep=1, totaltime=1000.0, interval=1.0,ntrials=10)
    Nstate = num_rates(transitions, R, S, insertstep)
    # traces = simulator(r, transitions, G, R, S, insertstep, traceinterval=interval, nhist=0, totaltime=totaltime, reporterfn=sum, a_grid=StochasticGene.make_a_grid(1.0, 4))[1]
    traces=StochasticGene.simulate_trace_vector(r,transitions,G,R,S,insertstep,interval,totaltime,ntrials,a_grid=a_grid)
    components = StochasticGene.make_components_TRG(transitions, G, R, S, insertstep, "")
    reporters_per_state = StochasticGene.num_reporters_per_state(G, R, S, insertstep)
    StochasticGene.ll_hmm_grid(r, p, Nstate, Ngrid, components, 4, reporters_per_state, StochasticGene.prob_Gaussian_grid, interval, (traces, [], 0.0, 1000)), StochasticGene.ll_hmm(r, Nstate, components, 4, reporters_per_state, StochasticGene.prob_Gaussian, StochasticGene.off_states(G, R, S, insertstep), interval, (traces, [], 0.0, 1000))

    # trace = StochasticGene.simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    # data = StochasticGene.TraceData("tracegrid", "test", interval, (trace, [], weight, nframes))
    # # model = StochasticGene.load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple())
    # elongationtime = StochasticGene.mean_elongationtime(rtarget, transitions, R)
    # priormean = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, elongationtime)
    # model = StochasticGene.load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], 1.0, propcv, "", prob_Gaussian, noisepriors, tuple(), tuple(), 1)
    # options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 0., 1.0)
    # fits, stats, measures = run_mh(data, model, options)
    # StochasticGene.get_rates(fits.parml, model), rtarget, fits, stats, measures, model, data
end

function test_fit_trace_grid(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70])

    trace = simulator(rtarget, transitions, G, R, S, insertstep, onstates=onstates, traceinterval=interval, nhist=0, totaltime=totaltime, reporterfn=reporterfn, ap=make_ap(1.0, G))[1]
    components = StochasticGene.make_components_TRG(transitions, G, R, S, insertstep, splicetype)
    reporters_per_state = StochasticGene.num_reporters_per_state(G, R, S, insertstep)
    a, p0 = make_ap(r, interval, components)
    ll_hmm_grid
    # trace = StochasticGene.simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    # data = StochasticGene.TraceData("tracegrid", "test", interval, (trace, [], weight, nframes))
    # # model = StochasticGene.load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple())
    # elongationtime = StochasticGene.mean_elongationtime(rtarget, transitions, R)
    # priormean = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, elongationtime)
    # model = StochasticGene.load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], 1.0, propcv, "", prob_Gaussian, noisepriors, tuple(), tuple(), 1)
    # options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    # fits, stats, measures = run_mh(data, model, options)
    # StochasticGene.get_rates(fits.parml, model), rtarget, fits, stats, measures, model, data
end
# components = StochasticGene.make_components_TRG(transitions, G, R, S, insertstep, splicetype)
#  Qtr = make_mat_TRG(components, r)
# a, p0 = make_ap(r, interval, components)
# reporters_per_state = StochasticGene.num_reporters_per_state(G, R, S, insertstep)