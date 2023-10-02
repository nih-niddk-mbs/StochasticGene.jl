### trace.jl

"""
    simulate_trace_vector(r, par, transitions, G, R, onstates, interval, steps, ntrials)

return vector of traces
"""
function simulate_trace_vector(r, transitions, G, R, S, interval, totaltime, ntrials; insertstep=1, onstates=Int[], reporterfn=sum)
    trace = Array{Array{Float64}}(undef, ntrials)
    for i in eachindex(trace)
        trace[i] = simulator(r[1:end-4], transitions, G, R, S, 1, 1, insertstep=insertstep, onstates=onstates, traceinterval=interval, totaltime=totaltime, par=r[end-3:end])[1:end-1, 2]
    end
    trace
end

"""
    simulate_trace(r,transitions,G,R,interval,totaltime,onstates=[G])

simulate a trace
"""
simulate_trace(r, transitions, G, R, S, interval, totaltime; insertstep=1, onstates=Int[], reporterfn=sum) = simulator(r, transitions, G, R, S, 2, 1, insertstep=insertstep, onstates=onstates, traceinterval=interval, reporterfn=reporterfn, totaltime=totaltime, par=r[end-4:end])[1:end-1, :]

"""
    trace_data(trace, interval)

set trace data
"""
function trace_data(trace, interval,nascent = 0.)
    if nascent > 0
        return TraceNascentData("trace", "test", interval, trace, nascent)
    else
        return TraceData("trace", "test", interval, trace)     
    end
end

"""
    trace_model(r::Vector, transitions::Tuple, G, R, fittedparam; onstates=[G], propcv=0.05, f=Normal, cv=1.)

set models
"""
function trace_model(r::Vector, transitions::Tuple, G, R, S, insertstep, fittedparam; noiseparams=5, probfn=prob_GaussianMixture, weightind=5, propcv=0.05, priorprob=Normal, priormean=[fill(.1, num_rates(transitions, R, S, insertstep)); fill(100,noiseparams-1);.9], priorcv=[fill(10, num_rates(transitions, R, S, insertstep)); fill(.5, noiseparams)], fixedeffects=tuple(), onstates::Vector=[G],genetrap=false,nhist=20,nalleles=2)
    if length(r) != num_rates(transitions,R,S,insertstep) + noiseparams
        throw("rate wrong length")
    end
    d = trace_prior(priormean, priorcv, fittedparam, num_rates(transitions, R, S, insertstep)+weightind, priorprob)
    method = 1
    if S > 0
        S = R
    end
    if genetrap
        #components = make_components_MTAI(transitions, G, R, S, insertstep, on_states(G, R, S, insertstep), nhist, r[num_rates(transitions, R, S, insertstep)])
        components = make_components_MT(transitions, G, R, S, insertstep, nhist, rr[num_rates(transitions, R, S, insertstep)])
    else
        components = make_components_T(transitions, G, R, S, insertstep, "")
    end

    # println(reporters)
    if R > 0
        reporter = ReporterComponents(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions,R,S,insertstep) + weightind)
        if isempty(fixedeffects)
            return GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(G, R, S, insertstep, nalleles, "", r, d, propcv, fittedparam, method, transitions, components, reporter)
        else
            return GRSMfixedeffectsmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(G, R, S, nalleles, "", r, d, propcv, fittedparam, fixedeffects, method, transitions, components, reporter)
        end
    else
        return GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G, nalleles, r, d, propcv, fittedparam, method, transitions, components, onstates)
    end
end

"""
    trace_options(samplesteps::Int=100000, warmupsteps=0, annealsteps=0, maxtime=1000.0, temp=1.0, tempanneal=100.0)

set options
"""
function trace_options(; samplesteps::Int=1000, warmupsteps=0, annealsteps=0, maxtime=1000.0, temp=1.0, tempanneal=1.0)

    MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)

end

"""
    trace_prior(r,fittedparam,f=Normal)

return prior distribution
"""
function trace_prior(r, rcv, fittedparam, ind, f=Normal)
    if typeof(rcv) <: Real
        rcv = rcv * ones(length(r))
    end
    distribution_array([logv(r[1:ind-1]);logit(r[ind:end])][fittedparam]  , sigmalognormal(rcv[fittedparam]), f)
    # distribution_array(transform_array(r,ind,fittedparam,logv,logit) , sigmalognormal(rcv[fittedparam]), f)
    # distribution_array(mulognormal(r[fittedparam],rcv[fittedparam]), sigmalognormal(rcv[fittedparam]), f)
end

# """
#     read_tracefiles(path::String,cond::String,col=3)

# read tracefiles
# """
# function read_tracefiles(path::String, cond::String, delim::AbstractChar, col=3)
#     readfn = delim == ',' ? read_tracefile_csv : read_tracefile
#     read_tracefiles(path, cond, readfn, col)
# end

# function read_tracefiles(path::String, cond::String="", readfn::Function=read_tracefile, col=3)
#     traces = Vector[]
#     for (root, dirs, files) in walkdir(path)
#         for file in files
#             target = joinpath(root, file)
#             if occursin(cond, target)
#                 push!(traces, readfn(target, col))
#             end
#         end
#     end
#     set = sum.(traces)
#     traces[unique(i -> set[i], eachindex(set))]  # only return unique traces
# end

# """
#     read_tracefile(target::String,col=3)

# read single trace file
# """
# read_tracefile(target::String, col=3) = readdlm(target)[:, col]


# read_tracefile_csv(target::String, col=3) = readdlm(target, ',')[:, col]