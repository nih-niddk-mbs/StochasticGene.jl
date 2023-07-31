
using Distributions
using StatsBase
using Statistics
using DelimitedFiles
using Dates
using Distributed
using LinearAlgebra
using Plots
using SparseArrays
using DifferentialEquations
using LSODA
using DataFrames
using FFTW
using Downloads
using CSV
using MultivariateStats
using Optim

include("common.jl")
include("transition_rate_matrices_transpose.jl")
include("chemical_master_new.jl")
include("simulator.jl")
include("utilities.jl")
include("metropolis_hastings_trace.jl")
include("rna.jl")
include("hmm.jl")

"""
test(r,transitions,G,nhist,nalleles,onstates,range)

returns dwell time and mRNA histograms for simulations and chemical master equation solutions

G = 2
r = [.01,.02,.03,.04,.001,.001]
transitions = ([1,2],[2,1],[2,3],[3,1])
OFF,ON,mhist,modelOFF,modelON,histF = test(r,transitions,G,20,1,[2],collect(1.:200))

e.g.

julia> incluce("test.jl")

julia> G = 3;

julia> r = [.01,.02,.03,.04,.001,.001];

julia> transitions = ([1,2],[2,1],[2,3],[3,1]);

julia> onstates = [2,3]

julia> OFFsim,ONsim,mRNAsim,OFFchem,ONchem,mRNAchem = test(r,transitions,G,20,1,onstates,collect(1.:200));


"""

"""
eg fit model to trace


r = [0.05,.1,.001,.001];
transitions = ([1,2],[2,1]);
G = 2;
R = 0;
onstates = [2];
interval = 5/3;
par=[30, 14, 200, 75];
trace=make_trace(r,par,transitions,G,R,onstates,interval,50,10);
data = test_data(trace,interval);
model = test_model(r,par,transitions,G,onstates);
options = test_options(1000);
@time fit,waic = metropolis_hastings(data,model,options);
@time fit,stats,measures = run_mh(data,model,options);

"""

function simulate_trace_vector(r, par, transitions, G, R, onstates, interval, steps, ntrials)
    trace = Array{Array{Float64}}(undef, ntrials)
    for i in eachindex(trace)
        trace[i] = simulator(r, transitions, G, R, 0, 1, 1, onstates=onstates, traceinterval=interval, totalsteps=steps, par=par)[1:end-1, 2]
    end
    trace
end

simulate_trace(r,transitions,G,R,interval,totaltime,onstates=[G]) = simulator(r, transitions, G, R, 0, 1, 1, onstates=onstates, traceinterval=interval, totaltime=totaltime, par=r[end-3:end])[1:end-1, :]

function trace_data(trace, interval)
    TraceData("trace", "test", interval, trace)
end

function trace_model(r::Vector, par::Vector, transitions::Tuple, G, R; onstates=[G], propcv=0.05, f=Normal, cv=1.0)
    ntransitions = length(transitions)
	npars = length(par)
	fittedparam = [1:ntransitions+R+1; ntransitions+R+3:ntransitions+R+2+npars]
	r = vcat(r, par)
	test_model(r, transitions, G, R, fittedparam, onstates=onstates, propcv=propcv, f=f,cv=cv)
 end

function trace_model(r::Vector, transitions::Tuple, G, R, fittedparam; onstates=[G], propcv=0.05, f=Normal, cv=1.)
	d = trace_prior(r, fittedparam,f,cv)
	method = 1
	if R == 0
		components = make_components_T(transitions, G)
		return GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G, 1, r, d, propcv, fittedparam, method, transitions, components, onstates)
	else
		components = make_components_T(transitions, G, R)
		return GRMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G,R,1,"",r,d,propcv,fittedparam,method,transitions,components,on_states(G,R))
	end
end

function trace_options(samplesteps::Int=100000, warmupsteps=0, annealsteps=0, maxtime=1000.0, temp=1.0, tempanneal=100.0)

    MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)

end

"""
    trace_prior(r,fittedparam,f=Normal)

TBW
"""
function trace_prior(r,fittedparam,f=Normal,cv = 1.)
	rcv = cv * ones(length(r))
	distribution_array(log.(r[fittedparam]),sigmalognormal(rcv[fittedparam]),f)
end

function test(r, transitions, G, R, S,nhist, nalleles, onstates, range)
    OFF, ON, mhist = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates, range=range)
    modelOFF, modelON, histF = test_cm(r, transitions, G, R, S, nhist, nalleles, onstates, range)
    OFF, ON, mhist, modelOFF, modelON, histF
end

function test_cm(r, transitions, G, R, S, nhist, nalleles, onstates, range)
	if R == 0
    	components = make_components(transitions, G, r, nhist, set_indices(length(transitions)), onstates)
	elseif S==0
		components = make_components(transitions, G, R, r, nhist, set_indices(length(transitions),R))
	else
		components = make_components(transitions, G, R, r, nhist, "", set_indices(length(transitions),R,R))
	end
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    M = make_mat_M(components.mcomponents, r)
	if R == 0
		 modelOFF, modelON = offonPDF(T, TA, TI, range, r, G, transitions, onstates)
	else
		modelOFF, modelON = offonPDF(T, TA, TI, range, r, G, R)
	end
    histF = steady_state(M, components.mcomponents.nT, nalleles, nhist)
    modelOFF, modelON, histF
end

test_sim(r, transitions, G, R, S, nhist, nalleles, onstates, range) = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates, range=range)

function read_tracefiles(path::String,cond::String,col=3)
    traces = Vector[]
    for (root,dirs,files) in walkdir(path)
        for file in files
            target = joinpath(root, file)
            if occursin(cond,target)
                # println(target)
                push!(traces, read_tracefile(target,col))
            end
        end
    end
    return traces
end

read_tracefile(target::String,col=3) = readdlm(target)[:,col]



