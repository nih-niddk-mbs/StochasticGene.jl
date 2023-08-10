
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
include("trace.jl")

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

function test(r, transitions, G, R, S, nhist, nalleles, onstates, range)
    OFF, ON, mhist = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates, range=range)
    modelOFF, modelON, histF = test_chem(r, transitions, G, R, S, nhist, nalleles, onstates, range)
    OFF, ON, mhist, modelOFF, modelON, histF
end

function test_chem(r, transitions, G, R, S, nhist, nalleles, onstates, range)
    if R == 0
        components = make_components(transitions, G, r, nhist, set_indices(length(transitions)), onstates)
    elseif S == 0
        components = make_components(transitions, G, R, r, nhist, set_indices(length(transitions), R))
    else
        components = make_components(transitions, G, R, r, nhist, "", set_indices(length(transitions), R, R))
    end
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    M = make_mat_M(components.mcomponents, r)
    if R == 0
        modelOFF, modelON = offonPDF(TA, TI, range, r, G, transitions, onstates)
    else
        modelOFF, modelON = offonPDF(T, TA, TI, range, r, G, R)
    end
    histF = steady_state(M, components.mcomponents.nT, nalleles, nhist)
    modelOFF, modelON, histF
end

test_sim(r, transitions, G, R, S, nhist, nalleles, onstates, range) = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates, range=range)

function test_fit_rna(; gene = "CENPL",cell = "HCT116", fish = false,G = 2,nalleles = 2,nsets = 1,propcv = 0.05,fittedparam = [1, 2, 3],fixedeffects = (),transitions = ([1, 2], [2, 1]), ejectprior = 0.05, r = [0.01, 0.1, 1.0, 0.01006327034802035],decayrate = 0.01006327034802035,datacond = "MOCK",datafolder = "data/HCT116_testdata",label = "scRNA_test",root=".")
    data = data_rna(gene, datacond, datafolder, fish, label)
    model = model_rna(data, r, G, nalleles, nsets, propcv, fittedparam, fixedeffects, transitions, decayrate, ejectprior)
    options = MHOptions(100000, 0, 0, 1000.0, 1.0, 1.0)
    fit, stats, measures = run_mh(data, model, options, 1)
    return stats.medparam, fit.llml, model
end


function test_fit_histograms(; rtarget=[0.02, 0.1, 0.5, 0.2, 0.01], rinit=[ones(num_rates(transitions,R,S)-1);rtarget[end]], nsamples=1000, G=2, R=1, S=0, nhist=20, nalleles=2, onstates=[G], bins=collect(0:1.0:200.0), transitions=([1, 2], [2, 1]), fittedparam=collect(1:num_rates(transitions,R,S)-1), propcv=0.05, cv=10.0)
    OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    data = RNALiveCellData("test", "test", nhist, mhist, bins[2:end], OFF[1:end-1], ON[1:end-1])
    println(rinit)
    model = histogram_model(rinit, transitions, G, R, S, nalleles, data.nRNA, fittedparam, onstates, propcv, cv)
    options = MHOptions(nsamples, 0, 0, 1000.0, 1.0, 1.0)
    fit, stats, measures = run_mh(data, model, options)
    fit, stats, measures, data, model, options
end


function test_fit_trace(; G=2, R=1, S=0,transitions=([1, 2], [2, 1]),rtarget=[0.02, 0.1, 0.5, 0.2, 0.01], par=[50,20,200,65], rinit=[ones(num_rates(transitions,R,S)); [20, 5, 100, 10]], nsamples=1000,  onstates=[G],  totaltime=1000., ntrials=10, fittedparam=[collect(1:num_rates(transitions,R,S)-1);num_rates(transitions,R,S)+1:num_rates(transitions,R,S)+4], propcv=0.05, cv=10.0, interval=2.0)
    traces = simulate_trace_vector(rtarget, par, transitions, G, R, S, interval, totaltime, onstates, ntrials)
    data = trace_data(traces, interval)
    model = trace_model(rinit, transitions, G, R, S, fittedparam)
    options = trace_options(nsamples)
    fit, stats, measures = run_mh(data, model, options)
    fit, stats, measures, data, model, options
end


function histogram_model(r, transitions, G::Int, R::Int, S::Int, nalleles::Int, nhist::Int, fittedparam,  onstates, propcv, cv)
    ntransitions = length(transitions)
    if length(r) != num_rates(transitions,R,S)
        throw("r does not have",num_rates(transitions,R,S),"elements")
    end
    method = 1
    if typeof(cv) <: Real
        rcv = fill(cv, length(r))
    else
        rcv = cv
    end
     d = distribution_array(log.(r[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
    if S > 0
        components = make_components(transitions, G, R, S, r, nhist, set_indices(ntransitions, R, S))
		return  GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G, R, S, nalleles, "", r, d, propcv, fittedparam, method, transitions, components)
	elseif R > 0
        components = make_components(transitions, G, R, r, nhist, set_indices(ntransitions, R))
        return GRMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G, R, nalleles, "", r, d, propcv, fittedparam, method, transitions, components)
	else
        components = make_components(transitions, G, r, nhist, set_indices(ntransitions), onstates)
        return GMmodel{typeof(r),typeof(d),Float64,typeof(fittedparam),typeof(method),typeof(components)}(G, nalleles, r, d, proposal, fittedparam, method, transitions, components, onstates)
  
	end
end