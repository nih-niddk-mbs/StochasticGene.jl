
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

function fit_rna_test(root=".")
    G = 2
    gene = "CENPL"
    cell = "HCT116"
    fish = false
    nalleles = 2
    nsets = 1
    propcv = 0.05
    fittedparam = [1,2,3]
    fixedeffects = ()
    transitions = ([1,2],[2,1])
    ejectprior = 0.05
    r = [0.01, 0.1, 1.0, 0.01006327034802035]
    decayrate = 0.01006327034802035
    datacond = "MOCK"
    datafolder = "data/HCT116_testdata"
    label = "scRNA_test"
    data = data_rna(gene,datacond,datafolder,fish,label)
    model = model_rna(data,r,G,nalleles,nsets,propcv,fittedparam,fixedeffects,transitions,decayrate,ejectprior)
    options = MHOptions(100000,0,0,30.,1.,100.)
    fit,stats,measures = run_mh(data,model,options,1);
    return stats.meanparam, fit.llml, model
end


function fit_histograms_test()
	G = 2
	R = 1
	S = 0
	nhist = 20
	nalleles = 2
	onstates = [2]
	bins = collect(0:2:200.)
	r = [0.02, 0.1, 0.5, 0.2, 0.01]
	transitions =  ([1,2],[2,1])
	fittedparam = [1,2,3,4]
	propcv = 0.05
	cv = 10.
	OFF,ON,mhist = test_sim(r, transitions, G, R, S, nhist, nalleles, onstates, bins)
	data = RNALiveCellData("test","test",nhist,mhist,bins[2:end],OFF[1:end-1],ON[1:end-1])
	model = model_histogram([.1,.1,.1,.1,.01],transitions,G,R,S,nalleles,fittedparam,data,onstates,propcv,cv)
	options = MHOptions(10000, 0, 0, 0, 1., 1.)
	fit,stats,measures = run_mh(data,model,options);
end

function fit_trace_test()
	G = 2
	R = 1
	S = 0
	onstates = [2]
	r = [0.02, 0.1, 0.5, 0.2, .01]
	par = [30,10,200,65]
	transitions =  ([1,2],[2,1])
	steps = 1000
	ntrials = 10
	fittedparam = [1,2,3,4,6,7,8,9]
	propcv = 0.05
	cv = 10.
	interval = 1.
	traces = simulate_trace_vector(r, par, transitions, G, R, onstates, interval, steps, ntrials)
	data = trace_data(traces, interval)
	model = trace_model([.1,.1,.1,.1,.01,20,5,100,10], transitions, G, R)
	options = trace_options(1000);
	fit,stats,measures = run_mh(data,model,options);
end


function model_histogram(r,transitions,G::Int,R::Int,S::Int,nalleles::Int,fittedparam,data,onstates=[G],propcv = .05,cv=1.)
    ntransitions = length(transitions)
    method = 1
    if R == 0
        decayrate = r[2*G]
        d = prior_rna(r,G,1,fittedparam,decayrate,1.)
        components = make_components(transitions,G,r,data.nRNA,Indices(collect(1:ntransitions),[ntransitions+1],Int[],ntransitions + 2),onstates)
        return GMmodel{typeof(r),typeof(d),Float64,typeof(fittedparam),typeof(method),typeof(components)}(G,nalleles,r,d,proposal,fittedparam,method,transitions,components,onstates)
	elseif S == 0
        if typeof(cv) <: Real
            rcv = fill(cv,length(r))
        end
        d = distribution_array(log.(r[fittedparam]),sigmalognormal(rcv[fittedparam]),Normal)
        components = make_components(transitions,G,R,r,data.nRNA,set_indices(ntransitions,R))
        return GRMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(G,R,nalleles,"",r,d,propcv,fittedparam,method,transitions,components)
	else

    end
end