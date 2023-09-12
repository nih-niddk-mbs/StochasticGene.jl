
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
using Test

include("common.jl")
include("transition_rate_matrices_transpose.jl")
include("chemical_master_new.jl")
include("simulator.jl")
include("utilities.jl")
include("metropolis_hastings_trace.jl")
include("rna.jl")
include("hmm.jl")
include("trace.jl")
include("io.jl")
include("biowulf.jl")
include("fit.jl")
include("genetrap.jl")
include("analysis.jl")

"""
test(r,transitions,G,nhist,nalleles,onstates,bins)

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
data = test_data(trace,interval);
model = test_model(r,par,transitions,G,onstates);
options = test_options(1000);
@time fit,stats,measures = run_mh(data,model,options);

"""

function test(r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins, total = 100000000,tol=1e-6)
    OFF, ON, mhist = simulator(r, transitions, G, R, S, nhist, nalleles, insertstep=insertstep, onstates=onstates, bins=bins,totalsteps = total,tol=tol)
    modelOFF, modelON, histF= test_chem(r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins)
    OFF, ON, mhist, modelOFF, modelON, histF
end

function test_chem(r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins)
    if isempty(onstates)
        onstates = on_states(G, R, S, insertstep)
    end
    components = make_components_MTAI(transitions, G, R, S, insertstep, onstates, nhist, r[num_rates(transitions,R,S,insertstep)])
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    M = make_mat_M(components.mcomponents, r)
    histF = steady_state(M, components.mcomponents.nT, nalleles, nhist)
    pss = normalized_nullspace(T)
    nonzeros = nonzero_rows(TI)
    base = S > 0 ? 3 : 2
    nT = G * base^R
    SAinit = init_SA(r, onstates, components.tcomponents.elementsT, pss)
    SIinit = init_SI(r, onstates, components.tcomponents.elementsT, pss, nonzeros)
    offstates = off_states(nT, onstates)
    modelON = ontimePDF(bins, TA, offstates, SAinit)
    modelOFF = offtimePDF(bins, TI[nonzeros, nonzeros], nonzero_states(onstates, nonzeros), SIinit)
    modelOFF, modelON, histF
end

test_sim(r, transitions, G, R, S, nhist, nalleles, onstates, bins) = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates, bins=bins)

function test_fit_rna(; gene="CENPL", cell="HCT116", fish=false, G=2, nalleles=2, nsets=1, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=(), transitions=([1, 2], [2, 1]), ejectprior=0.05, r=[0.01, 0.1, 1.0, 0.01006327034802035], decayrate=0.01006327034802035, datacond="MOCK", datafolder="data/HCT116_testdata", label="scRNA_test", root=".")
    data = data_rna(gene, datacond, datafolder, fish, label)
    model = model_rna(data, r, G, nalleles, nsets, propcv, fittedparam, fixedeffects, transitions, decayrate, ejectprior)
    options = MHOptions(100000, 0, 0, 1000.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, 1)
    return stats.medparam, fits.llml, model
end


function test_fit_histograms(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1,rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=[fill(0.01, length(rtarget) - 1); rtarget[end]], nsamples=1000, nhist=20, nalleles=2, onstates=Int[], bins=collect(0:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0)
    OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    data = RNALiveCellData("test", "test", nhist, mhist, bins[2:end], OFF[1:end-1], ON[1:end-1])
    model = model_genetrap("", rinit, transitions, G, R, S, insertstep, fittedparam, nalleles, nhist, priorcv, propcv, onstates, rnatype)
    options = MHOptions(nsamples, 0, 0, 1000.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end


function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.0,50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S,insertstep)); [20, 5, 100, 10]], nsamples=1000, onstates=[G], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S,insertstep)-1); num_rates(transitions, R, S,insertstep)+1:num_rates(transitions, R, S,insertstep)+4], propcv=0.01, cv=100.0, interval=1.0)
    traces = simulate_trace_vector(rtarget,transitions,G,R,S,interval,totaltime,ntrials)
    data = trace_data(traces, interval)
    model = trace_model(rinit, transitions, G, R, S, fittedparam)
    options = trace_options(samplesteps=nsamples)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end


function test_init(r,transitions,G,R,S,insertstep)
    onstates=on_states(G,R,S,insertstep)
    indices = set_indices(length(transitions), R, S,insertstep)
    base = S > 0 ? 3 : 2
	nT = G * base^R
    components =make_components_MTAI(transitions, G, R, S, insertstep, onstates, 2, r[num_rates(transitions,R,S,insertstep)], "")
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    pss = normalized_nullspace(T)
    nonzeros = nonzero_rows(TI)
	SAinit = init_SA(r,onstates,components.tcomponents.elementsT,pss)
	SIinit = init_SI(r,onstates,components.tcomponents.elementsT,pss,nonzeros)
    
    return SIinit,SAinit,onstates,components.tcomponents,pss,nonzeros,T,TA,TI[nonzeros,nonzeros]
end


function histogram_model(r, transitions, G::Int, R::Int, S::Int, insertstep::Int,nalleles::Int, nhist::Int, fittedparam, onstates, propcv, cv)
    ntransitions = length(transitions)
    if length(r) != num_rates(transitions, R, S, insertstep)
        throw("r does not have", num_rates(transitions, R, S,insertstep), "elements")
    end
    method = 1
    if typeof(cv) <: Real
        rcv = fill(cv, length(r))
    else
        rcv = cv
    end
    d = distribution_array(log.(r[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
    if R > 0
        components = make_components_MTAI(transitions, G, R, S, insertstep, on_states(G, R, S, insertstep), nhist, r[num_rates(transitions,R,S,insertstep)])
        return GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),Int}(G, R, S, nalleles, "", r, d, propcv, fittedparam, method, transitions, components,0)
    else
        components = make_components(transitions, G, r, nhist, set_indices(ntransitions), onstates)
        return GMmodel{typeof(r),typeof(d),Float64,typeof(fittedparam),typeof(method),typeof(components)}(G, nalleles, r, d, proposal, fittedparam, method, transitions, components, onstates)
    end
end