
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
include("transition_rate_matrices.jl")
include("chemical_master.jl")
include("simulator.jl")
include("utilities.jl")
include("metropolis_hastings.jl")
include("trace.jl")
include("hmm.jl")
include("io.jl")
include("biowulf.jl")
include("fit.jl")
include("analysis.jl")
include("biowulf.jl")

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

function test(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.01:0.1:10), collect(0.01:0.1:10)], total=1000000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])
    hs = test_sim(r, transitions, G, R, S, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    l = 0
    for i in eachindex(h)
        l += norm(h[i] - hs[i])
    end
    return h, hs, l, isapprox(make_array(h),make_array(hs),rtol=0.01)
end

function test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    for i in eachindex(onstates)
        if isempty(onstates[i])
            # onstates[i] = on_states(onstates[i], G, R, S, insertstep)
            onstates[i] = on_states(G, R, S, insertstep)
        end
        onstates[i] = Int64.(onstates[i])
    end
    components = make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nRNA, r[num_rates(transitions, R, S, insertstep)], "")
    likelihoodarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)
end

function test_sim(r, transitions, G, R, S, nhist, nalleles, onstates, bins, total, tol)
    h = Vector{Vector}(undef, 5)
    h[3], h[2], h[1] = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates[1], bins=bins[1], totalsteps=total, tol=tol)
    h[5], h[4], _ = simulator(r, transitions, G, R, S, nhist, nalleles, onstates=onstates[2], bins=bins[2], totalsteps=total, tol=tol)
    h
end

function test_fit_rna(; gene="CENPL", cell="HCT116", fish=false, G=2, nalleles=2, nsets=1, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=(), transitions=([1, 2], [2, 1]), ejectprior=0.05, r=[0.01, 0.1, 1.0, 0.01006327034802035], decayrate=0.01006327034802035, datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".")
    data = data_rna(gene, datacond, datapath, fish, label)
    model = model_rna(data, r, G, nalleles, nsets, propcv, fittedparam, fixedeffects, transitions, decayrate, ejectprior)
    options = MHOptions(100000, 0, 0, 1000.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nworkers())
    return stats.medparam, fits.llml, model
end

function test_fit_histograms(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=[fill(0.01, length(rtarget) - 1); rtarget[end]], nsamples=1000, nhist=20, nalleles=2, onstates=Int[], bins=collect(0:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="")
    OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    h = test_sim(r, transitions, G, R, S, nhist, nalleles, onstates, bins, total, tol)
    data = RNAOnOffData("test", "test", nhist, mhist, bins[2:end], OFF[1:end-1], ON[1:end-1])
    model = model_genetrap("", rinit, transitions, G, R, S, insertstep, fittedparam, nalleles, nhist, priorcv, propcv, onstates, splicetype)
    options = MHOptions(nsamples, 0, 0, 1000.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nworkers())
    model = model_genetrap("", get_rates(fits.parml, model), transitions, G, R, S, insertstep, fittedparam, nalleles, nhist, priorcv, propcv, onstates, splicetype)
    fits, stats, measures, data, model, options
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.0, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep)); [20, 5, 100, 10]], nsamples=1000, onstates=[G], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4], propcv=0.01, cv=100.0, interval=1.0)
    traces = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = trace_data(traces, interval)
    model = trace_model(rinit, transitions, G, R, S, insertstep, fittedparam)
    options = trace_options(samplesteps=nsamples)
    fits, stats, measures = run_mh(data, model, options, nworkers())
    model = trace_model(get_rates(fits.parml, model), transitions, G, R, S, insertstep, fittedparam)
    fits, stats, measures, data, model, options
end

function test_init(r, transitions, G, R, S, insertstep)
    onstates = on_states(G, R, S, insertstep)
    indices = set_indices(length(transitions), R, S, insertstep)
    base = S > 0 ? 3 : 2
    nT = G * base^R
    components = make_components_MTAI(transitions, G, R, S, insertstep, onstates, 2, r[num_rates(transitions, R, S, insertstep)], "")
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    pss = normalized_nullspace(T)
    nonzeros = nonzero_rows(TI)
    SAinit = init_SA(r, onstates, components.tcomponents.elementsT, pss)
    SIinit = init_SI(r, onstates, components.tcomponents.elementsT, pss, nonzeros)

    return SIinit, SAinit, onstates, components.tcomponents, pss, nonzeros, T, TA, TI[nonzeros, nonzeros]
end

function histogram_model(r, transitions, G::Int, R::Int, S::Int, insertstep::Int, nalleles::Int, nhist::Int, fittedparam, onstates, propcv, cv)
    ntransitions = length(transitions)
    if length(r) != num_rates(transitions, R, S, insertstep)
        throw("r does not have", num_rates(transitions, R, S, insertstep), "elements")
    end
    method = 1
    if typeof(cv) <: Real
        rcv = fill(cv, length(r))
    else
        rcv = cv
    end
    d = distribution_array(log.(r[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
    if R > 0
        components = make_components_MTAI(transitions, G, R, S, insertstep, on_states(G, R, S, insertstep), nhist, r[num_rates(transitions, R, S, insertstep)])
        return GRSMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),Int}(G, R, S, nalleles, "", r, d, propcv, fittedparam, method, transitions, components, 0)
    else
        components = make_components(transitions, G, r, nhist, set_indices(ntransitions), onstates)
        return GMmodel{typeof(r),typeof(d),Float64,typeof(fittedparam),typeof(method),typeof(components)}(G, nalleles, r, d, proposal, fittedparam, method, transitions, components, onstates)
    end
end

# """
# 0.006249532442813658
# 0.01590032848543878
# 3.3790567344108426
# 0.7267835250522062
# 0.0147184704573748
# 0.008557372599505498
# 61.15683730665482
# 20.930539320417463
# 254.92857617350322
# 213.0382725610589
# 0.6643622910563245
# """