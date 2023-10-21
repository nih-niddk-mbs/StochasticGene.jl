
# using Distributions
# using StatsBase
# using Statistics
# using DelimitedFiles
# using Dates
# using Distributed
# using LinearAlgebra
# using Plots
# using SparseArrays
# using DifferentialEquations
# using LSODA
# using DataFrames
# using FFTW
# using Downloads
# using CSV
# using MultivariateStats
# using Optim
# using Test


# include("common.jl")
# include("transition_rate_matrices.jl")
# include("chemical_master.jl")
# include("simulator.jl")
# include("utilities.jl")
# include("metropolis_hastings.jl")
# # include("trace.jl")
# include("hmm.jl")
# include("io.jl")
# include("biowulf.jl")
# include("fit.jl")
# include("analysis.jl")
# include("biowulf.jl")

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
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    hs = StochasticGene.normalize_histogram.(hs)
    l = 0
    for i in eachindex(h)
        l += StochasticGene.norm(h[i] - hs[i])
    end
    return make_array(h), make_array(hs), l, isapprox(make_array(h), make_array(hs), rtol=0.01)
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

function test_sim(r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins, total, tol)
    h = Vector{Vector}(undef, 3)
    h[3], h[2], h[1] = simulator(r, transitions, G, R, S, nhist, nalleles, insertstep=insertstep, onstates=onstates[1], bins=bins[1], totalsteps=total, tol=tol)
    for i in eachindex(onstates)[begin+1:end]
        hoff, hon, _ = simulator(r, transitions, G, R, S, nhist, nalleles, insertstep=insertstep, onstates=onstates[i], bins=bins[i], totalsteps=total, tol=tol)
        push!(h, hon)
        push!(h, hoff)
    end
    h
end

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
    data = StochasticGene.data_rna(gene,datacond,datafolder,fish,label)
    model = StochasticGene.model_rna(data,r,G,nalleles,nsets,propcv,fittedparam,fixedeffects,transitions,decayrate,ejectprior)
    options = StochasticGene.MHOptions(100000,0,0,30.,1.,100.)
    fit,stats,measures = StochasticGene.run_mh(data,model,options,1);
    return stats.meanparam, fit.llml, model
end

function test_fit_rna(; gene="CENPL", cell="HCT116", fish=false, G=2, nalleles=2, nsamples=100000, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), ejectprior=0.05, r=[0.01, 0.1, 1.0, 0.01006327034802035], decayrate=0.01006327034802035, datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".")
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, 1.0, 1.0, 1.0)
    model = load_model(data, r, Float64[], fittedparam, fixedeffects, transitions, G, 0, 0, 1, nalleles, 10.0, Int[], r[end], propcv, "", prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(100000,0,0,60.,1.,1.)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h,normalize_histogram(data.histRNA)
end

function test_fit_rnaonoff(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=fill(0.01, num_rates(transitions, R, S, insertstep)), nsamples=100000, nhist=20, nalleles=2, onstates=Int[], bins=collect(0:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="")
    # OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    OFF, ON, hRNA = simulator(rtarget, transitions, G, R, S, nhist, nalleles, onstates=onstates, bins=bins, totalsteps=3000)
    hRNA = div.(hRNA, 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, ON, OFF)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h,StochasticGene.datapdf(data),fits, stats, measures, data, model, options
end


function test_fit_rnadwelltime(; rtarget=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.01:0.1:10), collect(0.01:0.1:10)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="")
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 30000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype)
    println(typeof(data))
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], 0, 0)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nsamples, 1000, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nworkers())
    h = likelihoodfn(fits.parml, data, model)
    h,StochasticGene.datapdf(data),fits, stats, measures, data, model, options
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70, 0.1], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep)); [20, 5, 100, 10, 0.9]], nsamples=1000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+5)], propcv=0.01, cv=100.0, interval=1.0)
    traces = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, traces)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
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