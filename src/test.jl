
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

julia> OFFsim,ONstateim,mRNAsim,OFFchem,ONchem,mRNAchem = test(r,transitions,G,20,1,onstates,collect(1.:200));


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

function make_test(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70, 0.9], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10, 0.9]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+5)], propcv=0.01, cv=100.0, interval=1.0)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0))
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_GaussianMixture, 5, 5, tuple())
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    return data, model, options
end

function test(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=1000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    hs = StochasticGene.normalize_histogram.(hs)
    l = 0
    for i in eachindex(h)
        l += StochasticGene.norm(h[i] - hs[i])
    end
    return StochasticGene.make_array(h), StochasticGene.make_array(hs), l, isapprox(StochasticGene.make_array(h), StochasticGene.make_array(hs), rtol=0.01)
end

function test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    for i in eachindex(onstates)
        if isempty(onstates[i])
            # onstates[i] = on_states(onstates[i], G, R, S, insertstep)
            onstates[i] = StochasticGene.on_states(G, R, S, insertstep)
        end
        onstates[i] = Int64.(onstates[i])
    end
    components = StochasticGene.make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nRNA, r[StochasticGene.num_rates(transitions, R, S, insertstep)], "")
    StochasticGene.likelihoodarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)
end

function test_sim(r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins, total, tol)
    StochasticGene.simulator(r, transitions, G, R, S, insertstep, nhist=nhist, nalleles=nalleles, onstates=onstates, bins=bins, totalsteps=total, tol=tol)
end

function fit_rna_test(root=".")
    G = 2
    gene = "CENPL"
    cell = "HCT116"
    fish = false
    nalleles = 2
    nsets = 1
    propcv = 0.05
    fittedparam = [1, 2, 3]
    fixedeffects = ()
    transitions = ([1, 2], [2, 1])
    ejectprior = 0.05
    r = [0.01, 0.1, 1.0, 0.01006327034802035]
    decayrate = 0.01006327034802035
    datacond = "MOCK"
    datafolder = "data/HCT116_testdata"
    label = "scRNA_test"
    data = StochasticGene.data_rna(gene, datacond, datafolder, fish, label)
    model = StochasticGene.model_rna(data, r, G, nalleles, nsets, propcv, fittedparam, fixedeffects, transitions, decayrate, ejectprior)
    options = StochasticGene.MHOptions(100000, 0, 0, 30.0, 1.0, 100.0)
    fit, stats, measures = StochasticGene.run_mh(data, model, options, 1)
    return stats.meanparam, fit.llml, model
end

function test_fit_rna(; gene="CENPL", cell="HCT116", fish=false, G=2, nalleles=2, nunits=100000, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), ejectprior=0.05, r=[0.01, 0.1, 1.0, 0.01006327034802035], decayrate=0.01006327034802035, datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".")
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, 1.0, 1.0, 1.0)
    model = load_model(data, r, Float64[], fittedparam, fixedeffects, transitions, G, 0, 0, 1, nalleles, 10.0, Int[], r[end], propcv, "", prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(100000, 0, 0, 60.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

function test_fit_rnaonoff(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=fill(0.01, num_rates(transitions, R, S, insertstep)), nunits=100000, nhist=20, nalleles=2, onstates=Int[], bins=collect(0:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="")
    # OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    h = simulator(rtarget, transitions, G, R, S, nhist, nalleles, onstates=onstates, bins=bins, totalsteps=3000)
    hRNA = div.(h[1], 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, h[2], h[3])
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nunits, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data), fits, stats, measures, data, model, options
end


function test_fit_rnadwelltime(; rtarget=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nunits=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.01:0.1:10), collect(0.01:0.1:10)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="")
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 30000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype)
    println(typeof(data))
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], 0, 0)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nunits, 1000, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nworkers())
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data), fits, stats, measures, data, model, options
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70])
    trace = StochasticGene.simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], weight, nframes))
    # model = StochasticGene.load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple())
    elongationtime = StochasticGene.mean_elongationtime(rtarget, transitions, R)
    priormean = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, elongationtime)
    model = StochasticGene.load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], 1.0, propcv, "", prob_Gaussian, noisepriors, tuple(), tuple(), 1)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget, fits, stats, measures, model, data
end


function trace_structures(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70])
    trace = StochasticGene.simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], weight, nframes))
    # model = StochasticGene.load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple())
    elongationtime = StochasticGene.mean_elongationtime(rtarget, transitions, R)
    priormean = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, elongationtime)
    model = StochasticGene.load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], 1.0, propcv, "", prob_Gaussian, noisepriors, tuple(), tuple(), 1)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    return data, model, options
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

# function test_fit_trace_hierarchical1(; G=2, R=3, S=3, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, .2,.2,.2,.2,.2, 1.0, 50, 15, 200, 70], rinit=[], nsamples=1000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=collect(1:num_rates(transitions, R, S, insertstep)-1), propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], hierarchical=(2, collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+ length(noisepriors)), tuple()), method=(1, true))
function test_fit_trace_hierarchical(; G=2, R=1, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 1.0, 50, 15, 200, 70], rinit=[], nsamples=100000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=collect(1:num_rates(transitions, R, S, insertstep)-1), propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], hierarchical=(2, collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+length(noisepriors)), tuple()), method=(1, true))
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, hierarchical[1], mean_elongationtime(rtarget, transitions, R))
    isempty(rinit) && (rinit = StochasticGene.set_rinit(rm, transitions, R, S, insertstep, noisepriors, length(data.trace[1])))
    println(rinit)
    model = load_model(data, rinit, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, hierarchical, tuple(), method)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    nrates = num_rates(transitions, R, S, insertstep)
    h1 = StochasticGene.get_rates(fits.parml, model)[1:nrates]
    h2 = rtarget[1:nrates]
    return h1, h2
end

function test_fit_tracejoint(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(2, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.0, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([100, 50, 200, 100], [100, 50, 200, 100]), maxtime=120.0)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], 0.0, 1))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    noiseparams = length.(noisepriors)
    fittedparam = StochasticGene.set_fittedparam(fittedparam, data.label, transitions, R, S, insertstep, noiseparams, coupling)
    model = load_model(data, rm, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple(), coupling, 1)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    println(nsamples)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget, fits, stats, measures, data, model, options
end

function test_mat2(r, transitions, G, R, S, insertstep, nhist=20)
    components = StochasticGene.make_components_MTRG(transitions, G, R, S, insertstep, nhist, 1.01, "")
    T = StochasticGene.make_mat_TRG(components.tcomponents, r)
    M = StochasticGene.make_mat_MRG(components.mcomponents, r)
    B = StochasticGene.make_mat_B2(components.mcomponents, r)
    return T, M, B, components
end

function test_mat(r, transitions, G, R, S, insertstep, nhist=20)
    components = make_components_MT(transitions, G, R, S, insertstep, nhist, 1.01, "")
    T = StochasticGene.make_mat_T(components.tcomponents, r)
    M = StochasticGene.make_mat_M(components.mcomponents, r)
    T, M, components
end

function test_mat_Tc(coupling, r, coupling_strength, transitions, G, R, S, insertstep)
    components = make_components_Tcoupled(coupling, transitions, G, R, S, insertstep, "")
    return make_mat_TC(components, r, coupling_strength)
end

function test_mat_Tc2(coupling, r, coupling_strength, transitions, G, R, S, insertstep)
    T = []
    for i in eachindex(G)
        c = make_components_TRG(transitions[i], G[i], R[i], S[i], insertstep[i], "")
        push!(T, make_mat_TRG(c, r[i]))
    end
    components = make_components_Tcoupled(coupling, transitions, G, R, S, insertstep, "")
    return make_mat_TC(components, r, coupling_strength), T, kron(T[1], sparse(I, size(T[2]))) + kron(sparse(I, size(T[1])), T[2])
end

function compare(p0, x, y, z)
    p = reshape(p0, 2, 2, 2, 2)
    pm = marg(p, 4)
    a = p[x, y, z, :] / pm[x, y, z]
    p1 = marg(p, 2)
    p1m = marg(p, (2, 4))
    b = p1[x, z, :] / p1m[x, z]
    return a, b, a ≈ b
end

marg(p, dims) = dropdims(sum(p, dims=dims), dims=dims)

function make_test_matrices(r::Vector{Vector{type}}, coupling, transitions, G, R, S, insertstep, coupling_strength=[1.0, 1]) where {type<:Number}
    components = StochasticGene.make_components_Tcoupled(coupling, transitions, G, R, S, insertstep)
    Tc = StochasticGene.make_mat_TC(components, r, coupling_strength)
    T, GGv, Gt, Gs, IG, IR, IT = StochasticGene.make_matvec_C(components, r)
    GG = StochasticGene.make_mat_TC(coupling_strength, GGv, Gs, Gt, IG, components.sources, components.model)
    return Tc, T, GG, components
end

function extract_components(components)
    nR_values = []
    nG_values = []
    for c in components.modelcomponents
        push!(nR_values, c.nR)
        push!(nG_values, c.nG)
    end


    result = []
    for i in 1:length(nR_values)
        push!(result, nG_values[end-i+1])
        push!(result, nR_values[end-i+1])
    end

    return tuple(result...)
end

function test_a(Tc, GG)
    ac = StochasticGene.kolmogorov_forward(Tc', 1.0)
    ag = StochasticGene.kolmogorov_forward(GG', 1.0)
    return ac, ag
end


function test_distributions(r::Vector{Vector{type}}, coupling, transitions, G, R, S, insertstep, coupling_strength=[1.0, 1]) where {type<:Number}
    components = StochasticGene.make_components_Tcoupled(coupling, transitions, G, R, S, insertstep)
    Tc = StochasticGene.make_mat_TC(components, r, coupling_strength)
    T, GG, Gt, Gs, IG, IR, IT = StochasticGene.make_matvec_C(components, r)
    ac = StochasticGene.kolmogorov_forward(Tc', 1.0)
    ag = StochasticGene.kolmogorov_forward(GG', 1.0)

    p0 = StochasticGene.normalized_nullspace(Tc)
    p01 = StochasticGene.normalized_nullspace(Tcr[1])
    p02 = StochasticGene.normalized_nullspace(Tcr[2])
    p0G = StochasticGene.normalized_nullspace(GG)

    p0R = reshape(p0, 2, 2, 2, 2)
    p01R = reshape(p01, 2, 2)
    p02R = reshape(p02, 2, 2, 2)
    pG = reshape(p0G, 2, 2)
    return p0R, p01R, p02R, pG
end

function recursive_sum(matrix, result, indices1, indices2, depth, max_depth)
    if depth == max_depth
        i, j, k, l = indices1
        i1, j1, k1, l1 = indices2
        idx1 = i + 2*j + 4*k + 8*l + 1
        idx2 = i1 + 2*j1 + 4*k1 + 8*l1 + 1
        result[i + 2*k + 1, i1 + 2*k1 + 1] += matrix[idx1, idx2]
    else
        for x in 0:1
            if depth < 4
                new_indices1 = copy(indices1)
                new_indices1[depth + 1] = x
                recursive_sum(matrix, result, new_indices1, indices2, depth + 1, max_depth)
            else
                new_indices2 = copy(indices2)
                new_indices2[depth - 4 + 1] = x
                recursive_sum(matrix, result, indices1, new_indices2, depth + 1, max_depth)
            end
        end
    end
end

function sum_overR2(matrix)
    result = zeros(4, 4)
    recursive_sum(matrix, result, [0, 0, 0, 0], [0, 0, 0, 0], 0, 8)
    return result / 4
end

function sum_overR(matrix)
    result = zeros(4, 4)

    for i in 0:1
        for j in 0:1
            for k in 0:1
                for l in 0:1
                    for i1 in 0:1
                        for j1 in 0:1
                            for k1 in 0:1
                                for l1 in 0:1
                                    result[i+2*k+1, i1+2*k1+1] += matrix[i+2j+4k+8l+1, i1+2j1+4k1+8l1+1]
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    result/4
end

function recursive_sum3(matrix, result, indices1, indices2, depth, max_depth, ranges)
    if depth == max_depth
        idx1 = sum(indices1[i] * prod(ranges[1:i-1]) for i in 1:length(indices1)) + 1
        idx2 = sum(indices2[i] * prod(ranges[1:i-1]) for i in 1:length(indices2)) + 1
        result[indices1[1] + 2 * indices1[3] + 1, indices2[1] + 2 * indices2[3] + 1] += matrix[idx1, idx2]
    else
        for x in 0:ranges[depth + 1] - 1
            if depth < length(indices1)
                indices1[depth + 1] = x
                recursive_sum3(matrix, result, indices1, indices2, depth + 1, max_depth, ranges)
            else
                indices2[depth - length(indices1) + 1] = x
                recursive_sum3(matrix, result, indices1, indices2, depth + 1, max_depth, ranges)
            end
        end
    end
end

function sum_overR3(matrix, ranges)
    half_length = div(length(ranges), 2)
    result_size = prod(ranges[1:half_length])
    result = zeros(eltype(matrix), result_size, result_size)
    indices1 = zeros(Int, half_length)
    indices2 = zeros(Int, half_length)
    recursive_sum3(matrix, result, indices1, indices2, 0, length(ranges), ranges)
    return result / result_size
end


function conditional_distributiona(joint_prob::Array, dist_index::Int)
    # Sum over the distribution index to get the marginal distribution
    marginal_prob = sum(joint_prob, dims=dist_index)

    # Ensure no zero probabilities to avoid division by zero
    marginal_prob[marginal_prob.==0] .= eps()

    # Expand the marginal to the shape of the joint distribution for broadcasting
    marginal_prob_expanded = repeat(marginal_prob, inner=(1, 1, 1, size(joint_prob, dist_index)))

    # Divide the joint distribution by the marginal to get the conditional distribution
    cond_prob = joint_prob ./ marginal_prob_expanded

    return cond_prob
end

function conditional_distribution(joint_prob::Array, dist_index::Int, cond_indices::Vector{Int})
    # Sum over the dimensions that are not part of the conditional variables
    dims_to_sum = setdiff(1:ndims(joint_prob), vcat(dist_index, cond_indices))
    marginal_prob = sum(joint_prob, dims=Tuple(dims_to_sum))  # Pass dims as a tuple

    # Ensure no zero probabilities to avoid division by zero
    marginal_prob[marginal_prob.==0] .= eps()

    # Expand the marginal to the shape of the joint distribution for broadcasting
    repeat_dims = ones(Int, ndims(joint_prob))
    for i in cond_indices
        repeat_dims[i] = size(joint_prob, i)
    end
    repeat_dims[dist_index] = size(joint_prob, dist_index)

    marginal_prob_expanded = repeat(marginal_prob, repeat_dims...)

    # Compute the conditional distribution by dividing the joint distribution by the expanded marginal
    cond_prob = joint_prob ./ marginal_prob_expanded

    return cond_prob
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