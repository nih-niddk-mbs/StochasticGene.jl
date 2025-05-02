# This file is part of StochasticGene.jl  

# test.jl


"""
    test_steadystatemodel(model::AbstractGMmodel, nhist)

Compares the chemical master solution to a Gillespie simulation for the steady-state mRNA distribution.

# Arguments
- `model`: An instance of `AbstractGMmodel`.
- `nhist`: The number of histogram bins for the simulation.

# Description
This function compares the steady-state mRNA distribution obtained from the chemical master equation solution to that obtained from a Gillespie simulation. It uses the rates and number of states from the provided model.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: The steady-state mRNA distributions from the chemical master solution and the Gillespie simulation.
"""
function test_steadystatemodel(model::AbstractGMmodel, nhist)
    G = model.G
    r = model.rates
    g1 = steady_state(r[1:2*G], G - 1, nhist, model.nalleles)
    g2 = simulatorGM(r[1:2*G], G - 1, nhist, model.nalleles)
    return g1, g2
end

"""
    test_model(data::RNAOnOffData, model::AbstractGRSMmodel)

Simulates the RNA on-off model with splicing and compares it to the provided data.

# Arguments
- `data`: An instance of `RNAOnOffData` containing the observed data.
- `model`: An instance of `AbstractGRSMmodel` representing the model to be tested.

# Description
This function simulates the RNA on-off model with splicing using the provided model parameters and compares the simulated results to the observed data. It uses the `telegraphsplice0` function to perform the simulation.

# Returns
- `Nothing`: The function performs the simulation and comparison but does not return a value.
"""
function test_model(data::RNAOnOffData, model::AbstractGRSMmodel)
    telegraphsplice0(data.bins, data.nRNA, model.G - 1, model.R, model.rates, 1000000000, 1e-5, model.nalleles)
end

simDT_convert(v) = [[s[1], s[3]] for s in v]

function test_DT(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)])
    reporter, components = make_reporter_components_DT(transitions, G, R, S, insertstep, "", onstates, dttype)
    predictedarray(r, components, bins, reporter, dttype)
end

function test_sim(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nhist=150, nalleles=2, onstates=[Int[], [2, 3]], bins=[collect(5/3:5/3:200), collect(0.1:0.1:20)], total=10000000, tol=1e-6)
    simulator(r, transitions, G, R, S, insertstep, nhist=nhist, nalleles=nalleles, onstates=onstates, bins=bins, totalsteps=total, tol=tol)
end

function test_cm(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)])
    reporter, tcomponents = make_reporter_components_DT(transitions, G, R, S, insertstep, "", onstates, dttype)
    mcomponents = MComponents(transitions, G, R, nRNA, r[num_rates(transitions, R, S, insertstep)], "", 1)
    components = MTComponents{typeof(mcomponents),typeof(tcomponents)}(mcomponents, tcomponents)
    predictedarray(r, components, bins, reporter, dttype, nalleles, nRNA)
end

function test_CDT(r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.045, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [2, 0], [0, 3], 1))
    reporter, components = make_reporter_components_DT(transitions, G, R, S, insertstep, "", onstates, dttype, coupling)
    couplingindices = coupling_indices(transitions, R, S, insertstep, [0, 0], coupling, nothing)
    rates, couplingStrength = prepare_rates_coupled(r, [num_rates(transitions[1], R[1], S[1], insertstep[1]), num_rates(transitions[2], R[2], S[2], insertstep[2])], couplingindices)
    predictedarray(rates, couplingStrength, components, bins, reporter, dttype)
end


function test_sim_coupled(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 0.0, 0.2, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 0.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [2], [2]], [Int[], Int[], [2], [2]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [2, 0], [0, 2], 1), total=10000000, tol=1e-6, verbose=false)
    simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol, verbose=verbose)
end

function sim_grid(; r=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], p=0.2, Ngrid=4, transitions=([1, 2], [2, 1]), G=2, R=1, S=1, insertstep=1, totaltime=1000.0, interval=1.0, ntrials=10)
    Nstate = num_rates(transitions, R, S, insertstep)
    a_grid = StochasticGene.make_a_grid(p, Ngrid)
    # traces = simulator(r, transitions, G, R, S, insertstep, traceinterval=interval, nhist=0, totaltime=totaltime, reporterfn=sum, a_grid=StochasticGene.make_a_grid(1.0, 4))[1]
    StochasticGene.simulate_trace_vector(r, transitions, G, R, S, insertstep, interval, totaltime, ntrials, a_grid=a_grid)
end

### functions used in runtest

function test_compare(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=10000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_cm(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, dttype, bins)
    hs = StochasticGene.normalize_histogram.(hs)
    return h, make_array(hs)
end

function test_compare_coupling(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.45, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, 2.9], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [3, 0], [0, 4], 1), total=10000000, tol=1e-6)
    hs = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol)
    h = test_CDT(r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling)
    for i in eachindex(hs)
        hs[i] = StochasticGene.normalize_histogram.(hs[i])
    end
    return make_array(vcat(h...)), make_array(vcat(hs...))
end

function test_fit_simrna(; rtarget=[0.33, 0.19, 20.5, 1.0], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(), rinit=[0.1, 0.1, 0.1, 1.0], totalsteps=100000, nchains=1)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles)[1]
    data = RNAData("", "", nRNA, h)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0), fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(1000000, 0, 0, 20.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_rna(; gene="CENPL", G=2, nalleles=2, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), rinit=[0.01, 0.1, 1.0, 0.01006327034802035], datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".", nchains=1)
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, (), 1.0)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, 1.0, [], 1.0), fittedparam, fixedeffects, transitions, 2, 0, 0, 1, "", nalleles, 10.0, Int[], rinit[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(100000, 0, 0, 60.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

function test_fit_rnaonoff(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=fill(0.01, num_rates(transitions, R, S, insertstep)), nsamples=100000, nhist=20, nalleles=2, onstates=Int[], bins=collect(1:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="", nchains=1)
    # OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    hs = simulator(rtarget, transitions, G, R, S, insertstep, nalleles=nalleles, nhist=nhist, onstates=onstates, bins=bins, totalsteps=10000000)
    hRNA = div.(hs[1], 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, hs[2], hs[3])
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

function test_fit_rnadwelltime(; rtarget=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="", maxtime=360.0, nchains=1)
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 10000000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype)
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R))
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=10000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], nchains=1)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], noisepriors, mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, Tsit5(), tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_trace_hierarchical(; G=2, R=1, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 1.0, 50, 15, 200, 70], rinit=[], nsamples=100000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4], propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], hierarchical=(2, [6], tuple()), method=(Tsit5(), true), maxtime=180.0, nchains=1)
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, interval, totaltime, ntrials, hierarchical=(6, rh))
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, tuple(), nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    fittedparam = set_fittedparam(fittedparam, "trace", transitions, R, S, insertstep, noisepriors, tuple(), tuple())
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    # nrates = num_rates(transitions, R, S, insertstep)
    # h1 = StochasticGene.get_rates(fits.parml, model)[1:nrates+4]
    # h2 = rtarget[1:nrates+4]
    h1 = stats.medparam
    h2 = [rtarget[fittedparam]; [50.0; 10 / 50.0]; rh]
    return h1, h2
end

function test_fit_tracejoint(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 1.0, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 1.0, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, method=Tsit5())
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1))
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), coupling, nothing)
    rinit = set_rinit(rinit, priormean)
    fittedparam = set_fittedparam(fittedparam, "trace", transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    rfit = StochasticGene.get_rates(fits.parml, model)
    return rfit[1:length(rtarget)], rtarget
end

function test_fit_tracejoint_hierarchical(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, decayrate=1.0)
    # trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials, hierarchical=(6, rh))
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1))
    priormean = set_priormean([], transitions, R, S, insertstep, decayrate, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, coupling, nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), coupling, nothing)
    fittedparam = set_fittedparam(fittedparam, "tracejoint", transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    rfit = StochasticGene.get_rates(fits.parml, model)
    return rfit[1:length(rtarget)], rtarget, fits, stats, measures, model, data, options
end

function test_fit_trace_grid(; grid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [50, 15, 200, 20]; 0.1], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 6, 7, 8, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=grid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
    priormean = set_priormean([], transitions, R, S, insertstep, rtarget[num_rates(transitions, R, S, insertstep)], noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), grid)
    rinit = set_rinit(rinit, priormean)
    fittedparam = set_fittedparam(fittedparam, "tracegrid", transitions, R, S, insertstep, noisepriors, tuple(), grid)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian_grid, noisepriors, Tsit5(), tuple(), tuple(), grid)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end


function test_fit_trace_grid_hierarchical(; grid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[], nsamples=5000, onstates=Int[], hierarchical=(2, [7], tuple()), totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=grid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
    priormean = set_priormean([], transitions, R, S, insertstep, rtarget[num_rates(transitions, R, S, insertstep)], noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, tuple(), grid)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), grid)
    fittedparam = set_fittedparam(fittedparam, "tracegrid", transitions, R, S, insertstep, noisepriors, tuple(), grid)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian_grid, noisepriors, (Tsit5(), true), hierarchical, tuple(), grid)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end

### end of functions used in runtest



test_fit0(; nchains=1, maxtime=60.0, propcv=0.01, zeromedian=false) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 1.0, [0.0, 175.0]), "test", "test", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 100.0, 1200.0, 175.0, 1.0, 1.0, 1.0, 0.1, 0.1, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 290.0, 1200.0, 175.0], (), "ml", 0.004, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fit1(; nchains=1, maxtime=60.0, propcv=0.01, zeromedian=false) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 1.0, [290, 175.0]), "test", "test", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 290.0, 100.0, 1200.0, 175.0], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [290.0, 290.0, 1200.0, 175.0], (), "ml", 0.004, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(), zeromedian)

test_fitt(; nchains=1, maxtime=60.0, propcv=0.01, zeromedian=false) = fit(nchains, "trace", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 1., [25000.0, 10000.0]), "genetogther", "genetogther", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6,8,9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 20000.0, 15000.0, 80000.0, 4000.0], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [20000.0, 15000.0, 80000.0, 4000.0], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5(),zeromedian)

test_fit_h(; nchains=1, maxtime=60.0, propcv=0.01) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 1.0, [290.0, 140.0]), "test", "test", "trace-HBEC-nstate-h_gene", "trace-HBEC-nstate-h_gene", [1, 2, 3, 4, 5, 6], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 290.0, 290.0, 1200.0, 175.0, 1.0, 1.0, 1.0, 0.1, 0.1, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, Int64[], 1.0, "", prob_Gaussian, [290.0, 290.0, 1200.0, 175.0], (2, [8], ()), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (Tsit5(), true))

test_fit_coupled(; nchains=1, maxtime=60.0, propcv=0.01) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [[25000, 10000], [25000, 10000]]), "test", "test", "tracejoint-HBEC-nstate_enhancer-geneR5", "tracejoint-HBEC-nstate_enhancer-geneR5", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (3, 3), (1, 0), (1, 1), ((1, 2), ((), (1,)), ([4, 5, 6], 0), (0, 5), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([25000.0, 15000.0, 80000.0, 10000.0], [25000.0, 15000.0, 80000.0, 10000.0]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5())

test_fit_coupled_h(; nchains=1, maxtime=60.0, propcv=0.01) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [[25000, 10000], [25000, 10000]]), "test", "test", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (3, 3), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([25000.0, 15000.0, 80000.0, 10000.0], [25000.0, 15000.0, 80000.0, 10000.0]), (2, [11, 10], ()), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, true, false, false, (Tsit5(), true))

test_fit_coupled_h2(; nchains=1, maxtime=60.0, propcv=0.01) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [[25000, 0], [25000, 0]]), "test", "test", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([25000.0, 15000.0, 80000.0, 10000.0], [25000.0, 15000.0, 80000.0, 10000.0]), (2, [], ()), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (Tsit5(), true))

test_rnacount(;nchains=1, gene="RPLP1", datapath="data/U3AS4/WT-UTR", maxtime=60.0, propcv=0.01) = fit(nchains=nchains, datatype="rnacount", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath=datapath, cell="U3A", gene=gene, datacond="WT-UTR", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv);

test_rna(; nchains=1, gene="RPLP1", datapath="data/U3AS4/histograms", maxtime=60.0, propcv=0.01) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath=datapath, cell="U3A", gene=gene, datacond="WT-UTR", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv);

test_fit_G(; nchains=1, maxtime=60.0, propcv=0.01) = fit(nchains, "trace", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.86, [25000.0, 0.0]), "genetogther", "genetogther", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3], (), ([1, 2], [2, 1]), 2, 0, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 1.0, 20000.0, 15000.0, 80000.0, 4000.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 1.0, 0.5, 0.5, 1., .1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, [2], 1.0, "", prob_Gaussian, [20000.0, 15000.0, 80000.0, 4000.0], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5())

fit_myc3(;nchains=1, maxtime=60.0, propcv=0.05) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=0, S=0, insertstep=1, datapath="data/FISH", cell="HBEC", gene="MYC", datacond="gene", infolder="test", resultfolder="test", fittedparam=[1, 2, 3,4,5], maxtime=maxtime, ratetype="ml", propcv=propcv, burst=true);

fit_myc2(;nchains=1, maxtime=60.0, propcv=0.05) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath="data/FISH", cell="HBEC", gene="MYC", datacond="gene", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv, burst=true);

# @time fit(16, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, [290.0, 140.0]), "control-2025-03-15", "control-2025-03-15", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene",[1, 2,3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", 43000.0, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 290.0, 140.0, 1200.0, 175.0], [2.0, 2.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [290.0, 140.0, 1200.0, 175.0], (), "ml", 0.02, 2000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, Tsit5())

function test_tracejoint_h(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, decayrate=1.0)
    # trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials, hierarchical=(6, rh))
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1))
    priormean = set_priormean([], transitions, R, S, insertstep, decayrate, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, coupling, nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), coupling, nothing)
    fittedparam = set_fittedparam(fittedparam, "tracejoint", transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    return data, model, options
end


function test_tracejoint(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([100, 50, 200, 100], [100, 50, 200, 100]), maxtime=300.0, method=Tsit5())
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials; verbose=false)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], 0.0, 1))
    # coupling = ((1, 2), (tuple(), tuple(1)), (["G2"], String[]), (0, 1), 1)
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    rinit = set_rinit(rinit, rm)
    fittedparam = StochasticGene.set_fittedparam(fittedparam, data.label, transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rinit, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    return data, model, options
end

function test_predicted_states(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.1, 0.3, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 20, 200, 20, 0.1, 0.2, 0.4, 0.2, 0.1, 50, 20, 200, 20, -0.5], rinit=Float64[], nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 500, 20], [50, 30, 500, 20]), maxtime=300.0, decayrate=1.0, totalsteps=20, verbose=false, hierarchical=(2, [8, 17], tuple()), method=(Tsit5(), true))
    rh = 50.0 .+ 10 * randn(ntrials)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials, hierarchical=(6, rh))
    tracesingle = Vector{Vector{Vector{Float64}}}(undef, length(R))
    df = Vector{DataFrame}(undef, length(R))
    k = 0
    for i in eachindex(R)
        tracesingle = [t[:, i] for t in trace]
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + length(noisepriors[i])
        runit = rtarget[k+1:k+n]
        k += n
        df[i] = StochasticGene.make_traces_dataframe(tracesingle, interval, [runit; ones(length(runit)); repeat(runit, ntrials)], transitions[i], G[i], R[i], S[i], insertstep[i], prob_Gaussian, 4, "", true, true, tuple())
    end
    dfjoint = StochasticGene.make_traces_dataframe(trace, interval, [rtarget; ones(length(rtarget)); repeat(rtarget, ntrials)], transitions, G, R, S, insertstep, prob_Gaussian, 4, "", true, true, coupling)
    return df, dfjoint
end

function test_load(; traceinfo=(1.0, 1.0, -1, 1.0), datapath="data/inhibition/control/", label="", gene="MYC", datacond="gene", datatype=:trace, zscoretrace=true)
    data = load_data_trace(datapath, label, gene, datacond, traceinfo, datatype, 3, zscoretrace)
    return data
end

marg(p, dims) = dropdims(sum(p, dims=dims), dims=dims)


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
    components = StochasticGene.make_components_TCoupled(coupling, transitions, G, R, S, insertstep)
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
        idx1 = i + 2 * j + 4 * k + 8 * l + 1
        idx2 = i1 + 2 * j1 + 4 * k1 + 8 * l1 + 1
        result[i+2*k+1, i1+2*k1+1] += matrix[idx1, idx2]
    else
        for x in 0:1
            if depth < 4
                new_indices1 = copy(indices1)
                new_indices1[depth+1] = x
                recursive_sum(matrix, result, new_indices1, indices2, depth + 1, max_depth)
            else
                new_indices2 = copy(indices2)
                new_indices2[depth-4+1] = x
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
    result / 4
end

function recursive_sum3(matrix, result, indices1, indices2, depth, max_depth, ranges)
    if depth == max_depth
        idx1 = sum(indices1[i] * prod(ranges[1:i-1]) for i in 1:length(indices1)) + 1
        idx2 = sum(indices2[i] * prod(ranges[1:i-1]) for i in 1:length(indices2)) + 1
        result[indices1[1]+2*indices1[3]+1, indices2[1]+2*indices2[3]+1] += matrix[idx1, idx2]
    else
        for x in 0:ranges[depth+1]-1
            if depth < length(indices1)
                indices1[depth+1] = x
                recursive_sum3(matrix, result, indices1, indices2, depth + 1, max_depth, ranges)
            else
                indices2[depth-length(indices1)+1] = x
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



function split_matrix(mat::Matrix{Float64}, idx1::Vector{Int}, idx2::Vector{Int})
    return hcat([mat[idx1, j] for j in 1:size(mat, 2)]...,
        [mat[idx2, j] for j in 1:size(mat, 2)]...)
end

function split_matrix2(mat::Matrix{Float64})
    N, M = size(mat)
    half = div(N, 2)  # Divide N into two equal parts (assuming even N)

    return [mat[1:half, j] for j in 1:M],  # First row: first half
    [mat[half+1:N, j] for j in 1:M] # Second row: second half
end


function convert_to_m_dim(mat::Matrix{Float64})
    row1, row2 = split_matrix2(mat)
    return [[row1[j], row2[j]] for j in 1:size(mat, 2)]  # Each element is a vector of two vectors
end
