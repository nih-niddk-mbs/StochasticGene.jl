# This file is part of StochasticGene.jl  

# test.jl


simDT_convert(v) = [[s[1], s[3]] for s in v]


### functions used in runtest

function test_DT(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)])
    reporter, components = make_reporter_components(transitions, G, R, S, insertstep, onstates, dttype, "")
    likelihoodarray(r, components, bins, reporter, dttype)
end


function test_sim(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nhist=150, nalleles=2, onstates=[Int[], [2, 3]], bins=[collect(5/3:5/3:200), collect(0.1:0.1:20)], total=10000000, tol=1e-6)
    simulator(r, transitions, G, R, S, insertstep, nhist=nhist, nalleles=nalleles, onstates=onstates, bins=bins, totalsteps=total, tol=tol)
end

function test_cm(r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)])
    reporter, tcomponents = make_reporter_components(transitions, G, R, S, insertstep, onstates, dttype, "")
    components = MTDComponents(MComponents(transitions, G, R, nRNA, r[num_rates(transitions, R, S, insertstep)], ""), tcomponents)
    likelihoodarray(r, components, bins, reporter, dttype, nalleles, nRNA)
end

function test_compare(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=10000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_cm(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, dttype, bins)
    hs = StochasticGene.normalize_histogram.(hs)
    return h, make_array(hs)
end


function test_fit_rnadwelltime2(; rtarget=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="", maxtime=360.0)
    hs = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 10000000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, hs[1], bins, hs[2:end], dttype)
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R))
    # model = load_model(data, rtarget, rtarget, fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, nothing, [], 1, tuple(), tuple(), nothing)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    return h, StochasticGene.datapdf(data), fits, stats, measures, data, model, options
end


function test_CDT(r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.045, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [2, 0], [0, 3], 1))
    r, coupling_strength, _ = prepare_rates(r, coupling[3], transitions, G, R, S, insertstep, [0, 0])
    reporter, components = make_reporter_components(transitions, G, R, S, insertstep, onstates, dttype, "", coupling)
    likelihoodarray(r, coupling_strength, components::TCoupledComponents{Vector{TDCoupledUnitComponents}}, bins, reporter, dttype)
end

function test_sim_coupled(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 0.0, 0.2, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 0.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [2], [2]], [Int[], Int[], [2], [2]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [2, 0], [0, 2], 1), total=10000000, tol=1e-6, verbose=false)
    simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol, verbose=verbose)
end

function test_compare_coupling(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.45, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, 2.9], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [3, 0], [0, 4], 1), total=10000000, tol=1e-6)
    hs = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol)
    h = test_CDT(r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling)
    for i in eachindex(hs)
        hs[i] = StochasticGene.normalize_histogram.(hs[i])
    end
    return make_array(vcat(h...)), make_array(vcat(hs...))
end

function test_fit_simrna(; rtarget=[0.33, 0.19, 20.5, 1.0], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(), rinit=[0.1, 0.1, 0.1, 1.0], totalsteps=100000)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles)[1]
    data = RNAData("", "", nRNA, h)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0), fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(1000000, 0, 0, 20.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_rna(; gene="CENPL", G=2, nalleles=2, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), rinit=[0.01, 0.1, 1.0, 0.01006327034802035], datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".")
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, (), 1.0)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, 1.0, [], 1.0), fittedparam, fixedeffects, transitions, 2, 0, 0, 1, "", nalleles, 10.0, Int[], rinit[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(100000, 0, 0, 60.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

function test_fit_rnaonoff(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=fill(0.01, num_rates(transitions, R, S, insertstep)), nsamples=100000, nhist=20, nalleles=2, onstates=Int[], bins=collect(1:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="")
    # OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    hs = simulator(rtarget, transitions, G, R, S, insertstep, nalleles=nalleles, nhist=nhist, onstates=onstates, bins=bins, totalsteps=10000000)
    hRNA = div.(hs[1], 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, hs[2], hs[3])
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

function test_fit_rnadwelltime(; rtarget=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="", maxtime=360.0)
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 10000000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype)
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R))
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, rtarget[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=10000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70])
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], noisepriors, mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_trace_hierarchical(; G=2, R=1, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 1.0, 50, 15, 200, 70], rinit=[], nsamples=100000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=collect(1:num_rates(transitions, R, S, insertstep)-1), propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], hierarchical=(2, collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+length(noisepriors)), tuple()), method=(Tsit5(), true), maxtime=180.0)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, hierarchical[1], mean_elongationtime(rtarget, transitions, R))
    isempty(rinit) && (rinit = StochasticGene.set_rinit(rm, transitions, R, S, insertstep, noisepriors, length(data.trace[1])))
    model = load_model(data, rinit, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    nrates = num_rates(transitions, R, S, insertstep)
    h1 = StochasticGene.get_rates(fits.parml, model)[1:nrates+4]
    h2 = rtarget[1:nrates+4]
    return h1, h2
end

function test_fit_tracejoint(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(2, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.0, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([100, 50, 200, 100], [100, 50, 200, 100]), maxtime=300.0)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], 0.0, 1))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    fittedparam = StochasticGene.set_fittedparam(fittedparam, data.label, transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rm, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), coupling, nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function sim_grid(; r=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], p=0.2, Ngrid=4, transitions=([1, 2], [2, 1]), G=2, R=1, S=1, insertstep=1, totaltime=1000.0, interval=1.0, ntrials=10)
    Nstate = num_rates(transitions, R, S, insertstep)
    a_grid = StochasticGene.make_a_grid(p, Ngrid)
    # traces = simulator(r, transitions, G, R, S, insertstep, traceinterval=interval, nhist=0, totaltime=totaltime, reporterfn=sum, a_grid=StochasticGene.make_a_grid(1.0, 4))[1]
    StochasticGene.simulate_trace_vector(r, transitions, G, R, S, insertstep, interval, totaltime, ntrials, a_grid=a_grid)
end

function test_fit_trace_grid(; Ngrid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [50, 15, 200, 20]; 0.1], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 6, 7, 8, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=Ngrid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
    model = load_model(data, [rinit; 0.1], [rinit; 0.1], fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), tuple(), Ngrid)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_trace_gridhierarchical(; Ngrid=4, G=2, R=1, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 1.0, 50, 15, 200, 70, .02], rinit=[], nsamples=100000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1);num_rates(transitions, R, S, insertstep)+5], propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], hierarchical=(2, collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+length(noisepriors)), tuple()), method=(Tsit5(), true), maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=Ngrid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (traces, [], 0.0, 1))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, hierarchical[1], mean_elongationtime(rtarget, transitions, R))
    isempty(rinit) && (rinit = StochasticGene.set_rinit(rm, transitions, R, S, insertstep, noisepriors, length(data.trace[1])))
    model = load_model(data, rinit, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    # nrates = num_rates(transitions, R, S, insertstep)
    # h1 = StochasticGene.get_rates(fits.parml, model)[1:nrates+4]
    # h2 = rtarget[1:nrates+4]
    # return h1, h2
end
# function test_fit_trace_gridhierarchical(; Ngrid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [50, 15, 200, 20]; 0.1], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 6, 7, 8, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
#     traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=Ngrid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
#     data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
#     model = load_model(data, [rinit; 0.1], [rinit; 0.1], fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), tuple(), Ngrid)
#     options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
#     fits, stats, measures = run_mh(data, model, options)
# end

### end of functions used in runtest

function test_grid_ll(; Ngrid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70, 0.2], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+5)], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    # traces = sim_grid(r=rtarget, p=rtarget[end], Ngrid=Ngrid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=Ngrid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
    model = load_model(data, [rinit; 0.3], [rinit; 0.3], fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), tuple(), Ngrid)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    r, noiseparams, pgrid = prepare_rates(get_param(model), model)
    Nstate = model.components.nT
    a_grid = make_a_grid(pgrid[1], model.Ngrid)
    a, p0 = make_ap(r, data.interval, model.components)
    d = model.reporter.probfn(noiseparams, model.reporter.per_state, model.components.nT, model.Ngrid)
    logpredictions = Array{Float64}(undef, 0)

    t = data.trace[1][1]
    T = size(t, 2)
    b = set_b_grid(t, d, Nstate, Ngrid)
    _, C = forward_grid(a, a_grid, b, p0, Nstate, Ngrid, T)
    C, a, a_grid, b, p0, Nstate, Ngrid, T, t

end

function test_hierarchical(; G=2, R=1, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 1.0, 50, 15, 200, 70], rinit=[], nsamples=100000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=collect(1:num_rates(transitions, R, S, insertstep)-1), propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70], hierarchical=(2, collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+length(noisepriors)), tuple()), method=(1, true), maxtime=120.0)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, hierarchical[1], mean_elongationtime(rtarget, transitions, R))
    isempty(rinit) && (rinit = StochasticGene.set_rinit(rm, transitions, R, S, insertstep, noisepriors, length(data.trace[1])))
    model = load_model(data, rinit, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 120.0, 1.0, 1.0)
    return data, model, options
    # fits, stats, measures = run_mh(data, model, options)
    # nrates = num_rates(transitions, R, S, insertstep)
    # h1 = StochasticGene.get_rates(fits.parml, model)[1:nrates]
    # h2 = rtarget[1:nrates]
    # return h1, h2
end

function test_fit_coupled(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(2, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.0, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([100, 50, 200, 100], [100, 50, 200, 100]), maxtime=120.0)
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.5, 0.5], 900))
    rm = StochasticGene.prior_ratemean(transitions, R, S, insertstep, 1.0, noisepriors, [5.0, 5.0], coupling)
    fittedparam = StochasticGene.set_fittedparam(fittedparam, data.label, transitions, R, S, insertstep, noisepriors, coupling, nothing)
    model = load_model(data, rtarget, rm, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), coupling, nothing)
    return data, model
    # r, couplingStrength, noiseparams = prepare_rates(get_param(model), model)
end


function loglikelihood_test(data::TraceData, model::GRSMcoupledmodel)
    r, couplingStrength, noiseparams = prepare_rates(get_param(model), model)
    ll_hmm_coupled_test(r, couplingStrength, noiseparams, model.components, model.reporter, data.interval, data.trace)
end

function test_Cparts(r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.045, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, -0.5], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [2, 0], [0, 3], 1))
    reporter, components = make_reporter_components(transitions, G, R, S, insertstep, onstates, dttype, "", coupling)
    r, coupling_strength, _ = prepare_rates(r, coupling[3], transitions, G, R, S, insertstep, [0, 0])
    sojourn, nonzeros = reporter
    sources = components.sources
    model = components.model
    TDdims = get_TDdims(components)
    T, TD, Gm, Gt, Gs, IG, IR, IT = make_matvec_C(components, r)
    TCD = Vector{Vector{SparseMatrixCSC}}(undef, length(Gm))
    TCDnew = Vector{Vector{SparseMatrixCSC}}(undef, length(Gm))
    p = nothing
    pnew = nothing
    for α in eachindex(components.modelcomponents)
        dt = occursin.("G", dttype[α])
        p = steady_state_dist(α, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, dt)
        TCD[α] = make_mat_TCD(α, TD[α], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn[α])
    end
    return TCD, p
end

function test_trace_background(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=10000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0, noisepriors=[50, 15, 200, 70])
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.5, 1))
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], noisepriors, mean_elongationtime(rtarget, transitions, R)), fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    ll, lv = ll_hmm(rtarget, model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    # fits, stats, measures = run_mh(data, model, options)
    # StochasticGene.get_rates(fits.parml, model), rtarget
    return ll, trace
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

function test_mat_T(r, transitions, G, R, S, insertstep)
    components1 = StochasticGene.make_components_T(transitions, G, R, S, insertstep, "")
    components2 = StochasticGene.make_components_TGRS(transitions, G, R, S, insertstep, "")
    T1 = StochasticGene.make_mat_T(components1, r)
    T2 = StochasticGene.make_mat_T(components2, r)
    return T1, T2
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

function test_mat_TD(r, transitions, G, R, S, insertstep, onstates, dttype)
    _, components = make_reporter_components(transitions, G, R, S, insertstep, onstates, dttype, "")
    make_mat_TD(components, r)
end

function test_C(r=[0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1], transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), G=(2, 2), R=(1, 1), S=(0, 0), insertstep=(1, 1), bins=[[collect(1:20), collect(1:20)], [collect(1:20), collect(1:20)]], coupling=((1, 2), (Int[], [1]), [1, 0], [0, 1], 1), onstates=[[[2], Int[]], [[2], Int[]]], dttype=[["ONG", "ON"], ["ONG", "ON"]])
    reporter, components = make_reporter_components(transitions, G, R, S, insertstep, onstates, dttype, "", coupling)
    r, coupling_strength, _ = prepare_rates(r, coupling[3], transitions, G, R, S, insertstep, [0, 0])
    sojourn, nonzeros = reporter
    sources = components.sources
    model = components.model
    TDdims = get_TDdims(components)
    T, TD, Gm, Gt, Gs, IG, IR, IT = make_matvec_C(components, r)
    hists = Vector{Vector}[]
    TCD = Vector{Vector{SparseMatrixCSC}}(undef, length(Gm))
    for α in eachindex(components.modelcomponents)
        dt = occursin.("G", dttype[α])
        p = steady_state_dist(α, TDdims[α], T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, dt)
        h = Vector[]
        TCD[α] = make_mat_TCD(α, TD[α], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
        for i in eachindex(sojourn[α])
            if dt[i]
                push!(h, init_S(sojourn[α][i], p.GC, p.pssG))
            else
                push!(h, init_S(sojourn[α][i], p.TC, p.pss, nonzeros[α][i]))
            end
        end
        push!(hists, h)
    end
    return hists, reporter, components, r, coupling_strength, TCD
end

function test_mat_Tc(coupling, r, coupling_strength, transitions, G, R, S, insertstep)
    components = make_components_TRGCoupled(coupling, transitions, G, R, S, insertstep, "")
    return make_mat_TC(components, r, coupling_strength)
end

function test_mat_Tc2(coupling, r, coupling_strength, transitions, G, R, S, insertstep)
    T = []
    for i in eachindex(G)
        c = make_components_TRG(transitions[i], G[i], R[i], S[i], insertstep[i], "")
        push!(T, make_mat_TRG(c, r[i]))
    end
    components = make_components_TCoupled(coupling, transitions, G, R, S, insertstep, "")
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
    components = StochasticGene.make_components_TCoupled(coupling, transitions, G, R, S, insertstep)
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



