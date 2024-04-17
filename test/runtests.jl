
using StochasticGene
using Test

function test_sim(r, transitions, G, R, S, insertstep, nhist, nalleles, onstates, bins, total, tol)
    simulator(r, transitions, G, R, S, insertstep, nhist=nhist, nalleles=nalleles, onstates=onstates, bins=bins, totalsteps=total, tol=tol)
end

function test_compare(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=1000000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    hs = StochasticGene.normalize_histogram.(hs)
    return make_array(h), make_array(hs)
end

function test_chem(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, bins, dttype)
    for i in eachindex(onstates)
        if isempty(onstates[i])
            onstates[i] = on_states(G, R, S, insertstep)
        end
        onstates[i] = Int64.(onstates[i])
    end
    components = make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nRNA, r[num_rates(transitions, R, S, insertstep)], "")
    likelihoodarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)
end

function test_fit_simrna(; rtarget = [.33, .19, 20.5, 1.], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(),rinit=[.1,.1,.1,1.],totalsteps=100000)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps,nalleles=nalleles)
    data = RNAData("", "", nRNA, h)
    model = load_model(data, rinit, Float64[], fittedparam, fixedeffects, transitions, G, 0, 0, 0, nalleles, 10.0, Int[], rinit[end], .02, "", prob_GaussianMixture, [], tuple())
    options = StochasticGene.MHOptions(1000000, 0, 0, 20.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_rna(; gene="CENPL", G=2, nalleles=2, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), rinit=[0.01, 0.1, 1.0, 0.01006327034802035], datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".")
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, 1.0, 1.0, [1,2])
    model = load_model(data, rinit, Float64[], fittedparam, fixedeffects, transitions, 2, 0, 0, 1, nalleles, 10.0, Int[], rinit[end], propcv, "", prob_GaussianMixture, [], tuple())
    options = StochasticGene.MHOptions(100000, 0, 0, 60.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

function test_fit_rnaonoff(; G=2, R=1, S=1, transitions=([1, 2], [2, 1]), insertstep=1, rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01], rinit=fill(0.01, num_rates(transitions, R, S, insertstep)), nsamples=100000, nhist=20, nalleles=2, onstates=Int[], bins=collect(1:1.0:200.0), fittedparam=collect(1:length(rtarget)-1), propcv=0.05, priorcv=10.0, splicetype="")
    # OFF, ON, mhist = test_sim(rtarget, transitions, G, R, S, nhist, nalleles, onstates, bins)
    h = simulator(rtarget, transitions, G, R, S, insertstep, nalleles=nalleles, nhist=nhist, onstates=onstates, bins=bins, totalsteps=3000)
    hRNA = div.(h[1], 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, h[2], h[3])
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, [], tuple())
    options = StochasticGene.MHOptions(nsamples, 0, 0, 120.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

function test_fit_rnadwelltime(; rtarget=[0.038, 2.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:10), collect(0.1:0.1:10)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="")
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 30000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype)
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [])
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, [], tuple())
    options = StochasticGene.MHOptions(nsamples, 1000, 0, 200.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data)
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4)], propcv=0.01, cv=100.0, interval=1.0,noisepriors=[50, 15, 200, 70])
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.))
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, tuple())
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

function test_fit_trace_hierarchical(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [20, 5, 100, 10]], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=collect(1:num_rates(transitions, R, S, insertstep)-1), propcv=0.01, cv=100.0, interval=1.0,noisepriors=[50, 15, 200, 70],hierarchical=(2,collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+4),tuple()))
    trace = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.))
    println(hiearchical)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_Gaussian, noisepriors, hierarchical)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml, model), rtarget
end


@testset "StochasticGene" begin

    h1, h2 = test_compare()
    @test isapprox(h1, h2, rtol=0.2)

    h1, h2 = test_fit_simrna()
    @test isapprox(h1, h2, rtol=0.05)

    h1, h2 = test_fit_rna()
    @test isapprox(h1, h2, rtol=0.1)

    h1, h2 = test_fit_rnaonoff()
    @test isapprox(h1, h2, rtol=0.3)

    h1, h2 = test_fit_rnadwelltime()
    @test isapprox(h1, h2, rtol=0.3)

    h1, h2 = test_fit_trace()
    @test isapprox(h1, h2, rtol=0.05)

end