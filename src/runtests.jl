
using StochasticGene
using Test



"""
    test(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.01:0.1:10), collect(0.01:0.1:10)], total=1000000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])


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

function test_fit_rna(; gene="CENPL", cell="HCT116", fish=false, G=2, nalleles=2, nsamples=100000, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), ejectprior=0.05, r=[0.01, 0.1, 1.0, 0.01006327034802035], decayrate=0.01006327034802035, datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".")
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, 1.0, 1.0, 1.0)
    model = load_model(data, r, Float64[], fittedparam, fixedeffects, transitions, G, 0, 0, 1, nalleles, 10.0, Int[], r[end], propcv, "", prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(100000, 0, 0, 60.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
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
    h, StochasticGene.datapdf(data), fits, stats, measures, data, model, options
end

function test_fit_rnadwelltime(; rtarget=[0.038, 2.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, nRNA=150, nsamples=100000, nalleles=2, onstates=[Int[], Int[], [2, 3], [2, 3]], bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.01:0.1:10), collect(0.01:0.1:10)], dttype=["ON", "OFF", "ONG", "OFFG"], fittedparam=collect(1:length(rtarget)-1), propcv=0.01, priorcv=10.0, splicetype="")
    h = test_sim(rtarget, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], 30000, 1e-6)
    data = RNADwellTimeData("test", "test", nRNA, h[1], bins, h[2:end], dttype)
    rinit = StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], 0, 0)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, nalleles, priorcv, onstates, rtarget[end], propcv, splicetype, prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nsamples, 1000, 0, 200.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    h = likelihoodfn(fits.parml, data, model)
    h, StochasticGene.datapdf(data), fits, stats, measures, data, model, options
end

function test_fit_trace(; G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 70, 0.9], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep)); [20, 5, 100, 10, 0.9]], nsamples=1000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[collect(1:num_rates(transitions, R, S, insertstep)-1); collect(num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+5)], propcv=0.01, cv=100.0, interval=1.0)
    traces = simulate_trace_vector(rtarget, transitions, G, R, S, interval, totaltime, ntrials)
    data = StochasticGene.TraceData("trace", "test", interval, traces)
    model = load_model(data, rinit, Float64[], fittedparam, tuple(), transitions, G, R, S, insertstep, 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, "", prob_GaussianMixture, 5, 5)
    options = StochasticGene.MHOptions(nsamples, 0, 0, 100.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    StochasticGene.get_rates(fits.parml,model),rtarget,fits, stats, measures, data, model, options
end
