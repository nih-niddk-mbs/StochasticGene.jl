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

"""
    test_compare(; r, transitions, G, R, S, insertstep, nRNA, nalleles, bins, total, tol, onstates, dttype)

Compare simulated and chemical master equation histograms for a given parameter set.

# Arguments
- `r`: Rate parameters.
- `transitions`, `G`, `R`, `S`, `insertstep`: Model structure.
- `nRNA`, `nalleles`: RNA and allele counts.
- `bins`: Histogram bins.
- `total`: Number of simulation steps.
- `tol`: Simulation tolerance.
- `onstates`, `dttype`: State and dwell time types.

# Returns
- Tuple of (chemical master histogram, array of simulated histograms).
"""
function test_compare(; r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231], transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=2, S=2, insertstep=1, nRNA=150, nalleles=2, bins=[collect(5/3:5/3:200), collect(5/3:5/3:200), collect(0.1:0.1:20), collect(0.1:0.1:20)], total=10000000, tol=1e-6, onstates=[Int[], Int[], [2, 3], [2, 3]], dttype=["ON", "OFF", "ONG", "OFFG"])
    hs = test_sim(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates[[1, 3]], bins[[1, 3]], total, tol)
    h = test_cm(r, transitions, G, R, S, insertstep, nRNA, nalleles, onstates, dttype, bins)
    hs = StochasticGene.normalize_histogram.(hs)
    return h, make_array(hs)
end

"""
    test_compare_coupling(; r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling, total, tol)

Compare simulated and chemical master equation histograms for coupled models.

# Arguments
- `r`: Rate parameters.
- `transitions`, `G`, `R`, `S`, `insertstep`: Model structure.
- `onstates`, `dttype`, `bins`: State and dwell time types, histogram bins.
- `coupling`: Coupling structure.
- `total`: Number of simulation steps.
- `tol`: Simulation tolerance.

# Returns
- Tuple of (chemical master histogram, array of simulated histograms).
"""
function test_compare_coupling(; r=[0.38, 0.1, 0.23, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2, 1.0, 0.45, 0.2, 0.43, 0.3, 0.52, 0.31, 0.3, 0.86, 0.5, 1.0, 2.9], transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), G=(3, 3), R=(2, 2), S=(2, 2), insertstep=(1, 1), onstates=[[Int[], Int[], [3], [3]], [Int[], Int[], [3], [3]]], dttype=[["ON", "OFF", "ONG", "OFFG"], ["ON", "OFF", "ONG", "OFFG"]], bins=[[collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)], [collect(1:30), collect(1:30), collect(1.0:30), collect(1.0:30)]], coupling=((1, 2), (Int[], [1]), [3, 0], [0, 4], 1), total=10000000, tol=1e-6)
    hs = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, nhist=0, noiseparams=0, onstates=simDT_convert(onstates), bins=simDT_convert(bins), totalsteps=total, tol=tol)
    h = test_CDT(r, transitions, G, R, S, insertstep, onstates, dttype, bins, coupling)
    for i in eachindex(hs)
        hs[i] = StochasticGene.normalize_histogram.(hs[i])
    end
    return make_array(vcat(h...)), make_array(vcat(hs...))
end

"""
    test_fit_simrna(; rtarget, transitions, G, nRNA, nalleles, fittedparam, fixedeffects, rinit, totalsteps, nchains)

Fit a simulated RNA histogram using the provided parameters and compare to the target.

# Arguments
- `rtarget`: Target rate parameters.
- `transitions`, `G`, `nRNA`, `nalleles`: Model structure and counts.
- `fittedparam`, `fixedeffects`, `rinit`: Fitting options.
- `totalsteps`: Number of simulation steps.
- `nchains`: Number of MCMC chains.

# Returns
- Tuple of (fitted rates, target rates).
"""
function test_fit_simrna(; rtarget=[0.33, 0.19, 20.5, 1.0], transitions=([1, 2], [2, 1]), G=2, nRNA=100, nalleles=2, fittedparam=[1, 2, 3], fixedeffects=tuple(), rinit=[0.1, 0.1, 0.1, 1.0], totalsteps=100000, nchains=1)
    h = simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=totalsteps, nalleles=nalleles)[1]
    data = RNAData("", "", nRNA, h)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0), fittedparam, fixedeffects, transitions, G, 0, 0, 0, "", nalleles, 10.0, Int[], rtarget[end], 0.02, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(1000000, 0, 0, 20.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    StochasticGene.get_rates(fits.parml, model), rtarget
end

"""
    test_fit_rna(; gene, G, nalleles, propcv, fittedparam, fixedeffects, transitions, rinit, datacond, datapath, label, root, nchains)

Fit a real RNA histogram using the provided parameters and compare to the data.

# Arguments
- `gene`, `G`, `nalleles`: Gene and model structure.
- `propcv`, `fittedparam`, `fixedeffects`, `transitions`, `rinit`: Fitting options.
- `datacond`, `datapath`, `label`, `root`: Data and file options.
- `nchains`: Number of MCMC chains.

# Returns
- Tuple of (predicted histogram, normalized data histogram).
"""
function test_fit_rna(; gene="CENPL", G=2, nalleles=2, propcv=0.05, fittedparam=[1, 2, 3], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), rinit=[0.01, 0.1, 1.0, 0.01006327034802035], datacond="MOCK", datapath="data/HCT116_testdata", label="scRNA_test", root=".", nchains=1)
    data = load_data("rna", [], folder_path(datapath, root, "data"), label, gene, datacond, (), 1.0)
    model = load_model(data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, 1.0, [], 1.0), fittedparam, fixedeffects, transitions, 2, 0, 0, 1, "", nalleles, 10.0, Int[], rinit[end], propcv, prob_Gaussian, [], 1, tuple(), tuple(), nothing)
    options = StochasticGene.MHOptions(100000, 0, 0, 60.0, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    h = predictedfn(fits.parml, data, model)
    h, normalize_histogram(data.histRNA)
end

"""
    test_fit_rnaonoff(; G, R, S, transitions, insertstep, rtarget, rinit, nsamples, nhist, nalleles, onstates, bins, fittedparam, propcv, priorcv, splicetype, nchains)

Fit a simulated RNA on-off histogram using the provided parameters and compare to the data.

# Arguments
- `G`, `R`, `S`, `transitions`, `insertstep`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `nhist`, `nalleles`, `onstates`, `bins`: Data and simulation options.
- `fittedparam`, `propcv`, `priorcv`, `splicetype`, `nchains`: Fitting options.

# Returns
- Tuple of (predicted histogram, data PDF).
"""
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

"""
    test_fit_rnadwelltime(; rtarget, transitions, G, R, S, insertstep, nRNA, nsamples, nalleles, onstates, bins, dttype, fittedparam, propcv, priorcv, splicetype, maxtime, nchains)

Fit a simulated RNA dwell time histogram using the provided parameters and compare to the data.

# Arguments
- `rtarget`, `transitions`, `G`, `R`, `S`, `insertstep`: Model structure and rates.
- `nRNA`, `nsamples`, `nalleles`, `onstates`, `bins`, `dttype`: Data and simulation options.
- `fittedparam`, `propcv`, `priorcv`, `splicetype`, `maxtime`, `nchains`: Fitting options.

# Returns
- Tuple of (predicted histogram, data PDF).
"""
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

"""
    test_fit_trace(; G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, nchains)

Fit a simulated trace dataset using the provided parameters and compare to the target.

# Arguments
- `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `nchains`: Fitting options.

# Returns
- Tuple of (fitted rates, target rates).
"""
function test_fit_trace(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], nsamples=10000000, onstates=Int[], totaltime=1000.0, ntrials=100, fittedparam=[1:num_rates(transitions, R, S, insertstep)-1; num_rates(transitions, R, S, insertstep)+1:num_rates(transitions, R, S, insertstep)+1], propcv=0.01, cv=100.0, noisepriors=[0.0, 0.1, 1.0, 0.1], nchains=1, zeromedian=true, maxtime=100.0, initprior=0.1)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials)
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale))
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), nothing)
    elongationtime = mean_elongationtime(rtarget, transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R); fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]]

    rinit = isempty(tuple()) ? set_rinit([], priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, priorcv, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, lsoda(), tuple(), tuple(), nothing, zeromedian)
    options = StochasticGene.MHOptions(nsamples, nsamples, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    # return fits, stats, measures, data, model, options
    nrates = num_rates(model)
    stats.medparam[1:nrates-1], rtarget[1:nrates-1]
end

"""
    test_fit_trace_hierarchical(; G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, hierarchical, method, maxtime, nchains)

Fit a simulated trace dataset using a hierarchical model and compare to the target.

# Arguments
- `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `hierarchical`, `method`, `maxtime`, `nchains`: Fitting options.

# Returns
- Tuple of (median fitted parameters, target parameters).
"""
function test_fit_trace_hierarchical(; traceinfo=(1.0, 1.0, -1, 1.0, 0.5), G=2, R=2, S=0, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.05, 0.2, 0.1, 0.15, 0.1, 1.0, 50, 5, 50, 5], rinit=[], nsamples=10000000, onstates=Int[], totaltime=1000.0, ntrials=100, fittedparam=[1, 2, 3, 4, 5, 6], propcv=0.01, noisepriors=[0.0, 0.1, 1.0, 0.1], hierarchical=(2, [7], tuple()), method=(lsoda(), true), maxtime=180.0, nchains=1, zeromedian=true)
    rh = 50.0 .+  1.0 * randn(ntrials)
    tracer = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, traceinfo[1], totaltime, ntrials, hierarchical=(6, rh))
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] != 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    data = TraceData{String,String,Tuple}("trace", "gene", traceinfo[1], (trace, background, weight, nframes, tracescale))
    # trace, tracescale = zero_median(tracer, zeromedian)
    # data = StochasticGene.TraceData("trace", "test", interval, (trace, [], 0.0, 1))(trace, background, weight, nframes, tracescale)
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, tuple(), nothing)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), nothing)
    fittedparam = set_fittedparam(fittedparam, "trace", transitions, R, S, insertstep, noisepriors, tuple(), tuple())
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 0.1, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian, noisepriors, method, hierarchical, tuple(), nothing, zeromedian)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options, nchains)
    nrates = num_rates(model)
    # h1 = StochasticGene.get_rates(fits.parml, model)[1:nrates]
    h2 = rtarget[1:nrates-1]
    h1 = stats.medparam[1:nrates-1]
    # h2 = [rtarget[fittedparam]; [50.0; 10 / 50.0]; rh]
    return h1, h2
end

"""
    test_fit_tracejoint(; coupling, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, maxtime, method)

Fit a simulated joint trace dataset for coupled models and compare to the target.

# Arguments
- `coupling`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure and coupling.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `maxtime`, `method`: Fitting options.

# Returns
- Tuple of (fitted rates, target rates).
"""
function test_fit_tracejoint(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 1.0, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 1.0, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, method=lsoda())
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


### end of functions used in runtest

### functions to be used in the future
"""
    test_fit_tracejoint_hierarchical(; coupling, G, R, S, insertstep, transitions, rtarget, rinit, hierarchical, method, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, noisepriors, maxtime, decayrate)

Fit a simulated joint trace dataset for coupled models using a hierarchical model and compare to the target.

# Arguments
- `coupling`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure and coupling.
- `rtarget`, `rinit`: Rate parameters.
- `hierarchical`, `method`: Hierarchical model and method options.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `noisepriors`, `maxtime`, `decayrate`: Fitting options.

# Returns
- Tuple of (fitted rates, target rates, fits, stats, measures, model, data, options).
"""
function test_fit_tracejoint_hierarchical(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(lsoda(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, decayrate=1.0)
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

"""
    test_fit_trace_grid(; grid, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, totaltime, ntrials, fittedparam, propcv, cv, interval, weight, nframes, noisepriors, maxtime)

Fit a simulated trace grid dataset using the provided parameters and compare to the target.

# Arguments
- `grid`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `weight`, `nframes`, `noisepriors`, `maxtime`: Fitting options.

# Returns
- Tuple of (fits, stats, measures, data, model, options).
"""
function test_fit_trace_grid(; grid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[fill(0.1, num_rates(transitions, R, S, insertstep) - 1); 0.01; [50, 15, 200, 20]; 0.1], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 6, 7, 8, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=grid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
    priormean = set_priormean([], transitions, R, S, insertstep, rtarget[num_rates(transitions, R, S, insertstep)], noisepriors, mean_elongationtime(rtarget, transitions, R), tuple(), tuple(), grid)
    rinit = set_rinit(rinit, priormean)
    fittedparam = set_fittedparam(fittedparam, "tracegrid", transitions, R, S, insertstep, noisepriors, tuple(), grid)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian_grid, noisepriors, lsoda(), tuple(), tuple(), grid)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end

"""
    test_fit_trace_grid_hierarchical(; grid, G, R, S, insertstep, transitions, rtarget, rinit, nsamples, onstates, hierarchical, totaltime, ntrials, fittedparam, propcv, cv, interval, weight, nframes, noisepriors, maxtime)

Fit a simulated trace grid dataset using a hierarchical model and compare to the target.

# Arguments
- `grid`, `G`, `R`, `S`, `insertstep`, `transitions`: Model structure.
- `rtarget`, `rinit`: Rate parameters.
- `hierarchical`: Hierarchical model options.
- `nsamples`, `onstates`, `totaltime`, `ntrials`: Data and simulation options.
- `fittedparam`, `propcv`, `cv`, `interval`, `weight`, `nframes`, `noisepriors`, `maxtime`: Fitting options.

# Returns
- Tuple of (fits, stats, measures, data, model, options).
"""
function test_fit_trace_grid_hierarchical(; grid=4, G=2, R=1, S=1, insertstep=1, transitions=([1, 2], [2, 1]), rtarget=[0.02, 0.1, 0.5, 0.2, 0.1, 0.01, 50, 15, 200, 20, 0.2], rinit=[], nsamples=5000, onstates=Int[], hierarchical=(2, [7], tuple()), totaltime=1000.0, ntrials=10, fittedparam=[1, 2, 3, 4, 5, 11], propcv=0.01, cv=100.0, interval=1.0, weight=0, nframes=1, noisepriors=[50, 15, 200, 70], maxtime=10.0)
    traces = sim_grid(r=rtarget[1:end-1], p=rtarget[end], Ngrid=grid, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, totaltime=totaltime, interval=interval, ntrials=ntrials)
    data = StochasticGene.TraceData("tracegrid", "test", interval, (traces, [], weight, nframes))
    priormean = set_priormean([], transitions, R, S, insertstep, rtarget[num_rates(transitions, R, S, insertstep)], noisepriors, mean_elongationtime(rtarget, transitions, R), hierarchical, tuple(), grid)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), tuple(), grid)
    fittedparam = set_fittedparam(fittedparam, "tracegrid", transitions, R, S, insertstep, noisepriors, tuple(), grid)
    model = load_model(data, rinit, priormean, fittedparam, tuple(), transitions, G, R, S, insertstep, "", 1, 10.0, Int[], rtarget[num_rates(transitions, R, S, insertstep)], propcv, prob_Gaussian_grid, noisepriors, (lsoda(), true), hierarchical, tuple(), grid)
    options = StochasticGene.MHOptions(nsamples, 0, 0, maxtime, 1.0, 1.0)
    fits, stats, measures = run_mh(data, model, options)
    fits, stats, measures, data, model, options
end

### functions to be used in the future

### development test functions

test_fit0(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, 0.0), "test", "test", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.2, 1.0, 0.2], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.01, 0.01], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.2, 1.0, 0.2], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian)

test_fit1(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, 0.5), "test1", "test1", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.5, 0.25, 1.0, 0.2], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [0.5, 0.2, 1.0, 0.2], (), "ml", propcv, 20000, 10000, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian)

test_fit2(; maxtime=6.0, propcv=0.02, zeromedian=true) = fit(1, "trace", String[], "data/inhibition/inhibition/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.78, 0.0), "inhibition-test", "inhibition-test", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 7, 9], (), ([1, 2], [2, 1]), 2, 4, 0, 3, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.8, 0.8, 0.8, 0.8, 1.0, 0.0, 0.25, 1.0, 0.25], [1.0, 1.0, 0.2, 0.1, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.25, 1.0, 0.25], (), "ml", propcv, 2000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian, 3, 1);

test_fitt(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "trace", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 1.0, 0.5), "genetogther", "genetogther", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.5, 0.2, 1.0, 0.2], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [0.5, 0.2, 1.0, 0.2], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian)

test_fit_h(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.92, 0.0), "testh", "testh", "trace-HBEC-nstate-h_gene", "trace-HBEC-nstate-h_gene", [1, 2, 3, 4, 5, 6], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.2, 1.0, 0.1, 1.0, 1.0, 1.0, 0.1, 0.1, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.2, 1.0, 0.1], (2, [8], ()), "ml", propcv, 100000, 100000, 0, 1.0, 100.0, 1.0, false, false, false, (lsoda(), true), zeromedian)

test_fit_h2(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "trace", String[], "data/inhibition/inhibition/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.78, 0.0), "inhibition-test", "inhibition-test", "trace-HBEC-nstate-h_gene", "trace-HBEC-nstate-h_gene", [1, 2, 3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 0.0, 0.25, 1.0, 0.25, 1.0, 1.0, 1.0, 0.0, 0.1, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.1, 0.1, 0.1, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, Int64[], 1.0, "", prob_Gaussian, [0.0, 0.25, 1.0, 0.25], (2, [8], ()), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (lsoda(), true), zeromedian)

test_fit_coupleda(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.0, 0.0]), "test", "test", "tracejoint-HBEC-nstate_enhancer-geneR5", "tracejoint-HBEC-nstate_enhancer-geneR5", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (3, 3), (1, 0), (1, 1), ((1, 2), ((), (1,)), ([4, 5, 6], 0), (0, 5), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.0, 0.2, 1.0, 0.1], [0.0, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian)

test_fit_coupled2a(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.5, 0.5]), "test1", "test1", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25, 28], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (3, 3), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.5, 0.2, 1.0, 0.1], [0.5, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, true, false, false, lsoda(), zeromedian)

test_fit_coupled(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.5, 0.5]), "test1", "test1", "tracejoint-HBEC-nstate_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.5, 0.2, 1.0, 0.1], [0.5, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian)
test_fit_coupled0(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.0, 0.0]), "test", "test", "tracejoint-HBEC-nstate_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.0, 0.2, 1.0, 0.1], [0.0, 0.2, 1.0, 0.1]), (), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda(), zeromedian)

test_fit_coupled_h(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=false) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.5, 0.5]), "test1", "test1", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.5, 0.2, 1.0, 0.1], [0.5, 0.2, 1.0, 0.1]), (2, [9, 20], ()), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (lsoda(), true), zeromedian)
test_fit_coupled_h0(; nchains=1, maxtime=6.0, propcv=0.01, zeromedian=true) = fit(nchains, "tracejoint", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", ["enhancer", "gene"], (1.0, 1.0, -1, [0.46, 0.86], [0.0, 0.0]), "test", "test", "tracejoint-HBEC-nstate-h_enhancer-gene11", "tracejoint-HBEC-nstate-h_enhancer-gene11", [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 24], (), (([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])), (3, 3), (1, 1), (1, 0), (1, 1), ((1, 2), ((), (1,)), (1, 0), (0, 1), 1), nothing, ".", maxtime, (6.5, 5.0), Float64[], 10.0, 1, Int64[], 1.0, "", prob_Gaussian, ([0.0, 0.2, 1.0, 0.1], [0.0, 0.2, 1.0, 0.1]), (2, [9, 20], ()), "ml", propcv, 1000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, (lsoda(), true), zeromedian)


test_rnacount(; nchains=1, gene="RPLP1", datapath="data/U3AS4/counts/WT-UTR", maxtime=60.0, propcv=0.01) = fit(nchains=nchains, datatype="rnacount", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath=datapath, cell="U3A", gene=gene, datacond="WT-UTR", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv);

test_rna(; nchains=1, gene="RPLP1", datapath="data/U3AS4/histograms/", maxtime=60.0, propcv=0.01) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath=datapath, cell="U3A", gene=gene, datacond="WT-UTR", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv);

test_fit_G(; nchains=1, maxtime=60.0, propcv=0.01) = fit(nchains, "trace", String[], "data/MYC_gene_MYC_enhancer_together/traces/", "MYC", "HBEC", "gene", (1.0, 1.0, -1, 0.86, [25000.0, 0.0]), "genetogther", "genetogther", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene", [1, 2, 3], (), ([1, 2], [2, 1]), 2, 0, 0, 1, (), nothing, ".", maxtime, 5.0, [0.01, 0.01, 0.1, 1.0, 20000.0, 15000.0, 80000.0, 4000.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25], [5.0, 5.0, 0.2, 1.0, 0.5, 0.5, 1.0, 0.1, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0], 1, [2], 1.0, "", prob_Gaussian, [20000.0, 15000.0, 80000.0, 4000.0], (), "ml", propcv, 200000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda())

fit_myc3(; nchains=1, maxtime=60.0, propcv=0.05) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1], [2, 3], [3, 2]), G=3, R=0, S=0, insertstep=1, datapath="data/FISH", cell="HBEC", gene="MYC", datacond="gene", infolder="test", resultfolder="test", fittedparam=[1, 2, 3, 4, 5], maxtime=maxtime, ratetype="ml", propcv=propcv, burst=true);

fit_myc2(; nchains=1, maxtime=60.0, propcv=0.05) = fit(nchains=nchains, datatype="rna", transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, datapath="data/FISH", cell="HBEC", gene="MYC", datacond="gene", infolder="test", resultfolder="test", fittedparam=[1, 2, 3], maxtime=maxtime, ratetype="ml", propcv=propcv, burst=true);

# @time fit(16, "trace", String[], "data/inhibition/control/", "MYC", "HBEC", "gene", (1.6666666666666667, 1.0, -1, 0.92, [290.0, 140.0]), "control-2025-03-15", "control-2025-03-15", "trace-HBEC-nstate_gene", "trace-HBEC-nstate_gene",[1, 2,3, 4, 5, 6, 8, 9], (), ([1, 2], [2, 1]), 2, 3, 0, 1, (), nothing, ".", 43000.0, 5.0, [0.01, 0.01, 0.1, 0.6, 0.6, 0.6, 1.0, 290.0, 140.0, 1200.0, 175.0], [2.0, 2.0, 0.2, 0.1, 0.1, 0.1, 1.0, 0.5, 0.5, 0.1, 0.1], 1, Int64[], 1.0, "", prob_Gaussian, [290.0, 140.0, 1200.0, 175.0], (), "ml", 0.02, 2000000, 0, 0, 1.0, 100.0, 1.0, false, false, false, lsoda())

function test_tracejoint_h(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.03, 0.1, 0.5, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], hierarchical=(2, [8, 17], tuple()), method=(lsoda(), true), nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=20, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 100, 20], [50, 30, 100, 20]), maxtime=300.0, decayrate=1.0)
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


function test_tracejoint(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20, 0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5], rinit=Float64[], nsamples=5000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([100, 50, 200, 100], [100, 50, 200, 100]), maxtime=300.0, method=lsoda())
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

function test_predicted_states(; coupling=((1, 2), (tuple(), tuple(1)), (2, 0), (0, 1), 1), G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1), transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])), rtarget=[0.1, 0.3, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 20, 200, 20, 0.1, 0.2, 0.4, 0.2, 0.1, 50, 20, 200, 20, -0.5], rinit=Float64[], nsamples=20000, onstates=Int[], totaltime=1000.0, ntrials=10, fittedparam=Int[1, 2, 3, 4, 5, 6, 12, 13, 14, 15, 21], propcv=0.01, cv=100.0, interval=1.0, noisepriors=([50, 30, 500, 20], [50, 30, 500, 20]), maxtime=300.0, decayrate=1.0, totalsteps=20, verbose=false, hierarchical=(2, [8, 17], tuple()), method=(lsoda(), true))
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





"""
    compute_psis_loo(fits::Fit)

Compute PSIS-LOO (Pareto Smoothed Importance Sampling Leave-One-Out Cross-Validation).
Returns a tuple of (elpd_loo, elpd_loo_se) where:
- elpd_loo is the expected log pointwise predictive density for a new dataset
- elpd_loo_se is the standard error of elpd_loo

This implementation reuses the pointwise log-likelihoods stored in fits.lppd
that were computed during MCMC sampling for WAIC.
"""
function compute_psis_loo(fits::Fit)
    # Get pointwise log-likelihoods from fits.lppd
    n_samples = size(fits.param, 2)
    n_obs = length(fits.lppd)

    # Compute importance weights for each observation
    r_eff = compute_relative_eff(fits.param)  # Relative effective sample size

    # Compute LOO log-likelihoods using importance sampling
    loo_ll = zeros(n_obs)
    for i in 1:n_obs
        # Get pointwise log-likelihoods for this observation
        pointwise_ll = fits.lppd[i]

        # Compute importance weights for this observation
        weights = exp.(fits.param[i, :] .- maximum(fits.param[i, :]))
        weights .*= r_eff[i]
        weights ./= sum(weights)

        # Smooth weights using Pareto smoothing
        smoothed_weights = pareto_smooth_weights(weights)

        # Compute LOO estimate using smoothed weights
        loo_ll[i] = sum(smoothed_weights .* pointwise_ll)
    end

    # Compute elpd_loo and its standard error
    elpd_loo = sum(loo_ll)
    elpd_loo_se = sqrt(n_obs * var(loo_ll))

    return (elpd_loo, elpd_loo_se)
end

"""
    compute_relative_eff(params)

Compute relative effective sample size for importance sampling.
"""
function compute_relative_eff(params)
    n_samples = size(params, 2)
    ess = compute_ess(params)
    return ess ./ n_samples
end

"""
    pareto_smooth_weights(weights)

Apply Pareto smoothing to importance weights.
"""
function pareto_smooth_weights(weights)
    # Sort weights
    sorted_idx = sortperm(weights, rev=true)
    sorted_weights = weights[sorted_idx]

    # Find Pareto tail
    tail_idx = find_tail(sorted_weights)

    if tail_idx > 0
        # Fit Pareto distribution to tail
        k = fit_pareto(sorted_weights[tail_idx:end])

        # Replace tail with smoothed values
        smoothed = copy(weights)
        smoothed[sorted_idx[tail_idx:end]] = smooth_tail(sorted_weights[tail_idx:end], k)
        return smoothed
    end

    return weights
end

"""
    find_tail(weights)

Find the index where the Pareto tail begins.
"""
function find_tail(weights)
    n = length(weights)
    for i in 2:n
        if weights[i] / weights[1] < 0.1
            return i
        end
    end
    return 0
end

"""
    fit_pareto(weights)

Fit a Pareto distribution to the tail of weights.
Returns the shape parameter k.
"""
function fit_pareto(weights)
    n = length(weights)
    log_weights = log.(weights)
    k = 1 / (mean(log_weights) - log_weights[1])
    return k
end

"""
    smooth_tail(weights, k)

Smooth the tail of weights using the fitted Pareto distribution.
"""
function smooth_tail(weights, k)
    n = length(weights)
    smoothed = zeros(n)
    for i in 1:n
        smoothed[i] = weights[1] * (i / n)^(-1 / k)
    end
    return smoothed
end
