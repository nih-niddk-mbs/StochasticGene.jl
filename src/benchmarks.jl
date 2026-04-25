# This file is part of StochasticGene.jl
#
# benchmarks.jl
#
# Performance and experimental comparison helpers. Loaded from `src/test.jl` so the
# existing REPL workflow still exposes these `benchmark_*` entry points from the package.

"""
    profile_trace_prediction(; kwargs...)

Profile the trace prediction path used by `make_traces_dataframe` / `predict_traces`.
This is an interactive developer helper and is not called from `runtests.jl`.
"""
function profile_trace_prediction(;
    coupling=((1, 2), [(1, 2, 2, 1)]),
    G=(2, 2), R=(2, 1), S=(1, 0), insertstep=(1, 1),
    transitions=(([1, 2], [2, 1]), ([1, 2], [2, 1])),
    rtarget=[0.03, 0.1, 0.5, 0.4, 0.4, 0.01, 0.01, 50, 30, 100, 20,
             0.02, 0.05, 0.2, 0.2, 0.1, 50, 30, 100, 20, -0.5],
    interval=1.0, totaltime=200.0, ntrials=3,
    trace_specs=nothing)

    trace = simulate_trace_vector(rtarget, transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; verbose=false)
    ts = trace_specs === nothing ? default_trace_specs_for_coupled((interval, 1.0, -1.0), [false, false], 2) : trace_specs
    data = StochasticGene.TraceData("tracejoint", "test", interval, (trace, [], [0.0, 0.0], 1), [spec.unit for spec in ts])

    p = @profile StochasticGene.make_traces_dataframe(data, rtarget, transitions, G, R, S, insertstep,
                                                      prob_Gaussian, 4, "", true, false, coupling)
    return p
end

# ════════════════════════════════════════════════════════════════════════════════════
# SECTION 5: MH / NUTS / ADVI BENCHMARK COMPARISON HELPERS
# ════════════════════════════════════════════════════════════════════════════════════

"""
    benchmark_inference_simrna_small(; seed=42, totalsteps=80_000, nhist=60,
                                     transitions=([1,2],[2,1]), G=2,
                                     rtarget=[0.33, 0.19, 20.5, 1.0],
                                     rinit=[0.1, 0.1, 0.1, 1.0],
                                     fittedparam=[1,2,3])

Construct a small synthetic RNA histogram inference problem suitable for side-by-side
benchmarking of MH, NUTS, and ADVI.

Returns a named tuple with `data`, `model`, and `meta`.
"""
function benchmark_inference_simrna_small(; seed::Int=42,
                                          totalsteps::Int=80_000,
                                          nhist::Int=60,
                                          transitions=([1,2],[2,1]),
                                          G::Int=2,
                                          rtarget::AbstractVector{<:Real}=[0.33, 0.19, 20.5, 1.0],
                                          rinit::AbstractVector{<:Real}=[0.1, 0.1, 0.1, 1.0],
                                          fittedparam::AbstractVector{<:Integer}=[1,2,3])
    Random.seed!(seed)
    h = simulator(collect(Float64, rtarget), transitions, G, 0, 0, 0;
                  nhist=nhist, totalsteps=totalsteps, nalleles=2)[1]
    data = RNAData{typeof(nhist),typeof(h)}("benchmark_simrna", "synthetic", nhist, h, 1.0, [])
    priormean = prior_ratemean(transitions, 0, 0, 1, Float64(rtarget[end]), [], 1.0)
    model = load_model(data,
                       collect(Float64, rinit),
                       priormean,
                       collect(Int, fittedparam),
                       tuple(),
                       transitions,
                       G,
                       0,
                       0,
                       0,
                       "",
                       2,
                       10.0,
                       Int[],
                       Float64(rtarget[end]),
                       0.02,
                       prob_Gaussian,
                       [],
                       1,
                       tuple(),
                       tuple(),
                       nothing)
    meta = (scenario=:simrna_small,
            true_rates=collect(Float64, rtarget),
            fittedparam=collect(Int, fittedparam))
    return (data=data, model=model, meta=meta)
end

"""
    benchmark_inference_trace_gr2r2(; seed=42, interval=5/3, totaltime=60.0, ntrials=2,
                                    transitions=([1,2],[2,1]), G=2, R=2, S=0,
                                    insertstep=1,
                                    rtarget=[0.3, 0.15, 0.6, 0.3, 0.0, 0.1, 1.0, 0.1],
                                    fittedparam=[1,2,3])

Construct a small synthetic single-trace inference problem (`TraceData`) for comparing
MH / NUTS / ADVI on the trace likelihood path.

Returns a named tuple with `data`, `model`, and `meta`.
"""
function benchmark_inference_trace_gr2r2(; seed::Int=42,
                                         interval::Float64=5/3,
                                         totaltime::Float64=60.0,
                                         ntrials::Int=2,
                                         transitions=([1,2],[2,1]),
                                         G::Int=2,
                                         R::Int=2,
                                         S::Int=0,
                                         insertstep::Int=1,
                                         rtarget::AbstractVector{<:Real}=[0.3, 0.15, 0.6, 0.3, 0.0, 0.1, 1.0, 0.1],
                                         fittedparam::AbstractVector{<:Integer}=[1,2,3],
                                         propcv::Float64=0.01,
                                         noisepriors=[0.0, 0.1, 1.0, 0.1],
                                         zeromedian::Bool=true,
                                         initprior::Float64=0.1)
    Random.seed!(seed)
    tracer = simulate_trace_vector(collect(Float64, rtarget), transitions, G, R, S, insertstep,
                                   interval, totaltime, ntrials)
    trace, tracescale = zero_median(tracer, zeromedian)
    nframes = round(Int, mean(size.(trace, 1)))
    data = TraceData{String,String,Tuple}("trace", "gene", interval,
                                          (trace, 0.0, 0.0, nframes, tracescale), Int[])
    elongationtime = mean_elongationtime(collect(Float64, rtarget), transitions, R)
    priormean = [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R);
                 fill(0.05, max(0, S - insertstep + 1)); 1.0; noisepriors]
    priorcv = [fill(1.0, length(transitions)); 0.1; fill(0.1, R);
               fill(0.1, max(0, S - insertstep + 1)); 1.0; [0.5, 0.5, 0.1, 0.1]]
    rinit = set_rinit(collect(Float64, rtarget), priormean)
    model = load_model(data,
                       rinit,
                       priormean,
                       collect(Int, fittedparam),
                       tuple(),
                       transitions,
                       G,
                       R,
                       S,
                       insertstep,
                       "",
                       1,
                       priorcv,
                       Int[],
                       rtarget[num_rates(transitions, R, S, insertstep)],
                       propcv,
                       prob_Gaussian,
                       noisepriors,
                       Tsit5(),
                       tuple(),
                       tuple(),
                       nothing,
                       zeromedian)
    meta = (scenario=:trace_gr2r2,
            true_rates=collect(Float64, rtarget),
            fittedparam=collect(Int, fittedparam),
            zeromedian=zeromedian)
    return (data=data, model=model, meta=meta)
end

"""
    benchmark_inference_trace_coupled_3x3(; seed=42, interval=5/3, totaltime=60.0, ntrials=2,
                                          fittedparam=[1,2,3,14], coupled_stack=:full)

Construct a small synthetic coupled-trace inference problem with two 3-state promoter units
and 3 reporter positions each, suitable for comparing MH / NUTS / ADVI on the coupled trace path.

Returns a named tuple with `data`, `model`, and `meta`.
"""
function benchmark_inference_trace_coupled_3x3(; seed::Int=42,
                                               interval::Float64=1.0,
                                               totaltime::Float64=800.0,
                                               ntrials::Int=4,
                                               coupling=((1, 2), [(1, 2, 2, 1)], [:free]),
                                               G=(3, 3),
                                               R=(3, 3),
                                               S=(0, 0),
                                               insertstep=(1, 1),
                                               transitions=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
                                               rtarget::AbstractVector{<:Real}=Float64[
                                                   0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2,
                                                   0.0, 0.1, 0.5, 0.15,
                                                   0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2,
                                                   0.0, 0.1, 0.9, 0.2,
                                                   -0.4,
                                               ],
                                               noisepriors=([0., .1, 1., .1], [0., .1, 1., .1]),
                                               trace_specs=nothing,
                                               fittedparam::AbstractVector{<:Integer}=Int[19],
                                               propcv::Float64=0.2,
                                               method=Tsit5(),
                                               coupled_stack::Symbol=:full,
                                               zeromedian::Bool=true)
    Random.seed!(seed)
    trace_specs_eff = trace_specs === nothing ? StochasticGene.default_trace_specs_for_coupled((interval, 1.0, -1.0), [false, false], 2) : trace_specs
    units = [spec.unit for spec in trace_specs_eff]
    trace = simulate_trace_vector(collect(Float64, rtarget), transitions, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; noiseparams=[4, 4])
    data = TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units)), 1), units)
    priormean = set_priormean([], transitions, R, S, insertstep, 1.0, noisepriors,
                              mean_elongationtime(collect(Float64, rtarget), transitions, R), tuple(), coupling, nothing)
    rinit = collect(Float64, rtarget)
    nr = num_rates(transitions, R, S, insertstep)
    model = load_model(data, rinit, priormean, collect(Int, fittedparam), tuple(), transitions, G, R, S, insertstep,
                       "", 1, 10.0, Int[], rtarget[nr], propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing,
                       zeromedian; coupled_stack=coupled_stack)
    meta = (scenario=:trace_coupled_3x3,
            true_rates=collect(Float64, rtarget),
            fittedparam=collect(Int, fittedparam),
            trace_specs=trace_specs_eff,
            coupling=coupling,
            coupled_stack=coupled_stack)
    return (data=data, model=model, meta=meta)
end

"""
    benchmark_inference_trace_coupled_3x3_g3r0(; seed=42, interval=5/3, totaltime=60.0, ntrials=2,
                                               fittedparam=[1,2,3,9,10,11], coupled_stack=:full)

Construct a small synthetic coupled-trace inference problem with two 3-state promoter units
and no reporter progression (`R=(0,0)`), suitable for comparing MH / NUTS / ADVI on the
coupled trace path while varying only promoter and noise parameters.

Returns a named tuple with `data`, `model`, and `meta`.
"""
function benchmark_inference_trace_coupled_3x3_g3r0(; seed::Int=42,
                                                    interval::Float64=1.0,
                                                    totaltime::Float64=800.0,
                                                    ntrials::Int=4,
                                                    coupling=((1, 2, 3), [(3, 1, 1, 2), (3, 3, 2, 2)], [:inhibit, :inhibit]),
                                                    G=(3, 3, 3),
                                                    R=(3, 3, 0),
                                                    S=(0, 0, 0),
                                                    insertstep=(1, 1, 1),
                                                    transitions=((([1, 2], [2, 1], [2, 3], [3, 2])), (([1, 2], [2, 1], [2, 3], [3, 2])), (([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]))),
                                                    rtarget::Union{Nothing,Vector{Float64}}=nothing,
                                                    units=[1, 2],
                                                    noisepriors=([0., .1, 1., .1], [0., .1, 1., .1], Float64[]),
                                                    trace_specs=nothing,
                                                    fixedeffects=nothing,
                                                    fittedparam=nothing,
                                                    propcv::Float64=0.2,
                                                    method=Tsit5(),
                                                    coupled_stack::Symbol=:full,
                                                    zeromedian::Bool=true)
    Random.seed!(seed)
    transitions_eff = (
        ([1, 2], [2, 1], [2, 3], [3, 2]),
        ([1, 2], [2, 1], [2, 3], [3, 2]),
        ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]),
    )
    nrates_per_unit = num_rates(transitions_eff, R, S, insertstep)
    ncpl = ncoupling(coupling)
    n_noise = sum(length.(noisepriors))
    n_total = sum(nrates_per_unit) + n_noise + ncpl
    unit3_rate_start = sum(nrates_per_unit[1:2]) + sum(length.(noisepriors[1:2])) + 1
    unit3_rate_end = unit3_rate_start + nrates_per_unit[3] - 1
    coupling_start = sum(nrates_per_unit) + n_noise + 1
    if rtarget === nothing
        u1r = Float64[0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2]
        u1n = Float64[0.0, 0.1, 0.5, 0.15]
        u2r = Float64[0.03, 0.1, 0.5, 0.2, 0.25, 0.17, 0.2, 0.6, 0.2]
        u2n = Float64[0.0, 0.1, 0.9, 0.2]
        u3r = Float64[0.12, 0.08, 0.15, 0.1, 0.1, 0.1, 0.1, 0.2]
        rtarget = vcat(u1r, u1n, u2r, u2n, u3r, fill(-0.2, ncpl))
    else
        length(rtarget) == n_total || throw(ArgumentError("rtarget must have length $n_total (sum(num_rates)+noise+coupling)"))
    end
    if fixedeffects === nothing
        fixedeffects = (collect(unit3_rate_start:unit3_rate_end), collect(coupling_start:coupling_start + ncpl - 1))
    end
    if fittedparam === nothing
        fittedparam = Int[unit3_rate_start, coupling_start]
    end
    trace_specs_eff = if trace_specs === nothing
        zm = zeromedian isa Bool ? fill(zeromedian, length(units)) : zeromedian
        StochasticGene.default_trace_specs_for_coupled((interval, 1.0, -1.0), zm, units)
    else
        trace_specs
    end
    units_eff = [spec.unit for spec in trace_specs_eff]
    trace = simulate_trace_vector(rtarget, transitions_eff, G, R, S, insertstep, coupling,
                                  interval, totaltime, ntrials; observed_units=units_eff, noiseparams=[4, 4, 0])
    data = TraceData("tracejoint", "test", interval, (trace, [], fill(0.0, length(units_eff)), 1), units_eff)
    priormean = set_priormean([], transitions_eff, R, S, insertstep, 1.0, noisepriors,
                              mean_elongationtime(rtarget, transitions_eff, R), tuple(), coupling, nothing)
    rinit = rtarget
    block_ends = cumsum([nrates_per_unit[i] + length(noisepriors[i]) for i in eachindex(R)])
    block_starts = vcat(1, block_ends[1:end-1] .+ 1)
    decayrate = tuple((rtarget[block_starts[i] + nrates_per_unit[i] - 1] for i in eachindex(R))...)
    model = load_model(data, rinit, priormean, fittedparam, fixedeffects, transitions_eff, G, R, S, insertstep,
                       "", 1, 10.0, Int[], decayrate, propcv, prob_Gaussian, noisepriors, method, tuple(), coupling, nothing,
                       zeromedian; coupled_stack=coupled_stack)
    meta = (scenario=:trace_coupled_3x3_g3r0,
            true_rates=rtarget,
            fittedparam=collect(Int, fittedparam),
            trace_specs=trace_specs_eff,
            coupling=coupling,
            coupled_stack=coupled_stack)
    return (data=data, model=model, meta=meta)
end

function _benchmark_ess_summary(measures)
    ess = first(measures.ess)
    return (minimum(ess), median(ess), maximum(ess))
end

function benchmark_inference_run_mh(scenario;
                                    samplesteps::Int=300,
                                    warmupsteps::Int=0,
                                    maxtime::Float64=60.0,
                                    temp::Float64=1.0,
                                    verbose::Bool=true)
    data, model = scenario.data, scenario.model
    options = MHOptions(samplesteps, warmupsteps, maxtime, temp)
    t = @elapsed fits, stats, measures = run_mh(data, model, options)
    summary = (method=:MH,
               walltime=t,
               meanparam=stats.meanparam,
               loglik_max=maximum(fits.ll),
               waic=measures.waic,
               ess=measures.ess,
               rhat=measures.rhat)
    verbose && println("MH: walltime=$(round(t; digits=3)) s, loglik_max=$(round(summary.loglik_max; digits=3))")
    return merge(summary, (fits=fits, stats=stats, measures=measures))
end

"""
    benchmark_combined_likelihood_stack(data, model; kwargs...)

Micro-benchmark the independent [`CombinedData`](@ref) likelihood stack. Returns timings for
the primal and AD-friendly paths, plus optional legacy `TraceRNAData` parity/timing when a
legacy payload is supplied.
"""
function benchmark_combined_likelihood_stack(
    data::CombinedData,
    model::AbstractGeneTransitionModel;
    param=get_param(model),
    mh_options::MHOptions=MHOptions(1, 0, 1.0, 1.0),
    ad_options::NUTSOptions=NUTSOptions(),
    nruns::Int=10,
    steady_state_solver::Symbol=:augmented,
    verbose::Bool=true,
)
    nruns >= 1 || throw(ArgumentError("nruns must be ≥ 1, got $nruns"))
    loglikelihood(param, data, model, mh_options; steady_state_solver=steady_state_solver)
    loglikelihood_ad(param, data, model, ad_options; steady_state_solver=steady_state_solver)
    t_primal = @elapsed for _ in 1:nruns
        loglikelihood(param, data, model, mh_options; steady_state_solver=steady_state_solver)
    end
    t_ad = @elapsed for _ in 1:nruns
        loglikelihood_ad(param, data, model, ad_options; steady_state_solver=steady_state_solver)
    end
    result = (
        modalities=combined_modalities(data),
        nruns=nruns,
        primal_sec=t_primal / nruns,
        ad_sec=t_ad / nruns,
        ad_over_primal=(t_ad / nruns) / (t_primal / nruns),
    )
    if verbose
        println("Combined likelihood modalities=$(result.modalities)")
        println("  primal_sec=$(round(result.primal_sec; digits=6))")
        println("  ad_sec=$(round(result.ad_sec; digits=6))")
        println("  ad_over_primal=$(round(result.ad_over_primal; digits=3))")
    end
    return result
end

function benchmark_combined_likelihood_stack(
    data::TraceRNAData,
    model::AbstractGRSMmodel;
    param=get_param(model),
    mh_options::MHOptions=MHOptions(1, 0, 1.0, 1.0),
    ad_options::NUTSOptions=NUTSOptions(),
    nruns::Int=10,
    steady_state_solver::Symbol=:augmented,
    verbose::Bool=true,
)
    combined = CombinedData(data)
    base = benchmark_combined_likelihood_stack(
        combined, model;
        param=param,
        mh_options=mh_options,
        ad_options=ad_options,
        nruns=nruns,
        steady_state_solver=steady_state_solver,
        verbose=false,
    )
    loglikelihood(param, data, model; steady_state_solver=steady_state_solver, hmm_stack=mh_options.likelihood_executor)
    t_legacy = @elapsed for _ in 1:nruns
        loglikelihood(param, data, model; steady_state_solver=steady_state_solver, hmm_stack=mh_options.likelihood_executor)
    end
    ll_combined, pred_combined = loglikelihood(param, combined, model, mh_options; steady_state_solver=steady_state_solver)
    ll_legacy, pred_legacy = loglikelihood(param, data, model; steady_state_solver=steady_state_solver, hmm_stack=mh_options.likelihood_executor)
    result = merge(base, (
        legacy_sec=t_legacy / nruns,
        combined_over_legacy=base.primal_sec / (t_legacy / nruns),
        ll_matches_legacy=isapprox(ll_combined, ll_legacy; rtol=1e-8, atol=1e-8),
        logprediction_length_matches=length(pred_combined) == length(pred_legacy),
    ))
    if verbose
        println("Combined likelihood modalities=$(result.modalities)")
        println("  primal_sec=$(round(result.primal_sec; digits=6))")
        println("  ad_sec=$(round(result.ad_sec; digits=6))")
        println("  legacy_sec=$(round(result.legacy_sec; digits=6))")
        println("  combined_over_legacy=$(round(result.combined_over_legacy; digits=3))")
        println("  ll_matches_legacy=$(result.ll_matches_legacy), logprediction_length_matches=$(result.logprediction_length_matches)")
    end
    return result
end

function benchmark_inference_run_nuts_parallel(
    scenario;
    nchains::Int=4,
    n_samples::Int=300,
    n_adapts::Int=300,
    δ::Float64=0.8,
    max_depth::Int=10,
    nuts_gradient::Symbol=:finite,
    nuts_fd_ε::Float64=1e-4,
    rng=MersenneTwister(42),
    verbose::Bool=true,
    progress::Bool=false,
    report_chains::Bool=true,
    parallel::Bool=true,
)
    data, model = scenario.data, scenario.model
    if nchains < 1
        throw(ArgumentError("nchains must be ≥ 1, got $nchains"))
    end

    base_seed = rand(rng, UInt)
    chain_seeds = [xor(base_seed, UInt(i)) for i in 1:nchains]

    function run_one_chain(chain_id::Int)
        chain_rng = MersenneTwister(chain_seeds[chain_id])
        chain_verbose = verbose && report_chains
        opts = NUTSOptions(
            ;
            n_samples=n_samples,
            n_adapts=n_adapts,
            δ=δ,
            max_depth=max_depth,
            gradient=nuts_gradient,
            fd_ε=nuts_fd_ε,
            verbose=chain_verbose,
            progress=progress,
        )
        wall = @elapsed begin
            fits, stats, measures, nuts_info = run_nuts_fit(data, model, opts; rng=chain_rng)
            return (fits=fits, stats=stats, measures=measures, nuts_info=nuts_info)
        end
        return (walltime=wall,)
    end

    results = [begin
        chain_rng = MersenneTwister(chain_seeds[chain_id])
        chain_verbose = verbose && report_chains
        opts = NUTSOptions(
            ;
            n_samples=n_samples,
            n_adapts=n_adapts,
            δ=δ,
            max_depth=max_depth,
            gradient=nuts_gradient,
            fd_ε=nuts_fd_ε,
            verbose=chain_verbose,
            progress=progress,
        )
        local fits, stats, measures, nuts_info
        wall = @elapsed begin
            fits, stats, measures, nuts_info = run_nuts_fit(data, model, opts; rng=chain_rng)
        end
        (fits=fits, stats=stats, measures=measures, nuts_info=nuts_info, walltime=wall)
    end for chain_id in 1:nchains]

    fits_list = getfield.(results, :fits)
    walltimes = getfield.(results, :walltime)
    fit_arrays = [res.fits.param for res in results]
    ll_arrays = [reshape(res.fits.ll, 1, :) for res in results]
    combined_param = hcat(fit_arrays...)
    combined_ll = vec(hcat(ll_arrays...))
    combined_fits = fits_list[1]
    combined_fits = (;
        combined_fits...,
        param=combined_param,
        ll=combined_ll,
        total=sum(getfield.(fits_list, :total)),
        time=sum(getfield.(fits_list, :time)),
    )

    meanparam = vec(mean(combined_param; dims=2))
    medianparam = vec(median(combined_param; dims=2))
    stdparam = vec(std(combined_param; dims=2, corrected=true))
    param_quantiles = mapslices(x -> quantile(x, [0.025, 0.25, 0.5, 0.75, 0.975]), combined_param; dims=2)
    stats = (meanparam=meanparam, medianparam=medianparam, stdparam=stdparam, param_quantiles=param_quantiles)
    waic, waic_se = measure_waic(combined_fits.ll)
    ess = measure_ess(combined_param)
    rhat = if nchains > 1
        nch = length(fit_arrays)
        ns = size(fit_arrays[1], 2)
        p = size(fit_arrays[1], 1)
        chains3 = Array{Float64}(undef, p, ns, nch)
        for (j, arr) in enumerate(fit_arrays)
            chains3[:, :, j] = arr
        end
        measure_rhat(chains3)
    else
        fill(NaN, size(combined_param, 1))
    end
    measures = (waic=(waic, waic_se), ess=(ess,), rhat=(rhat,))

    all_step_sizes = [res.nuts_info[:step_size] for res in results if haskey(res.nuts_info, :step_size)]
    all_accepts = vcat([res.nuts_info[:acceptance] for res in results if haskey(res.nuts_info, :acceptance)]...)
    all_tree_depths = vcat([res.nuts_info[:tree_depth] for res in results if haskey(res.nuts_info, :tree_depth)]...)
    divergences = sum(sum(res.nuts_info[:divergent]) for res in results if haskey(res.nuts_info, :divergent))
    nuts_info = Dict(
        :nchains => nchains,
        :walltimes => walltimes,
        :step_size => all_step_sizes,
        :acceptance => all_accepts,
        :tree_depth => all_tree_depths,
        :divergent_total => divergences,
        :chain_infos => [res.nuts_info for res in results],
    )

    wall_summary = (sum=sum(walltimes), max=maximum(walltimes), mean=mean(walltimes))
    loglik_max = maximum(combined_ll)
    verbose && println("NUTS ($nchains chains, parallel=$parallel): walltime(sum)=$(round(wall_summary.sum; digits=3)) s, walltime(max)=$(round(wall_summary.max; digits=3)) s, loglik_max=$(round(loglik_max; digits=3))")
    return (method=:NUTS,
            walltime=wall_summary,
            meanparam=meanparam,
            loglik_max=loglik_max,
            waic=measures.waic,
            ess=measures.ess,
            rhat=measures.rhat,
            fits=combined_fits,
            stats=stats,
            measures=measures,
            nuts_info=nuts_info)
end

function benchmark_inference_run_advi(scenario;
                                      advi_options::ADVIOptions=ADVIOptions(maxiter=400,
                                                                            n_mc=8,
                                                                            σ_floor=1e-4,
                                                                            init_s_raw=-4.0,
                                                                            verbose=false,
                                                                            gradient=:Zygote,
                                                                            time_limit=60.0),
                                      rng=MersenneTwister(42),
                                      verbose::Bool=true)
    data, model = scenario.data, scenario.model
    gradient = advi_options.gradient
    if data isa AbstractTraceData && gradient == :Zygote
        verbose && @info "Switching trace ADVI benchmark from gradient=:Zygote to :finite for robustness on the trace likelihood path."
        advi_options = ADVIOptions(;
            maxiter=advi_options.maxiter,
            n_mc=advi_options.n_mc,
            σ_floor=advi_options.σ_floor,
            init_s_raw=advi_options.init_s_raw,
            verbose=advi_options.verbose,
            gradient=:finite,
            time_limit=advi_options.time_limit,
            device=advi_options.device,
            parallel=advi_options.parallel,
        )
    end
    t = @elapsed fits, stats, measures, advi_info = run_advi_fit(data, model, advi_options; rng=rng)
    μ = advi_info[:μ]
    μ_native = inverse_transform_params(reshape(μ, :, 1), model)[:, 1]
    summary = (method=:ADVI,
               walltime=t,
               meanparam=stats.meanparam,
               loglik_max=maximum(fits.ll),
               waic=measures.waic,
               ess=measures.ess,
               rhat=measures.rhat,
               μ=μ_native,
               σ=advi_info[:σ],
               neg_elbo_min=minimum(advi_info[:vb_elbo]))
    verbose && println("ADVI: walltime=$(round(t; digits=3)) s, neg_elbo_min=$(round(summary.neg_elbo_min; digits=3))")
    return merge(summary, (fits=fits, stats=stats, measures=measures, advi_info=advi_info))
end

"""
    benchmark_inference_compare_mh_nuts_advi(; scenario=:simrna_small, kwargs...)

Run the same inference problem with MH, NUTS, and ADVI and return a consolidated report.

# Keyword `scenario`
- `:simrna_small` → synthetic RNA histogram benchmark
- `:trace_gr2r2` → synthetic single-trace benchmark
- `:trace_coupled_3x3` → synthetic two-unit coupled-trace benchmark
- `:trace_coupled_3x3_g3r0` → synthetic two-unit coupled-trace benchmark without reporter progression

Additional keywords are forwarded to the scenario builder and individual method runners.
"""
function benchmark_inference_compare_mh_nuts_advi(; scenario::Symbol=:simrna_small,
                                                  seed::Int=42,
                                                  verbose::Bool=true,
                                                  mh_kwargs::NamedTuple=NamedTuple(),
                                                  nuts_kwargs::NamedTuple=NamedTuple(),
                                                  advi_kwargs::NamedTuple=NamedTuple(),
                                                  scenario_kwargs...)
    scen = if scenario == :simrna_small
        benchmark_inference_simrna_small(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_gr2r2
        benchmark_inference_trace_gr2r2(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_coupled_3x3
        benchmark_inference_trace_coupled_3x3(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_coupled_3x3_g3r0
        benchmark_inference_trace_coupled_3x3_g3r0(; seed=seed, scenario_kwargs...)
    else
        throw(ArgumentError("Unknown benchmark scenario: $scenario"))
    end

    mh = benchmark_inference_run_mh(scen; verbose=verbose, mh_kwargs...)
    nuts = benchmark_inference_run_nuts_parallel(scen; verbose=verbose, nuts_kwargs...)
    advi = benchmark_inference_run_advi(scen; verbose=verbose, advi_kwargs...)

    true_rates = get(scen.meta, :true_rates, nothing)
    report = (scenario=scen.meta,
              mh=mh,
              nuts=nuts,
              advi=advi,
              true_rates=true_rates)
    verbose && benchmark_inference_print_summary(report)
    return report
end

function _benchmark_trace_longest_timesteps(data::TraceData)
    traces = first(data.data)
    longest = 0
    for trial in traces
        if trial isa Tuple
            for unit_trace in trial
                longest = max(longest, length(unit_trace))
            end
        else
            longest = max(longest, length(trial))
        end
    end
    return longest
end

"""
    benchmark_trace_forwarddiff_gradient(scenario; param_indices=collect(1:length(get_param(scenario.model))),
                                         nruns=3, warmup=true, steady_state_solver=:augmented,
                                         hmm_checkpoint_steps=32, fd_ε=1e-4)

Benchmark full-gradient evaluation of the trace log-likelihood using `ForwardDiff.gradient`
for a selected subset of transformed parameters. Returns wall-time / allocation summaries and
last gradient values.
"""
function benchmark_trace_forwarddiff_gradient(scenario;
                                              param_indices::AbstractVector{<:Integer}=collect(1:length(get_param(scenario.model))),
                                              nruns::Int=3,
                                              warmup::Bool=true,
                                              steady_state_solver::Symbol=:augmented,
                                              hmm_checkpoint_steps::Union{Nothing,Integer}=32,
                                              fd_ε::Float64=1e-4)
    data, model = scenario.data, scenario.model
    data isa TraceData || throw(ArgumentError("benchmark_trace_forwarddiff_gradient requires TraceData scenario"))
    θ0 = Vector{Float64}(get_param(model))
    idx = collect(Int, param_indices)
    f = let θbase=copy(θ0), data=data, model=model, idx=idx, steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps
        function (θsub::AbstractVector)
            θ = copy(θbase)
            θ[idx] .= θsub
            loglikelihood_ad(θ, data, model;
                             mode=:ForwardDiff,
                             steady_state_solver=steady_state_solver,
                             hmm_checkpoint_steps=hmm_checkpoint_steps)[1]
        end
    end
    x0 = θ0[idx]
    warmup && ForwardDiff.gradient(f, x0)
    times = Float64[]
    allocs = Int[]
    g_last = nothing
    for _ in 1:nruns
        push!(allocs, @allocated g_last = ForwardDiff.gradient(f, x0))
        push!(times, @elapsed g_last = ForwardDiff.gradient(f, x0))
    end
    return (method=:ForwardDiff,
            param_indices=idx,
            nruns=nruns,
            mean_time=mean(times),
            median_time=median(times),
            min_time=minimum(times),
            mean_alloc=mean(allocs),
            median_alloc=median(allocs),
            grad=g_last)
end

"""
    benchmark_trace_finitediff_gradient(scenario; param_indices=collect(1:length(get_param(scenario.model))),
                                        nruns=3, warmup=true, steady_state_solver=:augmented,
                                        hmm_checkpoint_steps=32, fd_ε=1e-4)

Benchmark finite-difference gradient evaluation of the trace log-likelihood for a selected subset
of transformed parameters. Returns wall-time / allocation summaries and the last gradient vector.
"""
function benchmark_trace_finitediff_gradient(scenario;
                                             param_indices::AbstractVector{<:Integer}=collect(1:length(get_param(scenario.model))),
                                             nruns::Int=3,
                                             warmup::Bool=true,
                                             steady_state_solver::Symbol=:augmented,
                                             hmm_checkpoint_steps::Union{Nothing,Integer}=32,
                                             fd_ε::Float64=1e-4)
    data, model = scenario.data, scenario.model
    data isa TraceData || throw(ArgumentError("benchmark_trace_finitediff_gradient requires TraceData scenario"))
    θ0 = Vector{Float64}(get_param(model))
    idx = collect(Int, param_indices)
    f = let θbase=copy(θ0), data=data, model=model, idx=idx, steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps
        function (θsub::AbstractVector)
            θ = copy(θbase)
            θ[idx] .= θsub
            loglikelihood(θ, data, model;
                          rates_fn=get_rates,
                          steady_state_solver=steady_state_solver,
                          hmm_checkpoint_steps=hmm_checkpoint_steps)[1]
        end
    end
    x0 = θ0[idx]
    warmup && FiniteDiff.finite_difference_gradient(f, x0; absstep=fd_ε)
    times = Float64[]
    allocs = Int[]
    g_last = nothing
    for _ in 1:nruns
        push!(allocs, @allocated g_last = FiniteDiff.finite_difference_gradient(f, x0; absstep=fd_ε))
        push!(times, @elapsed g_last = FiniteDiff.finite_difference_gradient(f, x0; absstep=fd_ε))
    end
    return (method=:FiniteDiff,
            param_indices=idx,
            nruns=nruns,
            fd_ε=fd_ε,
            mean_time=mean(times),
            median_time=median(times),
            min_time=minimum(times),
            mean_alloc=mean(allocs),
            median_alloc=median(allocs),
            grad=g_last)
end

"""
    benchmark_trace_zygote_subset_gradient(scenario;
        param_indices=collect(1:min(length(get_param(scenario.model)), 3)),
        nruns=3, warmup=true, steady_state_solver=:augmented,
        hmm_checkpoint_steps=8, fd_ε=1e-4)

Benchmark `Zygote.gradient` on a selected subset of transformed parameters for a trace scenario.
The default subset keeps the benchmark practical for the current coupled trace likelihood path.

Keyword `fd_ε` is accepted for API compatibility with the finite-difference helper; it is not used.
"""
function benchmark_trace_zygote_subset_gradient(scenario;
                                                param_indices::AbstractVector{<:Integer}=collect(1:min(length(get_param(scenario.model)), 3)),
                                                nruns::Int=3,
                                                warmup::Bool=true,
                                                steady_state_solver::Symbol=:augmented,
                                                hmm_checkpoint_steps::Union{Nothing,Integer}=8,
                                                fd_ε::Float64=1e-4)
    data, model = scenario.data, scenario.model
    data isa TraceData || throw(ArgumentError("benchmark_trace_zygote_subset_gradient requires TraceData scenario"))
    θ0 = Vector{Float64}(get_param(model))
    idx = collect(Int, param_indices)
    f = let θbase=copy(θ0), data=data, model=model, idx=idx, steady_state_solver=steady_state_solver, hmm_checkpoint_steps=hmm_checkpoint_steps
        function (θsub::AbstractVector)
            θ = copy(θbase)
            θ[idx] .= θsub
            loglikelihood_ad(θ, data, model;
                             mode=:Zygote,
                             steady_state_solver=steady_state_solver,
                             hmm_checkpoint_steps=hmm_checkpoint_steps)[1]
        end
    end
    x0 = θ0[idx]
    warmup && Zygote.gradient(f, x0)
    times = Float64[]
    allocs = Int[]
    g_last = nothing
    for _ in 1:nruns
        push!(allocs, @allocated g_last = Zygote.gradient(f, x0)[1])
        push!(times, @elapsed g_last = Zygote.gradient(f, x0)[1])
    end
    return (method=:Zygote,
            param_indices=idx,
            nruns=nruns,
            mean_time=mean(times),
            median_time=median(times),
            min_time=minimum(times),
            mean_alloc=mean(allocs),
            median_alloc=median(allocs),
            grad=g_last)
end

"""
    compare_trace_subset_gradient_benchmarks(; scenario=:trace_gr2r2, param_indices=[1], nruns=3,
                                             warmup=true, steady_state_solver=:augmented,
                                             hmm_checkpoint_steps=8, fd_ε=1e-4,
                                             seed=42, scenario_kwargs...)

Run ForwardDiff, finite-difference, and Zygote subset gradient benchmarks on the same trace
scenario and print a concise timing / allocation comparison table.
"""
function compare_trace_subset_gradient_benchmarks(; scenario::Symbol=:trace_gr2r2,
                                                  param_indices::AbstractVector{<:Integer}=Int[1],
                                                  nruns::Int=3,
                                                  warmup::Bool=true,
                                                  steady_state_solver::Symbol=:augmented,
                                                  hmm_checkpoint_steps::Union{Nothing,Integer}=8,
                                                  fd_ε::Float64=1e-4,
                                                  seed::Int=42,
                                                  scenario_kwargs...)
    scen = if scenario == :trace_gr2r2
        benchmark_inference_trace_gr2r2(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_coupled_3x3
        benchmark_inference_trace_coupled_3x3(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_coupled_3x3_g3r0
        benchmark_inference_trace_coupled_3x3_g3r0(; seed=seed, scenario_kwargs...)
    else
        throw(ArgumentError("Unknown trace benchmark scenario: $scenario"))
    end

    fd = benchmark_trace_forwarddiff_gradient(
        scen;
        param_indices=param_indices,
        nruns=nruns,
        warmup=warmup,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )
    finite = benchmark_trace_finitediff_gradient(
        scen;
        param_indices=param_indices,
        nruns=nruns,
        warmup=warmup,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
        fd_ε=fd_ε,
    )
    zyg = benchmark_trace_zygote_subset_gradient(
        scen;
        param_indices=param_indices,
        nruns=nruns,
        warmup=warmup,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )

    println("Trace subset gradient benchmark ($(scen.meta.scenario), params=$(collect(param_indices))):")
    println("  ForwardDiff : mean=$(round(fd.mean_time; digits=4)) s, alloc=$(round(fd.mean_alloc; digits=0)) bytes")
    println("  FiniteDiff  : mean=$(round(finite.mean_time; digits=4)) s, alloc=$(round(finite.mean_alloc; digits=0)) bytes")
    println("  Zygote      : mean=$(round(zyg.mean_time; digits=4)) s, alloc=$(round(zyg.mean_alloc; digits=0)) bytes")

    return (scenario=scen.meta, forwarddiff=fd, finitediff=finite, zygote=zyg)
end

"""
    benchmark_trace_compare_forwarddiff_vs_finitediff(; scenario=:trace_gr2r2, param_indices=[1],
                                                      nruns=3, warmup=true,
                                                      steady_state_solver=:augmented,
                                                      hmm_checkpoint_steps=8,
                                                      fd_ε=1e-4, seed=42, scenario_kwargs...)

Convenience wrapper that benchmarks ForwardDiff and finite-difference gradients on the same
trace scenario and reports their timing, allocation, and gradient agreement.
"""
function benchmark_trace_compare_forwarddiff_vs_finitediff(; scenario::Symbol=:trace_gr2r2,
                                                           param_indices::AbstractVector{<:Integer}=Int[1],
                                                           nruns::Int=3,
                                                           warmup::Bool=true,
                                                           steady_state_solver::Symbol=:augmented,
                                                           hmm_checkpoint_steps::Union{Nothing,Integer}=8,
                                                           fd_ε::Float64=1e-4,
                                                           seed::Int=42,
                                                           scenario_kwargs...)
    scen = if scenario == :trace_gr2r2
        benchmark_inference_trace_gr2r2(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_coupled_3x3
        benchmark_inference_trace_coupled_3x3(; seed=seed, scenario_kwargs...)
    elseif scenario == :trace_coupled_3x3_g3r0
        benchmark_inference_trace_coupled_3x3_g3r0(; seed=seed, scenario_kwargs...)
    else
        throw(ArgumentError("Unknown trace benchmark scenario: $scenario"))
    end

    fd = benchmark_trace_forwarddiff_gradient(
        scen;
        param_indices=param_indices,
        nruns=nruns,
        warmup=warmup,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )
    finite = benchmark_trace_finitediff_gradient(
        scen;
        param_indices=param_indices,
        nruns=nruns,
        warmup=warmup,
        steady_state_solver=steady_state_solver,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
        fd_ε=fd_ε,
    )

    grad_diff = fd.grad .- finite.grad
    println("ForwardDiff vs FiniteDiff ($(scen.meta.scenario), params=$(collect(param_indices))):")
    println("  ForwardDiff : mean=$(round(fd.mean_time; digits=4)) s, alloc=$(round(fd.mean_alloc; digits=0)) bytes")
    println("  FiniteDiff  : mean=$(round(finite.mean_time; digits=4)) s, alloc=$(round(finite.mean_alloc; digits=0)) bytes")
    println("  max|Δgrad|  = $(maximum(abs.(grad_diff)))")

    return (scenario=scen.meta, forwarddiff=fd, finitediff=finite, grad_diff=grad_diff)
end

function benchmark_inference_print_summary(report)
    println("\nInference benchmark summary ($(report.scenario.scenario))")
    println("  MH   : wall=$(round(report.mh.walltime; digits=3)) s, loglik_max=$(round(report.mh.loglik_max; digits=3))")
    println("         WAIC=$(round(report.mh.waic[1]; digits=3)) ± $(round(report.mh.waic[2]; digits=3))")
    mh_ess = _benchmark_ess_summary(report.mh.measures)
    println("         ESS(min/med/max)=$(round(mh_ess[1]; digits=2)) / $(round(mh_ess[2]; digits=2)) / $(round(mh_ess[3]; digits=2))")

    nuts_wall = report.nuts.walltime
    println("  NUTS : wall(sum)=$(round(nuts_wall.sum; digits=3)) s, wall(max)=$(round(nuts_wall.max; digits=3)) s")
    println("         loglik_max=$(round(report.nuts.loglik_max; digits=3)), divergences=$(report.nuts.nuts_info[:divergent_total])")
    println("         WAIC=$(round(report.nuts.waic[1]; digits=3)) ± $(round(report.nuts.waic[2]; digits=3))")
    nuts_ess = _benchmark_ess_summary(report.nuts.measures)
    println("         ESS(min/med/max)=$(round(nuts_ess[1]; digits=2)) / $(round(nuts_ess[2]; digits=2)) / $(round(nuts_ess[3]; digits=2))")

    println("  ADVI : wall=$(round(report.advi.walltime; digits=3)) s, loglik_max=$(round(report.advi.loglik_max; digits=3))")
    println("         neg_ELBO_min=$(round(report.advi.neg_elbo_min; digits=3))")
    println("         WAIC=$(round(report.advi.waic[1]; digits=3)) ± $(round(report.advi.waic[2]; digits=3))")

    if report.true_rates !== nothing
        println("  True rates        : $(report.true_rates)")
        println("  MH mean estimate  : $(round.(report.mh.stats.meanparam; digits=4))")
        println("  NUTS mean estimate: $(round.(report.nuts.stats.meanparam; digits=4))")
        println("  ADVI μ estimate   : $(round.(report.advi.μ; digits=4))")
    end
    nothing
end
