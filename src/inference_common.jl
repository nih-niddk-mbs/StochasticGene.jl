
# Unified inference entry point for StochasticGene.jl

"""
    run_inference(data, model, options; rng=Random.default_rng(), nchains::Integer=1)

Unified inference entry point. Dispatches on `options` type (`MHOptions`, `NUTSOptions`, `ADVIOptions`).

Returns `(fits, stats, measures)` in the same convention as [`run_mh`](@ref).

`nchains` is the number of independent chains (MH: uses existing distributed chain runner when `nchains > 1`;
NUTS/ADVI: uses `options.parallel` together with `nchains` to run chains in parallel when requested).
"""
function run_inference(data, model, options; rng=Random.default_rng(), nchains::Integer=1)
    return run_inference(data, model, options, rng, Int(nchains))
end

function run_inference(data, model, options, rng, nchains::Int)
    error("No run_inference method for model=$(typeof(model)), data=$(typeof(data)), options=$(typeof(options))")
end

function run_inference(data, model, options::MHOptions, rng, nchains::Int)
    return _run_mh_inference(data, model, options, rng, nchains)
end

function run_inference(data, model, options::NUTSOptions, rng, nchains::Int)
    return _run_nuts_inference(data, model, options, rng, nchains)
end

function run_inference(data, model, options::ADVIOptions, rng, nchains::Int)
    return _run_advi_inference(data, model, options, rng, nchains)
end

function _run_mh_inference(data, model, options::MHOptions, rng, nchains::Int)
    if options.device === :gpu
        return run_mh_gpu(data, model, options, nchains)
    elseif nchains > 1
        return run_mh(data, model, options, nchains)
    else
        return run_mh(data, model, options)
    end
end

function _merge_waic_chain(fitsv::Vector{Fit}, data)
    n = length(fitsv)
    n >= 1 || throw(ArgumentError("need at least one chain to merge"))
    chain = Vector{Tuple{Fit,Tuple{Float64,Float64}}}(undef, n)
    for i in 1:n
        f = fitsv[i]
        chain[i] = (f, compute_waic(f.lppd, f.pwaic, data))
    end
    return chain
end

function _merge_nuts_chains(data, model, chain_rows::AbstractVector)
    fitsv = Fit[cr[1] for cr in chain_rows]
    chain = _merge_waic_chain(fitsv, data)
    waic = pooled_waic(chain)
    fits = merge_fit(fitsv)
    stats = compute_stats(fits.param, model)
    rhat = vec(compute_rhat(fitsv))
    ess, geweke, mcse = compute_measures(fitsv)
    return fits, stats, Measures(waic, rhat, ess, geweke, mcse)
end

function _merge_advi_chains(data, model, chain_rows::AbstractVector)
    # Same collation as NUTS: each "chain" is an independent ADVI run (optional parallelism).
    return _merge_nuts_chains(data, model, chain_rows)
end

function _nuts_options_progress(options::NUTSOptions, progress::Bool)
    return NUTSOptions(;
        n_samples=options.n_samples,
        n_adapts=options.n_adapts,
        δ=options.δ,
        gradient=options.gradient,
        fd_ε=options.fd_ε,
        verbose=options.verbose,
        progress=progress,
        device=options.device,
        parallel=options.parallel,
        likelihood_executor=options.likelihood_executor,
        gradient_checkpoint_length=options.gradient_checkpoint_length,
    )
end

function _run_nuts_inference(data, model, options::NUTSOptions, rng, nchains::Int)
    if options.device === :gpu
        error("NUTS inference on GPU is not implemented.")
    end
    # One ProgressMeter bar per process; disable when multiple chains interleave output.
    opts = nchains > 1 ? _nuts_options_progress(options, false) : options
    if nchains > 1 && options.parallel === :distributed
        @info "Running NUTS inference with Distributed parallelism ($(nchains) chains)."
        chain_rows = Distributed.pmap(1:nchains) do _
            run_nuts_fit(data, model, opts; rng=Random.default_rng())
        end
        return _merge_nuts_chains(data, model, chain_rows)
    elseif nchains > 1 && options.parallel === :threaded
        @info "Running NUTS inference with multithreading ($(nchains) chains)."
        seeds = Tuple{UInt64,UInt64,UInt64,UInt64}[
            (rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64)) for _ in 1:nchains
        ]
        chain_rows = Vector{Any}(undef, nchains)
        Threads.@threads for i in 1:nchains
            s = seeds[i]
            chain_rows[i] = run_nuts_fit(data, model, opts; rng=Random.Xoshiro(s[1], s[2], s[3], s[4]))
        end
        return _merge_nuts_chains(data, model, chain_rows)
    elseif nchains > 1
        @info "Running NUTS inference serially ($(nchains) chains); set options.parallel=:threaded or :distributed for parallel chains."
        chain_rows = Vector{Any}(undef, nchains)
        for i in 1:nchains
            chain_rows[i] = run_nuts_fit(data, model, opts; rng=Random.Xoshiro(rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64)))
        end
        return _merge_nuts_chains(data, model, chain_rows)
    else
        fits, stats, measures, _ = run_nuts_fit(data, model, opts; rng=rng)
        return fits, stats, measures
    end
end

function _run_advi_inference(data, model, options::ADVIOptions, rng, nchains::Int)
    if options.device === :gpu
        error("ADVI inference on GPU is not implemented.")
    elseif nchains > 1 && options.parallel === :distributed
        @info "Running ADVI inference with Distributed parallelism ($(nchains) repeats)."
        chain_rows = Distributed.pmap(1:nchains) do _
            run_advi_fit(data, model, options; rng=Random.default_rng(), zygote_trace=options.zygote_trace)
        end
        return _merge_advi_chains(data, model, chain_rows)
    elseif nchains > 1 && options.parallel === :threaded
        @info "Running ADVI inference with multithreading ($(nchains) repeats)."
        seeds = Tuple{UInt64,UInt64,UInt64,UInt64}[
            (rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64)) for _ in 1:nchains
        ]
        chain_rows = Vector{Any}(undef, nchains)
        Threads.@threads for i in 1:nchains
            s = seeds[i]
            chain_rows[i] = run_advi_fit(data, model, options; rng=Random.Xoshiro(s[1], s[2], s[3], s[4]), zygote_trace=options.zygote_trace)
        end
        return _merge_advi_chains(data, model, chain_rows)
    elseif nchains > 1
        @info "Running ADVI inference serially ($(nchains) repeats); set options.parallel=:threaded or :distributed for parallel runs."
        chain_rows = Vector{Any}(undef, nchains)
        for i in 1:nchains
            chain_rows[i] = run_advi_fit(data, model, options; rng=Random.Xoshiro(rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64), rand(rng, UInt64)), zygote_trace=options.zygote_trace)
        end
        return _merge_advi_chains(data, model, chain_rows)
    else
        fits, stats, measures, _ = run_advi_fit(data, model, options; rng=rng, zygote_trace=options.zygote_trace)
        return fits, stats, measures
    end
end
