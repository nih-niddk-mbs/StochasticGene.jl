function _decayrate(halflife)
    halflife === nothing && return 1.0
    halflife > 0 || throw(ArgumentError("halflife must be positive (minutes)"))
    return log(2) / halflife
end

# Stage traces as numbered `frame time intensity` files, the layout the trace loader reads.
function _stage_traces(data, tracedir, name, interval)
    mkpath(tracedir)
    for f in readdir(tracedir; join=true)
        startswith(basename(f), "0") && endswith(f, "_$(name)_DATA.txt") && rm(f)
    end
    tag = "_$(name)_DATA.txt"
    n = 0
    function target()
        n += 1
        return joinpath(tracedir, lpad(n, 4, '0') * tag)
    end
    if data isa AbstractString
        isdir(data) || throw(ArgumentError("trace folder not found: $data"))
        files = [joinpath(r, f) for (r, _, fs) in walkdir(data) for f in fs if !startswith(f, ".")]
        isempty(files) && throw(ArgumentError("no trace files in $data"))
        foreach(f -> cp(f, target(); force=true), sort(files))
    elseif data isa AbstractVector && all(x -> x isa AbstractString, data)
        isempty(data) && throw(ArgumentError("no trace files given"))
        for f in data
            isfile(f) || throw(ArgumentError("trace file not found: $f"))
            cp(f, target(); force=true)
        end
    elseif data isa AbstractVector && all(x -> x isa AbstractVector{<:Real}, data)
        isempty(data) && throw(ArgumentError("no traces given"))
        for tr in data
            all(isfinite, tr) || throw(ArgumentError("traces must be finite"))
            length(tr) >= 10 || throw(ArgumentError("each trace needs at least 10 frames"))
            frames = 1:length(tr)
            writedlm(target(), [collect(Float64, frames) (frames .- 1) .* interval collect(Float64, tr)])
        end
    else
        throw(ArgumentError("data must be a folder, a vector of trace files, or a vector of intensity vectors"))
    end
    return n
end

"""
    fit_traces(data; G=2, R=1, S=0, insertstep=1, name="gene", interval=1.0, halflife=nothing,
               noisepriors=[0.0, 0.1, 1.0, 0.1], zeromedian=true, datacol=3,
               outdir="stochasticgene_output", samplesteps=40_000, nchains=2, maxtime=300.0, kwargs...)

Fit a gene-state model with a reporter (MS2-style) to live-cell intensity traces, with defaults.

`data` is one of
- a folder of trace files,
- a vector of trace file paths,
- a vector of intensity vectors (one per trace, equal frame spacing).

Trace files are text/CSV with one frame per row; intensity is read from column `datacol`
(default 3: frame, time, intensity). `interval` is the frame interval in minutes. `halflife` is
the mRNA half-life in minutes (omitted: decay rate 1.0 per minute). `R` must be at least 1.
With `zeromedian=true` each trace is median-subtracted and scaled, so the default `noisepriors`
(background mean and SD, signal mean and SD) apply to that normalized scale; set
`zeromedian=false` and `noisepriors` for raw intensities. `warmupsteps` defaults to `samplesteps`;
trace fits mix poorly without a long adaptive warmup.

Inputs are copied to `outdir/data/traces`, results go to `outdir/results/<name>`. Extra `kwargs`
go to [`fit`](@ref). Returns the `fit` outputs plus `rates`, `outdir` and `resultdir`;
see [`describe_fit`](@ref).
"""
function fit_traces(data; G::Int=2, R::Int=1, S::Int=0, insertstep::Int=1, name::AbstractString="gene",
    interval::Real=1.0, halflife=nothing, noisepriors=[0.0, 0.1, 1.0, 0.1], zeromedian::Bool=true,
    datacol::Int=3, outdir::AbstractString="stochasticgene_output", samplesteps::Int=40_000,
    warmupsteps::Int=samplesteps, nchains::Int=2, maxtime::Real=300.0, kwargs...)
    R >= 1 || throw(ArgumentError("trace fits need R >= 1 (a reporter); use fit_rna for histograms"))
    interval > 0 || throw(ArgumentError("interval must be positive (minutes)"))
    outdir = abspath(outdir)
    datadir = joinpath(outdir, "data", "traces")
    mkpath(joinpath(outdir, "results"))
    ntraces = _stage_traces(data, datadir, name, interval)
    halflife === nothing && @info "No halflife given; decay rate set to 1.0 per minute"
    cond = "$(name)_DATA"
    transitions = G == 1 ? tuple() : get_transitions(G, "")
    spec = [(unit=1, interval=Float64(interval), start=1.0, t_end=-1.0, zeromedian=zeromedian,
        active_fraction=1.0, background=0.0)]
    out = fit(; datatype="trace", datapath=datadir, gene=String(name), cell="cell", datacond=cond,
        resultfolder=String(name), root=outdir, G=G, R=R, S=S, insertstep=insertstep,
        transitions=transitions, decayrate=_decayrate(halflife), samplesteps=samplesteps,
        warmupsteps=warmupsteps, nchains=nchains, maxtime=Float64(maxtime), noisepriors=noisepriors,
        datacol=datacol, trace_specs=spec, kwargs...)
    fits, stats, measures, d, model, options = out
    resultdir = joinpath(outdir, "results", String(name))
    return (; fits, stats, measures, data=d, model, options, outdir, resultdir, ntraces, rates=fit_summary(resultdir))
end

"""
    demo_traces(; outdir=mktempdir(), ntraces=30, totaltime=300.0, samplesteps=10_000)

Simulate MS2-style traces from a known two-state model, fit them with [`fit_traces`](@ref) and
print fitted and true rates side by side. A short run is only a rough check of the workflow.
"""
function demo_traces(; outdir::AbstractString=mktempdir(), ntraces::Int=30, totaltime::Real=300.0,
    samplesteps::Int=30_000)
    true_rates = [0.1, 0.2, 0.5, 0.5, 0.05, 0.0, 0.1, 1.0, 0.1]  # k12, k21, initiate, eject, decay, noise(4)
    traces = simulate_trace_vector(true_rates, ([1, 2], [2, 1]), 2, 1, 0, 1, 1.0, totaltime, ntraces)
    result = fit_traces([vec(Float64.(t)) for t in traces]; G=2, R=1, name="demo", halflife=log(2) / true_rates[5],
        zeromedian=false, outdir=outdir, samplesteps=samplesteps)
    describe_fit(result)
    println("True rates (per minute): Rate_12=0.1, Rate_21=0.2, Initiate=0.5, Eject=0.5")
    println("A short run from few traces recovers Eject and the noise well; switching and initiation rates are weakly identified. Use more traces and longer chains for real data.")
    return result
end
