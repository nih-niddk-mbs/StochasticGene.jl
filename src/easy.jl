"""
    rna_histogram(counts) -> Vector{Int}

Convert per-cell RNA counts into a histogram; element `i` is the number of cells with `i-1` RNA.
"""
function rna_histogram(counts::AbstractVector{<:Integer})
    isempty(counts) && throw(ArgumentError("counts is empty"))
    any(<(0), counts) && throw(ArgumentError("counts must be non-negative"))
    h = zeros(Int, maximum(counts) + 1)
    for c in counts
        h[c+1] += 1
    end
    return h
end

"""
    fit_rna(data; G=2, R=0, S=0, insertstep=1, name="gene", halflife=nothing,
            outdir="stochasticgene_output", samplesteps=100_000, nchains=2, kwargs...)

Fit a gene-state (telegraph-like) model to steady-state RNA data with sensible defaults.

`data` is one of
- a histogram vector (element `i` = number of cells with `i-1` RNA); for raw per-cell counts
  pass `rna_histogram(counts)`,
- a path to a text/CSV file whose first column is the histogram.

`halflife` is the mRNA half-life in minutes; if omitted the decay rate is set to 1.0 per minute,
so rates are then only meaningful relative to decay.
Inputs are staged under `outdir/data` and results are written to `outdir/results/<name>`.
Extra `kwargs` are passed to [`fit`](@ref). Returns the `fit` outputs plus a `rates` table
and `outdir` as a `NamedTuple`; see [`describe_fit`](@ref).
"""
function fit_rna(data; G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, name::AbstractString="gene",
    halflife=nothing, outdir::AbstractString="stochasticgene_output", samplesteps::Int=100_000,
    nchains::Int=2, kwargs...)
    outdir = abspath(outdir)
    datadir = joinpath(outdir, "data")
    mkpath(datadir)
    mkpath(joinpath(outdir, "results"))
    cond = "DATA"
    target = joinpath(datadir, "$(name)_$(cond).txt")
    if data isa AbstractString
        isfile(data) || throw(ArgumentError("data file not found: $data"))
        abspath(data) == abspath(target) || cp(data, target; force=true)
    elseif data isa AbstractVector{<:Real}
        any(<(0), data) && throw(ArgumentError("histogram entries must be non-negative"))
        writedlm(target, data)
    else
        throw(ArgumentError("data must be a file path or a histogram vector"))
    end
    if halflife === nothing
        @info "No halflife given; decay rate set to 1.0 per minute"
        decayrate = 1.0
    else
        halflife > 0 || throw(ArgumentError("halflife must be positive (minutes)"))
        decayrate = log(2) / halflife
    end
    transitions = G == 1 ? tuple() : get_transitions(G, "")
    out = fit(; datatype="rna", datapath=datadir, gene=String(name), cell="cell", datacond=cond,
        resultfolder=String(name), root=outdir, G=G, R=R, S=S, insertstep=insertstep,
        transitions=transitions, decayrate=decayrate, samplesteps=samplesteps, nchains=nchains,
        burst=G > 1, kwargs...)
    fits, stats, measures, d, model, options = out
    resultdir = joinpath(outdir, "results", String(name))
    return (; fits, stats, measures, data=d, model, options, outdir, resultdir, rates=fit_summary(resultdir))
end

"""
    fit_summary(resultdir) -> DataFrame

Parameter table (maximum-likelihood sample, posterior mean and median; rates per minute)
from the most recent `rates_*.txt` file in `resultdir`.
"""
function fit_summary(resultdir::AbstractString)
    files = filter(f -> startswith(basename(f), "rates_") && endswith(f, ".txt"), readdir(resultdir; join=true))
    isempty(files) && throw(ArgumentError("no rates_*.txt files in $resultdir"))
    file = last(sort(files; by=mtime))
    m, header = readdlm(file, ','; header=true)
    return DataFrame(parameter=vec(string.(header)), ml=m[1, :], mean=m[2, :], median=m[3, :])
end

"""
    describe_fit(result)

Print the parameter table and output location for the result of [`fit_rna`](@ref)
or for a results directory.
"""
function describe_fit(result)
    resultdir = result isa AbstractString ? String(result) : result.resultdir
    df = fit_summary(resultdir)
    println("Results in ", resultdir)
    println("Rates are per minute; 'median' is the recommended summary.")
    show(stdout, MIME("text/plain"), df)
    println()
    bursts = filter(f -> startswith(basename(f), "burst_"), readdir(resultdir))
    isempty(bursts) || println("Burst size summary: ", joinpath(resultdir, first(bursts)))
    return df
end

"""
    demo(; outdir=mktempdir(), samplesteps=20_000)

Run a short two-state fit on the bundled MYC smFISH histogram and print the result.
Takes a minute or so, mostly compilation on first use.
"""
function demo(; outdir::AbstractString=mktempdir(), samplesteps::Int=20_000)
    file = joinpath(pkgdir(@__MODULE__), "data", "HCT116_testdata", "MYC_MOCK.txt")
    isfile(file) || throw(ArgumentError("bundled test data not found: $file"))
    result = fit_rna(file; G=2, name="MYC", halflife=600.0, outdir=outdir, samplesteps=samplesteps)
    describe_fit(result)
    return result
end
