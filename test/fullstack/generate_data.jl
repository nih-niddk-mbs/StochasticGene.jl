#!/usr/bin/env julia

# Developer-only fixture generator for full-stack fit checks.
#
# This script writes small deterministic data files for the current loaders. It is
# intentionally not included from test/runtests.jl because it creates files.

using DelimitedFiles
using Random

using StochasticGene

const FULLSTACK_GENES = ("GENE1", "GENE2")
const FULLSTACK_CELL = "CELL"
const FULLSTACK_COND = "COND"
const FULLSTACK_TRACEJOINT_CONDS = ["$(g)_$(FULLSTACK_COND)" for g in FULLSTACK_GENES]

const DEFAULT_ROOT = @__DIR__
const DEFAULT_DATA_ROOT = joinpath(DEFAULT_ROOT, "data")

const DEFAULT_TRANSITIONS = ([1, 2], [2, 1])
const DEFAULT_G = 2
const DEFAULT_R = 0
const DEFAULT_S = 0
const DEFAULT_INSERTSTEP = 0
const DEFAULT_NALLELES = 2
const DEFAULT_RATES = Float64[0.33, 0.19, 20.5, 1.0]

const COUPLED_TRANSITIONS = (([1, 2], [2, 1]), ([1, 2], [2, 1]))
const COUPLED_G = (2, 2)
const COUPLED_R = (1, 1)
const COUPLED_S = (0, 0)
const COUPLED_INSERTSTEP = (1, 1)
const COUPLED_COUPLING = ((1, 2), [(1, 2, 2, 1)], [:free])
const COUPLED_RATES = Float64[
    0.30, 0.15, 0.60, 1.0, 0.10, 0.0, 0.10, 1.0, 0.10,
    0.25, 0.18, 0.55, 1.0, 0.12, 0.0, 0.10, 1.0, 0.10,
    -0.20,
]
const COUPLED_UNIT_RNA_RATES = (
    COUPLED_RATES[1:5],
    COUPLED_RATES[10:14],
)

function fixture_path(parts...)
    return joinpath(DEFAULT_DATA_ROOT, parts...)
end

function reset_dir!(path::AbstractString)
    isdir(path) && rm(path; recursive=true, force=true)
    mkpath(path)
    return path
end

function write_vector(path::AbstractString, values)
    mkpath(dirname(path))
    writedlm(path, collect(values))
    return path
end

function write_matrix(path::AbstractString, values)
    mkpath(dirname(path))
    writedlm(path, values)
    return path
end

function _sample_counts(rng::AbstractRNG, hist::AbstractVector, ncells::Integer)
    p = Float64.(hist)
    s = sum(p)
    s > 0 || throw(ArgumentError("RNA histogram must have positive mass"))
    p ./= s
    cdf = cumsum(p)
    counts = Vector{Int}(undef, ncells)
    for i in 1:ncells
        counts[i] = searchsortedfirst(cdf, rand(rng)) - 1
    end
    return counts
end

"""
    generate_rna_fixture!(data_root=DEFAULT_DATA_ROOT; kwargs...)

Write `data/rna/GENE*_COND.txt`, one-column RNA histograms used by `read_rna`.
"""
function generate_rna_fixture!(
    data_root::AbstractString=DEFAULT_DATA_ROOT;
    seed::Integer=11,
    totalsteps::Integer=50_000,
    nRNA::Integer=80,
)
    Random.seed!(seed)
    outdir = reset_dir!(joinpath(data_root, "rna"))
    paths = String[]
    for (gene, rates) in zip(FULLSTACK_GENES, COUPLED_UNIT_RNA_RATES)
        hist = StochasticGene.simulator(
            collect(rates),
            DEFAULT_TRANSITIONS,
            DEFAULT_G,
            1,
            0,
            1;
            nhist=nRNA,
            totalsteps=totalsteps,
            nalleles=DEFAULT_NALLELES,
        )[1]
        hist = round.(Int, hist)
        push!(paths, write_vector(joinpath(outdir, "$(gene)_$(FULLSTACK_COND).txt"), hist))
    end
    return paths
end

"""
    generate_rnacount_fixture!(data_root=DEFAULT_DATA_ROOT; kwargs...)

Write `data/rnacount/GENE*_COND.txt`, two-column per-cell count/yield files used by
`read_rnacount`.
"""
function generate_rnacount_fixture!(
    data_root::AbstractString=DEFAULT_DATA_ROOT;
    seed::Integer=12,
    totalsteps::Integer=50_000,
    nRNA::Integer=80,
    ncells::Integer=1_000,
)
    rng = MersenneTwister(seed)
    outdir = reset_dir!(joinpath(data_root, "rnacount"))
    paths = String[]
    for (gene, rates) in zip(FULLSTACK_GENES, COUPLED_UNIT_RNA_RATES)
        hist = StochasticGene.simulator(
            collect(rates),
            DEFAULT_TRANSITIONS,
            DEFAULT_G,
            1,
            0,
            1;
            nhist=nRNA,
            totalsteps=totalsteps,
            nalleles=DEFAULT_NALLELES,
        )[1]
        counts = _sample_counts(rng, hist, ncells)
        yield = ones(Float64, ncells)
        push!(paths, write_matrix(joinpath(outdir, "$(gene)_$(FULLSTACK_COND).txt"), hcat(counts, yield)))
    end
    return paths
end

"""
    generate_trace_fixture!(data_root=DEFAULT_DATA_ROOT; kwargs...)

Write paired two-unit trace files into one trace folder. Single trace fits discover one
gene by `gene` + `COND`; tracejoint fits use labels `GENE1_COND` and `GENE2_COND`.
"""
function generate_trace_fixture!(
    data_root::AbstractString=DEFAULT_DATA_ROOT;
    seed::Integer=13,
    interval::Real=1.0,
    totaltime::Real=120.0,
    ntrials::Integer=12,
)
    Random.seed!(seed)
    traces = StochasticGene.simulate_trace_vector(
        COUPLED_RATES,
        COUPLED_TRANSITIONS,
        COUPLED_G,
        COUPLED_R,
        COUPLED_S,
        COUPLED_INSERTSTEP,
        COUPLED_COUPLING,
        Float64(interval),
        Float64(totaltime),
        Int(ntrials);
        noiseparams=[4, 4],
    )
    outdir = joinpath(data_root, "trace")
    reset_dir!(outdir)
    for i in 1:ntrials
        trace_matrix = traces[i]
        for (unit, gene) in enumerate(FULLSTACK_GENES)
            trace = trace_matrix[:, unit]
            frame = collect(1:length(trace))
            time = (frame .- 1) .* Float64(interval)
            write_matrix(joinpath(outdir, lpad(string(i), 3, "0") * "_$(gene)_$(FULLSTACK_COND).txt"), hcat(frame, time, trace))
        end
    end
    return outdir
end

"""
    generate_dwelltime_fixture!(data_root=DEFAULT_DATA_ROOT; kwargs...)

Write per-gene two-column dwell-time histogram files with bins and normalized values.
The current `DwellTimeData` loader receives a vector of paths and reads each as
`(bin, value)`.
"""
function generate_dwelltime_fixture!(
    data_root::AbstractString=DEFAULT_DATA_ROOT;
    seed::Integer=14,
    nbins::Integer=40,
    nsamples::Integer=2_000,
)
    rng = MersenneTwister(seed)
    outdir = reset_dir!(joinpath(data_root, "dwelltime"))
    bins = collect(1:nbins)
    paths = String[]
    for (i, gene) in enumerate(FULLSTACK_GENES)
        samples = clamp.(ceil.(Int, randexp(rng, nsamples) .* (7 + i)), 1, nbins)
        hist = [count(==(b), samples) for b in bins]
        values = Float64.(hist) ./ max(sum(hist), 1)
        push!(paths, write_matrix(joinpath(outdir, "$(gene)_$(FULLSTACK_COND)_OFF.txt"), hcat(bins, values)))
    end
    return paths
end

"""
    generate_fullstack_data!(data_root=DEFAULT_DATA_ROOT)

Regenerate all developer full-stack fixture data.
"""
function generate_fullstack_data!(data_root::AbstractString=DEFAULT_DATA_ROOT)
    reset_dir!(data_root)
    paths = (
        rna=generate_rna_fixture!(data_root),
        rnacount=generate_rnacount_fixture!(data_root),
        trace=generate_trace_fixture!(data_root),
        dwelltime=generate_dwelltime_fixture!(data_root),
    )
    return paths
end

if abspath(PROGRAM_FILE) == @__FILE__
    paths = generate_fullstack_data!()
    println("Generated full-stack fixture data for genes=$(FULLSTACK_GENES), cell=$(FULLSTACK_CELL), datacond=$(FULLSTACK_COND):")
    println("  tracejoint datacond labels: ", FULLSTACK_TRACEJOINT_CONDS)
    for (name, path) in pairs(paths)
        println("  ", name, ": ", path)
    end
end
