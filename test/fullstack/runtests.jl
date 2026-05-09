#!/usr/bin/env julia

# Developer full-stack tests. These are intentionally separate from
# test/runtests.jl because they regenerate fixture files and write fit outputs.

using Test

using StochasticGene

include(joinpath(@__DIR__, "run_fits.jl"))

const FULLSTACK_TEST_RESULTFOLDER = "runtests"
const FULLSTACK_EXTENDED = get(ENV, "SG_FULLSTACK_EXTENDED", "0") == "1"

function _all_files(paths)
    out = String[]
    for p in paths
        if p isa AbstractString
            if isdir(p)
                append!(out, joinpath(root, file) for (root, _, files) in walkdir(p) for file in files)
            else
                push!(out, p)
            end
        elseif p isa AbstractVector || p isa Tuple
            append!(out, _all_files(p))
        else
            append!(out, _all_files(collect(p)))
        end
    end
    return out
end

function _write_key_spec!(key::AbstractString, resultfolder::AbstractString; overrides...)
    outdir = StochasticGene.folder_path(resultfolder, DEFAULT_ROOT, "results", make=true)
    spec = Dict{Symbol,Any}(
        :root => DEFAULT_ROOT,
        :resultfolder => resultfolder,
        :cell => FULLSTACK_CELL,
        :datatype => (:rna,),
        :datapath => (rna="rna",),
        :gene => FULLSTACK_GENES[1],
        :datacond => FULLSTACK_COND,
        :label => "will-be-overridden-by-key",
        :writesamples => false,
        :optimize => false,
        :burst => false,
        :inference_method => :mh,
        :samplesteps => 2,
        :warmupsteps => 0,
        :maxtime => 30.0,
        :nchains => 1,
    )
    merge!(spec, Dict{Symbol,Any}(pairs(_rna_model_kwargs(1))))
    merge!(spec, Dict{Symbol,Any}(pairs(overrides)))
    StochasticGene.write_run_spec_jld2(joinpath(outdir, "info_" * key * ".toml"), spec)
    return spec
end

@testset "full-stack developer fixtures" begin
    paths = generate_fullstack_data!()
    files = _all_files(collect(paths))

    @test all(isfile, files)
    @test length(paths.rna) == length(FULLSTACK_GENES)
    @test length(paths.rnacount) == length(FULLSTACK_GENES)
    @test length(paths.dwelltime) == length(FULLSTACK_GENES)
    @test count(endswith(".txt"), readdir(paths.trace)) == length(FULLSTACK_GENES) * 12
    @test isfile(joinpath(DEFAULT_DATA_ROOT, "rna", "GENE1_COND.txt"))
    @test isfile(joinpath(DEFAULT_DATA_ROOT, "rnacount", "GENE1_COND.txt"))
    @test isfile(joinpath(DEFAULT_DATA_ROOT, "dwelltime", "GENE1_COND_OFF.txt"))
end

@testset "full-stack fit smoke" begin
    smoke = fullstack_smoke_defaults()
    outputs = run_fullstack_fits!(
        cases="rna,trace";
        regenerate=false,
        resultfolder=FULLSTACK_TEST_RESULTFOLDER,
        smoke...,
    )
    @test haskey(outputs, :rna)
    @test haskey(outputs, :trace)
    @test assert_fullstack_fit_output(:rna, outputs[:rna])
    @test assert_fullstack_fit_output(:trace, outputs[:trace])

    _, _, _, data_rna, _, options_rna = outputs[:rna]
    @test data_rna isa StochasticGene.CombinedData
    @test StochasticGene.combined_modalities(data_rna) == (:rna,)
    @test options_rna isa StochasticGene.MHOptions
end

@testset "key-based fit smoke" begin
    key = "fullstack-key-rna"
    resultfolder = joinpath(FULLSTACK_TEST_RESULTFOLDER, "key")
    _write_key_spec!(key, resultfolder)

    out = StochasticGene.fit(;
        key=key,
        root=DEFAULT_ROOT,
        resultfolder=resultfolder,
        samplesteps=2,
        maxtime=30.0,
        nchains=1,
    )
    @test assert_fullstack_fit_output(:rna, out)
    _, _, _, data, _, _ = out
    @test getfield(data.legs.rna, :label) == key
end

@testset "modular CombinedData fit stack" begin
    out = fit_fullstack_rnatracejoint(;
        resultfolder=FULLSTACK_TEST_RESULTFOLDER,
        fullstack_smoke_defaults()...,
    )
    @test assert_fullstack_fit_output(:rnatracejoint, out)
    _, _, _, data, _, _ = out
    @test data isa StochasticGene.CombinedData
    @test StochasticGene.combined_modalities(data) == (:rna, :tracejoint)
end

if FULLSTACK_EXTENDED
    @testset "extended full-stack matrix" begin
        outputs = run_fullstack_inference_matrix!(
            cases="rna,trace",
            specs="mh_fast,nuts_forwarddiff";
            regenerate=false,
            resultfolder=joinpath(FULLSTACK_TEST_RESULTFOLDER, "extended"),
            samplesteps=20,
            warmupsteps=5,
            maxtime=120.0,
            maxiter=20,
            n_mc=2,
        )
        @test !isempty(outputs)
        for ((_, case), out) in outputs
            @test assert_fullstack_fit_output(case, out)
        end
    end
end
