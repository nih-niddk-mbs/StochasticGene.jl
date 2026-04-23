# This file is part of StochasticGene.jl
#
# runtests.jl
#
# Role: **orchestrate** `Pkg.test()` and **assert** outcomes. Reusable scenarios and
# computations called from this file live in `src/test.jl` as exported `test_*` helpers,
# so each check can also be rerun interactively (e.g. `using StochasticGene; test_fit_simrna()`).
#
# Standalone files under `test/` are developer-focused and are not included here.
#
# Default `Pkg.test()` runs AD checks and other fast tests only. Long MCMC / fit
# integration tests are gated: set environment variable
#   STOCHASTICGENE_FULL_TESTS=1
# when validating the full stack (e.g. after end-to-end AD is confirmed).

using Random
using StochasticGene
using Test

const FULL_TESTS = get(ENV, "STOCHASTICGENE_FULL_TESTS", "0") == "1"

@testset "StochasticGene" begin

    @testset "maxtime_seconds" begin
        @test StochasticGene.maxtime_seconds(3600) == 3600.0
        @test StochasticGene.maxtime_seconds(3600.0) == 3600.0
        @test StochasticGene.maxtime_seconds("2h") == 7200.0
        @test StochasticGene.maxtime_seconds("90m") == 5400.0
        @test StochasticGene.maxtime_seconds("120") == 120.0
        @test StochasticGene.maxtime_seconds(" 1.5h ") == 5400.0
    end

    @testset "RNA + ON/OFF + dwell" begin
        h1, h2 = StochasticGene.test_fit_rnaonoff()
        @test isapprox(h1, h2, rtol=0.05)

        h1, h2 = StochasticGene.test_fit_rnadwelltime()
        @test isapprox(h1, h2, rtol=0.3)

        @test StochasticGene.test_load_model_keyword_compatibility()
    end

    @testset "get_rates_ad consistency" begin
        @test StochasticGene.test_get_rates_ad_consistency()
    end

    @testset "run_nuts_fit smoke" begin
        @test StochasticGene.test_run_nuts_fit_smoke()
    end

    if FULL_TESTS
        @testset "RNA histograms" begin
            h1, h2 = StochasticGene.test_compare()
            @test isapprox(h1, h2, rtol=0.05)

            h1, h2 = StochasticGene.test_fit_simrna()
            @test isapprox(h1, h2, rtol=0.05)

            h1, h2 = StochasticGene.test_fit_rna()
            @test isapprox(h1, h2, rtol=0.05)
        end

        @testset "RNA + ON/OFF + dwell" begin
            h1, h2 = StochasticGene.test_fit_rnaonoff()
            @test isapprox(h1, h2, rtol=0.05)

            h1, h2 = StochasticGene.test_fit_rnadwelltime()
            @test isapprox(h1, h2, rtol=0.3)
        end

        @testset "Traces (single-unit)" begin
            lower, target, upper = StochasticGene.test_fit_trace()
            @test all(lower .<= target .<= upper)
            # `test_fit_trace_hierarchical` fits a very high-dimensional hierarchical
            # state (slow, machine-dependent sample counts under maxtime). Run it
            # interactively when validating hierarchical traces; see `src/test.jl`.
        end

        @testset "Coupled traces and dwell" begin
            lower, target, upper = StochasticGene.test_fit_tracejoint()
            @test all(lower .<= target .<= upper)

            # Basic coupled dwell-time vs simulator comparison (short 3-unit example)
            cme_vec, sim_vec = StochasticGene.test_compare_3unit()
            @test length(cme_vec) == length(sim_vec)
        end

        @testset "MH vs NUTS posterior (GM, simulated RNA histogram)" begin
            # Same likelihood and priors; compare posterior mean rates (Stats.meanparam).
            # Large MH chain + moderate NUTS; finite-diff gradients for robust NUTS.
            res = StochasticGene.test_fit_simrna_mh_nuts(;
                totalsteps=80_000,
                sim_seed=42,
                mh_samplesteps=35_000,
                mh_warmup=8_000,
                mh_maxtime=240.0,
                mh_seed=101,
                nuts_n_samples=1_200,
                nuts_n_adapts=600,
                nuts_gradient=:finite,
                nuts_fd_ε=1e-4,
                rng_nuts=MersenneTwister(202),
            )
            @test length(res.mean_mh) == length(res.mean_nuts)
            # Monte Carlo noise: allow generous margin on rate-space means
            @test res.max_abs_diff < 1.25
        end
    end

    @testset "trace_specs utilities" begin
        @test StochasticGene.test_trace_specs_utilities()
    end

    @testset "ObservationBundle (tracerna split / round-trip)" begin
        len = 5
        h = fill(0.1, len)
        trace = [rand(4, 2)]
        tup = (trace, 0.0, 0.0, 4, 1.0)
        d = StochasticGene.TraceRNAData("lab", "gene", 1.0, tup, len, h, 1.0, Int[])
        b = StochasticGene.ObservationBundle(d)
        @test StochasticGene.observation_modalities(b) == (:rna, :trace)
        d2 = StochasticGene.reconstruct_tracerna(b)
        @test d2.label == d.label && d2.gene == d.gene && d2.interval == d.interval
        @test d2.nRNA == d.nRNA && d2.histRNA == d.histRNA && d2.yield == d.yield
        @test d2.trace == d.trace && d2.units == d.units
    end

    @testset "inference benchmark (smoke)" begin
        scen = StochasticGene.benchmark_inference_simrna_small(seed=7)
        @test hasproperty(scen, :data) && hasproperty(scen, :model)
        @test length(StochasticGene.get_param(scen.model)) > 0
        tr = StochasticGene.benchmark_inference_trace_gr2r2(seed=11, totaltime=120.0, ntrials=2)
        @test hasproperty(tr, :data) && hasproperty(tr, :model)
        @test length(StochasticGene.get_param(tr.model)) > 0
        cj = StochasticGene.benchmark_inference_trace_coupled_3x3(seed=11, totaltime=120.0, ntrials=2)
        @test hasproperty(cj, :data) && hasproperty(cj, :model)
        @test length(StochasticGene.get_param(cj.model)) > 0
        cjh = StochasticGene.benchmark_inference_trace_coupled_3x3_g3r0(seed=11, totaltime=120.0, ntrials=2)
        @test hasproperty(cjh, :data) && hasproperty(cjh, :model)
        @test length(StochasticGene.get_param(cjh.model)) > 0
    end

end
