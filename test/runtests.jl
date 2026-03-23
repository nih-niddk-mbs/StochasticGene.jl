# This file is part of StochasticGene.jl
#
# runtests.jl
#
# All test functions are in test.jl
#

using StochasticGene
using Test

@testset "StochasticGene" begin

    h1, h2 = StochasticGene.test_compare()
    @test isapprox(h1, h2, rtol=0.05)

    h1, h2 = StochasticGene.test_compare_coupling()
    @test isapprox(h1, h2, rtol=0.05)

    h1, h2 = StochasticGene.test_fit_simrna()
    @test isapprox(h1, h2, rtol=0.05)

    h1, h2 = StochasticGene.test_fit_rna()
    @test isapprox(h1, h2, rtol=0.05)

    h1, h2 = StochasticGene.test_fit_rnaonoff()
    @test isapprox(h1, h2, rtol=0.05)

    h1, h2 = StochasticGene.test_fit_rnadwelltime()
    @test isapprox(h1, h2, rtol=0.3)

    lower, target, upper = StochasticGene.test_fit_trace()
    @test lower <= target <= upper

    lower, target, upper = StochasticGene.test_fit_trace_hierarchical()
    @test lower <= target <= upper

    lower, target, upper = StochasticGene.test_fit_tracejoint()
    @test lower <= target <= upper

    @testset "trace_specs defaults for coupled tracejoint" begin
        c2 = ((1, 2), NTuple{4,Int}[], Symbol[])
        @test StochasticGene.n_observed_trace_units(c2) == 2
        c3 = ((1, 2, 3), NTuple{4,Int}[], Symbol[])
        @test StochasticGene.n_observed_trace_units(c3) == 2
        sp = StochasticGene.default_trace_specs_for_coupled((1.0, 1.0, -1.0), [true, true], 2)
        @test length(sp) == 2 && sp[1].unit == 1 && sp[2].unit == 2 && sp[1].interval == 1.0
    end

end