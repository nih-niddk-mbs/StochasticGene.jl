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

    @test StochasticGene.test_spec_conversion()
    @test StochasticGene.test_spec_io_roundtrip()

    @test StochasticGene.test_spec_trace_in_fit()
    @test StochasticGene.test_spec_dwell_in_fit()

end