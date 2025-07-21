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
    @test isapprox(h1, h2, rtol=0.2)

    h1, h2 = StochasticGene.test_fit_trace()
    @test isapprox(h1, h2, rtol=1.0)

    h1, h2 = StochasticGene.test_fit_trace_hierarchical()
    @test isapprox(h1, h2, rtol=1.)

    h1, h2 = StochasticGene.test_fit_tracejoint()
    @test isapprox(h1, h2, rtol=0.2)

end