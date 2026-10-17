using Test, Random, LinearAlgebra, StochasticGene

@testset "Experimental batched forward vs reference" begin
    rng = MersenneTwister(524)
    for n in (1, 4, 24), lengths in ([1], [1, 3, 7], fill(15, 19))
        a = rand(rng, n, n)
        a ./= sum(a; dims=2)
        p0 = rand(rng, n); p0 ./= sum(p0)
        bs = [rand(rng, n, t) for t in lengths]
        old = [-sum(log.(StochasticGene.forward(a, b, p0)[2])) for b in bs]
        new = StochasticGene._forward_batched_loglik(a, p0, bs)
        @test isapprox(old, new; atol=1e-10, rtol=1e-12)
    end
    a = [0.9 0.1; 0.2 0.8]; p0 = [0.5, 0.5]
    for b in (zeros(2, 4), fill(1e-300, 2, 4), fill(10.0, 2, 4))
        old = -sum(log.(StochasticGene.forward(a, b, p0)[2]))
        @test StochasticGene._forward_batched_loglik(a, p0, [b])[1] ≈ old
    end
    @test isempty(StochasticGene._forward_batched_loglik(a, p0, Matrix{Float64}[]))
    @test_throws ArgumentError StochasticGene._forward_batched_loglik(a, p0, [zeros(2, 0)])
end
