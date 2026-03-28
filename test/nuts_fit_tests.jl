# run_nuts_fit: same return shape as run_mh (Fit, Stats, Measures) + nuts_info
using Random
using StochasticGene
using Test

@testset "run_nuts_fit (GM, finite-diff gradients, short chain)" begin
    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 40
    h = StochasticGene.simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=8000, nalleles=2)[1]
    data = StochasticGene.RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    model = StochasticGene.load_model(
        data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        StochasticGene.prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    opts = StochasticGene.NUTSOptions(;
        n_samples=12,
        n_adapts=12,
        δ=0.8,
        gradient=:finite,
        fd_ε=1e-4,
        verbose=false,
        progress=false,
    )
    rng = MersenneTwister(42)
    fits, stats, measures, nuts_info = StochasticGene.run_nuts_fit(
        data, model, opts;
        rng=rng,
        steady_state_solver=:augmented,
    )
    @test size(fits.param, 2) == 12
    @test length(fits.ll) == 12
    @test fits.total == opts.n_adapts + 12
    @test length(stats.meanparam) == size(fits.param, 1)
    @test measures.waic isa Tuple
    @test length(measures.ess) == size(fits.param, 1)
    @test haskey(nuts_info, :nuts_stats)
    @test haskey(nuts_info, :initial_θ)
    # Alias
    out = StochasticGene.run_NUTS(data, model, opts; rng=rng, steady_state_solver=:augmented)
    @test size(out[1].param, 2) == 12
end
