using Test, Random, LinearAlgebra, StochasticGene

@testset "Non-hierarchical batched emissions" begin
    SG = StochasticGene
    a = [0.85 0.15; 0.3 0.7]
    p0 = [2/3, 1/3]
    for d in ([SG.Bernoulli(0.0), SG.Bernoulli(1.0)],
              [SG.Normal(-0.1, 0.3), SG.Normal(1.0, 0.5)])
        for traces in ([Float64[0]], [Float64[0,1,1], Float64[1,0], Float64[0,1,0,1]])
            reference = SG._ll_hmm(a, p0, d, traces; hmm_stack=:fast)
            cache = SG.BatchedHMMCache()
            for batchsize in (1, 2, 16)
                actual = SG._ll_hmm_batched(a, p0, d, traces; cache, batchsize)
                @test isapprox(actual[1], reference[1]; atol=1e-10, rtol=1e-12)
                @test isapprox(actual[2], reference[2]; atol=1e-10, rtol=1e-12)
            end
            routed = SG._ll_hmm(a, p0, d, traces; hmm_stack=:batched)
            @test routed[2] ≈ reference[2]
            traces[1][1] = 1 - traces[1][1]
            @test SG._ll_hmm_batched(a, p0, d, traces; cache)[2] ≈
                  SG._ll_hmm(a, p0, d, traces; hmm_stack=:fast)[2]
        end
    end
    @test isempty(SG._ll_hmm_batched(a, p0, [SG.Bernoulli(0.0), SG.Bernoulli(1.0)], Vector{Float64}[])[2])
    @test_throws ArgumentError SG._ll_hmm_batched(a, p0, [], []; batchsize=0)
end

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

@testset "Budgeted batch tuning" begin
    SG = StochasticGene
    cache = SG.BatchedHMMCache()
    reference = (3.0, [1.0, 2.0])
    @test SG._tune_hmm_batchsize!(() -> reference, cache, reference).status == :unsupported
    cache.workspace.previous = zeros(2, 2)
    clock_value = Ref(0.0)
    clock() = clock_value[]
    evaluate() = begin
        clock_value[] += cache.batchsize == 32 ? 1.0 : 2.0
        reference
    end
    report = SG._tune_hmm_batchsize!(evaluate, cache, reference; seconds=100.0, clock=clock)
    @test report.batchsize == 32
    @test report.rounds == 3
    @test report.status == :selected
    cache.batchsize = 16
    report = SG._tune_hmm_batchsize!(evaluate, cache, reference; seconds=1.0, clock=clock)
    @test report.status == :budget_exhausted
    @test cache.batchsize == 16
    @test_throws ArgumentError SG._tune_hmm_batchsize!(() -> (4.0, [2.0, 2.0]), cache, reference)
    @test cache.batchsize == 16
    @test_throws ErrorException SG._tune_hmm_batchsize!(() -> error("failure"), cache, reference)
    @test cache.batchsize == 16
    @test_throws ArgumentError SG.MHOptions(2, 0, 10.0, 1.0; hmm_batchsize=0)
    @test_throws ArgumentError SG.MHOptions(2, 0, 10.0, 1.0; tune_hmm_batchsize=true)
    @test_throws ArgumentError SG.MHOptions(2, 0, 10.0, 1.0; hmm_batch_tuning_seconds=Inf)
    options = SG.load_options(Dict{Symbol,Any}(:likelihood_executor => :batched,
        :tune_hmm_batchsize => true, :hmm_batchsize => 32, :hmm_batch_tuning_seconds => 12.0))
    @test options.tune_hmm_batchsize
    @test options.hmm_batchsize == 32
    @test options.hmm_batch_tuning_seconds == 12.0
    mktempdir() do folder
        spec = SG.fit_default_spec()
        merge!(spec, Dict(:likelihood_executor => :batched, :tune_hmm_batchsize => true,
            :hmm_batchsize => 32, :hmm_batch_tuning_seconds => 12.0))
        path = SG.write_fitfile_full(joinpath(folder, "fit.jl"), spec)
        script = read(path, String)
        @test occursin("tune_hmm_batchsize=true", script)
        @test occursin("hmm_batchsize=32", script)
        @test occursin("hmm_batch_tuning_seconds=12.0", script)
    end
end

@testset "Batched emission cache invalidation and ownership" begin
    SG = StochasticGene
    rng = MersenneTwister(922)
    a = rand(rng, 6, 6); a ./= sum(a; dims=2)
    p0 = fill(1/6, 6)
    reporter = [SG.HMMReporter(4, [0, 1, 0, 1, 0, 1], SG.prob_Gaussian, 0, Int[], Int[]),
        SG.HMMReporter(4, [0, 0, 1, 1, 0, 0], SG.prob_Gaussian, 0, Int[], Int[])]
    noise = [[[-0.1, 0.2, 0.8, 0.15], [0.05, 0.3, 0.6, 0.2]] for _ in 1:19]
    traces = [rand(rng, 3 + i, 2) for i in 1:19]
    cache = SG.BatchedHMMCache()
    function check()
        ref = SG._ll_hmm(noise, a, p0, reporter, traces; hmm_stack=:fast)
        actual = SG._ll_hmm_batched(noise, a, p0, reporter, traces; cache=cache)
        @test isapprox(ref[2], actual[2]; atol=1e-9, rtol=1e-12)
        actual
    end
    original = check()
    @test length(cache.representatives) == 4
    @test all(b -> size(b, 1) == 4, cache.emissions)
    stored = copy(cache.emissions)
    buffer = cache.workspace.previous
    check()
    @test all(cache.emissions[i] === stored[i] for i in eachindex(stored))
    @test cache.workspace.previous === buffer
    for size in (1, 8, 32, 64)
        cache.batchsize = size
        check()
    end
    cache.batchsize = 0
    @test_throws ArgumentError SG._ll_hmm_batched(noise, a, p0, reporter, traces; cache=cache)
    @test SG._ll_hmm_batched(noise, a, p0, reporter, traces; cache=cache, batchsize=16)[2] ≈ original[2]
    cache.batchsize = 16
    noise[2][1][2] *= 1.2
    check()
    @test cache.emissions[1] === stored[1]
    @test cache.emissions[2] !== stored[2]
    noise[2][1][2] /= 1.2
    check() # A rejected proposal returning to old noise must rebuild correctly.
    traces[3][1, 1] += 0.5
    check()
    @test cache.emissions[3] !== stored[3]
    reporter[1].per_state[1] = 2
    check()
    @test length(cache.representatives) == 5
    traces[1] = rand(rng, 1, 2)
    check()
    pop!(traces); pop!(noise)
    check()
    # Changing kinetics never changes the emission cache.
    stored = copy(cache.emissions)
    a .= 1/6
    check()
    @test all(cache.emissions[i] === stored[i] for i in eachindex(stored))
    @test original[2] != check()[2] # Previous pointwise results own their storage.
    custom(par, counts) = SG.prob_Gaussian(par, counts)
    reporter[1] = SG.HMMReporter(4, copy(reporter[1].per_state), custom, 0, Int[], Int[])
    check()
    @test isempty(cache.emissions) # Unknown emission functions are never cached.
    slot = :stochasticgene_batched_hmm_cache
    outer = get(task_local_storage(), slot, nothing)
    SG.with_batched_hmm_cache() do
        first = task_local_storage(slot)
        SG.with_batched_hmm_cache() do
            @test task_local_storage(slot) !== first
        end
        @test task_local_storage(slot) === first
    end
    @test get(task_local_storage(), slot, nothing) === outer
    @test_throws ErrorException SG.with_batched_hmm_cache(() -> error("scope cleanup"))
    @test get(task_local_storage(), slot, nothing) === outer
    tasks = [@async SG.with_batched_hmm_cache(() -> task_local_storage(slot)) for _ in 1:2]
    @test fetch(tasks[1]) !== fetch(tasks[2])
end
