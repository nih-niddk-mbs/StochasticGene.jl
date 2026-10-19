# Read-only startup-tuning experiment. No sampling or production-file writes.
include("benchmark_hidden_likelihood.jl")
using Statistics, Random

function tune_hidden_batchsize(key, folder, root)
    data, model = load_hidden_benchmark(key, folder, root)
    p = SG.get_param(model)
    candidates = [8, 16, 32, 64]
    rng = MersenneTwister(318)
    SG.with_batched_hmm_cache() do
        cache = task_local_storage(:stochasticgene_batched_hmm_cache)
        ref = SG.loglikelihood(p, data, model; hmm_stack=:batched)
        check(value) = begin
            @test isapprox(value[1], ref[1]; atol=1e-7, rtol=1e-11)
            @test isapprox(value[2], ref[2]; atol=1e-8, rtol=1e-11)
        end
        # Compile/warm each size before timing, with size-appropriate buffers.
        for size in candidates
            cache.batchsize = size
            cache.workspace = SG.BatchedForwardWorkspace()
            check(SG.loglikelihood(p, data, model; hmm_stack=:batched))
        end
        timings = Dict(size => Float64[] for size in candidates)
        println("TUNING ", key); flush(stdout)
        start = time()
        for repetition in 1:3, size in shuffle(rng, candidates)
            cache.batchsize = size
            cache.workspace = SG.BatchedForwardWorkspace()
            # Allocate buffers outside the timed call, as in a warmed chain.
            SG.loglikelihood(p, data, model; hmm_stack=:batched)
            GC.gc()
            result = @timed SG.loglikelihood(p, data, model; hmm_stack=:batched)
            check(result.value)
            push!(timings[size], result.time)
            println("TIMING key=", key, " batch=", size, " repeat=", repetition,
                " seconds=", result.time, " bytes=", result.bytes)
            flush(stdout)
        end
        medians = [median(timings[size]) for size in candidates]
        best = minimum(medians)
        # Prefer a smaller buffer if it is within 3% of the measured minimum.
        selected = first(size for (size, elapsed) in zip(candidates, medians)
            if elapsed <= 1.03 * best)
        println("SELECTED key=", key, " batch=", selected, " medians=", medians,
            " tuning_seconds=", time()-start)
        flush(stdout)
        # Holdout proposals were not used to choose the batch size.
        for repetition in 1:3
            pp = p .+ repetition * 0.001
            results = Dict()
            for size in shuffle(rng, unique([16, selected]))
                cache.batchsize = size
                cache.workspace = SG.BatchedForwardWorkspace()
                SG.loglikelihood(pp, data, model; hmm_stack=:batched)
                GC.gc()
                results[size] = @timed SG.loglikelihood(pp, data, model; hmm_stack=:batched)
            end
            baseline, tuned = results[16], results[selected]
            @test isapprox(baseline.value[1], tuned.value[1]; atol=1e-7, rtol=1e-11)
            @test isapprox(baseline.value[2], tuned.value[2]; atol=1e-8, rtol=1e-11)
            println("HOLDOUT key=", key, " repeat=", repetition,
                " default_s=", baseline.time, " tuned_s=", tuned.time,
                " speedup=", baseline.time/tuned.time,
                " trace_delta=", maximum(abs.(baseline.value[2]-tuned.value[2])))
            flush(stdout)
        end
    end
end

root = ARGS[1]
folder = joinpath(root, "results", "3Prime-hidden-multitarget-2026-10-16-restart")
for key in ARGS[2:end]
    tune_hidden_batchsize(key, folder, root)
end
