# Read-only comparison against the uncached batched executor in main at c9bae8a.
include("benchmark_hidden_likelihood.jl")

function forward_before_cache(a, p0, emissions)
    n = size(a, 1); k = length(emissions)
    lengths = [size(b, 2) for b in emissions]
    previous = zeros(n, k); next = similar(previous); ll = zeros(k)
    for j in 1:k
        @views previous[:, j] .= p0 .* emissions[j][:, 1]
        c = 1 / max(sum(@view previous[:, j]), eps(Float64))
        @views previous[:, j] .*= c
        ll[j] = -log(c)
    end
    for t in 2:maximum(lengths)
        mul!(next, transpose(a), previous)
        for j in 1:k
            t > lengths[j] && continue
            @views next[:, j] .*= emissions[j][:, t]
            c = 1 / max(sum(@view next[:, j]), eps(Float64))
            @views next[:, j] .*= c
            ll[j] -= log(c)
        end
        previous, next = next, previous
    end
    ll
end

function before_cache(param, data, model)
    r = SG.prepare_rates(param, model)
    rates, _, noise, individualnoise, _, _, coupling, _ = r
    components = SG.get_components(model, data)
    @assert components isa SG.TCoupledFullComponents && model.method[2]
    cr = [coupling[1][k] * rates[1][components.targets[k][1]][components.targets[k][2]]
        for k in eachindex(components.targets)]
    a, p0 = SG.make_ap(rates[1], cr, data.interval, components, model.method[1])
    reporter = SG.get_reporter(model, data)
    observed = isempty(data.units) ? nothing : data.units
    re, ne = SG._filter_observed(reporter, individualnoise, observed)
    traces = data.trace[1]; ll = zeros(length(traces))
    for first in 1:16:length(traces)
        ids = first:min(first+15, length(traces))
        bs = [SG.set_b(traces[i], SG.set_d(ne[i], re)) for i in ids]
        ll[ids] = forward_before_cache(a, p0, bs)
    end
    lb = SG.ll_off(data.trace, rates[1], noise[1], reporter, data.interval, components, model.method[1])
    sum(ll) + lb + SG.hierarchy_logdensity(r, model), ll
end

function benchmark_cache(key, folder, root)
    data, model = load_hidden_benchmark(key, folder, root)
    p = SG.get_param(model)
    println("CACHE BENCHMARK ", key); flush(stdout)
    SG.with_batched_hmm_cache() do
        ref = before_cache(p, data, model)
        cold = @timed SG.loglikelihood(p, data, model; hmm_stack=:batched)
        @test isapprox(ref[2], cold.value[2]; atol=1e-8, rtol=1e-11)
        @test isapprox(ref[1], cold.value[1]; atol=1e-7, rtol=1e-11)
        for repetition in 1:3
            # Different proposal every repetition, including changed coupling.
            pp = p .+ repetition * 0.001
            GC.gc()
            old = @timed before_cache(pp, data, model)
            GC.gc()
            new = @timed SG.loglikelihood(pp, data, model; hmm_stack=:batched)
            @test isapprox(old.value[2], new.value[2]; atol=1e-8, rtol=1e-11)
            @test isapprox(old.value[1], new.value[1]; atol=1e-7, rtol=1e-11)
            println("RESULT ", key, " repetition=", repetition, " old_s=", old.time,
                " cached_s=", new.time, " speedup=", old.time/new.time,
                " old_bytes=", old.bytes, " cached_bytes=", new.bytes,
                " trace_delta=", maximum(abs.(old.value[2]-new.value[2])))
            flush(stdout)
        end
        cache = task_local_storage(:stochasticgene_batched_hmm_cache)
        println("CACHE states=", length(cache.state_map), " emission_groups=", length(cache.representatives),
            " retained_bytes=", Base.summarysize(cache), " cold_s=", cold.time)
        flush(stdout)
    end
end

root = ARGS[1]
folder = joinpath(root, "results", "3Prime-hidden-multitarget-2026-10-16-restart")
for key in ARGS[2:end]
    benchmark_cache(key, folder, root)
end
