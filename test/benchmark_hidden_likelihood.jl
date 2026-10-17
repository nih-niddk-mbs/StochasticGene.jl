# Read-only benchmark: no sampling and no writes to the results folder.
using StochasticGene, LinearAlgebra, Test
const SG = StochasticGene
BLAS.set_num_threads(1)
include("batched_forward_tests.jl")

function benchmark_hidden(key, folder, root)
    s = SG.read_run_spec(joinpath(folder, "info_$(key).jld2"))
    rates = SG.readrates(joinpath(folder, "rates_$(key).txt"), SG.get_row("median"))
    data = SG.load_data_trace(joinpath(root, s[:datapath]), s[:label], s[:gene],
        s[:datacond], (5/3, 1.0, -1.0), :tracejoint, s[:datacol], s[:zeromedian], 1.0, s[:trace_specs])
    model = SG.load_model(data, rates, s[:priormean], s[:fittedparam], s[:fixedeffects],
        s[:transitions], s[:G], s[:R], s[:S], s[:insertstep], s[:splicetype],
        s[:nalleles], s[:priorcv], s[:onstates], s[:decayrate], 0.02, s[:probfn],
        s[:noisepriors], s[:method], s[:hierarchical], s[:coupling], s[:grid], s[:zeromedian])
    p = SG.get_param(model)
    println("MODEL ", key, " traces=", length(data.trace[1]), " params=", length(p)); flush(stdout)
    # Warm both dispatch paths before the measured likelihood evaluation.
    ref = SG.loglikelihood(p, data, model; hmm_stack=:fast)
    opt = SG.loglikelihood(p, data, model; hmm_stack=:batched)
    @test isapprox(ref[1], opt[1]; atol=1e-7, rtol=1e-11)
    @test isapprox(ref[2], opt[2]; atol=1e-8, rtol=1e-11)
    println("Warmup agreement: total delta=", opt[1]-ref[1],
        " max trace delta=", maximum(abs.(opt[2]-ref[2]))); flush(stdout)
    GC.gc()
    baseline = @timed SG.loglikelihood(p, data, model; hmm_stack=:fast)
    GC.gc()
    batched = @timed SG.loglikelihood(p, data, model; hmm_stack=:batched)
    @test isapprox(baseline.value[2], batched.value[2]; atol=1e-8, rtol=1e-11)
    println("RESULT ", key, " reference_s=", baseline.time, " batched_s=", batched.time,
        " speedup=", baseline.time/batched.time,
        " reference_bytes=", baseline.bytes, " batched_bytes=", batched.bytes)
    flush(stdout)
    # A second point exercises changed kinetic/coupling parameters, not just one seed.
    p2 = p .+ 0.01
    r2 = SG.loglikelihood(p2, data, model; hmm_stack=:fast)
    o2 = SG.loglikelihood(p2, data, model; hmm_stack=:batched)
    @test isapprox(r2[2], o2[2]; atol=1e-8, rtol=1e-11)
    println("Perturbed-parameter agreement passed: ", key); flush(stdout)
end

root = ARGS[1]
folder = joinpath(root, "results", "3Prime-hidden-multitarget-2026-10-16-restart")
for key in ARGS[2:end]
    benchmark_hidden(key, folder, root)
end
