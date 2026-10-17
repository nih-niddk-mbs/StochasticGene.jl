# Explicit, opt-in integration check; never modifies the source results.
using Distributed, StochasticGene, Test, LinearAlgebra
@everywhere using StochasticGene, LinearAlgebra
@everywhere BLAS.set_num_threads(1)

root = ARGS[1]
key = "H2-targets135-free-h"
source = joinpath(root, "results", "3Prime-hidden-multitarget-2026-10-16-restart")
ratefile = joinpath(source, "rates_$(key).txt")
source_bytes = read(ratefile)
spec = StochasticGene.read_run_spec(joinpath(source, "info_$(key).jld2"))
initial = StochasticGene.readrates(ratefile, StochasticGene.get_row("median"))
out = mktempdir(; prefix="batched-fit-smoke-", cleanup=false)
kwargs = filter(p -> haskey(StochasticGene.fit_default_spec(), first(p)), spec)
merge!(kwargs, Dict{Symbol,Any}(
    :key => key, :root => root, :resultfolder => out,
    :rinit => initial, :init_ratefile => nothing,
    :nchains => 2, :parallel => :distributed,
    :likelihood_executor => :batched, :samplesteps => 2,
    :warmupsteps => 0, :maxtime => 180.0, :sample_stride => 1,
    :propcv => 0.02, :writesamples => false,
))
println("SMOKE OUTPUT: ", out); flush(stdout)
fit, stats, measures, data, model, options = StochasticGene.fit(; kwargs...)
@testset "Batched end-to-end fit" begin
    @test options.likelihood_executor == :batched
    @test fit.total > 0
    @test isfinite(fit.llml)
    rates, _ = StochasticGene.read_rates_table(joinpath(out, "rates_$(key).txt"))
    @test size(rates) == (4, length(initial))
    width = div(length(initial), length(data.trace[1]) + 2)
    fixed = setdiff(1:width, spec[:fittedparam])
    for row in 1:4, b in 0:length(data.trace[1])+1
        @test rates[row, b*width .+ fixed] == initial[b*width .+ fixed]
    end
    @test rates[:, width+1:2width] == repeat(permutedims(initial[width+1:2width]), 4, 1)
    @test read(ratefile) == source_bytes
    @test isfile(joinpath(out, "info_$(key).jld2"))
    @test isfile(joinpath(out, "measures_$(key).txt"))
    @test isfile(joinpath(out, "param-stats_$(key).txt"))
end
