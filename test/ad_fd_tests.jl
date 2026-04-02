# Zygote (reverse-mode AD) vs central finite differences for augmented steady state.
using LinearAlgebra
using SparseArrays
using StochasticGene
using Test
using Zygote

"""2×2 column CTMC generator: columns sum to 0; steady π ∝ [θ1, θ2]."""
function twostate_generator(θ::AbstractVector{T}) where {T}
    length(θ) == 2 || throw(ArgumentError("expected length-2 rate vector"))
    θ1, θ2 = θ[1], θ[2]
    # COO construction: avoid matrix literals like [a b; c d] (hvcat uses in-place fill; Zygote errors).
    sparse(Int[1, 2, 1, 2], Int[1, 1, 2, 2], T[-θ2, θ2, θ1, -θ1], 2, 2)
end

function central_grad(f, x::AbstractVector{<:Real}; ε=1e-6)
    g = similar(x, Float64)
    δ = zeros(length(x))
    for i in eachindex(x)
        δ[i] = ε
        g[i] = (f(x .+ δ) - f(x .- δ)) / (2ε)
        δ[i] = 0
    end
    return g
end

@testset "Zygote vs FD: normalized_nullspace_augmented (sum)" begin
    θ0 = Float64[1.2, 0.7]
    f(θ) = sum(StochasticGene.normalized_nullspace_augmented(twostate_generator(θ)))
    z = Zygote.gradient(f, θ0)[1]
    fd = central_grad(f, θ0; ε=1e-6)
    @test isapprox(z, fd; rtol=0.02, atol=1e-3)
end

@testset "Zygote vs FD: normalized_nullspace_augmented (dot observable)" begin
    w = Float64[0.3, 0.7]
    θ0 = Float64[2.0, 1.5]
    g(θ) = dot(w, StochasticGene.normalized_nullspace_augmented(twostate_generator(θ)))
    z = Zygote.gradient(g, θ0)[1]
    fd = central_grad(g, θ0; ε=1e-6)
    @test isapprox(z, fd; rtol=0.02, atol=1e-3)
end

@testset "benchmark_trace_compare: ForwardDiff vs finite-diff gradient agreement" begin
    sc = StochasticGene.benchmark_inference_trace_gr2r2(seed=1, totaltime=80.0, ntrials=1)
    r = StochasticGene.benchmark_trace_compare_forwarddiff_vs_finitediff(
        sc; param_indices=[1, 2], nruns=1, warmup=true, fd_ε=1e-4,
    )
    @test r.gradient_max_abs_diff < 1e-5
    @test r.gradient_max_rel_diff < 0.05
end

@testset "check_ad_gradient_feasibility: coupled trace + ForwardDiff path" begin
    sc = StochasticGene.benchmark_inference_trace_coupled_3x3(seed=1, totaltime=10.0, ntrials=1)
    StochasticGene.check_ad_gradient_feasibility(
        sc.data, sc.model;
        forwarddiff_through_likelihood=true,
        zygote_through_likelihood=false,
    )
    @test StochasticGene.longest_trace_timesteps(sc.data) isa Union{Nothing,Int}
end

@testset "Zygote vs FD: steady_state_vector (solver=:augmented)" begin
    θ0 = Float64[1.5, 1.1]
    h(θ) = dot([0.2, 0.8], StochasticGene.steady_state_vector(twostate_generator(θ); solver=:augmented))
    z = Zygote.gradient(h, θ0)[1]
    fd = central_grad(h, θ0; ε=1e-6)
    @test isapprox(z, fd; rtol=0.02, atol=1e-3)
end
