# Equality: mutating rate stack vs AD-friendly stack; AD check on get_rates_ad + forward match on predictions.
using LinearAlgebra
using StochasticGene
using Test
using Zygote

@testset "fixed_rates vs fixed_rates_ad" begin
    r = Float64[0.1, 0.2, 0.3, 0.4, 0.5]
    @test StochasticGene.fixed_rates(r, tuple()) ≈ StochasticGene.fixed_rates_ad(r, tuple())
    # Each effect is a vector: [reference, slave1, slave2, ...]
    fe = ([1, 3, 4],)
    @test StochasticGene.fixed_rates(copy(r), fe) ≈ StochasticGene.fixed_rates_ad(r, fe)
    fe2 = ([2, 4], [1, 3, 5])
    @test StochasticGene.fixed_rates(copy(r), fe2) ≈ StochasticGene.fixed_rates_ad(r, fe2)
end

@testset "get_rates vs get_rates_ad (single-unit GM)" begin
    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 60
    h = StochasticGene.simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=20_000, nalleles=2)[1]
    data = StochasticGene.RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    model = StochasticGene.load_model(
        data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        StochasticGene.prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    θ = StochasticGene.get_param(model)
    r1 = StochasticGene.get_rates(θ, model)
    r2 = StochasticGene.get_rates_ad(θ, model)
    @test r1 ≈ r2
    @test maximum(abs.(r1 .- r2)) == 0.0
end

function _central_grad(f, x::AbstractVector{Float64}; ε::Float64=1e-5)
    g = similar(x)
    δ = zeros(length(x))
    for i in eachindex(x)
        δ[i] = ε
        g[i] = (f(x .+ δ) - f(x .- δ)) / (2ε)
        δ[i] = 0.0
    end
    return g
end

@testset "predictedfn: get_rates vs get_rates_ad (same output)" begin
    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 60
    h = StochasticGene.simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=20_000, nalleles=2)[1]
    data = StochasticGene.RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    model = StochasticGene.load_model(
        data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        StochasticGene.prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    θ = StochasticGene.get_param(model)
    p1 = StochasticGene.predictedfn(θ, data, model; steady_state_solver=:augmented, rates_fn=StochasticGene.get_rates)
    p2 = StochasticGene.predictedfn(θ, data, model; steady_state_solver=:augmented, rates_fn=StochasticGene.get_rates_ad)
    @test p1 ≈ p2
end

@testset "get_rates vs get_rates_ad (GRSM, R=1, S=1)" begin
    transitions = ([1, 2], [2, 1])
    G = 2
    R = 1
    S = 1
    insertstep = 1
    rtarget = [0.02, 0.1, 0.5, 0.2, 0.1, 0.01]
    nhist = 24
    bins = collect(1:1.0:200.0)
    hs = StochasticGene.simulator(
        rtarget, transitions, G, R, S, insertstep;
        nalleles=2,
        nhist=nhist,
        totalsteps=80_000,
        bins=bins,
    )
    hRNA = div.(hs[1], 30)
    data = StochasticGene.RNAOnOffData("test", "test", nhist, hRNA, bins, hs[2], hs[3], 1.0)
    rinit = fill(0.01, StochasticGene.num_rates(transitions, R, S, insertstep))
    fittedparam = collect(1:length(rtarget) - 1)
    model = StochasticGene.load_model(
        data,
        rinit,
        StochasticGene.prior_ratemean(transitions, R, S, insertstep, rtarget[end], [], StochasticGene.mean_elongationtime(rtarget, transitions, R)),
        fittedparam,
        tuple(),
        transitions,
        G,
        R,
        S,
        insertstep,
        "",
        2,
        10.0,
        Int[],
        rtarget[end],
        0.05,
        StochasticGene.prob_Gaussian,
        [],
        1,
        tuple(),
        tuple(),
        nothing,
    )
    θ = StochasticGene.get_param(model)
    @test StochasticGene.get_rates(θ, model) ≈ StochasticGene.get_rates_ad(θ, model)
    p1 = StochasticGene.predictedfn(θ, data, model; steady_state_solver=:augmented, rates_fn=StochasticGene.get_rates)
    p2 = StochasticGene.predictedfn(θ, data, model; steady_state_solver=:augmented, rates_fn=StochasticGene.get_rates_ad)
    @test p1 ≈ p2
end

@testset "get_rates_ad: per-unit rate storage matches flat (GM)" begin
    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 50
    h = StochasticGene.simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=15_000, nalleles=2)[1]
    data = StochasticGene.RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    m = StochasticGene.load_model(
        data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        StochasticGene.prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    θ = StochasticGene.get_param(m)
    r_flat = StochasticGene.get_rates_ad(θ, m)
    m2 = StochasticGene.GMmodel(
        [copy(m.rates), copy(m.rates)],
        m.Gtransitions,
        m.G,
        m.nalleles,
        m.rateprior,
        m.proposal,
        m.fittedparam,
        m.fixedeffects,
        m.method,
        m.components,
        m.reporter,
    )
    @test StochasticGene.get_rates_ad(θ, m2) ≈ r_flat
end

@testset "Zygote vs FD: sum(get_rates_ad(θ))" begin
    transitions = ([1, 2], [2, 1])
    G = 2
    rtarget = [0.33, 0.19, 20.5, 1.0]
    nRNA = 40
    h = StochasticGene.simulator(rtarget, transitions, G, 0, 0, 0, nhist=nRNA, totalsteps=12_000, nalleles=2)[1]
    data = StochasticGene.RNAData{typeof(nRNA),typeof(h)}("", "", nRNA, h, 1.0, [])
    rinit = [0.1, 0.1, 0.1, 1.0]
    fittedparam = [1, 2, 3]
    model = StochasticGene.load_model(
        data, rinit, StochasticGene.prior_ratemean(transitions, 0, 0, 1, rtarget[end], [], 1.0),
        fittedparam, tuple(), transitions, G, 0, 0, 0, "", 2, 10.0, Int[], rtarget[end], 0.02,
        StochasticGene.prob_Gaussian, [], 1, tuple(), tuple(), nothing,
    )
    θ = Vector{Float64}(StochasticGene.get_param(model))
    f(θv) = sum(StochasticGene.get_rates_ad(θv, model))
    z = Zygote.gradient(f, θ)[1]
    fd = _central_grad(f, θ; ε=1e-4)
    @test isapprox(z, fd; rtol=1e-3, atol=1e-3)
end
