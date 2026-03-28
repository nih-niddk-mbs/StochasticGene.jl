# gradient_inference.jl — NUTS (AdvancedHMC) and mean-field ADVI on the same
# transformed parameter space as Metropolis–Hastings (`get_param`, `logprior`).
#
# Gradients default to reverse-mode AD (Zygote via LogDensityProblemsAD / AdvancedHMC).
# Use `gradient=:finite` in `NUTSOptions` / `ADVIOptions` for central finite differences.

"""
    GenePosteriorLogDensity

`LogDensityProblems.jl` wrapper for the unnormalized log posterior
`logprior(θ, model) + loglikelihood(θ, data, model)` with `θ` in **transformed**
space (same as `get_param(model)`).

# Fields
- `data`, `model`: passed through to `loglikelihood` / `logprior`
- `steady_state_solver`: forwarded to likelihood (default `:augmented` is best for downstream use)
- `ad_likelihood`: if `true` and `data isa RNACountData`, uses `loglikelihood_ad`; otherwise `loglikelihood`
"""
struct GenePosteriorLogDensity{D,M}
    data::D
    model::M
    steady_state_solver::Symbol
    ad_likelihood::Bool
end

function _ll_first(θ, data, model, steady_state_solver::Symbol, ad_likelihood::Bool)
    if ad_likelihood && data isa RNACountData
        return loglikelihood_ad(θ, data, model; steady_state_solver=steady_state_solver)[1]
    end
    return loglikelihood(θ, data, model; steady_state_solver=steady_state_solver)[1]
end

"""
    logposterior(θ, data, model; steady_state_solver=:augmented, ad_likelihood=nothing)

Unnormalized log posterior in transformed coordinates. `ad_likelihood` defaults to
`true` for `RNACountData` (uses `loglikelihood_ad`), else `false`.
"""
function logposterior(
    θ::AbstractVector,
    data,
    model::AbstractGeneTransitionModel;
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
)
    ad = something(ad_likelihood, data isa RNACountData)
    return logprior(θ, model) + _ll_first(θ, data, model, steady_state_solver, ad)
end

function LogDensityProblems.dimension(p::GenePosteriorLogDensity)
    return length(get_param(p.model))
end

function LogDensityProblems.logdensity(p::GenePosteriorLogDensity, θ)
    return logprior(θ, p.model) + _ll_first(θ, p.data, p.model, p.steady_state_solver, p.ad_likelihood)
end

function LogDensityProblems.capabilities(::Type{<:GenePosteriorLogDensity})
    return LogDensityProblems.LogDensityOrder{0}()
end

"""Central finite-difference gradient of `LogDensityProblems.logdensity(ℓ, ·)`."""
function _finitediff_grad!(
    g::AbstractVector{Float64},
    ℓ,
    θ::AbstractVector{Float64},
    ε::Float64,
)
    @assert length(g) == length(θ)
    for i in eachindex(θ)
        θp = copy(θ)
        θm = copy(θ)
        θp[i] += ε
        θm[i] -= ε
        g[i] = (LogDensityProblems.logdensity(ℓ, θp) - LogDensityProblems.logdensity(ℓ, θm)) / (2ε)
    end
    return g
end

function _make_finitediff_∂ℓπ∂θ(ℓ, ε::Float64)
    function ∂ℓπ∂θ(θ::AbstractVector)
        v = LogDensityProblems.logdensity(ℓ, θ)
        g = Vector{Float64}(undef, length(θ))
        _finitediff_grad!(g, ℓ, θ, ε)
        return (v, g)
    end
end

function _make_ℓπ(ℓ)
    return function (θ::AbstractVector)
        return LogDensityProblems.logdensity(ℓ, θ)
    end
end

"""
    NUTSOptions

Options for `run_nuts`.

# Fields
- `n_samples`, `n_adapts`: post-warmup samples and adaptation steps
- `δ`: target acceptance (NUTS dual averaging)
- `gradient`: `:Zygote` (default) or `:finite` (central differences; uses `fd_ε`)
- `fd_ε`: finite-difference step when `gradient === :finite`
- `verbose`, `progress`: passed to `AdvancedHMC.sample`
"""
struct NUTSOptions
    n_samples::Int
    n_adapts::Int
    δ::Float64
    gradient::Symbol
    fd_ε::Float64
    verbose::Bool
    progress::Bool
end

NUTSOptions(; n_samples=1000, n_adapts=1000, δ=0.8, gradient=:Zygote, fd_ε=1e-5, verbose=true, progress=false) =
    NUTSOptions(n_samples, n_adapts, δ, gradient, fd_ε, verbose, progress)

"""
    run_nuts(data, model, rng, options=NUTSOptions(); kwargs...)

Hamiltonian Monte Carlo with the No-U-Turn sampler (NUTS) and diagonal mass matrix
adaptation (Stan-style). Gradients of the log posterior use **Zygote** by default;
set `options.gradient=:finite` for central differences (`options.fd_ε`).

`rng` must be an `AbstractRNG` (e.g. `using Random; Random.default_rng()`).

Keyword arguments:
- `steady_state_solver`: passed to likelihood (default `:augmented`)
- `ad_likelihood`: override automatic choice (`nothing` → use `loglikelihood_ad` only for `RNACountData`)

Returns a `NamedTuple` with `samples` (`d × n` matrix), `sample_vectors`, `nuts_stats`, and `initial_θ`.
"""
function run_nuts(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    rng,
    options::NUTSOptions=NUTSOptions();
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
)
    ad = something(ad_likelihood, data isa RNACountData)
    ℓ = GenePosteriorLogDensity(data, model, steady_state_solver, ad)
    θ0 = Vector{Float64}(get_param(model))
    D = length(θ0)
    D == LogDensityProblems.dimension(ℓ) || throw(DimensionMismatch("parameter dimension mismatch"))

    metric = DiagEuclideanMetric(D)
    ham = if options.gradient === :finite
        ℓπ = _make_ℓπ(ℓ)
        ∂ℓπ∂θ = _make_finitediff_∂ℓπ∂θ(ℓ, Float64(options.fd_ε))
        Hamiltonian(metric, ℓπ, ∂ℓπ∂θ)
    else
        Hamiltonian(metric, ℓ, :Zygote)
    end
    ϵ = find_good_stepsize(rng, ham, θ0)
    integrator = Leapfrog(ϵ)
    kernel = HMCKernel(Trajectory{MultinomialTS}(integrator, GeneralisedNoUTurn()))
    adaptor = StanHMCAdaptor(MassMatrixAdaptor(metric), StepSizeAdaptor(options.δ, integrator))

    samples_θ, st = AdvancedHMC.sample(
        rng,
        ham,
        kernel,
        θ0,
        options.n_samples,
        adaptor,
        options.n_adapts;
        verbose=options.verbose,
        progress=options.progress,
    )
    mat = reduce(hcat, samples_θ)
    return (
        samples=mat,
        sample_vectors=samples_θ,
        nuts_stats=st,
        initial_θ=θ0,
    )
end

@inline function _softplus(x::Real)
    return log1p(exp(x))
end

"""
    ADVIOptions

Black-box ADVI-style mean-field Gaussian VI: `q(θ) = ∏_i N(θ_i | μ_i, σ_i^2)` in the
same transformed space as MCMC, with `σ_i = softplus(s_i) + ε`.

# Fields
- `maxiter`: `Optim` iterations
- `n_mc`: fixed Monte Carlo draws `ε` for the reparameterization gradient (deterministic objective)
- `σ_floor`: lower bound on `σ_i`
- `gradient`: `:Zygote` (default) or `:finite` (`Optim` finite-difference gradients)
- `verbose`: `Optim` show trace
"""
struct ADVIOptions
    maxiter::Int
    n_mc::Int
    σ_floor::Float64
    verbose::Bool
    gradient::Symbol
end

ADVIOptions(; maxiter=500, n_mc=8, σ_floor=1e-4, verbose=false, gradient=:Zygote) =
    ADVIOptions(maxiter, n_mc, σ_floor, verbose, gradient)

"""
    run_advi(data, model, rng, options=ADVIOptions(); kwargs...)

Mean-field Gaussian variational inference: maximize the ELBO using `Optim.LBFGS` with
Zygote gradients by default; set `options.gradient=:finite` for finite differences.

`rng` must be an `AbstractRNG`.

Returns `(; μ, σ, optimization, initial_θ)`.
"""
function run_advi(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    rng,
    options::ADVIOptions=ADVIOptions();
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
)
    ad = something(ad_likelihood, data isa RNACountData)
    θ0 = Vector{Float64}(get_param(model))
    D = length(θ0)
    εs = [randn(rng, D) for _ in 1:options.n_mc]

    neg_elbo(x::AbstractVector{T}) where {T<:Real} = _neg_elbo_internal(
        x, εs, data, model, steady_state_solver, ad, D, options.σ_floor,
    )

    x0 = Vector{Float64}(undef, 2D)
    x0[1:D] .= θ0
    x0[D+1:2D] .= 0.0

    opts = Optim.Options(iterations=options.maxiter, show_trace=options.verbose)
    result = if options.gradient === :finite
        Optim.optimize(neg_elbo, x0, LBFGS(), opts; autodiff=:finite)
    else
        function neg_elbo_grad!(g, x)
            gx = Zygote.gradient(neg_elbo, x)[1]
            g .= gx
            return g
        end
        od = OnceDifferentiable(neg_elbo, neg_elbo_grad!, x0)
        Optim.optimize(od, x0, LBFGS(), opts)
    end

    xmin = Optim.minimizer(result)
    μ = xmin[1:D]
    s_raw = xmin[D+1:2D]
    σ = _softplus.(s_raw) .+ options.σ_floor
    return (μ=μ, σ=σ, optimization=result, initial_θ=θ0)
end

function _neg_elbo_internal(
    x::AbstractVector{T},
    εs::Vector{Vector{Float64}},
    data,
    model::AbstractGeneTransitionModel,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
    D::Int,
    σ_floor::Float64,
) where {T<:Real}
    μ = view(x, 1:D)
    s_raw = view(x, D+1:2D)
    σ = _softplus.(s_raw) .+ σ_floor
    L = zero(T)
    for ε in εs
        η = μ .+ σ .* ε
        l = try
            logprior(η, model) + _ll_first(η, data, model, steady_state_solver, ad_likelihood)
        catch
            T(-Inf)
        end
        L = L + l
    end
    L = L / length(εs)
    ent = sum(log.(σ)) + (T(D) / 2) * (log(2 * T(π)) + one(T))
    elbo = L + ent
    return -elbo
end

# --- NUTS posterior pipeline (same outputs as `run_mh`) ---

"""Log-likelihood scalar and pointwise log predictions for WAIC (matches `run_mh` / `sample`)."""
function _loglikelihood_predictions(
    θ::AbstractVector,
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
)
    if ad_likelihood && data isa RNACountData
        return loglikelihood_ad(θ, data, model; steady_state_solver=steady_state_solver)
    end
    return loglikelihood(θ, data, model; steady_state_solver=steady_state_solver)
end

"""
    _nuts_samples_to_fit(param, data, model, n_adapts, steady_state_solver, ad_likelihood)

Build a [`Fit`](@ref) from a `d × n` posterior sample matrix (transformed parameters, same layout as `run_mh`).
"""
function _nuts_samples_to_fit(
    param::AbstractMatrix{Float64},
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    n_adapts::Int,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
)
    n = size(param, 2)
    n >= 1 || throw(ArgumentError("need at least one posterior sample"))
    ll = Vector{Float64}(undef, n)
    _, logpred0 = _loglikelihood_predictions(param[:, 1], data, model, steady_state_solver, ad_likelihood)
    pwaic = (0, log.(max.(logpred0, eps(Float64))), zeros(length(logpred0)))
    lppd = fill(-Inf, length(logpred0))
    for s in 1:n
        ll[s], logpred = _loglikelihood_predictions(param[:, s], data, model, steady_state_solver, ad_likelihood)
        lppd, pwaic = update_waic(lppd, pwaic, logpred)
    end
    pwaic = n > 1 ? pwaic[3] / (n - 1) : pwaic[3]
    lppd = lppd .- log(n)
    imax = argmax(ll)
    parml = param[:, imax]
    llml = ll[imax]
    prior = logprior(param[:, end], model)
    accept = n
    total = n_adapts + n
    return Fit(param, ll, parml, llml, lppd, pwaic, prior, accept, total)
end

"""
    run_nuts_fit(data, model, options=NUTSOptions(); rng=Random.default_rng(), kwargs...)

Hamiltonian Monte Carlo (NUTS) with the same **return convention as [`run_mh`](@ref)**:
`fits::Fit`, `stats::Stats`, `measures::Measures`, plus a fourth value `nuts_info` with sampler diagnostics.

Samples are in **transformed** parameter space (`get_param` / `get_rates`).

# Arguments
- `data`, `model`: experimental data and model (same as `run_mh`).
- `options`: [`NUTSOptions`](@ref) (`n_samples`, `n_adapts`, `δ`, `gradient`, `fd_ε`, …).
- `rng`: random number generator (default `Random.default_rng()`).
- `steady_state_solver`, `ad_likelihood`: forwarded to the log posterior (see [`run_nuts`](@ref)).

# Returns
- `fits`: [`Fit`](@ref) — `param` is `d × n_samples`; `ll` are **log-likelihood** values (same convention as MH); `accept`/`total` summarize retained draws vs. adaptation + sampling steps.
- `stats`: [`compute_stats`](@ref) on posterior samples.
- `measures`: WAIC (from posterior predictive accumulation), R-hat, ESS, Geweke, MCSE (same structs as `run_mh`).
- `nuts_info`: `NamedTuple` with `nuts_stats` (AdvancedHMC return), `initial_θ`, and `samples_matrix` (`d × n`).

# Notes
- Single-chain R-hat uses the same split-half heuristic as [`compute_rhat`](@ref) (needs enough samples).
- For gradients, prefer `RNACountData` + default Zygote when the AD likelihood path is valid; otherwise use `NUTSOptions(gradient=:finite)`.

See also: [`run_nuts`](@ref) for raw samples without `Fit`/`Stats` construction.
"""
function run_nuts_fit(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    options::NUTSOptions=NUTSOptions();
    rng::Random.AbstractRNG=Random.default_rng(),
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
)
    ad = something(ad_likelihood, data isa RNACountData)
    nt = run_nuts(data, model, rng, options; steady_state_solver=steady_state_solver, ad_likelihood=ad_likelihood)
    param = Matrix{Float64}(nt.samples)
    fits = _nuts_samples_to_fit(param, data, model, options.n_adapts, steady_state_solver, ad)
    stats = compute_stats(fits.param, model)
    rhat = vec(compute_rhat([fits]))
    ess, geweke, mcse = compute_measures(fits)
    waic = compute_waic(fits.lppd, fits.pwaic, data)
    measures = Measures(waic, rhat, ess, geweke, mcse)
    nuts_info = (
        nuts_stats=nt.nuts_stats,
        initial_θ=nt.initial_θ,
        samples_matrix=param,
    )
    return fits, stats, measures, nuts_info
end

"""Alias for [`run_nuts_fit`](@ref) (same signature and return values)."""
function run_NUTS(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    options::NUTSOptions=NUTSOptions();
    rng::Random.AbstractRNG=Random.default_rng(),
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
)
    return run_nuts_fit(data, model, options; rng=rng, steady_state_solver=steady_state_solver, ad_likelihood=ad_likelihood)
end
