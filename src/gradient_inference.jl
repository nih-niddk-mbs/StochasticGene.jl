# gradient_inference.jl — NUTS (AdvancedHMC) and mean-field ADVI on the same
# transformed parameter space as Metropolis–Hastings (`get_param`, `logprior`).
#
# Gradients default to reverse-mode AD (Zygote via LogDensityProblemsAD / AdvancedHMC).
# Use `gradient=:finite` in `NUTSOptions` for central finite differences; use `gradient=:ForwardDiff`
# for forward-mode AD (often good for **few** fitted parameters; not FD).
# ADVI: `gradient=:finite` or `gradient=:ForwardDiff` both use `ForwardDiff` on `neg_elbo` (see `run_advi`).
#
# ADVI with `gradient=:finite` uses `ForwardDiff` for ∇(neg_elbo): Optim’s own finite
# differences would evaluate the objective ~O(dim x) times per gradient (each trace HMM
# is expensive), which looks “stuck” between printed LBFGS iterations.

using ForwardDiff

"""
    check_ad_gradient_feasibility(
        data, model::AbstractGeneTransitionModel;
        forwarddiff_through_likelihood::Bool,
        zygote_through_likelihood::Bool,
        zygote_long_trace_warn_threshold::Union{Nothing,Int}=2048,
    )

Validate gradient choices **before** building NUTS/ADVI objectives. Callers set the boolean flags
according to how gradients are computed (see [`run_nuts`](@ref) / [`run_advi`](@ref)).

- **`forwarddiff_through_likelihood`**: `ForwardDiff.Dual` propagates through [`loglikelihood_ad`](@ref) on
  [`AbstractTraceData`](@ref) (e.g. NUTS `gradient=:ForwardDiff`, or ADVI with `gradient=:finite` / `:ForwardDiff`,
  which differentiate `neg_elbo` with ForwardDiff).

- **`zygote_through_likelihood`**: Zygote reverse-mode through [`loglikelihood_ad`](@ref) on trace data
  (NUTS `gradient=:Zygote`, or ADVI with `gradient=:Zygote` and `zygote_trace=true`).

Coupled **trace** models use [`kolmogorov_forward_ad`](@ref) with a `Dual`-generic matrix exponential; large
coupled state spaces can still be **slow** with ForwardDiff. This check does **not** block ForwardDiff on coupled
traces. Coupled **histogram** likelihoods use [`predictedarray`](@ref) / dwell-time and ON–OFF PDFs — separate code
paths (see the AD feasibility summary in the `likelihoods.jl` source comments near “Model loglikelihoods”).

**Warns** when `zygote_through_likelihood` on coupled traces (memory risk) or on long traces (compiler / tape issues),
unless `zygote_long_trace_warn_threshold === nothing` (skips length check).

Histogram / count likelihoods ([`AbstractHistogramData`](@ref), [`RNACountData`](@ref)) do not use this check at call sites;
they typically work with Zygote or ForwardDiff via the steady-state / [`predictedarray`](@ref) path when `ad_likelihood`
selects [`loglikelihood_ad`](@ref).
"""
function check_ad_gradient_feasibility(
    data,
    model::AbstractGeneTransitionModel;
    forwarddiff_through_likelihood::Bool,
    zygote_through_likelihood::Bool,
    zygote_long_trace_warn_threshold::Union{Nothing,Int}=2048,
)
    if zygote_through_likelihood && data isa AbstractTraceData && hastrait(model, :coupling)
        @warn "Coupled trace HMM with Zygote reverse-mode can allocate very large memory; prefer gradient=:finite, " *
            "or set hmm_checkpoint_steps for checkpointed Zygote through the HMM forward pass."
    end
    if zygote_through_likelihood && data isa AbstractTraceData && zygote_long_trace_warn_threshold !== nothing
        L = longest_trace_timesteps(data)
        if L !== nothing && L > zygote_long_trace_warn_threshold
            @warn "Long trace ($L frames > threshold $(zygote_long_trace_warn_threshold)): Zygote reverse-mode on trace HMMs " *
                "often overflows the compiler or uses huge memory; use gradient=:finite, hmm_checkpoint_steps, shorter traces, " *
                "or pass zygote_long_trace_warn_threshold=nothing to skip this warning."
        end
    end
    return nothing
end

"""
    GenePosteriorLogDensity

`LogDensityProblems.jl` wrapper for the unnormalized log posterior
`logprior(θ, model) + loglikelihood(θ, data, model)` with `θ` in **transformed**
space (same as `get_param(model)`).

# Fields
- `data`, `model`: passed through to `loglikelihood` / `logprior`
- `steady_state_solver`: forwarded to likelihood (default `:augmented` is best for downstream use)
- `ad_likelihood`: if `true`, uses `loglikelihood_ad` for `AbstractHistogramData` and `AbstractTraceData` (Zygote-friendly likelihood); otherwise `loglikelihood`
- `options`: inference options; [`CombinedData`](@ref) likelihoods receive the full struct, while legacy paths read `likelihood_executor`.
"""
struct GenePosteriorLogDensity{D,M,O}
    data::D
    model::M
    steady_state_solver::Symbol
    ad_likelihood::Bool
    options::O
end

function _ll_first(
    θ,
    data,
    model,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
    options::Options,
)
    if data isa CombinedData
        if ad_likelihood && _default_ad_likelihood(data)
            return loglikelihood_ad(θ, data, model, options; steady_state_solver=steady_state_solver)[1]
        end
        return loglikelihood(θ, data, model, options; steady_state_solver=steady_state_solver)[1]
    end
    likelihood_executor = options.likelihood_executor
    if ad_likelihood && _default_ad_likelihood(data)
        kwargs = (; steady_state_solver=steady_state_solver)
        if _trace_or_combined_for_hmm(data)
            kwargs = (; kwargs..., hmm_stack=likelihood_executor)
        end
        return loglikelihood_ad(θ, data, model; kwargs...)[1]
    end
    return loglikelihood(θ, data, model; steady_state_solver=steady_state_solver)[1]
end

@inline _default_ad_likelihood(data) =
    data isa AbstractHistogramData || data isa AbstractTraceData || data isa CombinedData

@inline _trace_or_combined_for_hmm(data) = data isa AbstractTraceData || data isa CombinedData

"""
    logposterior(θ, data, model; steady_state_solver=:augmented, ad_likelihood=nothing)

Unnormalized log posterior in transformed coordinates. `ad_likelihood` defaults to
`true` for `AbstractHistogramData`, `AbstractTraceData`, or [`CombinedData`](@ref) (uses `loglikelihood_ad`), else `false`.
"""
function logposterior(
    θ::AbstractVector,
    data,
    model::AbstractGeneTransitionModel;
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    likelihood_executor::Union{Nothing,Symbol}=nothing,
)
    ad = something(ad_likelihood, _default_ad_likelihood(data))
    lik_ad = something(likelihood_executor, HMM_STACK_AD)
    options = NUTSOptions(; likelihood_executor=lik_ad)
    return logprior(θ, model) + _ll_first(θ, data, model, steady_state_solver, ad, options)
end

function LogDensityProblems.dimension(p::GenePosteriorLogDensity)
    return length(get_param(p.model))
end

function LogDensityProblems.logdensity(p::GenePosteriorLogDensity, θ)
    return logprior(θ, p.model) + _ll_first(θ, p.data, p.model, p.steady_state_solver, p.ad_likelihood, p.options)
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
    run_nuts(data, model, rng, options=NUTSOptions(); kwargs...)

Hamiltonian Monte Carlo with the No-U-Turn sampler (NUTS) and diagonal mass matrix
adaptation (Stan-style). Gradients of the log posterior follow `options.gradient` (see [`NUTSOptions`](@ref)); the struct default is **ForwardDiff** (forward-mode AD, efficient for few parameters).
For **trace** data, `:finite` (central differences) is preferred; [`fit`](@ref) with `estimation=INFERENCE_NUTS` auto-selects `:finite` for `AbstractTraceData` unless overridden with `nuts_gradient=:ForwardDiff` or `:Zygote`.

`rng` must be an `AbstractRNG` (e.g. `using Random; Random.default_rng()`).

When **`options.progress`** is `true` (the default on [`NUTSOptions`](@ref)), AdvancedHMC shows a **ProgressMeter** sampling bar for the total number of draws (`n_samples`); multi-chain [`run_inference`](@ref) disables progress automatically to avoid garbled output.

Keyword arguments:
- `steady_state_solver`: passed to likelihood (default `:augmented`)
- `ad_likelihood`: override automatic choice (`nothing` → use `loglikelihood_ad` for `AbstractHistogramData` and `AbstractTraceData`)
- `hmm_checkpoint_steps`: if a positive integer, Zygote reverse-mode through trace HMMs uses gradient checkpointing in chunks of this many frames (see [`set_hmm_zygote_checkpoint_steps!`](@ref)); `nothing` leaves the global setting unchanged

Returns a `NamedTuple` with `samples` (`d × n` matrix), `sample_vectors`, `nuts_stats`, and `initial_θ`.
"""
function run_nuts(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    rng,
    options::NUTSOptions=NUTSOptions();
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    ck_outer = coalesce(hmm_checkpoint_steps, options.gradient_checkpoint_length)
    with_hmm_zygote_checkpoint(ck_outer) do
        ad = something(ad_likelihood, _default_ad_likelihood(data))
        ℓ = GenePosteriorLogDensity(data, model, steady_state_solver, ad, options)
        θ0 = Vector{Float64}(get_param(model))
        D = length(θ0)
        D == LogDensityProblems.dimension(ℓ) || throw(DimensionMismatch("parameter dimension mismatch"))

        gopt = options.gradient
        gopt in (:Zygote, :finite, :ForwardDiff) ||
            throw(ArgumentError("NUTSOptions.gradient must be :Zygote, :finite, or :ForwardDiff, got $(repr(gopt))"))

        check_ad_gradient_feasibility(
            data, model;
            forwarddiff_through_likelihood=(gopt === :ForwardDiff && ad && _trace_or_combined_for_hmm(data)),
            zygote_through_likelihood=(gopt === :Zygote && ad && _trace_or_combined_for_hmm(data)),
        )

        metric = DiagEuclideanMetric(D)
        ham = if options.gradient === :finite
            ℓπ = _make_ℓπ(ℓ)
            ∂ℓπ∂θ = _make_finitediff_∂ℓπ∂θ(ℓ, Float64(options.fd_ε))
            Hamiltonian(metric, ℓπ, ∂ℓπ∂θ)
        elseif options.gradient === :ForwardDiff
            Hamiltonian(metric, ℓ, Val(:ForwardDiff))
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
end

@inline function _softplus(x::Real)
    return x > 0 ? x + log1p(exp(-x)) : log1p(exp(x))
end


"""
    run_advi(data, model, rng, options=ADVIOptions(); kwargs...)

Mean-field Gaussian variational inference: maximize the ELBO using `Optim.LBFGS` with
Zygote gradients by default; set `options.gradient=:finite` or `options.gradient=:ForwardDiff` for `ForwardDiff` on `neg_elbo` (same implementation).

**Trace likelihoods — two tracks:** [`loglikelihood`](@ref) uses [`HMM_STACK_MH`](@ref) (MH / MCMC; in-place `forward`, `kolmogorov_forward`).
[`loglikelihood_ad`](@ref) uses [`HMM_STACK_AD`](@ref) (`forward_ad`, `kolmogorov_forward_ad`, …) for NUTS/ADVI/Zygote.
The package loads `SciMLSensitivity` for adjoints through `solve`. Without it, use `gradient=:finite` or `gradient=:ForwardDiff`.

`rng` must be an `AbstractRNG`.

Returns `(; μ, σ, optimization, initial_θ)`.

The ELBO uses Zygote-differentiable `logprior` + likelihood (`loglikelihood_ad` on trace data when `ad_likelihood` is true) only (no `try/catch` around them),
so invalid variational draws must not throw; they should yield `-Inf`/`NaN` log-values if unsupported.

# Trace data and `gradient=:Zygote`

Reverse-mode AD through a long trace HMM (`loglikelihood_ad`) builds a huge Zygote tape and can **overflow the compiler stack** (`_pullback_generator` / huge `NTuple` types).
Use **`hmm_checkpoint_steps`** on this function, [`run_nuts`](@ref), [`ll_hmm_trace`](@ref), or [`loglikelihood_ad`](@ref) for trace data, or [`set_hmm_zygote_checkpoint_steps!`](@ref), to enable gradient checkpointing on the HMM forward pass (Zygote only). Otherwise prefer [`ADVIOptions`](@ref) field `gradient_checkpoint_length` (or run-spec `gradient_checkpoint_length` / legacy `hmm_checkpoint_steps`).
By default, **`AbstractTraceData` forces the ForwardDiff ELBO gradient path** (`gradient=:finite` after override; same as explicit `gradient=:ForwardDiff`) even if `ADVIOptions` requests Zygote. Set `zygote_trace=true` to attempt Zygote anyway (short traces only).

# Keywords
- `zygote_trace`: if `true` and `data isa AbstractTraceData`, use Zygote when `options.gradient === :Zygote` (default `false`).
- `hmm_checkpoint_steps`: optional per-call override for Zygote checkpoint chunk size; `nothing` uses `options.gradient_checkpoint_length` then the process-global setting.
"""
function run_advi(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    rng,
    options::ADVIOptions=ADVIOptions();
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    zygote_trace::Bool=false,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    ad = something(ad_likelihood, _default_ad_likelihood(data))
    options.gradient in (:Zygote, :finite, :ForwardDiff) ||
        throw(ArgumentError("ADVIOptions.gradient must be :Zygote, :finite, or :ForwardDiff, got $(repr(options.gradient))"))
    options.n_mc >= 1 || throw(ArgumentError("ADVIOptions.n_mc must be ≥ 1, got $(options.n_mc)"))
    options.maxiter >= 1 || throw(ArgumentError("ADVIOptions.maxiter must be ≥ 1, got $(options.maxiter)"))

    options_eff = if !zygote_trace && _trace_or_combined_for_hmm(data) && options.gradient === :Zygote
        @info "ADVI: using finite-difference gradients for trace data (Zygote through long HMM likelihoods is not supported; set `zygote_trace=true` to try)."
        ADVIOptions(;
            maxiter=options.maxiter,
            n_mc=options.n_mc,
            σ_floor=options.σ_floor,
            init_s_raw=options.init_s_raw,
            verbose=options.verbose,
            gradient=:finite,
            time_limit=options.time_limit,
            device=options.device,
            parallel=options.parallel,
            likelihood_executor=options.likelihood_executor,
            gradient_checkpoint_length=options.gradient_checkpoint_length,
        )
    else
        options
    end

    check_ad_gradient_feasibility(
        data, model;
        forwarddiff_through_likelihood=(
            ad && _trace_or_combined_for_hmm(data) &&
            (options_eff.gradient === :finite || options_eff.gradient === :ForwardDiff)
        ),
        zygote_through_likelihood=(ad && _trace_or_combined_for_hmm(data) && options_eff.gradient === :Zygote),
    )

    θ0 = Vector{Float64}(get_param(model))
    D = length(θ0)
    εs = [randn(rng, D) for _ in 1:options_eff.n_mc]
    ck_outer = coalesce(hmm_checkpoint_steps, options_eff.gradient_checkpoint_length)
    likelihood_options = options_eff

    return with_hmm_zygote_checkpoint(ck_outer) do
        neg_elbo(x::AbstractVector{T}) where {T<:Real} = _neg_elbo_internal(
            x, εs, data, model, steady_state_solver, ad, D, options_eff.σ_floor, likelihood_options,
        )

        x0 = Vector{Float64}(undef, 2D)
        x0[1:D] .= θ0
        x0[D+1:2D] .= options_eff.init_s_raw

        opts = if options_eff.time_limit === nothing
            Optim.Options(iterations=options_eff.maxiter, show_trace=options_eff.verbose)
        else
            Optim.Options(iterations=options_eff.maxiter, show_trace=options_eff.verbose, time_limit=options_eff.time_limit)
        end
        lbfgs_method = LBFGS(linesearch=Optim.LineSearches.BackTracking())
        result = if options_eff.gradient === :finite || options_eff.gradient === :ForwardDiff
            if options_eff.verbose && data isa AbstractTraceData
                @info "ADVI: gradients use ForwardDiff on neg_elbo (not Optim finite-differences); " *
                    "LBFGS may still spend a long time on early line searches — try smaller `n_mc` (e.g. 2) if too slow."
            end
            function neg_elbo_grad_fd!(g::AbstractVector{Float64}, x::AbstractVector{Float64})
                gx = ForwardDiff.gradient(neg_elbo, x)
                g .= gx
                return g
            end
            od = OnceDifferentiable(neg_elbo, neg_elbo_grad_fd!, x0)
            Optim.optimize(od, x0, lbfgs_method, opts)
        else
            # One Zygote tape over the whole MC loop retains every forward intermediate → huge RAM.
            # ∇(-mean_i ℓ_i - H) = -mean_i ∇ℓ_i - ∇H, so accumulate per-ε gradients separately.
            function neg_elbo_grad!(g, x)
                _neg_elbo_grad_zygote!(
                    g, x, εs, data, model, steady_state_solver, ad, D, options_eff.σ_floor, likelihood_options,
                )
                return g
            end
            od = OnceDifferentiable(neg_elbo, neg_elbo_grad!, x0)
            Optim.optimize(od, x0, lbfgs_method, opts)
        end

        xmin = Optim.minimizer(result)
        μ = xmin[1:D]
        s_raw = xmin[D+1:2D]
        σ = _softplus.(s_raw) .+ options_eff.σ_floor
        return (μ=μ, σ=σ, optimization=result, initial_θ=θ0)
    end
end

function _elbo_single_lp(
    x::AbstractVector{T},
    ε::AbstractVector{Float64},
    data,
    model::AbstractGeneTransitionModel,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
    D::Int,
    σ_floor::Float64,
    options::Options,
) where {T<:Real}
    μ = view(x, 1:D)
    s_raw = view(x, D+1:2D)
    σ = _softplus.(s_raw) .+ σ_floor
    η = μ .+ σ .* ε
    # No try/catch here: Zygote cannot differentiate through `catch` (see ADVI + `gradient=:Zygote`).
    lp = logprior(η, model)
    isfinite(lp) || return oftype(lp, -Inf)
    ll = _ll_first(η, data, model, steady_state_solver, ad_likelihood, options)
    isfinite(ll) || return oftype(lp + ll, -Inf)
    return lp + ll
end

function _elbo_entropy(
    x::AbstractVector{T},
    D::Int,
    σ_floor::Float64,
) where {T<:Real}
    s_raw = view(x, D+1:2D)
    σ = _softplus.(s_raw) .+ σ_floor
    sum(log.(σ)) + (T(D) / 2) * (log(2 * T(π)) + one(T))
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
    options::Options,
) where {T<:Real}
    L = zero(T)
    for ε in εs
        L += _elbo_single_lp(x, ε, data, model, steady_state_solver, ad_likelihood, D, σ_floor, options)
    end
    L = L / length(εs)
    ent = _elbo_entropy(x, D, σ_floor)
    elbo = L + ent
    return -elbo
end

function _neg_elbo_grad_zygote!(
    g::AbstractVector{Float64},
    x::AbstractVector{Float64},
    εs::Vector{Vector{Float64}},
    data,
    model::AbstractGeneTransitionModel,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
    D::Int,
    σ_floor::Float64,
    options::Options,
)
    fill!(g, 0.0)
    n = length(εs)
    n >= 1 || throw(ArgumentError("n_mc must be ≥ 1"))
    invn = 1.0 / n
    for ε in εs
        gs = Zygote.gradient(
            x′ -> _elbo_single_lp(x′, ε, data, model, steady_state_solver, ad_likelihood, D, σ_floor, options),
            x,
        )
        gx = gs[1]
        gx === nothing && error("Zygote returned no gradient for ELBO MC term; check logprior and likelihood.")
        g .+= gx
    end
    g .*= -invn
    ge = Zygote.gradient(x′ -> -_elbo_entropy(x′, D, σ_floor), x)[1]
    ge === nothing && error("Zygote returned no gradient for ELBO entropy term.")
    g .+= ge
    return g
end

# --- NUTS posterior pipeline (same outputs as `run_mh`) ---

"""Log-likelihood scalar and pointwise log predictions for WAIC (matches `run_mh` / `sample`)."""
function _loglikelihood_predictions(
    θ::AbstractVector,
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
    options::Options,
)
    if data isa CombinedData
        if ad_likelihood && _default_ad_likelihood(data)
            return loglikelihood_ad(θ, data, model, options; steady_state_solver=steady_state_solver)
        end
        return loglikelihood(θ, data, model, options; steady_state_solver=steady_state_solver)
    end
    likelihood_executor = options.likelihood_executor
    if ad_likelihood && _default_ad_likelihood(data)
        kwargs = (; steady_state_solver=steady_state_solver)
        if _trace_or_combined_for_hmm(data)
            kwargs = (; kwargs..., hmm_stack=likelihood_executor)
        end
        return loglikelihood_ad(θ, data, model; kwargs...)
    end
    return loglikelihood(θ, data, model; steady_state_solver=steady_state_solver)
end

"""
    _nuts_samples_to_fit(param, data, model, n_adapts, steady_state_solver, ad_likelihood, likelihood_executor)

Build a [`Fit`](@ref) from a `d × n` posterior sample matrix (transformed parameters, same layout as `run_mh`).
"""
function _nuts_samples_to_fit(
    param::AbstractMatrix{Float64},
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    n_adapts::Int,
    steady_state_solver::Symbol,
    ad_likelihood::Bool,
    options::Options,
)
    n = size(param, 2)
    n >= 1 || throw(ArgumentError("need at least one posterior sample"))
    ll = Vector{Float64}(undef, n)
    _, logpred0 = _loglikelihood_predictions(param[:, 1], data, model, steady_state_solver, ad_likelihood, options)
    pwaic = (0, log.(max.(logpred0, eps(Float64))), zeros(length(logpred0)))
    lppd = fill(-Inf, length(logpred0))
    for s in 1:n
        ll[s], logpred = _loglikelihood_predictions(param[:, s], data, model, steady_state_solver, ad_likelihood, options)
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
- `steady_state_solver`, `ad_likelihood`, `hmm_checkpoint_steps`: forwarded to [`run_nuts`](@ref).

# Returns
- `fits`: [`Fit`](@ref) — `param` is `d × n_samples`; `ll` are **log-likelihood** values (same convention as MH); `accept`/`total` summarize retained draws vs. adaptation + sampling steps.
- `stats`: [`compute_stats`](@ref) on posterior samples.
- `measures`: WAIC (from posterior predictive accumulation), R-hat, ESS, Geweke, MCSE (same structs as `run_mh`).
- `nuts_info`: `NamedTuple` with `nuts_stats` (AdvancedHMC return), `initial_θ`, and `samples_matrix` (`d × n`).

# Notes
- Single-chain R-hat uses the same split-half heuristic as [`compute_rhat`](@ref) (needs enough samples).
- For gradients: prefer forward-mode AD (`:ForwardDiff`, the new default) for typical gene switching models with **few** parameters. For trace HMMs, `:finite` (central differences) is auto-selected by [`fit`](@ref); override with `nuts_gradient=:Zygote` for many parameters or if memory-intensive Zygote reverse mode is desired.

See also: [`run_nuts`](@ref) for raw samples without `Fit`/`Stats` construction.
"""
function run_nuts_fit(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    options::NUTSOptions=NUTSOptions();
    rng::Random.AbstractRNG=Random.default_rng(),
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    ad = something(ad_likelihood, _default_ad_likelihood(data))
    nt = run_nuts(
        data, model, rng, options;
        steady_state_solver=steady_state_solver,
        ad_likelihood=ad_likelihood,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )
    param = Matrix{Float64}(nt.samples)
    fits = _nuts_samples_to_fit(param, data, model, options.n_adapts, steady_state_solver, ad, options)
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
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    return run_nuts_fit(
        data, model, options;
        rng=rng,
        steady_state_solver=steady_state_solver,
        ad_likelihood=ad_likelihood,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )
end

function _advi_posterior_stats(
    μ::AbstractVector{Float64},
    σ::AbstractVector{Float64},
    model::AbstractGeneTransitionModel,
    rng::Random.AbstractRNG,
    posterior_samples::Int,
)
    posterior_samples >= 2 || return compute_stats(reshape(μ, :, 1), model)
    D = length(μ)
    Z = randn(rng, D, posterior_samples)
    param_draws = reshape(μ, :, 1) .+ reshape(σ, :, 1) .* Z
    return compute_stats(Matrix{Float64}(param_draws), model)
end

"""
    run_advi_fit(data, model, options=ADVIOptions(); rng=Random.default_rng(), kwargs...)

Mean-field ADVI with the same **return convention as [`run_nuts_fit`](@ref) and [`run_mh`](@ref)**:
`fits::Fit`, `stats::Stats`, `measures::Measures`, plus a fourth value `advi_info` with variational diagnostics.

Posterior summary uses the **variational mean** (μ) as the point estimate and the variational standard deviation (σ)
for interval reporting. Since ADVI returns a single point estimate, `Fit.param` is 1 column, and `R-hat`, `ESS`,
`Geweke` are computed as single-sample proxies (values will indicate insufficient samples).

# Arguments
- `data`, `model`: experimental data and model (same as `run_mh`).
- `options`: [`ADVIOptions`](@ref) (`maxiter`, `n_mc`, `σ_floor`, `gradient`, `time_limit`, …).
- `rng`: random number generator (default `Random.default_rng()`).
- `steady_state_solver`, `ad_likelihood`, `zygote_trace`, `hmm_checkpoint_steps`: forwarded to [`run_advi`](@ref).

# Returns
- `fits`: [`Fit`](@ref) — `param` is `d × 1` (the variational mean μ); `ll` contains the single log-likelihood; `accept`/`total` are both 1 (pseudo-sample count).
- `stats`: [`compute_stats`](@ref) on the variational mean.
- `measures`: R-hat (single sample; expect values ~1.0 or NA), ESS (also single-sample proxy), WAIC from μ, etc.
- `advi_info`: `NamedTuple` with `μ` (mean), `σ` (std), `optimization` (Optim.OptimizationResults), `initial_θ`, and `vb_elbo` (neg_elbo × -1).

# Notes
- ADVI provides a point estimate (μ) and uncertainty (σ), but does **not** produce a sample-based posterior. Results are for **variational approximation only** — posteriors are assumed Gaussian.
- For interpretation and model comparison, use WAIC from the variational mean; R-hat and ESS are formal-only (not valid for single-sample output).

See also: [`run_advi`](@ref) for raw variational output, [`run_nuts_fit`](@ref) for sample-based NUTS results.
"""
function run_advi_fit(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    options::ADVIOptions=ADVIOptions();
    rng::Random.AbstractRNG=Random.default_rng(),
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    zygote_trace::Bool=false,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
    posterior_samples::Int=1000,
)
    ad = something(ad_likelihood, _default_ad_likelihood(data))
    vb_out = run_advi(
        data, model, rng, options;
        steady_state_solver=steady_state_solver,
        ad_likelihood=ad_likelihood,
        zygote_trace=zygote_trace,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )

    μ_matrix = reshape(vb_out.μ, :, 1)
    fits = _nuts_samples_to_fit(μ_matrix, data, model, 0, steady_state_solver, ad, options)
    stats = _advi_posterior_stats(vb_out.μ, vb_out.σ, model, rng, posterior_samples)
    nparam = size(μ_matrix, 1)
    rhat = fill(NaN, nparam)
    ess = fill(NaN, nparam)
    geweke = fill(NaN, nparam)
    mcse = fill(NaN, nparam)
    waic_val = compute_waic(fits.lppd, fits.pwaic, data)
    measures = Measures(waic_val, rhat, ess, geweke, mcse)
    
    # ELBO = -neg_elbo; store for diagnostic reporting.
    vb_elbo = -Optim.minimum(vb_out.optimization)
    
    advi_info = (
        μ=vb_out.μ,
        σ=vb_out.σ,
        optimization=vb_out.optimization,
        initial_θ=vb_out.initial_θ,
        vb_elbo=vb_elbo,
        posterior_samples=posterior_samples,
    )
    
    return fits, stats, measures, advi_info
end

"""Alias for [`run_advi_fit`](@ref) (same signature and return values)."""
function run_ADVI(
    data::AbstractExperimentalData,
    model::AbstractGeneTransitionModel,
    options::ADVIOptions=ADVIOptions();
    rng::Random.AbstractRNG=Random.default_rng(),
    steady_state_solver::Symbol=:augmented,
    ad_likelihood::Union{Nothing,Bool}=nothing,
    zygote_trace::Bool=false,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
    posterior_samples::Int=1000,
)
    return run_advi_fit(
        data, model, options;
        rng=rng,
        steady_state_solver=steady_state_solver,
        ad_likelihood=ad_likelihood,
        zygote_trace=zygote_trace,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
        posterior_samples=posterior_samples,
    )
end
