# This file is part of StochasticGene.jl

"""Abstract category for optional penalties, separate from priors and data likelihoods."""
abstract type AbstractRegularizer end

"""
    JointResidenceRegularizer(target; covariance, strength=1.0, units=(1, 2))

Penalize stationary reporter ON/OFF occupancy for two coupled units. `target`
and the rows/columns of `covariance` use order `(ON_ON, ON_OFF, OFF_ON, OFF_OFF)`.
ON means at least one reporter, not promoter activation. Unit order is explicit.
The penalty is `strength/2 * delta' * pinv(covariance) * delta`; covariance
may be singular because the fractions sum to one. Supply uncertainty of the
estimated fractions (e.g. paired-trace bootstrap covariance), not frame counts.

Pass one or more penalties as `fit(; regularizers=(regularizer, ...))`.
Initially supported only for coupled MH trace fits with one shared stationary
distribution. This is a composite objective when targets come from the fitted
traces; it is not additional independent data or a trace-centered statistic.
`strength=0` disables this term. No regularization is enabled by default.
"""
struct JointResidenceRegularizer <: AbstractRegularizer
    target::NTuple{4,Float64}
    covariance::Matrix{Float64}
    strength::Float64
    units::Tuple{Int,Int}
    function JointResidenceRegularizer(target; covariance, strength=1.0, units=(1, 2))
        length(target) == 4 || throw(ArgumentError("joint residence target needs four fractions"))
        t = Tuple(Float64.(target))
        all(x -> isfinite(x) && 0 <= x <= 1, t) && isapprox(sum(t), 1; atol=1e-10, rtol=0) ||
            throw(ArgumentError("joint residence fractions must be finite, nonnegative and sum to one"))
        c = Matrix{Float64}(covariance)
        size(c) == (4, 4) && all(isfinite, c) && isapprox(c, c'; atol=1e-14, rtol=1e-10) ||
            throw(ArgumentError("joint residence covariance must be a finite symmetric 4x4 matrix"))
        e = eigvals(Symmetric(c))
        maximum(e) > 0 && minimum(e) >= -1e-10 * maximum(e) ||
            throw(ArgumentError("joint residence covariance must be positive semidefinite and nonzero"))
        s = Float64(strength)
        isfinite(s) && s >= 0 || throw(ArgumentError("regularizer strength must be finite and nonnegative"))
        length(units) == 2 && all(u -> u isa Integer && u > 0, units) && units[1] != units[2] ||
            throw(ArgumentError("units must be two distinct positive unit indices"))
        new(t, c, s, (Int(units[1]), Int(units[2])))
    end
end

# Keep generated fitscripts reproducible using the validated public constructor.
function Base.show(io::IO, r::JointResidenceRegularizer)
    print(io, "StochasticGene.JointResidenceRegularizer(")
    show(io, r.target)
    print(io, "; covariance="); show(io, r.covariance)
    print(io, ", strength=", r.strength, ", units=", r.units, ")")
end

function _normalize_regularizers(regularizers)
    terms = regularizers === nothing ? () : regularizers isa AbstractRegularizer ? (regularizers,) : Tuple(regularizers)
    all(r -> r isa AbstractRegularizer, terms) || throw(ArgumentError("regularizers must contain AbstractRegularizer objects"))
    Tuple(r for r in terms if !(r isa JointResidenceRegularizer && r.strength == 0))
end

function _prepare_regularizer(r::JointResidenceRegularizer, data, model)
    data isa AbstractTraceData && model isa AbstractGRSMmodel && model.G isa Tuple &&
        model.components isa Union{TCoupledComponents,TCoupledFullComponents} && !hastrait(model, :grid) ||
        throw(ArgumentError("joint residence regularization requires a coupled trace model without a grid"))
    if hastrait(model, :hierarchical)
        model.method isa Tuple && model.method[2] === true ||
            throw(ArgumentError("joint residence regularization requires shared hierarchical kinetics"))
    end
    maximum(r.units) <= length(model.reporter) || throw(ArgumentError("regularizer unit index out of range"))
    all(u -> model.R[u] > 0, r.units) || throw(ArgumentError("joint residence ON requires reporter-bearing units (R > 0)"))
    a = model.reporter[r.units[1]].per_state .> 0
    b = model.reporter[r.units[2]].per_state .> 0
    length(a) == length(b) || throw(DimensionMismatch("joint reporter state counts"))
    masks = Float64.(permutedims(hcat(a .& b, a .& .!b, .!a .& b, .!a .& .!b)))
    e = eigen(Symmetric(r.covariance))
    threshold = 1e-10 * maximum(e.values)
    whitener = Diagonal([v > threshold ? inv(sqrt(v)) : 0.0 for v in e.values]) * e.vectors'
    (spec=r, masks=masks, whitener=whitener)
end

mutable struct RegularizerContext
    prepared::Tuple
    predictions::Vector{Vector{Float64}}
    penalties::Vector{Float64}
    seen::Int
    active::Bool
end

function _with_regularizers(f, terms, data, model)
    isempty(terms) && return f()
    prepared = map(r -> _prepare_regularizer(r, data, model), terms)
    context = RegularizerContext(prepared, [zeros(4) for _ in terms], zeros(length(terms)), 0, false)
    task_local_storage(f, :stochasticgene_regularizers, context)
end

function _record_regularizer_stationary!(p0)
    context = get(task_local_storage(), :stochasticgene_regularizers, nothing)
    context === nothing && return nothing
    context.active || return nothing
    context.seen += 1
    context.seen == 1 || throw(ArgumentError("regularizer expected one shared stationary distribution"))
    for i in eachindex(context.prepared)
        term = context.prepared[i]
        q = context.predictions[i]
        mul!(q, term.masks, p0)
        delta = q .- collect(term.spec.target)
        context.penalties[i] = term.spec.strength / 2 * sum(abs2, term.whitener * delta)
    end
    nothing
end

function _evaluate_with_regularizers(f, terms)
    isempty(terms) && return (f()..., 0.0)
    context = get(task_local_storage(), :stochasticgene_regularizers, nothing)
    context === nothing && throw(ArgumentError("regularized evaluation requires a regularizer context"))
    context.seen = 0
    fill!(context.penalties, 0)
    context.active = true
    ll, pointwise = try
        f()
    finally
        context.active = false
    end
    context.seen == 1 || throw(ArgumentError("likelihood path did not supply a shared stationary distribution"))
    penalty = sum(context.penalties)
    isfinite(penalty) && penalty >= 0 || throw(ArgumentError("invalid regularization penalty"))
    ll, pointwise, penalty
end

"""
    regularization_report(param, data, model; regularizers=(), likelihood_executor=:fast)

Evaluate raw likelihood, prior, regularization penalties and stationary fractions
without sampling. `objective = loglikelihood + logprior - penalty`. Pointwise
likelihoods remain unpenalized. Does not refit or modify rates or source data.
"""
function regularization_report(param, data, model; regularizers=(), likelihood_executor=HMM_STACK_MH)
    terms = _normalize_regularizers(regularizers)
    _with_regularizers(terms, data, model) do
        ll, pointwise, penalty = _evaluate_with_regularizers(terms) do
            _mh_loglikelihood(param, data, model, likelihood_executor)
        end
        context = get(task_local_storage(), :stochasticgene_regularizers, nothing)
        prior = logprior(param, model)
        (loglikelihood=ll, pointwise=pointwise, logprior=prior, penalty=penalty,
            objective=ll + prior - penalty,
            fractions=isempty(terms) ? [] : deepcopy(context.predictions),
            penalties=isempty(terms) ? Float64[] : copy(context.penalties))
    end
end

function _write_regularization_report(folder, data, model, fits, stats, options)
    rows = DataFrame()
    for (label, param) in (("ml", fits.parml), ("median", transform_rates(vec(stats.medparam), model)))
        report = regularization_report(param, data, model;
            regularizers=options.regularizers, likelihood_executor=options.likelihood_executor)
        for i in eachindex(options.regularizers)
            r = options.regularizers[i]
            q = report.fractions[i]
            push!(rows, (point=label, term=i, kind=string(nameof(typeof(r))),
                unit1=r.units[1], unit2=r.units[2], strength=r.strength,
                ON_ON=q[1], ON_OFF=q[2], OFF_ON=q[3], OFF_OFF=q[4],
                penalty=report.penalties[i], total_penalty=report.penalty,
                loglikelihood=report.loglikelihood, logprior=report.logprior, objective=report.objective))
        end
    end
    name = _current_name_override[] === nothing ? filename(data, model) : _current_name_override[]
    file = splitext("regularization" * name)[1] * ".csv"
    CSV.write(joinpath(folder, file), rows)
    @info "Regularized MH fit: reported likelihood/AIC are unpenalized; WAIC uses trace likelihoods under a regularized sampling distribution" file=file
end
