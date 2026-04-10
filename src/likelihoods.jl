# Model and Data dependent functions

"""
datahistogram(data)
Return the RNA histogram data as one vector

# Arguments
- `data::RNAData`: RNA histogram data

# Methods
- `datahistogram(data::RNAData)`
- `datahistogram(data::RNAOnOffData)`
- `datahistogram(data::RNADwellTimeData)`
- `datahistogram(data::AbstractTraceHistogramData)`
- `datahistogram(data::AbstractRNAData{Array{Float64,1}})`

# Returns
- `data.histRNA::Vector`: RNA histogram data as one vector
"""


"""
    data_histogram(data)

Returns the RNA histogram data as one vector.

# Arguments
- `data::AbstractRNAData{Array{Float64,1}}`: RNA histogram data.
- `data::RNAOnOffData`: RNA On/Off data.
- `data::AbstractTraceHistogramData`: Trace histogram data.
- `data::AbstractRNAData{Array{Array,1}}`: RNA histogram data with nested arrays.
- `data::RNADwellTimeData`: RNA dwell time data.

# Description
This function returns the RNA histogram data as one vector. It supports various types of RNA data, including RNA histogram data, RNA On/Off data, trace histogram data, and RNA dwell time data.

# Methods
- `data_histogram(data::AbstractRNAData{Array{Float64,1}})`: Returns RNA histogram data.
- `data_histogram(data::RNAOnOffData)`: Returns RNA On/Off data.
- `data_histogram(data::AbstractTraceHistogramData)`: Returns trace histogram data.
- `data_histogram(data::AbstractRNAData{Array{Array,1}})`: Returns RNA histogram data with nested arrays.
- `data_histogram(data::RNADwellTimeData)`: Returns RNA dwell time data.

# Returns
- `Vector{Float64}`: RNA histogram data as one vector.
"""
# datahistogram(data::RNAData) = data.histRNA
function datahistogram(data::AbstractRNAData{Array{Array,1}})
    v = data.histRNA[1]
    for i in 2:length(data.histRNA)
        v = vcat(v, data.histRNA[i])
    end
    return v
end
datahistogram(data::AbstractRNAData{Array{Float64,1}}) = data.histRNA
datahistogram(data::RNAOnOffData) = [data.OFF; data.ON; data.histRNA]
datahistogram(data::AbstractTraceHistogramData) = data.histRNA
datahistogram(data::RNACountData) = data.countsRNA

function datahistogram(data::RNADwellTimeData)
    v = data.histRNA
    for d in data.DwellTimes
        v = vcat(v, d)
    end
    return v
end

function datahistogram(data::DwellTimeData)
    v = data.DwellTimes[1]
    for d in data.DwellTimes[2:end]
        v = vcat(v, d)
    end
    return v
end


"""
    datapdf(data)

Returns the normalized RNA histogram data as one vector.

# Arguments
- `data::AbstractRNAData{Array{Float64,1}}`: RNA histogram data.
- `data::RNAOnOffData`: RNA On/Off data.
- `data::AbstractTraceHistogramData`: Trace histogram data.
- `data::AbstractRNAData{Array{Array,1}}`: RNA histogram data with nested arrays.
- `data::RNADwellTimeData`: RNA dwell time data.

# Description
This function normalizes the RNA histogram data and returns it as one vector. It supports various types of RNA data, including RNA histogram data, RNA On/Off data, trace histogram data, and RNA dwell time data.

# Methods
- `datapdf(data::AbstractRNAData{Array{Float64,1}})`: Normalizes and returns RNA histogram data.
- `datapdf(data::RNAOnOffData)`: Normalizes and returns RNA On/Off data.
- `datapdf(data::AbstractTraceHistogramData)`: Normalizes and returns trace histogram data.
- `datapdf(data::AbstractRNAData{Array{Array,1}})`: Normalizes and returns RNA histogram data with nested arrays.
- `datapdf(data::RNADwellTimeData)`: Normalizes and returns RNA dwell time data.

# Returns
- `Vector{Float64}`: Normalized RNA histogram data as one vector.
"""
datapdf(data::AbstractRNAData{Array{Float64,1}}) = normalize_histogram(data.histRNA)
datapdf(data::RNAOnOffData) = [normalize_histogram(data.OFF); normalize_histogram(data.ON); normalize_histogram(data.histRNA)]
datapdf(data::AbstractTraceHistogramData) = normalize_histogram(data.histRNA)

function datapdf(data::AbstractRNAData{Array{Array,1}})
    v = normalize_histogram(data.histRNA[1])
    for i in 2:length(data.histRNA)
        v = vcat(v, normalize_histogram(data.histRNA[i]))
    end
    return v
end

function datapdf(data::RNADwellTimeData)
    # Normalize each component separately to match predictedfn structure
    # (which returns [RNA steady state; dwell time distributions...])
    v = normalize_histogram(data.histRNA)
    for d in data.DwellTimes
        v = vcat(v, normalize_histogram(d))
    end
    return v
end



####### GRMSmodel

"""
1. if hierarchical trait is present then, reshape rates into matrix with rows pertaining to rate/parameter and columns pertaining to shared, hyper, or individuals
2. if coupled trait is present then divide rows into each unit, extract out coupling strength row
3. if reporter type is HMMReprter then divide unit rates into rates and noise parameters
4. if grid trait is present then extract out grid rows
"""

"""
    prepare_hyper(r, param, hierarchy::HierarchicalTrait)
"""
function prepare_hyper(r, param, hierarchy::HierarchicalTrait)
    pindividual = collect(eachcol(reshape(param[hierarchy.paramstart:end], hierarchy.nindividualparams, hierarchy.nindividuals)))
    rhyper = [r[i] for i in hierarchy.hyperindices]
    return pindividual, rhyper
end

"""
    prepare_rates(r, hierarchy::HierarchicalTrait)
"""

function prepare_rates(r, hierarchy::HierarchicalTrait)
    # rates reshaped from a vector into a vector of vectors pertaining to shared params, hyper params and individual params 
    # (shared parameters are considered to be hyper parameters without other hyper parameters (e.g. mean without variance))
    nallparams = hierarchy.nrates
    rshared = reshape(r[1:hierarchy.individualstart-1], nallparams, hierarchy.nhypersets)
    rindividual = reshape(r[hierarchy.individualstart:end], nallparams, hierarchy.nindividuals)
    rindividual[hierarchy.fittedshared, :] .= rshared[hierarchy.fittedshared, 1]
    return collect(eachcol(rshared)), collect(eachcol(rindividual))
end

"""
    prepare_rates_ad(r, hierarchy::HierarchicalTrait)

Non-mutating counterpart to [`prepare_rates(r, hierarchy)`](@ref) for hierarchical
traits: returns shared and individual rate columns without mutating `r`. Use with AD
when splitting hierarchical rate vectors.
"""
function prepare_rates_ad(r, hierarchy::HierarchicalTrait)
    nallparams = hierarchy.nrates
    rshared = reshape(r[1:hierarchy.individualstart-1], nallparams, hierarchy.nhypersets)
    rindividual = reshape(r[hierarchy.individualstart:end], nallparams, hierarchy.nindividuals)
    # Build a new matrix, column by column, without mutation
    rindividual_new = [i in hierarchy.fittedshared ? rshared[i, 1] : rindividual[i, j]
                       for i in 1:nallparams, j in 1:hierarchy.nindividuals]
    return collect(eachcol(rshared)), collect(eachcol(rindividual_new))
end

"""
    prepare_rates(r, param, hierarchy::HierarchicalTrait)
"""
function prepare_rates(r, param, hierarchy::HierarchicalTrait)
    rshared, rindividual = prepare_rates(r, hierarchy)
    pindividual, rhyper = prepare_hyper(r, param, hierarchy)
    return rshared, rindividual, pindividual, rhyper
end

"""
prepare_rates_noiseparams(rates, nrates, reporter)
"""
function prepare_rates_noiseparams(rates::AbstractArray, nrates::Int, reporter::HMMReporter)
    rates[1:nrates], rates[nrates+1:nrates+reporter.n]
end

function prepare_rates_noiseparams(rates::Vector{T}, nrates::Int, reporter::HMMReporter) where {T<:AbstractArray}
    parts = [prepare_rates_noiseparams(ri, nrates, reporter) for ri in rates]
    r = [p[1] for p in parts]
    noiseparams = [p[2] for p in parts]
    return r, noiseparams
end

function prepare_rates_noiseparams(rates::Vector, nrates::Vector{Int}, reporter::Vector{HMMReporter})
    lens = nrates .+ getproperty.(reporter, :n)
    starts = cumsum(vcat(1, lens[1:end-1]))
    r = [rates[k:k+nrates[i]-1] for (i, k) in enumerate(starts)]
    noiseparams = [rates[k+nrates[i]:k+nrates[i]+reporter[i].n-1] for (i, k) in enumerate(starts)]
    return r, noiseparams
end

function prepare_rates_noiseparams(rates::Vector{T}, nrates::Vector{Int}, reporter::Vector{HMMReporter}) where {T<:AbstractArray}
    lens = nrates .+ getproperty.(reporter, :n)
    starts = cumsum(vcat(1, lens[1:end-1]))
    r = [
        [rates[i][k:k+nrates[j]-1] for (j, k) in enumerate(starts)]
        for i in eachindex(rates)
    ]
    noiseparams = [
        [rates[i][k+nrates[j]:k+nrates[j]+reporter[j].n-1] for (j, k) in enumerate(starts)]
        for i in eachindex(rates)
    ]
    return r, noiseparams
end

"""
    prepare_coupling(rates, couplingindices)
    
"""

function prepare_coupling(rates::AbstractVector, couplingindices)
    _typed_rate_copy(rates[couplingindices])
end

function prepare_coupling_inplace(rates::Vector{T}, couplingindices) where {T<:AbstractArray}
    coupling = Vector{Float64}[]
    for i in eachindex(rates)
        push!(coupling, rates[i][couplingindices])
    end
    return coupling
end

function prepare_coupling(rates::Vector{T}, couplingindices) where {T<:AbstractArray}
    [rates[i][couplingindices] for i in eachindex(rates)]
end

"""
    prepare_grid(rates, ngrid)
"""
function prepare_grid(rates, ngrid)
    rates[ngrid]
end

function prepare_grid_inplace(rates::Vector{T}, ngrid) where {T<:AbstractArray}
    grid = Vector{Float64}[]
    for i in eachindex(rates)
        push!(grid, rates[i][ngrid])
    end
    return grid
end

function prepare_grid(rates::Vector{T}, ngrid) where {T<:AbstractArray}
    [rates[i][ngrid] for i in eachindex(rates)]
end

"""
    prepare_rates_coupled(r, nrates, reporter, couplingindices)
"""
function prepare_rates_coupled(r, nrates, reporter, couplingindices)
    rates, noiseparams = prepare_rates_noiseparams(r, nrates, reporter)
    couplingStrength = prepare_coupling(r, couplingindices)
    return rates, noiseparams, couplingStrength
end

function prepare_coupling_ad(rates::AbstractVector, couplingindices)
    rates[couplingindices]
end

function prepare_rates_coupled_ad(r, nrates, reporter, couplingindices)
    rates, noiseparams = prepare_rates_noiseparams(r, nrates, reporter)
    couplingStrength = prepare_coupling_ad(r, couplingindices)
    return rates, noiseparams, couplingStrength
end

function prepare_rates_coupled(rates::Vector, nrates::Vector{Int}, couplingindices)
    starts = cumsum(vcat(1, nrates[1:end-1]))
    r = [rates[k:k+nrates[i]-1] for (i, k) in enumerate(starts)]
    return r, rates[couplingindices]
end

"""
    prepare_rates_coupled_full(rates::Vector, nrates::Vector{Int},
                               couplingindices::Vector{Int},
                               target_flat_indices::Vector{Int})

Split a flat parameter vector into per-unit base rates and precomputed coupling
rates for the full coupled stack.

- `rates`: Flat parameter vector (base rates followed by coupling strengths).
- `nrates`: Number of base rates per unit (same convention as `prepare_rates_coupled`).
- `couplingindices`: Flat indices into `rates` selecting the coupling strengths γₖ.
- `target_flat_indices`: Flat indices into the base-rate part of `rates` giving,
   for each coupling k, the associated target base rate index.

Returns `(unit_rates, coupling_rates)` where:

- `unit_rates::Vector{Vector{Float64}}`: Per-unit base rates, as in `prepare_rates_coupled`.
- `coupling_rates::Vector{Float64}`: Length `length(couplingindices)`, with
  `coupling_rates[k] = rates[couplingindices[k]] * rates[target_flat_indices[k]]`.
"""
function prepare_rates_coupled_full(rates::Vector, nrates::Vector{Int},
                                    couplingindices, targets)
    # Per-unit base rates (same layout as prepare_rates_coupled)
    model_rates, couplingStrengths = prepare_rates_coupled(rates, nrates, couplingindices)
    coupling_rates = [
        couplingStrengths[k] * model_rates[targets[k][1]][targets[k][2]]
        for k in eachindex(couplingStrengths)
    ]

    return model_rates, coupling_rates
end

"""
    prepare_rates_coupled_full(r, nrates, reporter, couplingindices, targets)

Split a flat parameter vector into per-unit base rates, per-unit noise params, and
pre-multiplied coupling rates for the full coupled stack (TCoupledFullComponents).

- `r`: Flat parameter vector (base rates, noise params, then coupling strengths).
- `nrates`: Number of base rates per unit.
- `reporter`: Reporter structure (used to determine noise-param layout).
- `couplingindices`: Flat indices selecting the coupling strengths γₖ.
- `targets`: `Vector{Tuple{Int,Int}}` — for coupling k: `(model, localindex)` of the
  target base rate used to compute `coupling_rates[k] = γₖ * rates[model][localindex]`.

Returns `(unit_rates, noiseparams, coupling_rates)`.
"""
function prepare_rates_coupled_full(r, nrates, reporter, couplingindices, targets)
    rates, noiseparams, couplingStrengths = prepare_rates_coupled(r, nrates, reporter, couplingindices)
    coupling_rates = [couplingStrengths[k] * rates[targets[k][1]][targets[k][2]]
                      for k in eachindex(couplingStrengths)]
    return rates, noiseparams, coupling_rates
end

function prepare_rates_coupled(rates, sourceStates, transitions, R::Tuple, S, insertstep, n_noise)
    r = Vector{Float64}[]
    noiseparams = Vector{Float64}[]
    couplingStrength = Float64[]
    j = 1
    for i in eachindex(R)
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + n_noise[i]
        push!(r, rates[j:j+n-1])
        j += n
    end
    for i in eachindex(R)
        s = sourceStates[i]
        if (s isa Integer && s > 0) || (s isa Vector && !isempty(s))
            push!(couplingStrength, rates[j])
            j += 1
        else
            push!(couplingStrength, 0.0)
        end
    end
    for i in eachindex(r)
        push!(noiseparams, r[i][end-n_noise[i]+1:end])
    end
    return r, couplingStrength, noiseparams
end

# Overload using full coupling (unit_model, connections); coupling strength order matches connection list.
function prepare_rates_coupled(rates, coupling::Tuple, transitions, R::Tuple, S, insertstep, n_noise)
    r = Vector{Float64}[]
    couplingStrength = Float64[]
    j = 1
    for i in eachindex(R)
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + n_noise[i]
        push!(r, rates[j:j+n-1])
        j += n
    end
    conns = (isempty(coupling) || length(coupling) < 2) ? [] : coupling[2]
    for k in 1:length(conns)
        push!(couplingStrength, rates[j])
        j += 1
    end
    noiseparams = [r[i][end-n_noise[i]+1:end] for i in eachindex(r)]
    return r, couplingStrength, noiseparams
end

"""
    prepare_rates_hierarchical(r, param, nrates, hierarchy, reporter)
"""
function prepare_rates_hierarchical(r, param, nrates, hierarchy, reporter)
    rshared, rindividual, pindividual, rhyper = prepare_rates(r, param, hierarchy)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, nrates, reporter)
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, nrates, reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper
end

"""
    prepare_rates_coupled_hierarchical(r, param, nrates, coupling, hierarchy, reporter)
"""
function prepare_rates_coupled_hierarchical(r, param, nrates, coupling, hierarchy, reporter)
    rshared, rindividual, pindividual, rhyper = prepare_rates(r, param, hierarchy)
    couplingshared = prepare_coupling(rshared, coupling.couplingindices)
    couplingindividual = prepare_coupling(rindividual, coupling.couplingindices)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, nrates, reporter)
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, nrates, reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, couplingshared, couplingindividual
end

"""
    prepare_rates_grid(r, nrates, reporter, gridindices)
"""
function prepare_rates_grid(r, nrates, reporter, gridindices)
    rates, noiseparams = prepare_rates_noiseparams(r, nrates, reporter)
    pgrid = prepare_grid(r, gridindices)
    return rates, noiseparams, pgrid
end

"""
    prepare_rates_coupled_grid(r, nrates, reporter, couplingindices, gridindices)
"""
function prepare_rates_coupled_grid(r, nrates, reporter, couplingindices, gridindices)
    rates, noiseparams = prepare_rates_noiseparams(r, nrates, reporter)
    couplingStrength = prepare_coupling(r, couplingindices)
    pgrid = prepare_grid(r, gridindices)
    return rates, noiseparams, couplingStrength, pgrid
end

"""
    prepare_rates_hierarchical_grid(r, param, nrates, hierarchy, grid, reporter)
"""
function prepare_rates_hierarchical_grid(r, param, nrates, hierarchy, grid, reporter)
    rshared, rindividual, pindividual, rhyper = prepare_rates(r, param, hierarchy)
    pgridshared = prepare_grid(rshared, grid)
    pgridindividual = prepare_grid(rindividual, grid)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, nrates, reporter)
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, nrates, reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, pgridshared, pgridindividual
end

"""
    prepare_rates_coupled_hierarchical_grid(r, param, nrates, coupling, hierarchy, grid, reporter)
"""
function prepare_rates_coupled_hierarchical_grid(r, param, nrates, coupling, hierarchy, grid, reporter)
    rshared, rindividual, pindividual, rhyper = prepare_rates(r, param, hierarchy)
    couplingshared = prepare_coupling(rshared, coupling.couplingindices)
    couplingindividual = prepare_coupling(rindividual, coupling.couplingindices)
    pgridshared = prepare_grid(rshared, grid)
    pgridindividual = prepare_grid(rindividual, grid)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, nrates, reporter)
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, nrates, reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, couplingshared, couplingindividual, pgridshared, pgridindividual
end

"""
    prepare_rates(param, model::AbstractGRSMmodel)
"""
function prepare_rates(param, model::AbstractGRSMmodel)
    r = get_rates(param, model)
    prepare_rates_noiseparams(r, model.nrates, model.reporter)
end

function prepare_rates_ad(param, model::AbstractGRSMmodel)
    r = get_rates_ad(param, model)
    prepare_rates_noiseparams(r, model.nrates, model.reporter)
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling,)}}
    r = get_rates(param, model)
    prepare_rates_coupled(r, model.nrates, model.reporter, model.trait.coupling.couplingindices)
end

function prepare_rates_ad(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling,)}}
    r = get_rates_ad(param, model)
    prepare_rates_coupled_ad(r, model.nrates, model.reporter, model.trait.coupling.couplingindices)
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical,)}}
    r = get_rates(param, model)
    prepare_rates_hierarchical(r, param, model.nrates, model.trait.hierarchical, model.reporter)
end

function prepare_rates_ad(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical,)}}
    r = get_rates_ad(param, model)
    rshared, rindividual = prepare_rates_ad(r, model.trait.hierarchical)
    pindividual, rhyper = prepare_hyper(r, param, model.trait.hierarchical)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, model.nrates, model.reporter)
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, model.nrates, model.reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical)}}
    r = get_rates(param, model)
    prepare_rates_coupled_hierarchical(r, param, model.nrates, model.trait.coupling, model.trait.hierarchical, model.reporter)
end

function prepare_rates_ad(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical)}}
    r = get_rates_ad(param, model)
    rshared, rindividual = prepare_rates_ad(r, model.trait.hierarchical)
    pindividual, rhyper = prepare_hyper(r, param, model.trait.hierarchical)
    couplingshared = prepare_coupling_ad(rshared, model.trait.coupling.couplingindices)
    couplingindividual = prepare_coupling_ad(rindividual, model.trait.coupling.couplingindices)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, model.nrates, model.reporter)
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, model.nrates, model.reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, rhyper, couplingshared, couplingindividual
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:grid,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:grid,)}}
    r = get_rates(param, model)
    prepare_rates_grid(r, model.nrates, model.reporter, model.trait.grid.gridindices)   
end


"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :grid,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :grid,)}}
    r = get_rates(param, model)
    prepare_rates_coupled_grid(r, model.nrates, model.reporter, model.trait.coupling.couplingindices, model.trait.grid.gridindices)
end



"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical, :grid,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical, :grid,)}}
    r = get_rates(param, model)
    prepare_rates_hierarchical_grid(r, param, model.nrates, model.trait.hierarchical, model.trait.grid.gridindices, model.reporter)
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical, :grid,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical, :grid,)}}
    r = get_rates(param, model)
    prepare_rates_coupled_hierarchical_grid(r, param, model.nrates, model.trait.coupling, model.trait.hierarchical, model.trait.grid.gridindices, model.reporter)
end

# Model loglikelihoods
#
# ## Automatic differentiation feasibility (summary)
#
# These notes refer to gradients of the **log-likelihood** (and log-posterior) with respect to parameters.
# **Finite differences** (`NUTSOptions(gradient=:finite)`, or benchmarks using central differences) evaluate
# the likelihood at **Float64** parameters only — they avoid pushing `Dual` or Zygote through the whole graph
# and are the most robust option when other modes fail.
#
# ### Histogram / count likelihoods (no trace HMM forward–backward)
#
# These use [`predictedfn`](@ref) → [`predictedarray`](@ref) / [`make_array`](@ref), built from promoter / reporter
# transition matrices, [`steady_state_vector`](@ref) or [`dwelltimePDF`](@ref) / [`offtimePDF`](@ref) /
# [`ontimePDF`](@ref) / [`offonPDF`](@ref) — **not** [`ll_hmm_trace`](@ref) and **not** [`kolmogorov_forward_ad`](@ref).
#
# | Data (all via [`loglikelihood_ad`](@ref) when `AbstractHistogramData` / count) | Zygote | ForwardDiff | Finite diff |
# |:--|:--|:--|:--|
# | [`RNAData`](@ref), [`RNACountData`](@ref) | Usually OK with `steady_state_solver=:augmented` | Usually OK if dimension moderate | Always OK |
# | [`RNAOnOffData`](@ref) (RNA + ON/OFF **time** histograms): [`predictedarray`](@ref) → [`offonPDF`](@ref) etc. | Usually OK (same steady-state / matrix stack as other histograms) | Usually OK | Always OK |
# | [`RNADwellTimeData`](@ref) (RNA + **dwell-time** histograms per type): [`MTComponents`](@ref) / [`TDComponents`](@ref), [`dwelltimePDF`](@ref) | Usually OK uncoupled; **coupled** units → large `make_mat_TC` / joint steady state — watch memory | Usually OK uncoupled; coupled cost scales with generator size (different failure mode than trace matrix-exp) | Always OK |
# | [`DwellTimeData`](@ref) (dwell-time **only**): same PDF machinery; if [`hastrait`](@ref)`(model, :coupling)` uses [`predictedarray`](@ref) with [`TCoupledComponents`](@ref) / [`TDCoupledFullComponents`](@ref) — joint generators, not the trace HMM | Coupled multi-unit dwell can be **heavy** for Zygote (large graphs) | Not blocked by the trace [`kolmogorov_forward_ad`](@ref) / `Dual` issue; may still be slow for huge state | Always OK |
#
# Specialized non-mutating helpers ([`predictedarray_ad`](@ref)) exist for some shapes (e.g. [`TDComponents`](@ref),
# [`RNAData`](@ref)); generic [`predictedfn`](@ref) for [`AbstractHistogramData`](@ref) calls [`predictedarray`](@ref)
# — use `steady_state_solver=:augmented` for AD through the nullspace.
#
# ### Trace likelihoods ([`AbstractTraceData`](@ref), [`TraceRNAData`](@ref))
#
# | Data / method | Zygote | ForwardDiff | Finite diff |
# |:--|:--|:--|:--|
# | [`AbstractTraceData`](@ref) (**uncoupled**) via [`ll_hmm_trace`](@ref) | Often OK on **short** traces; long traces → compiler / tape issues — [`hmm_checkpoint_steps`](@ref) or `gradient=:finite` | Often OK if HMM state dimension modest | Always OK |
# | [`AbstractTraceData`](@ref) + **coupling** | Very high memory; not for routine use | OK (Kolmogorov [`kolmogorov_forward_ad`](@ref) uses `Dual`-generic matrix exp); large state → slow | Always OK |
# | [`TraceRNAData`](@ref) | Same as trace for HMM + histogram part as for [`RNAData`](@ref) | Same as trace HMM part | Always OK |
#
# [`check_ad_gradient_feasibility`](@ref) does **not** block ForwardDiff on coupled traces (Kolmogorov path uses a
# `Dual`-generic matrix exponential). Coupled **dwell / ON–OFF histogram** models use separate code paths;
# use `gradient=:finite` if any AD mode misbehaves.
#
# Low-level calls (e.g. `ForwardDiff.gradient` on [`loglikelihood_ad`](@ref) alone) are not guarded — use this table
# and [`benchmark_trace_finitediff_gradient`](@ref) for coupled **trace** benchmarks.
#
# For gradients (e.g. Zygote), pass `steady_state_solver=:augmented` into `loglikelihood`,
# `predictedfn`, `predictedarray`, `ll_hmm_trace`, etc. That uses [`steady_state_vector`](@ref)
# with the augmented linear solve ([`normalized_nullspace_augmented`](@ref)). The default
# `steady_state_solver=:default` keeps the fast QR-based [`normalized_nullspace`](@ref).
# Trace HMMs use two parallel tracks: `HMM_STACK_MH` (default for `loglikelihood` / MCMC) vs `HMM_STACK_AD`
# (`loglikelihood_ad` / NUTS / ADVI). Pass `hmm_stack=HMM_STACK_MH` or `HMM_STACK_AD` to `ll_hmm_trace` if needed.

"""
    loglikelihood(param, data, model)

Calculates the log-likelihood for various types of data and models.

# Arguments
- `param`: The model parameters.
- `data`: The data, which can be of various types (e.g., `AbstractHistogramData`, `AbstractTraceData`, `TraceData`, `TraceRNAData`).
- `model`: The model, which can be of various types (e.g., `AbstractGmodel`, `GRSMcoupledmodel`, `AbstractGRSMmodel`, `GRSMhierarchicalmodel`).

# Description
This function calculates the log-likelihood for different types of data and models. It supports histogram data, trace data, coupled models, and hierarchical models. The specific calculation method depends on the types of `data` and `model` provided.

# Methods
- `loglikelihood(param, data::AbstractHistogramData, model::AbstractGmodel)`: Returns the log-likelihood of all data and a vector of the prediction histogram log-likelihood.
- `loglikelihood(param, data::AbstractTraceData, model::AbstractGmodel)`: Returns the log-likelihood of combined time series traces and each trace.
- `loglikelihood(param, data::TraceData, model::GRSMcoupledmodel)`: Returns the log-likelihood for a coupled model.
- `loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel)`: Returns the log-likelihood of time series traces and mRNA FISH steady state histogram.
- `loglikelihood(param, data::AbstractTraceData, model::GRSMhierarchicalmodel)`: Returns the log-likelihood for a hierarchical model.

# Returns
- `Float64`: The log-likelihood for the combined time series traces and each trace.
- `Tuple{Float64, Vector{Float64}}`: The log-likelihood of all data and a vector of the prediction histogram log-likelihood.

# Convention
- All loglikelihood functions in this codebase return the log-likelihood (higher is better), not the negative log-likelihood.
"""
function loglikelihood(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:default)
    predictions = predictedfn(param, data, model; steady_state_solver=steady_state_solver)
    hist = datahistogram(data)
    logpredictions = log.(max.(predictions, eps()))
    return sum(hist .* logpredictions), hist .* logpredictions  # Convention: return log-likelihoods
end

function loglikelihood(param, data::RNACountData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:default)
    predictions = predictedfn(param, data, model; steady_state_solver=steady_state_solver)
    logpredictions = Array{Float64,1}(undef, length(data.countsRNA))
    for k in eachindex(data.countsRNA)
        # Use per-cell yieldfactor (BUG FIX: was hardcoded to 1.0)
        p = technical_loss_at_k(data.countsRNA[k], predictions, data.yieldfactor[k], data.nRNA)
        logpredictions[k] = log(max(p, eps()))
    end
    return sum(logpredictions), logpredictions  # Convention: return log-likelihoods
end

# AD-friendly version (no in-place mutation); defaults to augmented steady-state for AD
function loglikelihood_ad(param, data::RNACountData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:augmented)
    predictions = predictedfn(param, data, model; steady_state_solver=steady_state_solver, rates_fn=get_rates_ad)
    logpredictions = [
        log(max(technical_loss_at_k(data.countsRNA[k], predictions, data.yieldfactor[k], data.nRNA), eps()))
        for k in eachindex(data.countsRNA)
    ]
    return sum(logpredictions), logpredictions
end

"""
    loglikelihood_ad(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel; steady_state_solver=:augmented)

Same as [`loglikelihood`](@ref) for histogram data but uses [`get_rates_ad`](@ref) in [`predictedfn`](@ref) (Zygote-friendly).
[`RNACountData`](@ref) and trace types use more specific methods.
"""
function loglikelihood_ad(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:augmented)
    predictions = predictedfn(param, data, model; steady_state_solver=steady_state_solver, rates_fn=get_rates_ad)
    hist = datahistogram(data)
    logpredictions = log.(max.(predictions, eps()))
    return sum(hist .* logpredictions), hist .* logpredictions
end

function loglikelihood(param, data::AbstractTraceData, model::AbstractGRSMmodel; steady_state_solver::Symbol=:default, hmm_stack::Symbol=nothing)
    stack = isnothing(hmm_stack) ? getfield(model, :hmm_stack) : hmm_stack
    ll_hmm_trace(param, data, model; steady_state_solver=steady_state_solver, hmm_stack=stack)
end

"""
    loglikelihood_ad(param, data::AbstractTraceData, model::AbstractGRSMmodel; ...)

Same as [`loglikelihood`](@ref) but uses the **AD track** ([`HMM_STACK_AD`](@ref)): `forward_ad`, [`kolmogorov_forward_ad`](@ref), …
Use for NUTS/ADVI/Zygote; use [`loglikelihood`](@ref) (MH track, [`HMM_STACK_MH`](@ref)) for Metropolis–Hastings.
"""
function loglikelihood_ad(
    param,
    data::AbstractTraceData,
    model::AbstractGRSMmodel;
    steady_state_solver::Symbol=:augmented,
    hmm_stack::Symbol=nothing,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    stack = isnothing(hmm_stack) ? getfield(model, :hmm_stack) : hmm_stack
    ll_hmm_trace(
        param, data, model;
        steady_state_solver=steady_state_solver,
        hmm_stack=stack,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )
end

function loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel; steady_state_solver::Symbol=:default, hmm_stack::Symbol=nothing)
    stack = isnothing(hmm_stack) ? getfield(model, :hmm_stack) : hmm_stack
    llg, llgp = ll_hmm_trace(param, data, model; steady_state_solver=steady_state_solver, hmm_stack=stack)
    r = get_rates(param, model)
    
    # Get nRNA_true from data structure (extract from yield tuple)
    nRNA_true = get_nRNA_true(data.yield, data.nRNA)
    p_true = predictedRNA(r[1:num_rates(model)], model.components.mcomponents, model.nalleles, nRNA_true; steady_state_solver=steady_state_solver)
    
    # Apply loss matrix if yield < 1.0 (observation noise for RNA histogram)
    yield_val = get_yield_value(data.yield)
    if yield_val < 1.0
        L = make_loss_matrix(data.nRNA, nRNA_true, yield_val)
        predictions = L * p_true
    else
        predictions = p_true
    end
    
    logpredictions = log.(max.(predictions, eps()))
    hist = datahistogram(data)
    return sum(hist .* logpredictions) + llg, vcat(hist .* logpredictions, llgp)  # Convention: all log-likelihoods
end

function loglikelihood_ad(
    param,
    data::TraceRNAData,
    model::AbstractGRSMmodel;
    steady_state_solver::Symbol=:augmented,
    hmm_stack::Symbol=nothing,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    stack = isnothing(hmm_stack) ? getfield(model, :hmm_stack) : hmm_stack
    llg, llgp = ll_hmm_trace(
        param, data, model;
        steady_state_solver=steady_state_solver,
        hmm_stack=stack,
        hmm_checkpoint_steps=hmm_checkpoint_steps,
    )
    r = get_rates(param, model)
    nRNA_true = get_nRNA_true(data.yield, data.nRNA)
    p_true = predictedRNA(r[1:num_rates(model)], model.components.mcomponents, model.nalleles, nRNA_true; steady_state_solver=steady_state_solver)
    yield_val = get_yield_value(data.yield)
    if yield_val < 1.0
        L = make_loss_matrix(data.nRNA, nRNA_true, yield_val)
        predictions = L * p_true
    else
        predictions = p_true
    end
    logpredictions = log.(max.(predictions, eps()))
    hist = datahistogram(data)
    return sum(hist .* logpredictions) + llg, vcat(hist .* logpredictions, llgp)
end

# Helper to get the right component
get_components(model::AbstractGRSMmodel, data::AbstractTraceHistogramData) = model.components.tcomponents
get_components(model::AbstractGRSMmodel, data::AbstractTraceData) = model.components
# (Add more as needed)

function ll_hmm_trace(
    param,
    data,
    model::AbstractGRSMmodel;
    steady_state_solver::Symbol=:default,
    hmm_stack::Symbol=HMM_STACK_MH,
    hmm_checkpoint_steps::Union{Nothing,Integer}=nothing,
)
    with_hmm_zygote_checkpoint(hmm_checkpoint_steps) do
        r = hmm_stack === HMM_STACK_AD ? prepare_rates_ad(param, model) : prepare_rates(param, model)
        components = get_components(model, data)
        observed_units = isempty(data.units) ? nothing : data.units
        if !isnothing(model.trait) && haskey(model.trait, :grid)
            ll_hmm(
                r, model.trait.grid.ngrid, components, model.reporter, data.interval, data.trace, model.method;
                steady_state_solver=steady_state_solver, hmm_stack=hmm_stack,
            )
        else
            # Only pass observed_units for CoupledFull stack; legacy Coupled and single-unit don't support it
            if components isa TCoupledFullComponents
                ll_hmm(
                    r, components, model.reporter, data.interval, data.trace, model.method;
                    observed_units=observed_units, steady_state_solver=steady_state_solver, hmm_stack=hmm_stack,
                )
            else
                ll_hmm(
                    r, components, model.reporter, data.interval, data.trace, model.method;
                    steady_state_solver=steady_state_solver, hmm_stack=hmm_stack,
                )
            end
        end
    end
end


# Predicted histogram functions

function predictedRNA(r, mcomponents, nalleles, nRNA; steady_state_solver::Symbol=:default)
    M = make_mat_M(mcomponents, r)
    steady_state(M, mcomponents.nT, nalleles, nRNA; steady_state_solver=steady_state_solver)
end


"""
    predictedfn(param, data::RNAData, model::AbstractGeneTransitionModel; steady_state_solver, rates_fn=get_rates)

Calculates the likelihood for a single RNA histogram.

# Arguments
- `param`: The model parameters.
- `data::RNAData`: The RNA data.
- `model::AbstractGeneTransitionModel`: The model.
- `rates_fn`: Maps `(param, model)` to the full rate vector; use [`get_rates_ad`](@ref) for AD.

# Returns
- `Vector{Float64}`: The steady-state probabilities for the RNA histogram.
"""
function predictedfn(
    param,
    data::AbstractRNAData,
    model::AbstractGeneTransitionModel;
    steady_state_solver::Symbol=:default,
    rates_fn::Function=get_rates,
)
    r = rates_fn(param, model)
    M = make_mat_M(model.components, r)
    
    # Get nRNA_true from data structure (already computed in load_data)
    # For RNACountData: data.nRNA is already nRNA_true
    # For other types: extract from yield tuple using get_nRNA_true
    if typeof(data) <: RNACountData
        nRNA_true = data.nRNA
    else
        nRNA_true = get_nRNA_true(data.yield, data.nRNA)
    end
    p_true = steady_state(M, model.components.nT, model.nalleles, nRNA_true; steady_state_solver=steady_state_solver)
    
    # Apply loss matrix if yield < 1.0 (observation noise)
    # Note: RNACountData has per-cell yields (Vector), handled separately in loglikelihood()
    # For histogram data (RNAData, RNAOnOffData, etc.), use global yield (Union{Float64, Tuple})
    if typeof(data) <: RNACountData
        # RNACountData uses per-cell yields in loglikelihood(), no loss matrix here
        # Return full distribution at nRNA_true (used by technical_loss_at_k)
        return p_true
    else
        # Extract yield value from data structure
        yield_val = get_yield_value(data.yield)
        if yield_val < 1.0
            # Apply loss matrix: observed = L * true
            L = make_loss_matrix(data.nRNA, nRNA_true, yield_val)
            p_observed = L * p_true
            return p_observed
        else
            # No loss: return distribution at observed size (truncate if needed)
            return p_true[1:data.nRNA]
        end
    end
end




"""
    predictedfn(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel; steady_state_solver, rates_fn=get_rates)

Calculates the likelihood for multiple histograms.

# Arguments
- `param`: The model parameters.
- `data::AbstractHistogramData`: The histogram data.
- `model::AbstractGeneTransitionModel`: The model.
- `rates_fn`: Maps `(param, model)` to the full rate vector; use [`get_rates_ad`](@ref) for AD.

# Returns
- `Array{Float64,2}`: An array of likelihoods for the histograms.
"""
function predictedfn(
    param,
    data::AbstractHistogramData,
    model::AbstractGeneTransitionModel;
    steady_state_solver::Symbol=:default,
    rates_fn::Function=get_rates,
)
    h = predictedarray(rates_fn(param, model), data, model; steady_state_solver=steady_state_solver)
    make_array(h)
end


"""
    predictedarray(r, data::RNAData{T1,T2}, model::AbstractGMmodel) where {T1<:Array, T2<:Array}

Calculates the likelihood for multiple histograms and returns an array of PDFs.

# Arguments
- `r`: The rates.
- `data::RNAData{T1,T2}`: The RNA data.
- `model::AbstractGMmodel`: The model.

# Returns
- `Array{Array{Float64,1},1}`: An array of PDFs for the histograms.
"""
function predictedarray(r, data::RNAData{T1,T2}, model::AbstractGMmodel; steady_state_solver::Symbol=:default) where {T1<:Array,T2<:Array}
    h = Array{Array{Float64,1},1}(undef, length(data.nRNA))
    for i in eachindex(data.nRNA)
        M = make_mat_M(model.components[i], r[(i-1)*2*model.G+1:i*2*model.G])
        h[i] = steady_state(M, model.G, model.nalleles, data.nRNA[i]; steady_state_solver=steady_state_solver)
    end
    trim_hist(h, data.nRNA)
end

# AD-friendly version (no in-place mutation); defaults to augmented steady-state for AD
function predictedarray_ad(r, data::RNAData{T1,T2}, model::AbstractGMmodel; steady_state_solver::Symbol=:augmented) where {T1<:Array,T2<:Array}
    h = [
        steady_state(
            make_mat_M(model.components[i], r[(i-1)*2*model.G+1:i*2*model.G]),
            model.G, model.nalleles, data.nRNA[i];
            steady_state_solver=steady_state_solver,
        )
        for i in eachindex(data.nRNA)
    ]
    trim_hist(h, data.nRNA)
end


"""
    predictedarray(r, data::RNAOnOffData, model::AbstractGeneTransitionModel)

Calculates the likelihood for RNA On/Off data.

# Arguments
- `r`: The rates.
- `data::RNAOnOffData`: The RNA On/Off data.
- `model::AbstractGeneTransitionModel`: The model.

# Returns
- `Array{Float64,1}`: The likelihoods for the RNA On/Off data.
"""
function predictedarray(r, data::RNAOnOffData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:default)
    #     if model.splicetype == "offdecay"
    #         # r[end-1] *= survival_fraction(nu, eta, model.R)
    #     end
    components = model.components
    onstates = model.reporter
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    M = make_mat_M(components.mcomponents, r)
    histF = steady_state(M, components.mcomponents.nT, model.nalleles, data.nRNA; steady_state_solver=steady_state_solver)
    modelOFF, modelON = offonPDF(data.bins, r, T, TA, TI, components.tcomponents.nT, components.tcomponents.elementsT, onstates; steady_state_solver=steady_state_solver)
    return [modelOFF, modelON, histF]
end


"""
    predictedarray(r, data::RNADwellTimeData, model::AbstractGeneTransitionModel)

Calculates the likelihood array for RNA dwell time data.

# Arguments
- `r`: The rates.
- `data::RNADwellTimeData`: The RNA dwell time data.
- `model::AbstractGeneTransitionModel`: The model.

# Description
This function calculates the likelihood array for RNA dwell time data using the specified model and rates. It constructs the transition matrices and calculates the probability density functions for the dwell times.

# Returns
- `Vector{Vector{Float64}}`: A vector of histograms representing the likelihoods for the dwell times.

"""
function predictedarray(r, data::RNADwellTimeData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:default)
    predictedarray(r, model.components, data.bins, model.reporter, data.DTtypes, model.nalleles, data.nRNA; steady_state_solver=steady_state_solver)
end

function predictedarray(r, data::DwellTimeData, model::AbstractGeneTransitionModel; steady_state_solver::Symbol=:default)
    predictedarray(r, model.components, data.bins, model.reporter, data.DTtypes; steady_state_solver=steady_state_solver)
end

function predictedarray(r, data::DwellTimeData, model::AbstractGRSMmodel{T}; steady_state_solver::Symbol=:default) where {T<:NamedTuple{(:coupling,)}}
    r, coupling_strength = prepare_rates_coupled(r, model.nrates, model.trait.coupling.couplingindices)
    predictedarray(r, coupling_strength, model.components, data.bins, model.reporter, data.DTtypes; steady_state_solver=steady_state_solver)
end

function predictedarray(r, G::Int, tcomponents, bins, onstates, dttype; steady_state_solver::Symbol=:default)
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = steady_state_vector(T; solver=steady_state_solver)
    elementsG = tcomponents.elementsG
    TG = make_mat(elementsG, r, G)
    pssG = steady_state_vector(TG; solver=steady_state_solver)
    hists = Vector[]
    for (i, Dtype) in enumerate(dttype)
        if Dtype == "OFF"
            TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
            nonzeros = nonzero_rows(TD)
            h = offtimePDF(bins[i], TD[nonzeros, nonzeros], nonzero_states(onstates[i], nonzeros), init_SI(r, onstates[i], elementsT, pss, nonzeros))
        elseif Dtype == "ON"
            TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
            h = ontimePDF(bins[i], TD, off_states(tcomponents.nT, onstates[i]), init_SA(r, onstates[i], elementsT, pss))
        elseif Dtype == "OFFG"
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            h = offtimePDF(bins[i], TD, onstates[i], init_SI(r, onstates[i], elementsG, pssG, collect(1:G)))
        elseif Dtype == "ONG"
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            h = ontimePDF(bins[i], TD, off_states(G, onstates[i]), init_SA(r, onstates[i], elementsG, pssG))
        end
        push!(hists, h)
    end
    return hists
end

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, G::Int, tcomponents, bins, onstates, dttype; steady_state_solver::Symbol=:augmented)
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = steady_state_vector(T; solver=steady_state_solver)
    elementsG = tcomponents.elementsG
    TG = make_mat(elementsG, r, G)
    pssG = steady_state_vector(TG; solver=steady_state_solver)
    [
        Dtype == "OFF" ? begin
            TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
            nonzeros = nonzero_rows(TD)
            offtimePDF(bins[i], TD[nonzeros, nonzeros], nonzero_states(onstates[i], nonzeros), init_SI(r, onstates[i], elementsT, pss, nonzeros))
        end :
        Dtype == "ON" ? begin
            TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
            ontimePDF(bins[i], TD, off_states(tcomponents.nT, onstates[i]), init_SA(r, onstates[i], elementsT, pss))
        end :
        Dtype == "OFFG" ? begin
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            offtimePDF(bins[i], TD, onstates[i], init_SI(r, onstates[i], elementsG, pssG, collect(1:G)))
        end :
        Dtype == "ONG" ? begin
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            ontimePDF(bins[i], TD, off_states(G, onstates[i]), init_SA(r, onstates[i], elementsG, pssG))
        end :
        error("Unknown Dtype: $Dtype")
        for (i, Dtype) in enumerate(dttype)
    ]
end

function steady_state_dist(r, components, dt; steady_state_solver::Symbol=:default)
    pss = nothing
    pssG = nothing
    for b in union(dt)
        if b
            pssG = steady_state_vector(make_mat_G(components, r); solver=steady_state_solver)
        else
            pss = steady_state_vector(make_mat_T(components, r); solver=steady_state_solver)
        end
    end
    return (pss=pss, pssG=pssG, dt=dt)
end

function predictedarray(r, components::TDComponents, bins, reporter, dttype; steady_state_solver::Symbol=:default)
    sojourn, nonzeros = reporter
    dt = occursin.("G", dttype)
    p = steady_state_dist(r, components, dt; steady_state_solver=steady_state_solver)
    hists = Vector[]
    for i in eachindex(sojourn)
        if dt[i]
            TD = make_mat(components.elementsTD[i], r, components.TDdims[i])
            push!(hists, dwelltimePDF(bins[i], TD, sojourn[i], init_S(r, sojourn[i], components.elementsG, p.pssG)))
        else
            TD = make_mat(components.elementsTD[i], r, components.TDdims[i])
            push!(hists, dwelltimePDF(bins[i], TD[nonzeros[i], nonzeros[i]], nonzero_states(sojourn[i], nonzeros[i]), init_S(r, sojourn[i], components.elementsT, p.pss, nonzeros[i])))
        end
    end
    hists
end

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, components::TDComponents, bins, reporter, dttype; steady_state_solver::Symbol=:augmented)
    sojourn, nonzeros = reporter
    dt = occursin.("G", dttype)
    p = steady_state_dist(r, components, dt; steady_state_solver=steady_state_solver)
    [
        dt[i] ?
            dwelltimePDF(bins[i], make_mat(components.elementsTD[i], r, components.TDdims[i]), sojourn[i], init_S(r, sojourn[i], components.elementsG, p.pssG)) :
            dwelltimePDF(bins[i], make_mat(components.elementsTD[i], r, components.TDdims[i])[nonzeros[i], nonzeros[i]], nonzero_states(sojourn[i], nonzeros[i]), init_S(r, sojourn[i], components.elementsT, p.pss, nonzeros[i]))
        for i in eachindex(sojourn)
    ]
end

function predictedarray(r, components::MTComponents, bins, reporter, dttype, nalleles, nRNA; steady_state_solver::Symbol=:default)
    M = make_mat_M(components.mcomponents, r)
    [steady_state(M, components.mcomponents.nT, nalleles, nRNA; steady_state_solver=steady_state_solver); predictedarray(r, components.tcomponents, bins, reporter, dttype; steady_state_solver=steady_state_solver)...]
end

function steady_state_dist(unit::Int, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, dt; steady_state_solver::Symbol=:default)
    pssG = nothing
    pss = nothing
    TC = nothing
    GC = nothing
    for b in union(dt)
        if b
            GC = make_mat_TC(coupling_strength, Gm, Gs, Gt, IG, sources, model)
            pssG = steady_state_vector(GC; solver=steady_state_solver)
        else
            TC = make_mat_TC(unit, T[unit], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
            pss = steady_state_vector(TC; solver=steady_state_solver)
        end
    end
    return (pss=pss, pssG=pssG, TC=TC, GC=GC)
end

function compute_dwelltime!(cache::CoupledDTCache, unit, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn, dt, nonzeros, bins)
    if dt
        if !haskey(cache.TC, dt)
            cache.TC[dt] = make_mat_TC(coupling_strength, Gm, Gs, Gt, IG, sources, model)
            cache.pss[dt] = steady_state_vector(cache.TC[dt]; solver=cache.steady_state_solver)
        end
        TCD = make_mat_TCD(cache.TC[dt], sojourn)
        return dwelltimePDF(bins, TCD, sojourn, init_S(sojourn, cache.TC[dt], cache.pss[dt]))
    else
        if !haskey(cache.TC, dt)
            cache.TC[dt] = make_mat_TC(unit, T[unit], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
            cache.pss[dt] = steady_state_vector(cache.TC[dt]; solver=cache.steady_state_solver)
        end
        TCD = make_mat_TCD(cache.TC[dt], sojourn)
        return dwelltimePDF(bins, TCD[nonzeros, nonzeros], nonzero_states(sojourn, nonzeros), init_S(sojourn, cache.TC[dt], cache.pss[dt], nonzeros))
    end
end

function empty_cache!(cache)
    if haskey(cache.TC, false)
        delete!(cache.TC, false)
        delete!(cache.pss, false)
    end
end

function predictedarray(r, coupling_strength, components::TCoupledComponents{Vector{TDCoupledUnitComponents}}, bins, reporter, dttype; steady_state_solver::Symbol=:default)
    sojourn, nonzeros = reporter
    sources = components.sources
    model = components.model
    cache = CoupledDTCache(Dict(), Dict(), steady_state_solver)
    T, TD, Gm, Gt, Gs, IG, IR, IT = make_matvec_C(components, r)
    hists = Vector{Vector}[]
    for α in eachindex(components.modelcomponents)
        dt = occursin.("G", dttype[α])
        h = Vector[]
        for i in eachindex(sojourn[α])
            hdt = compute_dwelltime!(cache, α, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn[α][i], dt[i], nonzeros[α][i], bins[α][i])
            push!(h, hdt)
        end
        push!(hists, h)
        empty_cache!(cache)
    end
    return hists
end

# Full stack dwell-time prediction.
# r::Vector{Vector{Float64}} (per-unit rates), coupling_strength::Vector{Float64} (γ values).
# Sojourn sets come from reporter (same pattern as TDComponents); components store filtered elements.
function predictedarray(r, coupling_strength, components::TDCoupledFullComponents, bins, reporter, dttype; steady_state_solver::Symbol=:default)
    sojourn_full, _ = reporter
    coupling_rates = [coupling_strength[k] * r[components.targets[k][1]][components.targets[k][2]]
                      for k in eachindex(components.targets)]
    T = make_mat_TC(components, r, coupling_rates)
    pss = steady_state_vector(T; solver=steady_state_solver)
    hists = Vector{Vector}[]
    for α in eachindex(sojourn_full)
        h = Vector[]
        for i in eachindex(sojourn_full[α])
            soj = sojourn_full[α][i]
            TCD = make_mat_TCD(components, α, i, r, coupling_rates)
            init = init_S(soj, T, pss, dttype[α][i])
            init = init[soj]
            s = sum(init)
            (iszero(s) || !isfinite(s)) && continue
            init = init / s
            push!(h, dwelltimePDF(bins[α][i], TCD[soj, soj], 1:length(soj), init))
        end
        push!(hists, h)
    end
    return hists
end

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, coupling_strength, components::TCoupledComponents{Vector{TDCoupledUnitComponents}}, bins, reporter, dttype; steady_state_solver::Symbol=:augmented)
    sojourn, nonzeros = reporter
    sources = components.sources
    model = components.model
    T, TD, Gm, Gt, Gs, IG, IR, IT = make_matvec_C(components, r)
    [
        [
            compute_dwelltime!(
                CoupledDTCache(Dict(), Dict(), steady_state_solver), α, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model,
                sojourn[α][i], occursin("G", dttype[α][i]), nonzeros[α][i], bins[α][i]
            )
            for i in eachindex(sojourn[α])
        ]
        for α in eachindex(components.modelcomponents)
    ]
end

"""
    transform_array(v, index, f1, f2)
    transform_array(v, index, mask, f1, f2)

Applies transformation functions to an array or vector at specified indices.

# Description
This function applies the transformation functions `f1` and `f2` to an array or vector `v` at specified indices. If a `mask` is provided, it adjusts the indexing accordingly before applying the transformations.

# Arguments
- `v`: The array or vector to be transformed.
- `index`: The index up to which `f1` is applied and after which `f2` is applied.
- `mask`: A vector used to adjust the indexing for the transformation.
- `f1`: The function to apply to the elements up to `index`.
- `f2`: The function to apply to the elements after `index`.

# Methods
- `transform_array(v::Array, index::Int, f1::Function, f2::Function)`: Applies `f1` to the array `v` up to `index` and `f2` after `index`.
- `transform_array(v::Vector, index::Int, f1::Function, f2::Function)`: Applies `f1` to the vector `v` up to `index` and `f2` after `index`.
- `transform_array(v::Array, index::Int, mask::Vector, f1::Function, f2::Function)`: Applies `f1` to the array `v` up to `index` and `f2` after `index`, adjusting for the `mask`.
- `transform_array(v::Vector, index::Int, mask::Vector, f1::Function, f2::Function)`: Applies `f1` to the vector `v` up to `index` and `f2` after `index`, adjusting for the `mask`.

# Returns
- `Array` or `Vector`: The transformed array or vector.
"""
function transform_array(v::Array, index, f1::Function, f2::Function)
    vcat(f1(v[1:index-1, :]), f2(v[index, :]))
end

function transform_array(v::Vector, index, f1::Function, f2::Function)
    vcat(f1(v[1:index-1]), f2(v[index:end]))
end

function transform_array(v::Vector, index, mask::Vector, f1::Function, f2::Function)
    # if index > 0 && index ∈ mask
    if !isdisjoint(index, mask)
        n = findfirst(index .== mask)
        return vcat(f1(v[1:n-1]), f2(v[n]), f1(v[n+1:end]))
    else
        return f1(v)
    end
end

function transform_array(v::Array, index, mask::Vector, f1::Function, f2::Function)
    # if index > 0 && index ∈ mask
    if !isdisjoint(index, mask)
        n = findfirst(index .== mask)
        # return vcat(f1(v[1:n-1, :]), f2(v[n:n:end, :]))
        return vcat(f1(v[1:n-1, :]), f2(v[n:n, :]), f1(v[n+1:end, :]))
    else
        return f1(v)
    end
end

# function transform_array(v::Vector, index, mask::Vector, f1::Function, f2::Function)
#     if all(x -> x > 0, index) && index ∈ mask
#         f = f1(v)
#         f[index] = f2(v[index])
#         return f[mask]
#         # n = findfirst(index .== mask)
#         # return vcat(f1(v[1:n-1]), f2(v[n:end]))
#     else
#         return f1(v)
#     end
# end


function apply_transform(p::AbstractVector, fs)
    length(p) == length(fs) || throw(DimensionMismatch("apply_transform: length(p) != length(fs)"))
    return [fs[i](p[i]) for i in eachindex(p, fs)]
end

function apply_transform(p::AbstractMatrix, fs)
    size(p, 1) == length(fs) || throw(DimensionMismatch("apply_transform: size(p,1) != length(fs)"))
    size(p, 2) == 0 && return Matrix{eltype(p)}(undef, size(p, 1), 0)
    cols = [apply_transform(@view(p[:, j]), fs) for j in axes(p, 2)]
    return reduce(hcat, cols)
end

"""
    transform_rates(r, model::AbstractGeneTransitionModel)

Transform rates to log space.

# Arguments
- `r`: Vector of rates.
- `model`: The model, which can be of various types (e.g., `AbstractGeneTransitionModel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Returns
- Transformed rates.
"""
transform_rates(r, model::AbstractGeneTransitionModel) = log.(r)

function transform_rates(r, model::AbstractGRSMmodel)
    apply_transform(r, model.transforms.f[model.fittedparam])
    # rtransformed = map((f, x) -> f(x), model.transforms.f[model.fittedparam], r)
    # rtransformed[model.fittedparam]
end

# function transform_rates(r::AbstractVector, f::Vector{Function})
#     map((f, x) -> f(x), f, r)
# end

# function transform_rates(r::AbstractMatrix, f::Vector{Function})
#     hcat([map((f, x) -> f(x), f, r[:, i]) for i in 1:size(r, 2)]...)
# end

function transform_rates2(r, model::AbstractGRSMmodel)
    if hastrait(model, :coupling)
        return transform_array(r, model.trait.coupling.couplingindices, model.fittedparam, logv, log_shift1)
    else
        # future fix: include reporter.weightind
        return transform_array(r, 0, model.fittedparam, logv, logit)
    end
end

# transform_rates(r, model::AbstractGRSMmodel{Vector{Float64},HMMReporter}) = transform_array(r, model.reporter.weightind, model.fittedparam, logv, logit)

# transform_rates(r, model::GRSMcoupledmodel) = transform_array(r, length(model.rates), model.fittedparam, logv, log_shift1)

# transform_rates(r, model::GRSMcoupledhierarchicalmodel) = transform_array(r, model.coupling[2], model.fittedparam, logv, log_shift1)

"""
    inverse_transform_params(x, model::AbstractGeneTransitionModel)

Transform rates from log space back to original space.

# Arguments
- `x`: Vector of log-transformed rates.
- `model`: The model, which can be of various types (e.g., `AbstractGeneTransitionModel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Returns
- Transformed rates.
"""
inverse_transform_params(p, model::AbstractGeneTransitionModel) = exp.(p)

function inverse_transform_params(p, model::AbstractGRSMmodel)
    apply_transform(p, model.transforms.f_inv[model.fittedparam])
end

# function inverse_transform_params(p::AbstractMatrix, model::AbstractGRSMmodel)
#     # Apply inverse transforms to each column of p
#     hcat([map((f, x) -> f(x), model.transforms.f_inv[model.fittedparam], p[:, i]) for i in 1:size(p, 2)]...)
# end

# function inverse_transform_params(p::AbstractVector, f_inv::Vector{Function})
#     map((f, x) -> f(x), f_inv, p)
# end

# function inverse_transform_params(p::AbstractMatrix, f_inv::Vector{Function})
#     hcat([map((f, x) -> f(x), f_inv, p[:, i]) for i in 1:size(p, 2)]...)
# end



function inverse_transform_params2(p, model::AbstractGRSMmodel)
    if hastrait(model, :coupling)
        return transform_array(p, model.trait.coupling.couplingindices, model.fittedparam, expv, invlog_shift1)
    else
        return exp.(p)
    end
end

# inverse_transform_params(p, model::AbstractGRSMmodel{Vector{Float64},HMMReporter}) = transform_array(p, model.reporter.weightind, model.fittedparam, expv, invlogit)

# inverse_transform_params(p, model::GRSMcoupledmodel) = transform_array(p, length(model.rates), model.fittedparam, expv, invlog_shift1)

# inverse_transform_params(p, model::GRSMcoupledhierarchicalmodel) = transform_array(p, model.coupling[2], model.fittedparam, expv, invlog_shift1)


"""
    get_param(model::AbstractGeneTransitionModel)

Return the fitted parameters in the same (transformed) space used for MCMC (e.g. log-scale for positive rates, logit for probabilities).

# Arguments
- `model`: The model (e.g. `GMmodel`, `GRSMmodel`).

# Returns
- `Vector{Float64}`: Transformed fitted parameters. For simple models, log.(rates); for GRSM, model-specific transforms (e.g. log, logit) are applied.
"""
get_param(model::AbstractGeneTransitionModel) = log.(model.rates[model.fittedparam])

function get_param(model::AbstractGRSMmodel) 
    transform_rates(model.rates[model.fittedparam], model)
end

"""Write inverse-transformed fitted parameters into rate vector `r` (mutates `r`)."""
function get_rates!(r, param, model, inverse)
    if inverse
        r[model.fittedparam] = inverse_transform_params(param, model)
    else
        r[model.fittedparam] = param
    end
end

"""True if `rates` is stored as one vector per unit (e.g. `Vector{Vector{Float64}}`)."""
_rates_is_per_unit_vectors(rates) =
    rates isa AbstractVector && !isempty(rates) && eltype(rates) <: AbstractVector{<:Real}

"""Concrete element type for a possibly abstract-typed rate container (e.g. `Vector{Real}` holding `Dual`s)."""
function _concrete_rate_eltype(vals)
    isempty(vals) && return Float64
    T = typeof(first(vals))
    for x in Iterators.drop(vals, 1)
        T = promote_type(T, typeof(x))
    end
    return T
end

function _typed_rate_copy(vals::AbstractVector)
    map(identity, vals)
end

"""Merge `fittedparam` slots from `u` into base rates `r0` without mutating (AD-friendly)."""
function _merge_param_into_rates(r0::AbstractVector, u::AbstractVector, fittedparam)
    isempty(r0) && return _typed_rate_copy(r0)
    seed = isempty(u) ? first(r0) : first(u)
    return [
        begin
            i = findfirst(==(k), fittedparam)
            i === nothing ? (r0[k] + zero(seed)) : u[i]
        end for k in eachindex(r0)
    ]
end

"""
    fixed_rates_ad(r, fixedeffects)

Same result as [`fixed_rates`](@ref), but implemented without in-place broadcast
mutation so reverse-mode AD (e.g. Zygote) can differentiate through it.
"""
function fixed_rates_ad(r::AbstractVector, fixedeffects)
    out = r
    isempty(fixedeffects) && return out
    for effect in fixedeffects
        ref = out[effect[1]]
        slaves = effect[2:end]
        out = [i ∈ slaves ? ref : out[i] for i in eachindex(out)]
    end
    return out
end

"""
    get_rates_ad(param, model::AbstractGeneTransitionModel, inverse=true)

AD-friendly counterpart to [`get_rates`](@ref): builds the full rate vector from
transformed `param` without mutating buffers (`get_rates!` or shared model storage).
Use this inside `predictedfn` / `loglikelihood_ad` when computing gradients.

For `param::Vector{Vector}` with per-unit `model.rates`, see the specialized
[`get_rates_ad(param::Vector{Vector}, model)`](@ref) method.

When `model.rates` is a vector of per-unit rate vectors, this method merges `param`
into **each** unit (same rule as [`get_rates`](@ref) with `param::Vector{Vector}`) and
applies [`fixed_rates_ad`](@ref) to the **last** unit only, matching legacy [`fixed_rates`](@ref)
on that layout.
"""
function get_rates_ad(param, model::AbstractGeneTransitionModel, inverse::Bool=true)
    r0 = copy_r(model)
    if _rates_is_per_unit_vectors(r0)
        rv = [copy(ri) for ri in r0]
        u = inverse ? inverse_transform_params(param, model) : param
        for i in eachindex(rv)
            rv[i] = collect(_merge_param_into_rates(rv[i], u, model.fittedparam))
        end
        return fixed_rates_ad(rv[end], model.fixedeffects)
    end
    u = inverse ? inverse_transform_params(param, model) : param
    r = _merge_param_into_rates(r0, u, model.fittedparam)
    fixed_rates_ad(r, model.fixedeffects)
end

"""
    get_rates_ad(param::Vector{Vector}, model::AbstractGeneTransitionModel, inverse=true)

Same role as [`get_rates`](@ref) for the `param::Vector{Vector}` layout when
`model.rates` stores one rate vector per unit (`Vector{<:AbstractVector{<:Real}}`).
Uses deep copies of each unit vector so AD and callers never mutate `model.rates`
(unlike [`get_rates`](@ref), which uses `copy` and can alias inner vectors).

If `model.rates` is a flat `Vector{<:Real}`, use a single vector `param` instead;
this method throws `ArgumentError` in that case.
"""
function get_rates_ad(param::Vector{Vector}, model::AbstractGeneTransitionModel, inverse::Bool=true)
    rates = model.rates
    if !(rates isa AbstractVector) || !(eltype(rates) <: AbstractVector)
        throw(ArgumentError(
            "get_rates_ad(param::Vector{Vector}, model): model.rates must be a Vector of per-unit rate vectors. " *
            "For a flat rate vector, use get_rates_ad(param::AbstractVector, model).",
        ))
    end
    rv = [copy(ri) for ri in rates]
    u = inverse ? inverse_transform_params(param, model) : param
    for i in eachindex(rv)
        rv[i] = collect(_merge_param_into_rates(rv[i], u, model.fittedparam))
    end
    return fixed_rates_ad(rv[end], model.fixedeffects)
end

"""
    get_rates(param, model::AbstractGeneTransitionModel, inverse=true)

Map (transformed) parameters back to rate space and return the full rate vector.

For a flat `param` vector, this delegates to [`get_rates_ad`](@ref) (same numerics as the
legacy `get_rates!` + [`fixed_rates`](@ref) path, but type-generic so `ForwardDiff.Dual`
and other reals propagate). For `param::Vector{Vector}`, the implementation still uses
[`get_rates!`](@ref) on copied buffers.

# Arguments
- `param`: Vector of parameters in transformed space (e.g. from `get_param(model)` or MCMC samples).
- `model`: The model defining the transform and fitted indices.
- `inverse`: If `true` (default), apply inverse transform to `param`; if `false`, treat `param` as already in rate space for the fitted indices.

# Returns
- Full rate vector (or vector of vectors for coupled/hierarchical) with fitted indices set from `param` and fixed effects applied.
"""

"""
    get_rates(param, model::AbstractGeneTransitionModel, inverse=true)

Classic (non-AD) parameter-to-rate logic. Copies model rates, inserts (optionally inverse-transformed) params, applies fixed effects. Does NOT call or delegate to get_rates_ad. For AD/gradient code, use get_rates_ad instead.

# Arguments
- `param`: Vector of parameters in transformed space (e.g. from `get_param(model)` or MCMC samples).
- `model`: The model defining the transform and fitted indices.
- `inverse`: If `true` (default), apply inverse transform to `param`; if `false`, treat `param` as already in rate space for the fitted indices.

# Returns
- Full rate vector (or vector of vectors for coupled/hierarchical) with fitted indices set from `param` and fixed effects applied.
"""
function get_rates(param, model::AbstractGeneTransitionModel, inverse=true)
    r = copy_r(model)
    get_rates!(r, param, model, inverse)
    fixed_rates(r, model.fixedeffects)
end

function get_rates(param::Vector{Vector}, model::AbstractGeneTransitionModel, inverse=true)
    rv = copy_r(model)
    for (i, r) in enumerate(rv)
        get_rates!(r, param[i], model, inverse)
        rv[i] = fixed_rates(r, model.fixedeffects)
    end
    rv
end


"""
    get_rates(fits, stats, model, ratetype)

Return rate vector from MCMC results using the requested summary type.

# Arguments
- `fits`: `Fit` struct from a completed MCMC run.
- `stats`: `Stats` struct from the same run.
- `model`: The model used for the fit.
- `ratetype`: One of `"ml"` (maximum likelihood), `"median"` (posterior median), or `"mean"` (posterior mean).

# Returns
- `Vector` or `Vector{Vector}` of rates in rate space (not transformed). For single-unit models, a single vector; for coupled/hierarchical, may be a vector of vectors.
"""
function get_rates(fits, stats, model, ratetype)
    if ratetype == "ml"
        return get_rates(fits.parml, model)
    elseif ratetype == "median"
        return get_rates(stats.medparam, model, false)
    elseif ratetype == "mean"
        return get_rates(stats.meanparam, model, false)
    else
        println("ratetype unknown")
    end
end


# function get_rates(param, model::GRSMcoupledmodel, inverse=true)
#     r = copy_r(model)
#     get_rates!(r, param, model, inverse)
#     fixed_rates(r, model.fixedeffects)
# end

# function get_rates_hierarchical(param, model::AbstractGRSMhierarchicalmodel, inverse=true)
#     r = get_rates(param, model, inverse)
#     rshared, rindividual, pindividual, rhyper = prepare_rates(r, param, model.hierarchy)
#     vcat(vec(rshared), vec(rindividual))
# end

"""
    fixed_rates(r, fixedeffects)

Applies fixed effects to the rates.

# Arguments
- `r`: The rates to be modified.
- `fixedeffects`: A vector of fixed effects, where each fixed effect is a vector of indices. The first index is the reference rate, and the remaining indices are the rates to be set equal to the reference rate.

# Description
This function applies fixed effects to the rates by setting specified rates equal to a reference rate. If the `fixedeffects` vector is not empty, it iterates through each fixed effect and sets the rates at the specified indices equal to the reference rate.

# Returns
- `Vector{Float64}`: The modified rates with fixed effects applied.
"""
function fixed_rates(r, fixedeffects)
    rnew = copy(r)
    if ~isempty(fixedeffects)
        for effect in fixedeffects
            rnew[effect[2:end]] .= rnew[effect[1]]
        end
    end
    return rnew
end

function fixed_rates_inplace(r, fixedeffects)
    if ~isempty(fixedeffects)
        for effect in fixedeffects
            r[effect[2:end]] .= r[effect[1]]
        end
    end
    return r
end
"""
    copy_r(model)

Copies the rates from the model structure.

# Arguments
- `model`: The model from which to copy the rates.

# Description
This function creates a copy of the rates from the provided model structure. It is useful for preserving the original rates while making modifications.

# Returns
- `Vector{Float64}`: A copy of the rates from the model.
"""
copy_r(model) = copy(model.rates)

# """
#     resetr(r, model::AbstractGRSMmodel)

# """
# function resetr(r, model::AbstractGRSMmodel)
#     n = model.G - 1
#     nr = model.R
#     eta = get_eta(r, n, nr)
#     r[2*n+1+nr+1:2*n+1+nr+nr] = eta
#     r
# end
"""
    logprior(param, model::AbstractGeneTransitionModel)

Compute log prior of parameters.

# Arguments
- `param`: Vector of parameters.
- `model::AbstractGeneTransitionModel`: The model, which includes the prior distributions for the rates.

# Returns
- Log prior of parameters.
"""
function logprior(param, model::AbstractGeneTransitionModel)
    d = model.rateprior
    p = 0
    for i in eachindex(d)
        p += logpdf(d[i], param[i])
    end
    return p
end

"""
    num_rates(transitions, R, S, insertstep)
    num_rates(transitions, R::Tuple, S::Tuple, insertstep::Tuple)

compute number of transition rates (not counting noise parameters)
"""
function num_rates(transitions::Tuple, R, S, insertstep)
    if R > 0
        if insertstep > R
            throw(DomainError("insertstep>R"))
        end
        if S == 0
            insertstep = 1
        elseif S > R - insertstep + 1
            S = R - insertstep + 1
        end
        return length(transitions) + 1 + R + S + 1
    else
        return length(transitions) + 2
    end
end


function num_rates(transitions::Tuple, R::Tuple, S::Tuple, insertstep::Tuple)
    n = Vector{Int}(undef, length(R))
    for i in eachindex(R)
        n[i] = num_rates(transitions[i], R[i], S[i], insertstep[i])
    end
    n
end

function num_rates_summed(transitions, R::Tuple, S::Tuple, insertstep::Tuple)
    sum(num_rates(transitions, R, S, insertstep))
end




"""
    num_rates(model::AbstractGeneTransitionModel)

compute number of transition rates

# Arguments
- `model`: The model, which can be of various types (e.g., `AbstractGRSMmodel`, `GRSMcoupledmodel`) or a model string.

# Description
This function computes the number of transition rates for the provided GRSM model. It calculates the number of rates based on the model's G transitions, R steps, S splicing indicator, and insert step.

#Methods
- `num_rates(model::AbstractGRSMmodel)`: Returns the number of transition rates for the GRSM model.
- `num_rates(model::String)`: Computes the number of transition rates given a model string.

# Returns
- `Int`: The number of transition rates.

"""

"""
    num_rates(model::AbstractGRSMmodel)

# Description
This function computes the number of transition rates for the provided GRSM model. It calculates the number of rates based on the model's G transitions, R steps, S splicing indicator, and insert step.

"""
num_rates(model::AbstractGRSMmodel) = num_rates(model.Gtransitions, model.R, model.S, model.insertstep)

num_rates(model::AbstractGMmodel) = length(model.Gtransitions) + 2


"""
    num_rates(model::String)

compute number of transition rates given model string
"""
function num_rates(model::String)
    m = digits(parse(Int, model))
    if length(m) == 1
        return 2 * m[1]
    else
        return 2 * (m[4] - 1) + m[3] + m[2] - m[1] + 2
    end
end


"""
    num_all_parameters(transitions, R, S, insertstep, reporter, coupling=tuple(), grid=nothing)

Compute the total number of parameters including rates, noise parameters, coupling parameters, and grid parameters.

# Arguments
- `transitions`: Model transitions
- `R`: Number of RNA steps (Int or Tuple)
- `S`: Number of splice states (Int or Tuple) 
- `insertstep`: Insertion step (Int or Tuple)
- `reporter`: Reporter structure (HMMReporter, Vector, or Int)
- `coupling`: Coupling structure (default: empty tuple)
- `grid`: Grid parameter (default: nothing)

# Returns
- `Int`: Total number of parameters

# Notes
- Calculates base rates using num_rates()
- Adds noise parameters based on reporter type
- Adds coupling parameters if coupling is not empty
- Adds grid parameter if grid is not nothing
- Used for parameter counting in model setup
"""
function num_all_parameters(transitions, R::Int, S, insertstep, reporter, coupling=tuple(), grid=nothing)
    if typeof(reporter) <: HMMReporter
        n = reporter.n
    elseif typeof(reporter) <: Vector
        n = length(reporter)
    elseif typeof(reporter) <: Int
        n = reporter
    else
        n = 0
    end
    c = isempty(coupling) ? 0 : length(coupling[2])
    g = isnothing(grid) ? 0 : 1
    num_rates(transitions, R, S, insertstep) + n + c + g
end

function num_all_parameters(transitions, R::Tuple, S::Tuple, insertstep::Tuple, reporter, coupling=tuple(), grid=nothing)
    n = 0
    for i in eachindex(R)
        n += num_all_parameters(transitions[i], R[i], S[i], insertstep[i], reporter[i])
    end
    c = isempty(coupling) ? 0 : length(coupling[2])
    g = isnothing(grid) ? 0 : 1
    n + c + g
end

function num_all_parameters(model::AbstractGeneTransitionModel)
    num_all_parameters(model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
end

# GRSM models: include coupling and grid from trait so total param count matches length(model.rates)
function num_all_parameters(model::AbstractGRSMmodel)
    n = num_all_parameters(model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
    c = hastrait(model, :coupling) ? model.trait.coupling.ncoupling : 0
    g = hastrait(model, :grid) ? 1 : 0
    n + c + g
end

function num_fitted_core_params(model::AbstractGeneTransitionModel)
    count(x -> x <= num_all_parameters(model), model.fittedparam)
end
