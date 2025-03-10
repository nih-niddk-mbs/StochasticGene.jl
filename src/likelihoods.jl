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
    v = normalize_histogram(data.histRNA)
    for d in data.DwellTimes
        v = vcat(v, normalize_histogram(d))
    end
    return v
end



####### trait system

"""
1. if hierarchical trait is present then, reshape rates into matrix with rows pertaining to rate/parameter and columns pertaining to shared, hyper, or individuals
2. if coupled trait is present then divide rows into each unit, extract out coupling strength row
3. if reporter type is HMMReprter then divide unit rates into rates and noise parameters
4. if grid trait is present then extract out grid rows
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

function prepare_hyper(r, param, hierarchy::HierarchicalTrait)
    pindividual = collect(eachcol(reshape(param[hierarchy.paramstart:end], hierarchy.nparams, hierarchy.nindividuals)))
    phyper = Vector{Float64}[]
    for i in hierarchy.hyperindices
        push!(phyper, r[i])
    end
    return pindividual, phyper
end

function prepare_rates(r, param, hierarchy::HierarchicalTrait)
    rshared, rindividual = prepare_rates(r, hierarchy)
    pindividual, phyper = prepare_hyper(r, param, hierarchy)
    return rshared, rindividual, pindividual, phyper
end

function prepare_rates_noiseparams(rates::Vector{T}, nrates::Vector{Int}, reporter::Vector{HMMReporter}) where T <: AbstractArray
    r = Vector{Vector{Vector{Float64}}}(undef, length(rates))
    noiseparams = Vector{Vector{Vector{Float64}}}(undef, length(rates))
    for i in eachindex(r)
        r[i] = Vector{Vector{Float64}}(undef, length(nrates))
        noiseparams[i] = Vector{Vector{Float64}}(undef, length(nrates))
        k = 1
        for j in eachindex(nrates)
            n = nrates[j]
            r[i][j] = rates[i][k:k+n-1]
            k += n
            noiseparams[i][j] = rates[i][k:k+reporter[j].n-1]
            k += reporter[j].n
        end
    end
    return r, noiseparams
end

# function prepare_rates_noiseparams(rates::Vector{T}, nrates::Vector{Int}, reporter::Vector{HMMReporter}) where T <: AbstractArray 
#     r = Vector{Vector{Float64}}(undef, length(nrates))
#     noiseparams = Vector{Vector{Float64}}(undef, length(nrates))
#     k = 1
#     for i in eachindex(r)
#         n = nrates[i]
#         r[i] = rates[i][k:k+n-1]
#         k += n
#         noiseparams[i] = rates[i][k:k+reporter[i].n-1]
#         k += reporter[i].n
#     end
#     return r, noiseparams
# end

function prepare_rates_noiseparams(rates::Vector{T}, nrates::Int, reporter::HMMReporter) where T <: AbstractArray
    r = Vector{Vector{Float64}}(undef, length(rates))
    noiseparams = Vector{Vector{Float64}}(undef, length(rates))
    k = 1
    for i in eachindex(r)
        r[i], noiseparams[i] = prepare_rates_noiseparams(rates[i], nrates, reporter)
    end
    return r, noiseparams
end

function prepare_rates_noiseparams(rates::AbstractArray, nrates::Int, reporter::HMMReporter)
    rates[1:nrates], rates[nrates+1:nrates+reporter.n]
end

function prepare_coupling(rates, couplingindices)
    rates[couplingindices]
end

function prepare_grid(rates, ngrid)
    rates[ngrid]
end

function prepare_rates_hierarchical(r, param, nrates, hierarchy, reporter)
    rshared, rindividual, pindividual, phyper = prepare_rates(r, param, hierarchy)
    rshared, noiseshared = prepare_rates_noiseparams(rshared, nrates, reporter)
    println(typeof(rindividual))
    println(length(rindividual))
    rindividual, noiseindividual = prepare_rates_noiseparams(rindividual, nrates, reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper
end

function prepare_rates(param, model::AbstractGRSMtraitmodel)
    r = get_rates(param, model)
    prepare_rates_noiseparams(r, model.nrates, model.reporter)
end

function prepare_rates(param, model::AbstractGRSMtraitmodel{T}) where {T<:NamedTuple{(:coupling,)}}
    r = get_rates(param, model)
    rates, noiseparams = prepare_rates_noiseparams(r, model.nrates, model.reporter)
    couplingStrength = prepare_coupling(r, model.trait.coupling.couplingindices)
    return rates, noiseparams, couplingStrength
end

function prepare_rates(param, model::AbstractGRSMtraitmodel{T}) where {T<:NamedTuple{(:hierarchical,)}}
    r = get_rates(param, model)
    prepare_rates_hierarchical(r, param, model.nrates, model.trait.hierarchical, model.reporter)
end

function prepare_rates(param, model::AbstractGRSMtraitmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical)}}
    r = get_rates(param, model)
    rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper = prepare_rates_hierarchical(r, param, model.nrates, model.trait.hierarchical, model.reporter)
    couplingStrength = prepare_coupling(r, model.trait.coupling.couplingindices)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper, couplingStrength
end

function prepare_rates(param, model::AbstractGRSMtraitmodel{T}) where {T<:NamedTuple{(:grid,)}}
    r = get_rates(param, model)
    rates, noiseparams = prepare_rates_noiseparams(r, model.nrates, model.reporter)
    pgrid = prepare_rates_grid(r, model.trait.grid.ngrid, model.reporter)
    return rates, noiseparams, pgrid
end

function prepare_rates(param, model::AbstractGRSMtraitmodel{T}) where {T<:NamedTuple{(:hierarchical, :grid,)}}
    r = get_rates(param, model)
    rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper = prepare_rates_hierarchical(r, param, model.nrates, model.trait.hierarchical, model.reporter)
    pgrid = prepare_rates_grid(r, model.trait.grid.ngrid, model.reporter)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper, pgrid
end

#####

"""
    prepare_coupling(rates, sourceStates::Vector, transitions, R, S, insertstep, reporter)

Prepare the coupling strength for the coupled model.

# Arguments
- `rates`: The rates.
- `sourceStates`: The source states.
- `transitions`: The transitions.
"""
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


function prepare_coupling(rates, sourceStates::Vector, transitions, R, S, insertstep, reporter)
    couplingStrength = Float64[]
    j = num_all_parameters(transitions, R, S, insertstep, reporter) + 1
    for s in sourceStates
        if (s isa Integer && s > 0) || (s isa Vector && !isempty(s))
            push!(couplingStrength, rates[j])
            j += 1
        else
            push!(couplingStrength, 0.0)
        end
    end
    couplingStrength
end

function prepare_coupled_rates(rates, transitions, R::Tuple, S, insertstep, reporter)
    r = Matrix{Float64}[]
    noiseparams = Matrix{Float64}[]
    j = 1
    for i in eachindex(R)
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + reporter[i].n
        push!(r, rates[j:j+n-1, :])
        j += n
    end
    for i in eachindex(r)
        push!(noiseparams, r[i][end-reporter[i].n+1:end, :])
    end
    return r, noiseparams
end


function prepare_hyper(r, param, hierarchy::Hierarchy)
    pindividual = collect(eachcol(reshape(param[hierarchy.paramstart:end], hierarchy.nparams, hierarchy.nindividuals)))
    phyper = Vector{Float64}[]
    for i in hierarchy.hyperindices
        push!(phyper, r[i])
    end
    return pindividual, phyper
end

function prepare_rates(r, hierarchy::Hierarchy)
    # rates reshaped from a vector into a vector of vectors pertaining to shared params, hyper params and individual params 
    # (shared parameters are considered to be hyper parameters without other hyper parameters (e.g. mean without variance))

    rshared = reshape(r[1:hierarchy.individualstart-1], hierarchy.nrates, hierarchy.nhypersets)

    rindividual = reshape(r[hierarchy.individualstart:end], hierarchy.nrates, hierarchy.nindividuals)
    rindividual[hierarchy.fittedshared, :] .= rshared[hierarchy.fittedshared, 1]

    return rshared, rindividual
end

function prepare_rates(r, param, hierarchy::Hierarchy)
    # rates reshaped from a vector into a vector of vectors pertaining to shared params, hyper params and individual params 
    # (shared parameters are considered to be hyper parameters without other hyper parameters (e.g. mean without variance))\

    rshared, rindividual = prepare_rates(r, hierarchy)
    pindividual, phyper = prepare_hyper(r, param, hierarchy)

    return rshared, rindividual, pindividual, phyper
end

"""
    prepare_rates(param, model::GRSMcoupledmodel)

convert MCMC params into form to compute likelihood for coupled model
"""
function prepare_rates(param, model::GRSMcoupledmodel)
    rates = get_rates(param, model)
    n_noise = [r.n for r in model.reporter]
    sourceStates = [c.sourceState for c in model.components.modelcomponents]
    prepare_rates_coupled(rates, sourceStates, model.Gtransitions, model.R, model.S, model.insertstep, n_noise)
end
function prepare_rates_new(param, model::GRSMhierarchicalmodel)
    r = get_rates(param, model)
    rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper = prepare_rates(r, param, model.hierarchy)
    return rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper
end
"""
    prepare_rates(param, model::AbstractGRSMhierarchicalmodel)

TBW
"""
function prepare_rates(param, model::AbstractGRSMhierarchicalmodel)
    # rates reshaped from a vector into a matrix with columns pertaining to shared params, hyper params and individual params 
    # (shared parameters are considered to be hyper parameters without other hyper parameters (e.g. mean without variance))
    r = get_rates(param, model)
    rshared = reshape(r[1:model.hierarchy.individualstart-1], model.hierarchy.nrates, model.hierarchy.nhypersets)
    phyper = Vector{Float64}[]
    for i in model.hierarchy.hyperindices
        push!(phyper, r[i])
    end
    rindividual = reshape(r[model.hierarchy.individualstart:end], model.hierarchy.nrates, model.hierarchy.nindividuals)
    rindividual[model.hierarchy.fittedshared, :] .= rshared[model.hierarchy.fittedshared, 1]
    pindividual = reshape(param[model.hierarchy.paramstart:end], model.hierarchy.nparams, model.hierarchy.nindividuals)
    return rshared, rindividual, pindividual, phyper
end

function prepare_rates(param, model::GRSMcoupledhierarchicalmodel)
    r = get_rates(param, model)
    rshared, rindividual, pindividual, phyper = prepare_rates(r, param, model.hierarchy)
    sourceStates = [c.sourceState for c in model.components.modelcomponents]
    rshared, noiseshared = prepare_coupled_rates(rshared, model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
    rindividual, noiseindividual = prepare_coupled_rates(rindividual, model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
    couplingStrength = prepare_coupling(r, sourceStates, model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
    rshared = [[p[:, j] for p in rshared] for j in 1:size(rshared[1], 2)]
    noiseshared = [[p[:, j] for p in noiseshared] for j in 1:size(noiseshared[1], 2)]
    rindividual = [[p[:, j] for p in rindividual] for j in 1:size(rindividual[1], 2)]
    noiseindividual = [[p[:, j] for p in noiseindividual] for j in 1:size(noiseindividual[1], 2)]
    rshared, rindividual, noiseshared, noiseindividual, couplingStrength, pindividual, phyper
end


"""
    prepare_rates(param, model::GRSMgridmodel)

TBW
"""
function prepare_rates(param, model::GRSMgridmodel)
    r = get_rates(param, model)
    r[model.raterange], r[model.noiserange], r[model.gridrange]
end



# Model loglikelihoods


"""
    ll_hierarchy(pindividual, phyper)

Loglikelihood for hierarchical model individual parameters.
"""
function ll_hierarchy(pindividual, phyper)
    d = distribution_array(mulognormal(phyper[1], phyper[2]), sigmalognormal(phyper[2]))
    lhp = Float64[]
    for pc in eachcol(pindividual)
        lhpc = 0
        for i in eachindex(pc)
            lhpc -= logpdf(d[i], pc[i])
        end
        push!(lhp, lhpc)
    end
    lhp
end

"""
    ll_hierarchy_c(pindividual, phyper)

Loglikelihood for coupled hierarchical model individual parameters.
"""
function ll_hierarchy_c(pindividual, phyper)
    d = distribution_array(mulognormal(phyper[1], phyper[2]), sigmalognormal(phyper[2]))
    lhp = Float64[]
    for pc in pindividual
        lhpc = 0
        println(pc)
        for i in eachindex(pc)
            lhpc -= logpdf(d[i], pc[i])
        end
        push!(lhp, lhpc)
    end
    lhp
end


"""
    loglikelihood(param, data, model)

Calculates the negative loglikelihood for various types of data and models.

# Arguments
- `param`: The model parameters.
- `data`: The data, which can be of various types (e.g., `AbstractHistogramData`, `AbstractTraceData`, `TraceData`, `TraceRNAData`).
- `model`: The model, which can be of various types (e.g., `AbstractGmodel`, `GRSMcoupledmodel`, `AbstractGRSMmodel`, `GRSMhierarchicalmodel`).

# Description
This function calculates the negative loglikelihood for different types of data and models. It supports histogram data, trace data, coupled models, and hierarchical models. The specific calculation method depends on the types of `data` and `model` provided.

# Methods
- `loglikelihood(param, data::AbstractHistogramData, model::AbstractGmodel)`: Returns the negative loglikelihood of all data and a vector of the prediction histogram negative loglikelihood.
- `loglikelihood(param, data::AbstractTraceData, model::AbstractGmodel)`: Returns the negative loglikelihood of combined time series traces and each trace.
- `loglikelihood(param, data::TraceData, model::GRSMcoupledmodel)`: Returns the negative loglikelihood for a coupled model.
- `loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel)`: Returns the negative loglikelihood of time series traces and mRNA FISH steady state histogram.
- `loglikelihood(param, data::AbstractTraceData, model::GRSMhierarchicalmodel)`: Returns the negative loglikelihood for a hierarchical model.

# Returns
- `Float64`: The negative loglikelihood for the combined time series traces and each trace.
- `Tuple{Float64, Vector{Float64}}`: The negative loglikelihood of all data and a vector of the prediction histogram negative loglikelihood.
"""
function loglikelihood(param, data::AbstractHistogramData, model::AbstractGmodel)
    predictions = predictedfn(param, data, model)
    hist = datahistogram(data)
    logpredictions = log.(max.(predictions, eps()))
    return crossentropy(logpredictions, hist), -logpredictions
end

function loglikelihood(param, data::AbstractTraceData, model::AbstractGmodel)
    ll_hmm(get_rates(param, model), model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
end

function loglikelihood(param, data::TraceData, model::GRSMcoupledmodel)
    r, couplingStrength, noiseparams = prepare_rates(param, model)
    ll_hmm_coupled(r, couplingStrength, noiseparams, model.components, model.reporter, data.interval, data.trace)
end

function loglikelihood(param, data::AbstractTraceData, model::GRSMhierarchicalmodel)
    rshared, rindividual, pindividual, phyper = prepare_rates(param, model)
    if model.method[2]
        llg, llgp = ll_hmm_hierarchical_rateshared(rshared, rindividual, model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    else
        llg, llgp = ll_hmm_hierarchical(rshared, rindividual, model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    end
    lhp = ll_hierarchy(pindividual, phyper)
    return llg + sum(lhp), vcat(llgp, lhp)
end

function loglikelihoodnew(param, data::AbstractTraceData, model::GRSMhierarchicalmodel)
    rshared, rindividual, noiseshared, noiseindividual, pindividual, phyper = prepare_rates(param, model)
    ll_hmm_hierarchical((rshared, rindividual, noiseshared, noiseindividual), model.components, model.reporter, data.interval, data.trace, model.method)
    lhp = ll_hierarchy(pindividual, phyper)
    return llg + sum(lhp), vcat(llgp, lhp)
end

function loglikelihood(param, data::TraceData, model::GRSMcoupledhierarchicalmodel)
    rshared, rindividual, noiseshared, noiseindividual, couplingStrength, pindividual, phyper = prepare_rates(param, model)
    llg, llgp = ll_hmm_coupled_hierarchical((rshared, rindividual, noiseshared, noiseindividual, couplingStrength), model.components, model.reporter, data.interval, data.trace, model.method)
    lhp = ll_hierarchy_c(pindividual, phyper)
    return llg + sum(lhp), vcat(llgp, lhp)
end

function loglikelihood(param, data::TraceData, model::GRSMgridmodel)
    r, noiseparams, pgrid = prepare_rates(param, model)
    ll_hmm_grid(r, noiseparams, pgrid[1], model.components.nT, model.Ngrid, model.components, model.reporter, data.interval, data.trace)
end

function loglikelihood(param, data::TraceData, model::GRSMgridhierarchicalmodel)
    r, noiseparams, pgrid = prepare_rates(param, model)
    ll_hmm_grid_hierarchical((r, noiseparams, pgrid[1]), model.components.nT, model.Ngrid, model.components, model.reporter, data.interval, data.trace)
end

function loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel)
    r = get_rates(param, model)
    llg, llgp = ll_hmm(r, model.components.tcomponents.nT, model.components.tcomponents, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    predictions = predictedRNA(r[1:num_rates(model)], model.components.mcomponents, model.nalleles, data.nRNA)
    logpredictions = log.(max.(predictions, eps()))
    return crossentropy(logpredictions, datahistogram(data)) + llg, vcat(-logpredictions, llgp)  # concatenate logpdf of histogram data with loglikelihood of traces
end

#### Trait model

# function extract_rates(r, rindices)
#     extracted = Dict{Symbol,Any}()

#     # Possible keys
#     possible_keys = (:rates, :noiseparams, :coupling, :grid)

#     for key in possible_keys
#         if hasfield(typeof(rindices), key)
#             extracted[key] = r[rindices[key]]
#         end
#     end
#     return extracted
# end


# function prepare_rates(rates, coupling, rindices)
#     r = Vector{Float64}(undef, length(rindices))
#     noiseparams = similar(r)
#     couplingStrength = similar(r)
#     for i in eachindex(rindices)
#         r[i] = rates[rindices[i].rates]
#         couplingStrength[i] = rates[rindices[i].coupling]
#     end
#     rateset = (rate=r, coupling=couplingStrength)
#     if :noise ∈ keys(rindices)
#         noiseparams = similar(r)
#         for i in eachindex(rindices)
#             noiseparams[i] = rates[rindices[i].noise]
#         end
#         rateset = merge(rateset, (noise=noiseparams,))
#     end

# end

# function reshape_rates(r, htrait::HierarchicalTrait)
#     reshape(r[1:htrait.sharedindices], h.nrates, htrait.nhypersets)
# end

# function prepare_rates(r, htrait::HierarchicalTrait)
#     # rates reshaped from a vector into a matrix with columns pertaining to shared params, hyper params and individual params 
#     # (shared parameters are considered to be hyper parameters without other hyper parameters (e.g. mean without variance))
#     rshared = reshape(r[1:htrait.sharedindices], h.nrates, htrait.nhypersets)
#     rindividual = reshape(r[htrait.individualindices], htrait.nrates, htrait.nindividuals)
#     rindividual[htrait.fittedshared, :] .= rshared[htrait.fittedshared, 1]
#     pindividual = reshape(param[htrait.paramindices], htrait.nparams, htrait.nindividuals)
#     phyper = Vector{Float64}[]
#     for i in htrait.hyperindices
#         push!(phyper, r[i])
#     end
#     return rshared, rindividual, pindividual, phyper
# end

# function prepare_rates(allparams, ctrait::CouplingTrait)
#     r = Vector{Float64}[]
#     noiseparams = Vector{Float64}[]
#     couplingStrength = Float64[]
#     for i in ctrait.raterange
#         push!(r, allparams[i, :])
#     end
#     for i in ctrait.noiserange
#         push!(noiseparams, allparams[i, :])
#     end
#     for i in ctrait.couplingrange
#         push!(couplingStrength, allparams[i, :])
#     end
# end


# function prepare_rates(param, model::AbstractGRSMtraitmodel)
#     r = get_rates(param, model)
#     traits = keys(model.traits)
#     if :hierarchical ∈ traits
#         rshared, rindividual, pindividual, phyper = prepare_rates(r, model.hierarchical)
#         if :coupled ∈ traits
#             rshared, rindividual, noiseparams, couplingStrength = prepare_rates(rshared, rindividual, model.coupling[1])
#             rateset = (rateshared=rshared, rindividual=rindividual, pindividual=pindividual, phyper=phyper, coupling=couplingStrength, noiseparams=noiseparams)
#         else
#             rateset = (rateshared=rshared, rateindividual=rindividual, pindividual=pindividual, phyper=phyper)
#         end
#     else
#         if :coupled ∈ traits
#             r, noiseparams, couplingStrength = prepare_rates(r, model.coupling[1])
#             rateset = (rates=r, noiseparams=noiseparams, coupling=couplingStrength)
#         else
#             r, noiseparams = prepare_rates(r, model)
#         end
#     end
#     if :grid ∈ traits
#         merge(rateset, (grid=r[model.traits.gridindices],))
#     end
#     return rateset
# end

# function prepare_nstates(tcomponents, traits)
#     if :grid ∈ keys(traits)
#         return tcomponents.nT
#     else
#         return tcomponents.nT, traits.grid.Ngrid
#     end
# end


# function loglikelihood(param, data::TraceData, model::AbstractGRSMtraitmodel)
#     r = prepare_rates(param, model)
#     N = prepare_nstates(model.component, model.traits)
#     ll_hmm_trait(r, N, model.components, model.reporter, data.interval, data.trace, model.method)
# end

# function loglikelihood(param, data::TraceRNAData, model::AbstractGRSMtraitmodel)
#     r = prepare_rates(param, model)
#     N = prepare_nstates(model.components.tcomponents, model.traits)
#     llg, llgp = ll_hmm_trai(r, N, model.components.tcomponents, model.reporter, data.interval, data.trace, model.method)
#     predictions = predictedRNA(r[1:num_rates(model)], model.components.mcomponents, model.nalleles, data.nRNA)
#     logpredictions = log.(max.(predictions, eps()))
#     return crossentropy(logpredictions, datahistogram(data)) + llg, vcat(-logpredictions, llgp)  # concatenate logpdf of histogram data with loglikelihood of traces
# end
#####


# Predicted histogram functions

function predictedRNA(r, mcomponents, nalleles, nRNA)
    M = make_mat_M(mcomponents, r)
    steady_state(M, mcomponents.nT, nalleles, nRNA)
end


"""
    predictedfn(param, data::RNAData, model::AbstractGMmodel)

Calculates the likelihood for a single RNA histogram.

# Arguments
- `param`: The model parameters.
- `data::RNAData`: The RNA data.
- `model::AbstractGMmodel`: The model.

# Returns
- `Vector{Float64}`: The steady-state probabilities for the RNA histogram.
"""
function predictedfn(param, data::RNAData, model::AbstractGmodel)
    r = get_rates(param, model)
    M = make_mat_M(model.components, r)
    steady_state(M, model.components.nT, model.nalleles, data.nRNA)
end




"""
    predictedfn(param, data::AbstractHistogramData, model::AbstractGmodel)

Calculates the likelihood for multiple histograms.

# Arguments
- `param`: The model parameters.
- `data::AbstractHistogramData`: The histogram data.
- `model::AbstractGmodel`: The model.

# Returns
- `Array{Float64,2}`: An array of likelihoods for the histograms.
"""
function predictedfn(param, data::AbstractHistogramData, model::AbstractGmodel)
    h = predictedarray(get_rates(param, model), data, model)
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
function predictedarray(r, data::RNAData{T1,T2}, model::AbstractGMmodel) where {T1<:Array,T2<:Array}
    h = Array{Array{Float64,1},1}(undef, length(data.nRNA))
    for i in eachindex(data.nRNA)
        M = make_mat_M(model.components[i], r[(i-1)*2*model.G+1:i*2*model.G])
        h[i] = steady_state(M, model.G, model.nalleles, data.nRNA[i])
    end
    trim_hist(h, data.nRNA)
end


"""
    predictedarray(r, data::RNAOnOffData, model::AbstractGmodel)

Calculates the likelihood for RNA On/Off data.

# Arguments
- `r`: The rates.
- `data::RNAOnOffData`: The RNA On/Off data.
- `model::AbstractGmodel`: The model.

# Returns
- `Array{Float64,1}`: The likelihoods for the RNA On/Off data.
"""
function predictedarray(r, data::RNAOnOffData, model::AbstractGmodel)
    #     if model.splicetype == "offdecay"
    #         # r[end-1] *= survival_fraction(nu, eta, model.R)
    #     end
    components = model.components
    onstates = model.reporter
    T = make_mat_T(components.tcomponents, r)
    TA = make_mat_TA(components.tcomponents, r)
    TI = make_mat_TI(components.tcomponents, r)
    M = make_mat_M(components.mcomponents, r)
    histF = steady_state(M, components.mcomponents.nT, model.nalleles, data.nRNA)
    modelOFF, modelON = offonPDF(data.bins, r, T, TA, TI, components.tcomponents.nT, components.tcomponents.elementsT, onstates)
    return [modelOFF, modelON, histF]
end


"""
    predictedarray(r, data::RNADwellTimeData, model::AbstractGmodel)

Calculates the likelihood array for RNA dwell time data.

# Arguments
- `r`: The rates.
- `data::RNADwellTimeData`: The RNA dwell time data.
- `model::AbstractGmodel`: The model.

# Description
This function calculates the likelihood array for RNA dwell time data using the specified model and rates. It constructs the transition matrices and calculates the probability density functions for the dwell times.

# Returns
- `Vector{Vector{Float64}}`: A vector of histograms representing the likelihoods for the dwell times.

"""
function predictedarray(r, data::RNADwellTimeData, model::AbstractGmodel)
    predictedarray(r, model.components, data.bins, model.reporter, data.DTtypes, model.nalleles, data.nRNA)
end

function predictedarray(r, data::DwellTimeData, model::AbstractGmodel)
    predictedarray(r, model.components, data.bins, model.reporter, data.DTtypes)
end

function predictedarray(r, data::DwellTimeData, model::GRSMcoupledmodel)
    r, coupling_strength, _ = prepare_rates(r, model)
    predictedarray(r, coupling_strength, model.components, data.bins, model.reporter, data.DTtypes)
end
"""
    predictedarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)

Calculates the likelihood array

# Arguments
- `r`: The rates.
- `G`: The number of gene states.
- `components`: The model components.
- `bins`: The bins for the data.
- `onstates`: The on states for the data.
- `dttype`: The data types (e.g., "OFF", "ON", "OFFG", "ONG").
- `nalleles`: The number of alleles.
- `nRNA`: The RNA data.

# Description
This function calculates the likelihood array for various types of data, including "OFF", "ON", "OFFG", and "ONG". It constructs the transition matrices and calculates the probability density functions for each data type. The function returns a vector of histograms representing the likelihoods.

# Returns
- `Vector{Vector{Float64}}`: A vector of histograms representing the likelihoods.
"""


function predictedarray(r, G::Int, tcomponents, bins, onstates, dttype)
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = normalized_nullspace(T)
    elementsG = tcomponents.elementsG
    TG = make_mat(elementsG, r, G)
    pssG = normalized_nullspace(TG)
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

function steady_state_dist(r, components, dt)
    pss = nothing
    pssG = nothing
    for b in union(dt)
        if b
            pssG = normalized_nullspace(make_mat_G(components, r))
        else
            pss = normalized_nullspace(make_mat_T(components, r))
        end
    end
    return (pss=pss, pssG=pssG, dt=dt)
end

function predictedarray(r, components::TDComponents, bins, reporter, dttype)
    sojourn, nonzeros = reporter
    dt = occursin.("G", dttype)
    p = steady_state_dist(r, components, dt)
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

function predictedarray(r, components::MTComponents, bins, reporter, dttype, nalleles, nRNA)
    M = make_mat_M(components.mcomponents, r)
    [steady_state(M, components.mcomponents.nT, nalleles, nRNA); predictedarray(r, components.tcomponents, bins, reporter, dttype)...]
end

function steady_state_dist(unit::Int, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, dt)
    pssG = nothing
    pss = nothing
    TC = nothing
    GC = nothing
    for b in union(dt)
        if b
            GC = make_mat_TC(coupling_strength, Gm, Gs, Gt, IG, sources, model)
            pssG = normalized_nullspace(GC)
        else
            TC = make_mat_TC(unit, T[unit], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
            pss = normalized_nullspace(TC)
        end
    end
    return (pss=pss, pssG=pssG, TC=TC, GC=GC)
end

function compute_dwelltime!(cache::CoupledDTCache, unit, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn, dt, nonzeros, bins)
    if dt
        if !haskey(cache.TC, dt)
            cache.TC[dt] = make_mat_TC(coupling_strength, Gm, Gs, Gt, IG, sources, model)
            cache.pss[dt] = normalized_nullspace(cache.TC[dt])
        end
        TCD = make_mat_TCD(cache.TC[dt], sojourn)
        return dwelltimePDF(bins, TCD, sojourn, init_S(sojourn, cache.TC[dt], cache.pss[dt]))
    else
        if !haskey(cache.TC, dt)
            cache.TC[dt] = make_mat_TC(unit, T[unit], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
            cache.pss[dt] = normalized_nullspace(cache.TC[dt])
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

function predictedarray(r, coupling_strength, components::TCoupledComponents{Vector{TDCoupledUnitComponents}}, bins, reporter, dttype)
    sojourn, nonzeros = reporter
    sources = components.sources
    model = components.model
    cache = CoupledDTCache(Dict(), Dict())
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

"""
    transform_rates(r, model)

Transforms rates to the real domain using log transformation.

# Arguments
- `r`: The rates to be transformed.
- `model`: The model, which can be of various types (e.g., `AbstractGmodel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Description
This function applies a log transformation to the rates to map them to the real domain. It supports various types of models, including `AbstractGmodel`, `AbstractGRSMmodel`, and `GRSMcoupledmodel`. The specific transformation method depends on the type of `model` provided.

# Methods
- `transform_rates(r, model::AbstractGmodel)`: Applies a log transformation to the rates.
- `transform_rates(r, model::AbstractGRSMmodel{Vector{Float64},HMMReporter})`: Applies a log transformation to the rates, with specific handling for HMM reporters.
- `transform_rates(r, model::GRSMcoupledmodel)`: Applies a log transformation to the rates, with specific handling for coupled models.

# Returns
- `Vector{Float64}`: The transformed rates.
"""
transform_rates(r, model::AbstractGmodel) = log.(r)

transform_rates(r, model::AbstractGRSMmodel{Vector{Float64},HMMReporter}) = transform_array(r, model.reporter.weightind, model.fittedparam, logv, logit)

transform_rates(r, model::GRSMcoupledmodel) = transform_array(r, length(model.rates), model.fittedparam, logv, log_shift1)

transform_rates(r, model::GRSMcoupledhierarchicalmodel) = transform_array(r, model.coupling[2], model.fittedparam, logv, log_shift1)


# function transform_rates(pin, model::GRSMcoupledmodel)
#     p = copy(pin)
#     for f in model.transformations
#         p[i] = f(p[i])
#     end
#     return p
# end

"""
    inverse_transform_rates(x,model)

Inverse transforms MH parameters on the real domain back to the rate domain.

# Arguments
- `x`: The parameters to be inverse transformed.
- `model`: The model, which can be of various types (e.g., `AbstractGmodel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Description
This function applies an inverse transformation to the parameters to map them back to the rate domain. It supports various types of models, including `AbstractGmodel`, `AbstractGRSMmodel`, and `GRSMcoupledmodel`. The specific inverse transformation method depends on the type of `model` provided.

# Methods
- `inverse_transform_rates(x, model::AbstractGmodel)`: Applies an exponential transformation to the parameters.
- `inverse_transform_rates(x, model::AbstractGRSMmodel{Vector{Float64},HMMReporter})`: Applies an inverse transformation to the parameters, with specific handling for HMM reporters.
- `inverse_transform_rates(x, model::GRSMcoupledmodel)`: Applies an inverse transformation to the parameters, with specific handling for coupled models.

# Returns
- `Vector{Float64}`: The inverse transformed parameters.

"""
inverse_transform_rates(p, model::AbstractGmodel) = exp.(p)

inverse_transform_rates(p, model::AbstractGRSMmodel{Vector{Float64},HMMReporter}) = transform_array(p, model.reporter.weightind, model.fittedparam, expv, invlogit)

inverse_transform_rates(p, model::GRSMcoupledmodel) = transform_array(p, length(model.rates), model.fittedparam, expv, invlog_shift1)

inverse_transform_rates(p, model::GRSMcoupledhierarchicalmodel) = transform_array(p, model.coupling[2], model.fittedparam, expv, invlog_shift1)


"""
    get_param(model)

Retrieves the fitted parameters from the model.

# Arguments
- `model`: The model, which can be of various types (e.g., `AbstractGmodel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Description
This function retrieves the fitted parameters from the model. It supports various types of models, including `AbstractGmodel`, `AbstractGRSMmodel`, and `GRSMcoupledmodel`. The specific retrieval method depends on the type of `model` provided.

# Methods
- `get_param(model::AbstractGmodel)`: Retrieves the log-transformed fitted parameters from the model.
- `get_param(model::AbstractGRSMmodel)`: Retrieves the transformed fitted parameters from the model.
- `get_param(model::GRSMcoupledmodel)`: Retrieves the transformed fitted parameters from the model.

# Returns
- `Vector{Float64}`: The fitted parameters.

"""
get_param(model::AbstractGmodel) = log.(model.rates[model.fittedparam])

get_param(model::AbstractGRSMmodel) = transform_rates(model.rates[model.fittedparam], model)

get_param(model::GRSMcoupledmodel) = transform_rates(model.rates[model.fittedparam], model)

get_param(model::GRSMcoupledhierarchicalmodel) = transform_rates(model.rates[model.fittedparam], model)

"""
    get_rates(param, model)

Replaces fitted rates with new values and returns the updated rates.

# Arguments
- `param`: The new parameter values.
- `model`: The model, which can be of various types (e.g., `AbstractGmodel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Description
This function replaces the fitted rates in the model with the provided new parameter values and returns the updated rates. It supports various types of models, including `AbstractGmodel`, `AbstractGRSMmodel`, and `GRSMcoupledmodel`. The specific replacement method depends on the type of `model` provided.

# Methods
- `get_rates(param, model::AbstractGmodel)`: Replaces and returns the updated rates for an `AbstractGmodel`.
- `get_rates(param, model::AbstractGRSMmodel)`: Replaces and returns the updated rates for an `AbstractGRSMmodel`.
- `get_rates(param, model::GRSMcoupledmodel)`: Replaces and returns the updated rates for a `GRSMcoupledmodel`.

# Returns
- `Vector{Float64}`: The updated rates.
"""
function get_rates(param, model::AbstractGmodel, inverse=true)
    r = copy_r(model)
    get_rates!(r, param, model, inverse)
    fixed_rates(r, model.fixedeffects)
end

function get_rates!(r, param, model, inverse)
    if inverse
        r[model.fittedparam] = inverse_transform_rates(param, model)
    else
        r[model.fittedparam] = param
    end
end

function get_rates(param::Vector{Vector}, model::AbstractGmodel, inverse=true)
    rv = copy_r(model)
    for r in rv
        get_rates!(r, param, model, inverse)
    end
    fixed_rates(r, model.fixedeffects)
end

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

function get_rates(param, model::GRSMcoupledmodel, inverse=true)
    r = copy_r(model)
    get_rates!(r, param, model, inverse)
    fixed_rates(r, model.fixedeffects)
end

function get_rates_hierarchical(param, model::AbstractGRSMhierarchicalmodel, inverse=true)
    r = get_rates(param, model, inverse)
    rshared, rindividual, pindividual, phyper = prepare_rates(r, param, model.hierarchy)
    vcat(vec(rshared), vec(rindividual))
end

# function get_rates(param, model::GRSMModel, inverse=true)
#     r = copy(model.rates)
#     if inverse
#         transformed_param = inverse_transform_rates(param, model)
#     else
#         transformed_param = param
#     end
#     r[model.fittedparam] = transformed_param
#     # Apply fixed effects if any
#     r = fixed_rates(r, model.fixedeffects)
#     return r
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
    logprior(param, model::AbstractGmodel)

Computes the log of the prior distribution for the model parameters.

# Arguments
- `param`: The model parameters.
- `model::AbstractGmodel`: The model, which includes the prior distributions for the rates.

# Description
This function computes the log of the prior distribution for the provided model parameters. It iterates through each prior distribution specified in the model and calculates the log probability density function (logpdf) for each parameter. The sum of these logpdf values is returned as the log prior.

# Returns
- `Float64`: The log of the prior distribution for the model parameters.
"""
function logprior(param, model::AbstractGmodel)
    d = model.rateprior
    p = 0
    for i in eachindex(d)
        p -= logpdf(d[i], param[i])
    end
    return p
end

"""
    num_rates(transitions, R, S, insertstep)
    num_rates(transitions, R::Tuple, S::Tuple, insertstep::Tuple)

compute number of transition rates (not counting noise parameters)
"""
function num_rates(transitions, R, S, insertstep)
    if R > 0
        if insertstep > R
            throw("insertstep>R")
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


function num_rates(transitions, R::Tuple, S::Tuple, insertstep::Tuple)
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
    num_rates(model)

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
    num_parameters(model::AbstractGmodel)

total number of parameters

# Arguments
- `model`: The model, which can be of various types (e.g., `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Description
This function calculates the total number of parameters for the provided model. It includes the number of transition rates, R steps, S splicing indicator, insert step, and the number of reporters

# Methods
- `total_parameters(model::AbstractGmodel)`: Returns the total number of parameters for the model.

# Returns
- `Int`: The total number of parameters.

"""

"""
    num_all_parameters(transitions, R::Int, S, insertstep, reporter, coupling=tuple(), grid=nothing)
    n = typeof(model.reporter) <: HMMReporterReporter ? model.reporter.n : 0
    num_rates(transitions, R, S, insertstep) + n

    `reporter` can either be a vector of noisepriors or of type HMMReporter
    

TBW
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
    c = isempty(coupling) ? 0 : coupling[5]
    g = isnothing(grid) ? 0 : grid
    # n = typeof(reporter) <: HMMReporter ? reporter.n : 0
    # num_rates(transitions, R, S, insertstep) + n
    num_rates(transitions, R, S, insertstep) + n + c + g
end

function num_all_parameters(transitions, R::Tuple, S::Tuple, insertstep::Tuple, reporter, coupling=tuple(), grid=nothing)
    n = 0
    for i in eachindex(R)
        n += num_all_parameters(transitions[i], R[i], S[i], insertstep[i], reporter[i])
    end
    c = isempty(coupling) ? 0 : coupling[5]
    g = isnothing(grid) ? 0 : grid
    n + c + g
end

function num_all_parameters(model::AbstractGmodel)
    n = typeof(model.reporter) <: HMMReporterReporter ? model.reporter.n : 0
    num_all_parameters(model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
end


