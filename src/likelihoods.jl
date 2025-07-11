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

function prepare_rates_coupled(rates::Vector, nrates::Vector{Int}, couplingindices)
    r = Vector{Vector{Float64}}(undef, length(nrates))
    k = 1
    for i in eachindex(r)
        n = nrates[i]
        r[i] = rates[k:k+n-1]
        k += n
    end
    return r, rates[couplingindices]
end
"""
prepare_rates_noiseparams(rates, nrates, reporter)
"""
function prepare_rates_noiseparams(rates::AbstractArray, nrates::Int, reporter::HMMReporter)
    rates[1:nrates], rates[nrates+1:nrates+reporter.n]
end

function prepare_rates_noiseparams(rates::Vector{T}, nrates::Int, reporter::HMMReporter) where {T<:AbstractArray}
    r = Vector{Vector{Float64}}(undef, length(rates))
    noiseparams = Vector{Vector{Float64}}(undef, length(rates))
    k = 1
    for i in eachindex(r)
        r[i], noiseparams[i] = prepare_rates_noiseparams(rates[i], nrates, reporter)
    end
    return r, noiseparams
end

function prepare_rates_noiseparams(rates::Vector, nrates::Vector{Int}, reporter::Vector{HMMReporter})
    r = Vector{Vector{Float64}}(undef, length(nrates))
    noiseparams = Vector{Vector{Float64}}(undef, length(nrates))
    k = 1
    for i in eachindex(r)
        n = nrates[i]
        r[i] = rates[k:k+n-1]
        k += n
        noiseparams[i] = rates[k:k+reporter[i].n-1]
        k += reporter[i].n
    end
    return r, noiseparams
end

function prepare_rates_noiseparams(rates::Vector{T}, nrates::Vector{Int}, reporter::Vector{HMMReporter}) where {T<:AbstractArray}
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

"""
    prepare_coupling(rates, couplingindices)
    
"""

function prepare_coupling(rates::Vector{Float64}, couplingindices)
    rates[couplingindices]
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

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling,)}}
    r = get_rates(param, model)
    prepare_rates_coupled(r, model.nrates, model.reporter, model.trait.coupling.couplingindices)
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical,)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:hierarchical,)}}
    r = get_rates(param, model)
    prepare_rates_hierarchical(r, param, model.nrates, model.trait.hierarchical, model.reporter)
end

"""
    prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical)}}
"""
function prepare_rates(param, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling, :hierarchical)}}
    r = get_rates(param, model)
    prepare_rates_coupled_hierarchical(r, param, model.nrates, model.trait.coupling, model.trait.hierarchical, model.reporter)
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
function loglikelihood(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel)
    predictions = predictedfn(param, data, model)
    hist = datahistogram(data)
    logpredictions = log.(max.(predictions, eps()))
    return sum(hist .* logpredictions), hist .* logpredictions  # Convention: return log-likelihoods
end

function loglikelihood(param, data::RNACountData, model::AbstractGeneTransitionModel)
    predictions = predictedfn(param, data, model)
    logpredictions = Array{Float64,1}(undef, length(data.countsRNA))
    for k in eachindex(data.countsRNA)
        p = technical_loss_at_k(data.countsRNA[k], predictions, 1., data.nRNA)
        logpredictions[k] = log(max(p, eps()))
    end
    return sum(logpredictions), logpredictions  # Convention: return log-likelihoods
end

# AD-friendly version (no in-place mutation)
function loglikelihood_ad(param, data::RNACountData, model::AbstractGeneTransitionModel)
    predictions = predictedfn(param, data, model)
    logpredictions = [
        log(max(technical_loss_at_k(data.countsRNA[k], predictions, 1., data.nRNA), eps()))
        for k in eachindex(data.countsRNA)
    ]
    return sum(logpredictions), logpredictions
end

# Helper to get the right component
get_components(model::AbstractGRSMmodel, data::AbstractTraceHistogramData) = model.components.tcomponents
get_components(model::AbstractGRSMmodel, data::AbstractTraceData) = model.components
# (Add more as needed)

function ll_hmm_trace(param, data, model::AbstractGRSMmodel)
    r = prepare_rates(param, model)
    components = get_components(model, data)
    if !isnothing(model.trait) && haskey(model.trait, :grid)
        ll_hmm(r, model.trait.grid.ngrid, components, model.reporter, data.interval, data.trace, model.method)
    else
        ll_hmm(r, components, model.reporter, data.interval, data.trace, model.method)
    end
end

function loglikelihood(param, data::AbstractTraceData, model::AbstractGRSMmodel)
    ll_hmm_trace(param, data, model)
end

function loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel)
    llg, llgp = ll_hmm_trace(param, data, model)
    r = get_rates(param, model)
    predictions = predictedRNA(r[1:num_rates(model)], model.components.mcomponents, model.nalleles, data.nRNA)
    logpredictions = log.(max.(predictions, eps()))
    hist = datahistogram(data)
    return sum(hist .* logpredictions) + llg, vcat(hist .* logpredictions, llgp)  # Convention: all log-likelihoods
end


# Predicted histogram functions

function predictedRNA(r, mcomponents, nalleles, nRNA)
    M = make_mat_M(mcomponents, r)
    steady_state(M, mcomponents.nT, nalleles, nRNA)
end


"""
    predictedfn(param, data::RNAData, model::AbstractGeneTransitionModel)

Calculates the likelihood for a single RNA histogram.

# Arguments
- `param`: The model parameters.
- `data::RNAData`: The RNA data.
- `model::AbstractGeneTransitionModel`: The model.

# Returns
- `Vector{Float64}`: The steady-state probabilities for the RNA histogram.
"""
function predictedfn(param, data::AbstractRNAData, model::AbstractGeneTransitionModel)
    r = get_rates(param, model)
    M = make_mat_M(model.components, r)
    steady_state(M, model.components.nT, model.nalleles, data.nRNA)
end




"""
    predictedfn(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel)

Calculates the likelihood for multiple histograms.

# Arguments
- `param`: The model parameters.
- `data::AbstractHistogramData`: The histogram data.
- `model::AbstractGeneTransitionModel`: The model.

# Returns
- `Array{Float64,2}`: An array of likelihoods for the histograms.
"""
function predictedfn(param, data::AbstractHistogramData, model::AbstractGeneTransitionModel)
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

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, data::RNAData{T1,T2}, model::AbstractGMmodel) where {T1<:Array,T2<:Array}
    h = [
        steady_state(
            make_mat_M(model.components[i], r[(i-1)*2*model.G+1:i*2*model.G]),
            model.G, model.nalleles, data.nRNA[i]
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
function predictedarray(r, data::RNAOnOffData, model::AbstractGeneTransitionModel)
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
function predictedarray(r, data::RNADwellTimeData, model::AbstractGeneTransitionModel)
    predictedarray(r, model.components, data.bins, model.reporter, data.DTtypes, model.nalleles, data.nRNA)
end

function predictedarray(r, data::DwellTimeData, model::AbstractGeneTransitionModel)
    predictedarray(r, model.components, data.bins, model.reporter, data.DTtypes)
end

function predictedarray(r, data::DwellTimeData, model::AbstractGRSMmodel{T}) where {T<:NamedTuple{(:coupling,)}}
    r, coupling_strength = prepare_rates_coupled(r, model.nrates, model.trait.coupling.couplingindices)
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

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, G::Int, tcomponents, bins, onstates, dttype)
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = normalized_nullspace(T)
    elementsG = tcomponents.elementsG
    TG = make_mat(elementsG, r, G)
    pssG = normalized_nullspace(TG)
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

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, components::TDComponents, bins, reporter, dttype)
    sojourn, nonzeros = reporter
    dt = occursin.("G", dttype)
    p = steady_state_dist(r, components, dt)
    [
        dt[i] ?
            dwelltimePDF(bins[i], make_mat(components.elementsTD[i], r, components.TDdims[i]), sojourn[i], init_S(r, sojourn[i], components.elementsG, p.pssG)) :
            dwelltimePDF(bins[i], make_mat(components.elementsTD[i], r, components.TDdims[i])[nonzeros[i], nonzeros[i]], nonzero_states(sojourn[i], nonzeros[i]), init_S(r, sojourn[i], components.elementsT, p.pss, nonzeros[i]))
        for i in eachindex(sojourn)
    ]
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

# AD-friendly version (no in-place mutation)
function predictedarray_ad(r, coupling_strength, components::TCoupledComponents{Vector{TDCoupledUnitComponents}}, bins, reporter, dttype)
    sojourn, nonzeros = reporter
    sources = components.sources
    model = components.model
    T, TD, Gm, Gt, Gs, IG, IR, IT = make_matvec_C(components, r)
    [
        [
            compute_dwelltime!(
                CoupledDTCache(Dict(), Dict()), α, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model,
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


function apply_transform(p::AbstractVector, f)
    map((f, x) -> f(x), f, p)
end

function apply_transform(p::AbstractMatrix, f)
    hcat([map((f, x) -> f(x), f, p[:, i]) for i in 1:size(p, 2)]...)
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

Get fitted parameters from model.

# Arguments
- `model`: The model, which can be of various types (e.g., `AbstractGeneTransitionModel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Returns
- Log-transformed fitted parameters.
"""
get_param(model::AbstractGeneTransitionModel) = log.(model.rates[model.fittedparam])

function get_param(model::AbstractGRSMmodel) 
    transform_rates(model.rates[model.fittedparam], model)
end

"""
    get_rates(param, model::AbstractGeneTransitionModel)

Get rates from parameters.

# Arguments
- `param`: Vector of parameters.
- `model`: The model, which can be of various types (e.g., `AbstractGeneTransitionModel`, `AbstractGRSMmodel`, `GRSMcoupledmodel`).

# Returns
- Updated rates.
"""

function get_rates!(r, param, model, inverse)
    if inverse
        r[model.fittedparam] = inverse_transform_params(param, model)
    else
        r[model.fittedparam] = param
    end
end

function get_rates(param, model::AbstractGeneTransitionModel, inverse=true)
    r = copy_r(model)
    get_rates!(r, param, model, inverse)
    fixed_rates(r, model.fixedeffects)
end

function get_rates(param::Vector{Vector}, model::AbstractGeneTransitionModel, inverse=true)
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
    c = isempty(coupling) ? 0 : coupling[5]
    g = isnothing(grid) ? 0 : 1
    num_rates(transitions, R, S, insertstep) + n + c + g
end

function num_all_parameters(transitions, R::Tuple, S::Tuple, insertstep::Tuple, reporter, coupling=tuple(), grid=nothing)
    n = 0
    for i in eachindex(R)
        n += num_all_parameters(transitions[i], R[i], S[i], insertstep[i], reporter[i])
    end
    c = isempty(coupling) ? 0 : coupling[5]
    g = isnothing(grid) ? 0 : 1
    n + c + g
end

function num_all_parameters(model::AbstractGeneTransitionModel)
    num_all_parameters(model.Gtransitions, model.R, model.S, model.insertstep, model.reporter)
end

function num_fitted_core_params(model::AbstractGeneTransitionModel)
    count(x -> x <= num_all_parameters(model.Gtransitions, model.R, model.S, model.insertstep, model.reporter), model.fittedparam)
end
