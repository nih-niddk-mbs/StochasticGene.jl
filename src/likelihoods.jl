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

"""
    hyper_distribution(p)

hyper parameter distribution for hierarchical model

# Arguments
- `p::Vector`: Vector of hyperparameters

#Description
This function returns a distribution for the hyper parameters of a hierarchical model.

# Returns
- `Distribution`: Distribution for the hierarchical model.
"""
function hyper_distribution(p)
    distribution_array(p[1], sigmalognormal(p[2]))
end

"""
    prepare_rates(param, model)

Extracts and reassembles parameters for use in likelihood calculations.

# Arguments
- `param`: The model parameters.
- `model`: The model, which can be of various types (e.g., `GRSMhierarchicalmodel`, `GRSMcoupledmodel`).

# Description
This function extracts and reassembles parameters from the provided model parameters for use in likelihood calculations. It supports hierarchical models and coupled models. The specific extraction and reassembly process depends on the type of `model` provided.

# Methods
- `prepare_rates(param, model::GRSMhierarchicalmodel)`: Extracts and reassembles parameters for a hierarchical model.
- `prepare_rates(param, model::GRSMcoupledmodel)`: Converts MCMC parameters into a form to compute likelihood for a coupled model.

# Returns
- `Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Int}}`: For hierarchical models, returns reshaped rates, parameters, and hyperparameters.
- `Tuple{Vector{Float64}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}, Vector{Int}}`: For coupled models, returns prepared rates and other necessary components for likelihood calculations.
"""
function prepare_rates(param, model::GRSMhierarchicalmodel)
    # rates reshaped from a vector into a matrix with columns pertaining to hyperparams and individuals 
    # (shared parameters are considered to be hyper parameters without other hyper parameters (e.g. mean without variance))
    h = Vector{Int}[]
    for i in model.pool.hyperindices
        push!(h, i)
    end
    r = reshape(get_rates(param, model)[model.pool.ratestart:end], model.pool.nrates, model.pool.nindividuals)
    p = reshape(param[model.pool.paramstart:end], model.pool.nparams, model.pool.nindividuals)
    return r, p, h
end

"""
    prepare_rates(param, model::GRSMcoupledmodel)

convert MCMC params into form to compute likelihood for coupled model
"""
function prepare_rates(param, model::GRSMcoupledmodel)
    rates = get_rates(param, model)
    n_noise = [r.n for r in model.reporter]
    sourceStates = [c.sourceState for c in model.components.modelcomponents]
    prepare_rates(rates, sourceStates, model.Gtransitions, model.G, model.R, model.S, model.insertstep, n_noise)
end

"""
    prepare_rates(rates, sourceStates, transitions, G::Tuple, R, S, insertstep, n_noise)

convert MCMC params into form to compute likelihood for coupled model

# Arguments
- `rates`: The model rates.
- `sourceStates`: The source states.
- `transitions`: The transitions.
- `G::Tuple`: The G steps.
- `R`: The R steps.
- `S`: The splicing indicator.
- `insertstep`: The R step where the reporter is inserted.
- `n_noise`: The number of noise parameters.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}, Vector{Vector{Float64}}}`: Prepared rates, coupling strength, and noise parameters.

"""
function prepare_rates(rates, sourceStates, transitions, G, R, S, insertstep, n_noise)
    r = Vector{Float64}[]
    noiseparams = Vector{Float64}[]
    couplingStrength = Float64[]
    j = 1
    for i in eachindex(G)
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + n_noise[i]
        push!(r, rates[j:j+n-1])
        j += n
    end
    for i in eachindex(G)
        if sourceStates[i] > 0
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

function prepare_rates(param, model::GRSMgridmodel)
    r = get_rates(param, model)
    r[model.raterange], r[model.noiserange], r[model.gridrange]
end
# Model loglikelihoods

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
    predictions = likelihoodfn(param, data, model)
    hist = datahistogram(data)
    logpredictions = log.(max.(predictions, eps()))
    return crossentropy(logpredictions, hist), -logpredictions
end

function loglikelihood(param, data::AbstractTraceData, model::AbstractGmodel)
    ll_hmm_2(get_rates(param, model), model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, model.reporter.offstates, data.interval, data.trace)
end

function loglikelihood(param, data::TraceData, model::GRSMcoupledmodel)
    r, couplingStrength, noiseparams = prepare_rates(param, model)
    ll_hmm_coupled(r, couplingStrength, noiseparams, model.components, model.reporter, data.interval, data.trace)
end

function loglikelihood(param, data::TraceData, model::GRSMgridmodel)
    r, noiseparams, pgrid = prepare_rates(param, model)
    ll_hmm_grid(r, noiseparams, pgrid[1], model.components.nT, model.Ngrid, model.components, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
end

function loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel)
    r = get_rates(param, model)
    # llg, llgp = ll_hmm(r, model.components.tcomponents.nT, model.components.tcomponents.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn, model.reporter.offstates, data.interval, data.trace)
    llg, llgp = ll_hmm(get_rates(param, model), model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, model.reporter.offstates, data.interval, data.trace)
    # M = make_mat_M(model.components.mcomponents, r[1:num_rates(model)])
    M = make_mat_MRG(model.components.mcomponents, r[1:num_rates(model)])
    logpredictions = log.(max.(steady_state(M, model.components.mcomponents.nT, model.nalleles, data.nRNA), eps()))
    return crossentropy(logpredictions, datahistogram(data)) + llg, vcat(-logpredictions, llgp)  # concatenate logpdf of histogram data with loglikelihood of traces
end

function loglikelihood(param, data::AbstractTraceData, model::GRSMhierarchicalmodel)
    r, p, hyper = prepare_rates(param, model)
    if model.method[2]
        # llg, llgp = ll_hmm_hierarchical_rateshared_background(r, model.components.nT, model.components.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn, model.reporter.offstates, data.interval, data.trace)
        base = model.S > 0 ? 3 : 2
        nT = model.G * base^model.R
        llg, llgp = ll_hmm_hierarchical_rateshared_background(r, model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, model.reporter.offstates, data.interval, data.trace)
    else
        llg, llgp = ll_hmm_hierarchical(r, model.components.nT, model.components, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    end
    # d = distribution_array(pm, psig)
    d = hyper_distribution(hyper)
    lhp = Float64[]
    for pc in eachcol(p)
        lhpc = 0
        for i in eachindex(pc)
            lhpc -= logpdf(d[i], pc[i])
        end
        push!(lhp, lhpc)
    end
    return llg + sum(lhp), vcat(llgp, lhp)
end

# Likelihood functions

"""
    likelihoodfn(param, data::RNAData, model::AbstractGMmodel)

Calculates the likelihood for a single RNA histogram.

# Arguments
- `param`: The model parameters.
- `data::RNAData`: The RNA data.
- `model::AbstractGMmodel`: The model.

# Returns
- `Vector{Float64}`: The steady-state probabilities for the RNA histogram.
"""
function likelihoodfn(param, data::RNAData, model::AbstractGMmodel)
    r = get_rates(param, model)
    M = make_mat_M(model.components, r)
    steady_state(M, model.G, model.nalleles, data.nRNA)
end
"""
    likelihoodfn(param, data::RNAData, model::AbstractGRSMmodel)

Calculates the likelihood for a single RNA histogram using a GRSM model.

# Arguments
- `param`: The model parameters.
- `data::RNAData`: The RNA data.
- `model::AbstractGRSMmodel`: The GRSM model.

# Returns
- `Vector{Float64}`: The steady-state probabilities for the RNA histogram.
"""
function likelihoodfn(param, data::RNAData, model::AbstractGRSMmodel)
    r = get_rates(param, model)
    M = make_mat_M(components.mcomponents, r)
    steady_state(M, components.mcomponents.nT, model.nalleles, data.nRNA)
end

"""
    likelihoodfn(param, data::AbstractHistogramData, model::AbstractGmodel)

Calculates the likelihood for multiple histograms.

# Arguments
- `param`: The model parameters.
- `data::AbstractHistogramData`: The histogram data.
- `model::AbstractGmodel`: The model.

# Returns
- `Array{Float64,2}`: An array of likelihoods for the histograms.
"""
function likelihoodfn(param, data::AbstractHistogramData, model::AbstractGmodel)
    h = likelihoodarray(get_rates(param, model), data, model)
    make_array(h)
end

"""
    likelihoodarray(r, data::RNAData{T1,T2}, model::AbstractGMmodel) where {T1<:Array, T2<:Array}

Calculates the likelihood for multiple histograms and returns an array of PDFs.

# Arguments
- `r`: The rates.
- `data::RNAData{T1,T2}`: The RNA data.
- `model::AbstractGMmodel`: The model.

# Returns
- `Array{Array{Float64,1},1}`: An array of PDFs for the histograms.
"""
function likelihoodarray(r, data::RNAData{T1,T2}, model::AbstractGMmodel) where {T1<:Array,T2<:Array}
    h = Array{Array{Float64,1},1}(undef, length(data.nRNA))
    for i in eachindex(data.nRNA)
        M = make_mat_M(model.components[i], r[(i-1)*2*model.G+1:i*2*model.G])
        h[i] = steady_state(M, model.G, model.nalleles, data.nRNA[i])
    end
    trim_hist(h, data.nRNA)
end


"""
    likelihoodarray(r, data::RNAOnOffData, model::AbstractGmodel)

Calculates the likelihood for RNA On/Off data.

# Arguments
- `r`: The rates.
- `data::RNAOnOffData`: The RNA On/Off data.
- `model::AbstractGmodel`: The model.

# Returns
- `Array{Float64,1}`: The likelihoods for the RNA On/Off data.
"""
function likelihoodarray(r, data::RNAOnOffData, model::AbstractGmodel)
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
    likelihoodarray(rin,data::RNAOnOffData,model::AbstractGRSMmodel)

"""

# function likelihoodarray(rin, data::RNAOnOffData, model::AbstractGRSMmodel)
#     r = copy(rin)
#     if model.splicetype == "offdecay"
#         r[end-1] *= survival_fraction(nu, eta, model.R)
#     end
#     likelihoodarray(r, data, model)
# end

"""
    likelihoodarray(r, data::RNADwellTimeData, model::AbstractGmodel)

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

function likelihoodarray(r, data::DwellTimeData, model::GRSMcoupledmodel)
    components = model.components
    for i in eachindex(data.DTtypes)
        if data.DTtypes[i] == "OFF"
            TD = make_mat(components.tcomponents.elementsTD[i], r, components.tcomponents.nT)
            nonzeros = nonzero_rows(TD)
            h = offtimePDF(data.bins[i], TD[nonzeros, nonzeros], nonzero_states(model.reporter.offstates[i], nonzeros), init_SI(r, model.reporter.offstates[i], components.tcomponents.elementsT, normalized_nullspace(TD), nonzeros))
        elseif data.DTtypes[i] == "ON"
            TD = make_mat(components.tcomponents.elementsTD[i], r, components.tcomponents.nT)
            h = ontimePDF(data.bins[i], TD, off_states(components.tcomponents.nT, model.reporter.offstates[i]), init_SA(r, model.reporter.offstates[i], components.tcomponents.elementsT, normalized_nullspace(TD)))
        end
        push!(hists, h)
    end

end

function likelihoodarray(r, data::RNADwellTimeData, model::AbstractGmodel)
    likelihoodarray(r, model.G, model.components, data.bins, model.reporter, data.DTtypes, model.nalleles, data.nRNA)
end

function likelihoodarray(r, data::DwellTimeData, model::AbstractGmodel)
    likelihoodarray(r, model.G, model.components, data.bins, model.reporter, data.DTtypes)
end


"""
    likelihoodarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)

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
function likelihoodarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)
    mcomponents = components.mcomponents
    M = make_mat_M(components.mcomponents, r)
    [steady_state(M, mcomponents.nT, nalleles, nRNA); likelihoodarray(r, G, components.tcomponents, bins, onstates, dttype)]
end

function likelihoodarray(r, G, tcomponents, bins, onstates, dttype)
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = normalized_nullspace(T)
    elementsTG = tcomponents.elementsTG
    TG = make_mat(elementsTG, r, G)
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
            h = offtimePDF(bins[i], TD, onstates[i], init_SI(r, onstates[i], elementsTG, pssG, collect(1:G)))
        elseif Dtype == "ONG"
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            h = ontimePDF(bins[i], TD, off_states(G, onstates[i]), init_SA(r, onstates[i], elementsTG, pssG))
        end
        push!(hists, h)
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
function transform_array(v::Array, index::Int, f1::Function, f2::Function)
    vcat(f1(v[1:index-1, :]), f2(v[index:end, :]))
end

function transform_array(v::Vector, index::Int, f1::Function, f2::Function)
    vcat(f1(v[1:index-1]), f2(v[index:end]))
end

function transform_array(v::Array, index::Int, mask::Vector, f1::Function, f2::Function)
    if index > 0 && index ∈ mask
        n = findfirst(index .== mask)
        return vcat(f1(v[1:n-1, :]), f2(v[n:n:end, :]))
    else
        return f1(v)
    end
end

function transform_array(v::Vector, index::Int, mask::Vector, f1::Function, f2::Function)
    if index > 0 && index ∈ mask
        n = findfirst(index .== mask)
        return vcat(f1(v[1:n-1]), f2(v[n:end]))
    else
        return f1(v)
    end
end

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
    n = 0
    for i in eachindex(R)
        n += num_rates(transitions[i], R[i], S[i], insertstep[i])
    end
    n
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
    total_parameters(model::AbstractGmodel)

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
function num_total_parameters(model::AbstractGmodel)
    n = typeof(model.reporter) <: HMMReporterReporter ? model.reporter.n : 0
    num_rates(model.Gtransitions, model.R, model.S, model.insertstep) + n
end