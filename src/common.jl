# This file is part of StochasticGene.jl   

###  common.jl

# Data types
"""
    AbstractExperimentalData

abstract type for experimental data
"""
abstract type AbstractExperimentalData end
"""
    AbstractSampleData

abstract type for data in the form of samples
"""
abstract type AbstractSampleData <: AbstractExperimentalData end
"""
    AbstractHistogramData

abstract type for data in the form of a histogram (probability distribution)
"""
abstract type AbstractHistogramData <: AbstractExperimentalData end
"""
    AbstractTraceData
    
abstract type for intensity time series data
"""
abstract type AbstractTraceData <: AbstractExperimentalData end

"""
  AbstractRNAData{hType}

abstract type for steady state RNA histogram data
"""
abstract type AbstractRNAData{hType} <: AbstractHistogramData end


"""
    AbstractTraceHistogramData
    
abstract type for intensity time series data with RNA histogram data
"""
abstract type AbstractTraceHistogramData <: AbstractExperimentalData end

# Data structures 
#
# Do not use underscore "_" in label

"""
    Data structures

arguments:
label: label for the data set
gene: gene name (case sensitive)
nRNA: length of histogram
histRNA: RNA histograms
bins: number of live cell recording time bins
OFF:  OFF time probability density
ON:: ON time probability density
"""
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    label::String
    gene::String
    nRNA::nType
    histRNA::hType
end
struct RNAOnOffData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Vector
    bins::Vector
    ON::Vector
    OFF::Vector
end
struct RNADwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Array
    bins::Vector{Vector}
    DwellTimes::Vector{Vector}
    DTtypes::Vector
end
struct TraceData{traceType} <: AbstractTraceData
    label::String
    gene::String
    interval::Float64
    trace::traceType
end
struct TraceNascentData{traceType} <: AbstractTraceData
    label::String
    gene::String
    interval::Float64
    trace::traceType
    nascent::Float64
end
struct TraceRNAData{traceType,hType} <: AbstractTraceHistogramData
    label::String
    gene::String
    interval::Float64
    trace::traceType
    nRNA::Int
    histRNA::hType
end

# Model structures

"""
    Abstract model types
"""
abstract type AbstractModel end
abstract type AbstractGmodel <: AbstractModel end
abstract type AbstractGMmodel <: AbstractGmodel end
abstract type AbstractGRSMmodel{ReporterType} <: AbstractGmodel end

"""
    Model structures

fields:
- `rates`: transition rates
- `Gtransitions`:  tuple of vectors of G state transitions
- `G`: number of G steps
- `R`: number of R steps
- `S`: indicator for splicing, 0 no splicing, > 1 splicing
- `insertstep`: R step where reporter is inserted (first step where reporter is visible)
- `nalleles`: number of alleles producing RNA
- `splicetype`: choices are "", "offeject", "offdecay"
- `rateprior`: prior distribution for rates
- `proposal`: MCMC proposal distribution
- `fittedparam`: indices of rates to be fitted
- `fixedeffects`: indices of rates that are fixed to each other, in the form of a 2 tuple of vectors
    with index 1 the tied index vector and 2 the corresponding fitted index vector
- `fixedeffects`: tuple of vectors of rates that are locked together
- `method`: numerical method for solving Master equation
- `components`: rate marix components
-` reporter`: vector of reporters or sojorn states (onstates) or vectors of vectors depending on model and data

"""
struct GMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGMmodel
    rates::RateType
    Gtransitions::Tuple
    G::Int
    nalleles::Int
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end
struct GMfixedeffectsmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGMmodel
    rates::RateType
    Gtransitions::Tuple
    G::Int
    nalleles::Int
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end

struct GRSMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{ReporterType}
    rates::RateType
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
    nalleles::Int
    splicetype::String
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end

struct GRSMfixedeffectsmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{ReporterType}
    rates::RateType
    Gtransitions::Tuple
    G::Int
    R::Int
    S::Int
    insertstep::Int
    nalleles::Int
    splicetype::String
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end

"""
hmmReporter

structure for reporters

- `n`: number of noise parameters
- `per_state`: number of reporters per state
- `probfn`: noise distribution e.g. prob_GaussianMixture
- `weightind`: index for mixture model bias parameter (restricted to range [0,1])
"""
struct hmmReporter
    n::Int
    per_state::Vector{Int}
    probfn::Function
    weightind::Int
end

"""
    print_model(model::AbstractModel)

print all fields of model
"""
function print_model(model::AbstractModel)
    for fname in fieldnames(model)
        println("$fname =", getfield(model, fname))
    end
end

"""
Abstract Option types for fitting methods
"""
abstract type Options end
abstract type Results end

# Model and Data dependent functions

"""
datahistogram(data)
Return the RNA histogram data as one vector
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

datapdf(data::AbstractRNAData{Array{Float64,1}}) = normalize_histogram(data.histRNA)
datapdf(data::RNAOnOffData) = [normalize_histogram(data.OFF); normalize_histogram(data.ON); normalize_histogram(data.histRNA)]
datapdf(data::AbstractTraceHistogramData) = normalize_histogram(data.histRNA)
datapdf(data::RNADwellTimeData) =[data.histRNA;(make_array(data.DwellTimes))]

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

# Model loglikelihoods

"""
loglikelihood(param,data::AbstractHistogramData,model)

returns negative loglikelihood of all data and vector of the prediction histogram negative loglikelihood
Calls likelihoodfn and datahistogram
provided for each data and model type
"""
function loglikelihood(param, data::AbstractHistogramData, model::AbstractGmodel)
    predictions = likelihoodfn(param, data, model)
    hist = datahistogram(data)
    logpredictions = log.(max.(predictions, eps()))
    return crossentropy(logpredictions, hist), -logpredictions
end
"""
    loglikelihood(param,data::AbstractTraceData,model::GMmodel)

negative loglikelihood of combined time series traces and each trace
"""
function loglikelihood(param, data::AbstractTraceData, model::AbstractGmodel)
    ll_hmm(get_rates(param, model), model.components.nT, model.components.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
end
"""
    loglikelihood(param, data::TraceRNAData{Float64}, model::AbstractGmodel)

negative loglikelihood of trace data with nascent RNA FISH active fraction (stored in data.histRNA field)
"""
function loglikelihood(param, data::TraceNascentData, model::AbstractGmodel)
    ll_hmm(get_rates(param, model), model.components.nT, model.components.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace, data.nascent)
end

"""
loglikelihood(param, data::TraceRNAData{Vector{Float64}}, model::AbstractGRSMmodel)

negative loglikelihood of time series traces and mRNA FISH steady state histogram
"""
function loglikelihood(param, data::TraceRNAData, model::AbstractGRSMmodel)
    r = get_rates(param, model)
    llg, llgp = ll_hmm(r, model.components.tcomponents.nT, model.components.tcomponents.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    M = make_mat_M(model.components.mcomponents, r[1:num_rates(model)])
    logpredictions = log.(max.(steady_state(M, model.components.mcomponents.nT, model.nalleles, data.nRNA), eps()))
    return crossentropy(logpredictions, datahistogram(data)) + llg, vcat(-logpredictions, llgp)  # concatenate logpdf of histogram data with loglikelihood of traces
end


# Likelihood functions

"""
    likelihoodfn(param,data::RNAData,model::AbstractGMmodel)

likelihood for single RNA histogram
"""
function likelihoodfn(param, data::RNAData, model::AbstractGMmodel)
    r = get_rates(param, model)
    M = make_mat_M(model.components, r)
    steady_state(M, model.G, model.nalleles, data.nRNA)
end
"""
    likelihoodfn(param,data::AbstractHistogramArrayData,model::AbstractGmodel)

likelihood for multiple histograms
"""
function likelihoodfn(param, data::AbstractHistogramData, model::AbstractGmodel)
    h = likelihoodarray(get_rates(param, model), data, model)
    make_array(h)
end

"""
    likelihoodarray(r, data::RNAData{T1,T2}, model::AbstractGMmodel) where {T1<:Array,T2<:Array}

likelihood for multiple histograms, returns an array of pdfs
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

    #     if model.splicetype == "offdecay"
    #         # r[end-1] *= survival_fraction(nu, eta, model.R)
    #     end
"""
function likelihoodarray(r, data::RNAOnOffData, model::AbstractGmodel)
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

function likelihoodarray(rin, data::RNAOnOffData, model::AbstractGRSMmodel)
    r = copy(rin)
    if model.splicetype == "offdecay"
        r[end-1] *= survival_fraction(nu, eta, model.R)
    end
    likelihoodarray(r,data,model)
end

"""
    likelihoodarray(rin,data::RNADwellTimeData,model::AbstractGRSMmodel)

likelihood of an array of dwell time histograms
"""
function likelihoodarray(rin, data::RNADwellTimeData, model::AbstractGRSMmodel)
    r = copy(rin)
    tcomponents = model.components.tcomponents
    onstates = model.reporter
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = normalized_nullspace(T)
    hists = Vector[]
    M = make_mat_M(model.components.mcomponents, r)
    histF = steady_state(M, model.components.mcomponents.nT, model.nalleles, data.nRNA)
    push!(hists, histF)
    for (i, Dtype) in enumerate(data.DTtypes)
        TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
        if Dtype == "OFF"
            nonzeros = nonzero_rows(TD)
            h = offtimePDF(data.bins[i], TD[nonzeros, nonzeros], nonzero_states(onstates[i], nonzeros), init_SI(r, onstates[i], elementsT, pss, nonzeros))
        elseif Dtype == "ON"
            h = ontimePDF(data.bins[i], TD, off_states(tcomponents.nT, onstates[i]), init_SA(r, onstates[i], elementsT, pss))
        else
            # h = ontimePDF(data.bins[i], TD, [1,0], [0,1])
            h = normalize_histogram(data.DwellTimes[3])
        end
        push!(hists, h)
    end
    return hists
end

"""
    transform_array(v::Array, index::Int, f1::Function, f2::Function)

apply function f1 to actionrray v at indices up to index and f2 after index
"""
function transform_array(v::Array, index::Int, f1::Function, f2::Function)
    vcat(f1(v[1:index-1, :]), f2(v[index:end, :]))
end

"""
    transform_array(v::Array, index, mask, f1, f2)

apply function f1 to array v at indices up to index and f2 after index accounting for application of mask (which changes indexing)
"""
function transform_array(v::Array, index::Int, mask::Vector, f1::Function, f2::Function)
    if index ∈ mask
        n = findfirst(index .== mask)
        if typeof(v) <: Vector
            return vcat(f1(v[1:n-1]), f2(v[n:end]))
        else
            return vcat(f1(v[1:n-1, :]), f2(v[n:end, :]))
        end
    else
        return f1(v)
    end
end



"""
    transform_rates(r,model::AbstractGmodel)

log transform rates to real domain
"""
transform_rates(r, model::AbstractGmodel) = log.(r)

# function transform_rates(r, model::AbstractGRSMmodel{hmmReporter}) 
#     n = num_rates(model)
#     [log.(r[1:n + model.reporter.weightind-1]); logit(r[n + model.reporter.weightind:end])]
# end

transform_rates(r, model::AbstractGRSMmodel{hmmReporter}) = transform_array(r, model.reporter.weightind, model.fittedparam, logv, logit)


"""
    inverse_transform_rates(x,model::AbstractGmodel)

inverse transform MH parameters on real domain back to rate domain

"""
inverse_transform_rates(p, model::AbstractGmodel) = exp.(p)

inverse_transform_rates(p, model::AbstractGRSMmodel{hmmReporter}) = transform_array(p, model.reporter.weightind, model.fittedparam, expv, invlogit)


# function inverse_transform_rates(p::Vector, model::AbstractGRSMmodel{hmmReporter})
#     n = num_rates(model)
#     vcat(exp.(p[1:n + model.reporter.weightind-1]), invlogit(p[n + model.reporter.weightind:end]))
# end
# function inverse_transform_rates(p::Vector, model::AbstractGRSMmodel{hmmReporter})
#     wind = num_rates(model) + model.reporter.weightind
#     if wind ∈ model.fittedparam
#         n = findfirst(wind .== model.fittedparam)
#         return vcat(exp.(p[1:n-1]), invlogit(p[n:end]))
#     else
#        return exp.(p)
#     end
# end

# function inverse_transform_rates(p::Matrix, model::AbstractGRSMmodel{hmmReporter})
#     wind = num_rates(model) + model.reporter.weightind
#     if wind ∈ model.fittedparam
#         n = findfirst(wind .== model.fittedparam)
#         return vcat(exp.(p[1:n-1,:]), invlogit(p[n:end,:]))
#     else
#        return exp.(p)
#     end
# end

"""
    get_param(model)

get fitted parameters from model
"""
get_param(model::AbstractGmodel) = log.(model.rates[model.fittedparam])

get_param(model::AbstractGRSMmodel) = transform_rates(model.rates[model.fittedparam], model)

# get_param(r, model::AbstractGmodel) = transform_rates(r[model.fittedparam], model)

# get_param(r, model::AbstractGRSMmodel) = transform_rates(r, model)[model.fittedparam]


"""
    get_rates(param,model)

replace fitted rates with new values and return
"""
function get_rates(param, model::AbstractGmodel, inverse=true)
    r = copy_r(model)
    if inverse
        r[model.fittedparam] = inverse_transform_rates(param, model)
    else
        r[model.fittedparam] = param
    end
    return r
end


# function get_rates(param, model::AbstractGRSMmodel{hmmReporter}, inverse=true)
#     if inverse
#         p = transform_rates(model.rates,model)
#         p[model.fittedparam] = param
#         r = inverse_transform_rates(p,model)
#     else
#         r = copy_r(model)
#         r[model.fittedparam] = param
#     end
#     return r
# end

get_rates(param, model::GRSMfixedeffectsmodel; inverse=true) = fixed_rates(param, model, inverse)

"""
    fixed_rates(param,model,inverse)

apply fixed effects on rates
"""
function fixed_rates(param, model::GRSMfixedeffectsmodel, inverse)
    r = get_rates(param, model, inverse)
    for effect in model.fixedeffects
        r[effect[2:end]] .= r[effect[1]]
    end
    return r
end


"""
    copy_r(model)

copy rates from model structure
"""
copy_r(model) = copy(model.rates)

"""
setr(r,model)

"""
function setr(r, model::AbstractGRSMmodel)
    n = model.G - 1
    nr = model.R
    eta = get_eta(r, n, nr)
    r[2*n+1+nr+1:2*n+1+nr+nr] = eta
    r
end
"""
    logprior(param,model::AbstractGMmodel)

compute log of the prior
"""
function logprior(param, model::AbstractGmodel)
    d = model.rateprior
    p = 0
    for i in eachindex(param)
        p -= logpdf(d[i], param[i])
    end
    return p
end

"""
    num_rates(transitions, R, S, insertstep)

compute number of transition rates 
"""
function num_rates(transitions, R, S, insertstep)
    if R > 0
        if insertstep > R
            throw("insertstep>R")
        end
        if S == 0
            insertstep = 1
        else
            S = R
        end
        return length(transitions) + 1 + R + S + 1 - insertstep + 1
    else
        return length(transitions) + 2
    end
end

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
