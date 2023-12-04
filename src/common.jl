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
abstract type AbstractGRSMmodel{RateType,ReporterType} <: AbstractGmodel end

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
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end

struct GRSMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
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

struct GRSMhierarchicalmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{RateType,ReporterType}
    rates::RateType
    pool::Pool
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

struct Pool
    npoolsets::Int
    npoolparams::Int
    npoolrates::Int
    nindividualparams::Int
    nindividuals::Int
end

"""
HMMReporter

structure for reporters

- `n`: number of noise parameters
- `per_state`: number of reporters per state
- `probfn`: noise distribution e.g. prob_GaussianMixture
- `weightind`: index for mixture model bias parameter (restricted to range [0,1])
"""
struct HMMReporter
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

function loglikelihood(param, data::AbstractTraceData, model::GRSMhierarchicalmodel)
    r, p, pm, psig = prepare_params(param, model)
    llg, llgp = ll_hmm_hierarchical(r, model.components.nT, model.components.elementsT, model.reporter.n, model.reporter.per_state, model.reporter.probfn, data.interval, data.trace)
    d = distribution_array(pm, psig)
    lhp = Float64[]
    for pc in eachcol(p)
        lhpc = 0
        for i in eachindex(pc)
            lhpc += logpdf(d[i], pc[i])
        end
        push!(lhp,lhpc)
    end
    return llg + sum(lhp), vcat(llgp, lhp)
end

function prepare_params(param, model)
    r = reshape(get_rates(param, model), model.nrates, model.nsets + model.nindividuals)
    pm = param[1:model.npars]
    if model.nsets == 2
        psig = sigmalognormal(param[model.nrates+1:model.nrates+model.npars])
    end
    p = reshape(param[end-model.npars*model.nindividuals+1:end],model.npars,model.nindividuals)
    return r[:, model.nsets+1:end], p, pm, psig
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

# function likelihoodarray(rin, data::RNAOnOffData, model::AbstractGRSMmodel)
#     r = copy(rin)
#     if model.splicetype == "offdecay"
#         r[end-1] *= survival_fraction(nu, eta, model.R)
#     end
#     likelihoodarray(r, data, model)
# end

"""
    likelihoodarray(rin,data::RNADwellTimeData,model::AbstractGRSMmodel)

likelihood of an array of dwell time histograms
"""
function likelihoodarray(r, data::RNADwellTimeData, model::AbstractGRSMmodel)
    # likelihoodarray(r,model.G,model.components,data.bins,model.reporter,data.DTtypes,model.nalleles,data.nRNA)
    G = model.G
    tcomponents = model.components.tcomponents
    onstates = model.reporter
    elementsT = tcomponents.elementsT
    elementsTG = tcomponents.elementsTG
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = normalized_nullspace(T)
    TG = make_mat(elementsTG, r, G)
    pssG = normalized_nullspace(TG)
    hists = Vector[]
    M = make_mat_M(model.components.mcomponents, r)
    histF = steady_state(M, model.components.mcomponents.nT, model.nalleles, data.nRNA)
    push!(hists, histF)
    for (i, Dtype) in enumerate(data.DTtypes)
        if Dtype == "OFF"
            TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
            nonzeros = nonzero_rows(TD)
            h = offtimePDF(data.bins[i], TD[nonzeros, nonzeros], nonzero_states(onstates[i], nonzeros), init_SI(r, onstates[i], elementsT, pss, nonzeros))
        elseif Dtype == "ON"
            TD = make_mat(tcomponents.elementsTD[i], r, tcomponents.nT)
            h = ontimePDF(data.bins[i], TD, off_states(tcomponents.nT, onstates[i]), init_SA(r, onstates[i], elementsT, pss))
        elseif Dtype == "OFFG"
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            h = offtimePDF(data.bins[i], TD, onstates[i], init_SI(r, onstates[i], elementsTG, pssG, collect(1:G)))
        elseif Dtype == "ONG"
            TD = make_mat(tcomponents.elementsTD[i], r, G)
            h = ontimePDF(data.bins[i], TD, off_states(G, onstates[i]), init_SA(r, onstates[i], elementsTG, pssG))
        end
        push!(hists, h)
    end
    return hists
end

"""
    likelihoodarray(r,G,components,bins,onstates,dttype,nalleles,nRNA)


"""
function likelihoodarray(r, G, components, bins, onstates, dttype, nalleles, nRNA)
    tcomponents = components.tcomponents
    mcomponents = components.mcomponents
    elementsT = tcomponents.elementsT
    T = make_mat(elementsT, r, tcomponents.nT)
    pss = normalized_nullspace(T)
    elementsTG = tcomponents.elementsTG
    TG = make_mat(elementsTG, r, G)
    pssG = normalized_nullspace(TG)
    hists = Vector[]
    M = make_mat_M(mcomponents, r)
    histF = steady_state(M, mcomponents.nT, nalleles, nRNA)
    push!(hists, histF)
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
            return vcat(f1(v[1:n-1]), f2(v[n]), f1(v[n+1:end]))
        else
            return vcat(f1(v[1:n-1, :]), f2(v[n, :]), f2(v[n+1:end,:]))
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

transform_rates(r, model::AbstractGRSMmodel{Vector{Float64},HMMReporter}) = transform_array(r, model.reporter.weightind, model.fittedparam, logv, logit)


"""
    inverse_transform_rates(x,model::AbstractGmodel)

inverse transform MH parameters on real domain back to rate domain

"""
inverse_transform_rates(p, model::AbstractGmodel) = exp.(p)

inverse_transform_rates(p, model::AbstractGRSMmodel{Vector{Float64},HMMReporter}) = transform_array(p, model.reporter.weightind, model.fittedparam, expv, invlogit)


"""
    get_param(model)

get fitted parameters from model
"""
get_param(model::AbstractGmodel) = log.(model.rates[model.fittedparam])

get_param(model::AbstractGRSMmodel) = transform_rates(model.rates[model.fittedparam], model)


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
    fixed_rates(r, model.fixedeffects)
end

function get_rates(param::Vector{Vector}, model::AbstractGmodel, inverse=true)
    rv = copy_r(model)
    for r in rv
        if inverse
            r[model.fittedparam] = inverse_transform_rates(param, model)
        else
            r[model.fittedparam] = param
        end
    end
    fixed_rates(r, model.fixedeffects)
end



"""
    fixed_rates(r, fixedeffects)

TBW
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

function logprior(param, model::GRSMhierarchicalmodel)


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
