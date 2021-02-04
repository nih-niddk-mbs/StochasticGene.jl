"""
Abstract Experimental Data types
"""
abstract type ExperimentalData end
"""
Base for data in the form of samples
"""
abstract type SampleData <: ExperimentalData end
"""
Base for data in the form of a distribution
"""
abstract type HistogramData <: ExperimentalData end

abstract type AbstractRNAData{hType} <: HistogramData end

"""
Data structures

Do not use underscore "_" in name
"""
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    name::String
    gene::String
    nRNA::nType
    histRNA::hType
end
struct TransientRNAData{nType,tType,hType} <: AbstractRNAData{hType}
    name::String
    gene::String
    nRNA::nType
    time::tType
    histRNA::hType
end
struct RNAMixedData{hType} <: AbstractRNAData{hType}
    name::String
    gene::String
    nRNA::Array
    fish::Array{Bool,1}
    histRNA::hType
end
struct LiveCellData <: HistogramData
    name::String
    gene::String
    bins::Array
    OFF::Array
    ON::Array
end
struct RNALiveCellData <: HistogramData
    name::String
    gene::String
    bins::Array
    OFF::Array
    ON::Array
    nRNA::Int
    histRNA::Array
end
struct MultiRNALiveCellData <: HistogramData
    name::String
    gene::String
    nRNA::Array
    histRNA::Tuple
    bins::Tuple
    ON::Tuple
    OFF::Tuple
end

"""
Abstract model types
"""
abstract type Model end
abstract type StochasticGRmodel <: Model end
abstract type AbstractGMmodel <: StochasticGRmodel end
abstract type AbstractGRMmodel <: StochasticGRmodel end


"""
Model structures
"""
struct GMmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMrescaledmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMmultimodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMdelaymodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMlossmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GMmixedmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGMmodel
    G::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GRMmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGRMmodel
    G::Int
    R::Int
    nalleles::Int
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

struct GRSMmodel{RateType,PriorType,ProposalType,ParamType,MethodType} <: AbstractGRMmodel
    G::Int
    R::Int
    nalleles::Int
    type::String
    rates::RateType
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    method::MethodType
end

function write_model(model)


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
datahistogram(data::RNAData) = data.histRNA
function datahistogram(data::AbstractRNAData{Array{Array,1}})
# function datahistogram(data::TransientRNAData)
    v = data.histRNA[1]
    for i in 2:length(data.histRNA)
        v = vcat(v,data.histRNA[i])
    end
    return v
end
datahistogram(data::AbstractRNAData{Array{Float64,1}}) = data.histRNA


function datahistogram(data::RNALiveCellData)
    return [data.OFF;data.ON;data.histRNA]
end



#### Model likelihoods   ####

"""
likelihoodfn(param,data,model)
model likelihoodfn called by loglikelihood
"""
function likelihoodfn(param,data::RNAData,model::GMmodel)
    r = get_rates(param,model)
    n = model.G-1
    steady_state(r[1:2*n+2],n,data.nRNA,model.nalleles)
end
function likelihoodfn(param,data::RNAData,model::GMlossmodel)
    r = get_rates(param,model)
    yieldfactor = r[end]
    n = model.G-1
    steady_state(r[1:2*n+2],yieldfactor,n,data.nRNA,model.nalleles)
end
function likelihoodfn(param,data::AbstractRNAData{Array{Array,1}},model::AbstractGMmodel)
    h = likelihoodarray(param,data,model)
    hconcat = Array{Float64,1}(undef,0)
    for h in h
        hconcat = vcat(hconcat,h)
    end
    return hconcat
end
function likelihoodfn(param::Vector,data::RNALiveCellData,model::GRSMmodel,method)
    pdftuple = likelihoodtuple(param,data,model)
    pdfvector = Array{Float64,1}(undef,0)
    for p in pdftuple
        vcat(pdfvector,p)
    end
    return pdfvector
end
function likelihoodfn(param::Vector,data::RNALiveCellData,model::GRSMmodel)
    modelOFF, modelON, histF = likelihoodtuple(param,data,model)
    return [modelOFF;modelON;histF]
end

"""
likelihoodarray(param,data,model)

Compute time dependent model likelihoods
first set of parameters gives the initial histogram
2nd set gives the new parameters at time 0
data.histRNA holds array of histograms for time points given by data.time
transient computes the time evolution of the histogram
model.method=1 specifies finite difference solution otherwise use eigendecomposition solution,
"""
function likelihoodarray(param,data::TransientRNAData,model::GMlossmodel)
    yieldfactor = get_rates(param,model)[end]
    h = likelihoodarray(param,data::TransientRNAData,model,maximum(data.nRNA))
    technical_loss!(h,yieldfactor)
    trim(h,data.nRNA)
end
function likelihoodarray(param,data::TransientRNAData,model,maxdata)
    r = get_rates(param,model)
    G = model.G
    h0 = initial_distribution(param,r,G,model,maxdata)
    transient(r,G,data.time,model,h0)
end
function likelihoodarray(param,data::TransientRNAData,model::AbstractGMmodel)
    h=likelihoodarray(param,data,model,maximum(data.nRNA))
    # r = get_rates(param,model)
    # G = model.G
    # h0 = initial_distribution(param,r,G,model,maximum(data.nRNA))
    # h = transient(r,G,data.time,model,h0)
    trim(h,data.nRNA)
end
function likelihoodarray(param,data::RNAData,model::GMmultimodel)
    r = get_rates(param,model)
    G = model.G
    h = Array{Array{Float64,1},1}(undef,length(data.nRNA))
    for i in eachindex(data.nRNA)
        g = steady_state(r[1:2*G],G-1,data.nRNA[i],model.nalleles)
        h[i] = threshold_noise(g,r[2*G+1],r[2*G+1+i],data.nRNA[i])
    end
    return h
end
function likelihoodarray(param,data::RNAMixedData,model::AbstractGMmodel)
    r = get_rates(param,model)
    G = model.G
    h = Array{Array{Float64,1},1}(undef,length(data.nRNA))
    j = 1
    for i in eachindex(data.fish)
        g = steady_state(r[1:2*G],G-1,maximum(data.nRNA),model.nalleles)
        if data.fish[i]
            h[i] = threshold_noise(g,r[2*G+j],r[2*G+j+1],data.nRNA[i])
            j += 2
        else
            h[i] = technical_loss(g,r[2*G+j],data.nRNA[i])
            j += 1
        end
    end
    return h
end

function likelihoodtuple(param,data::RNALiveCellData,model::GRSMmodel)
    r = get_rates(param,model)
    if model.method < 2
        modelOFF, modelON = offonPDF(data.bins,r,model.G-1,model.R,model.method)
        if model.type == "off"
            histF = steady_state_offpath(r,model.G-1,model.R,data.nRNA,model.nalleles)
        else
            histF = steady_state(r,model.G-1,model.R,data.nRNA,model.nalleles)
        end
    else
        if model.type == "off"
            modelOFF,modelON,histF = telegraphoff(data.bins,data.nRNA,r,model.G-1,model.R,model.nalleles)
        else
            modelOFF,modelON,histF = telegraph(data.bins,data.nRNA,r,model.G-1,model.R,model.nalleles)
        end
    end
    return modelOFF, modelON, histF
end

function likelihoodtuple(r,data::RNALiveCellData,model,method)
    telegraph(data.bins,model.nRNA,n,nr,r,model.nalleles)
end

function likelihoodtuples(bins,nRNA,r,n,nr,nalleles,type)
    if type == "off"
        modelOFFt, modelONt, histFt = telegraphoff(bins,nRNA,r,n,nr,nalleles)
        histF = steady_state_offpath(r,n,nr,nRNA,nalleles)
    else
        modelOFFt, modelONt, histFt = telegraph(bins,nRNA,r,n,nr,nalleles)
        histF = steady_state(r,n,nr,nRNA,nalleles)
    end
    modelOFF, modelON = offonPDF(bins,r,n,nr)
    return modelOFF,modelON,histF,modelOFFt,modelONt,histFt
end


# Functions for transient chemical master solutions

function transient(r::Vector,G::Int,times::Vector,model::GMmodel,h0::Vector)
    transient(times,r[2*G+1:4*G],G-1,model.nalleles,h0,model.method)
end
function transient(r::Vector,G::Int,times::Vector,model::GMdelaymodel,h0::Vector)
    transient_delay(times,r[1:2*G],r[2*G+1:4*G],r[end],G-1,model.nalleles,h0)
end

function initial_distribution(param,r,G::Int,model::AbstractGMmodel,nRNAmax)
    steady_state_full(r[1:2*G],G-1,nRNAmax)
end


"""
get_rates(param,model)
replace fitted rates with new values and return
"""
function get_rates(param,model::AbstractGMmodel)
    r = copy(model.rates)
    r[model.fittedparam] = param
    return r
end


function get_rates(params,model::GRSMmodel)
    r = copy(model.rates)
    r[model.fittedparam] = params
    return r
end

"""
get_rates(param,model::GMrescaledmodel)

gammas are scaled by nu
"""
function get_rates(param,model::GMrescaledmodel)
    r = copy(model.rates)
    n = 2*model.G - 1
    nu = n in model.fittedparam ? param[findfirst(model.fittedparam .== n)] : r[n]
    r[1:n-1] /= r[n]
    r[model.fittedparam] = param
    r[1:n-1] *= nu
    if r[2*model.G + 3] > 1
        r[2*model.G + 3] = 1
    end
    return r
end

"""
get_param(model)

get fitted parameters from model
"""
get_param(model::StochasticGRmodel) = model.rates[model.fittedparam]

function get_param(model::GMrescaledmodel)
    r = copy(model.rates)
    n = 2*model.G - 1
    r[1:n-1] /= r[n]
    r[model.fittedparam]
end


"""
logprior(param,model::AbstractGMmodel)

compute log of the prior
"""
function logprior(param,model::AbstractGMmodel)
    d = model.rateprior
    p=0
    for i in eachindex(d)
        p -= logpdf(d[i],param[i])
    end
    return p
end

"""
logprior(x,model::GRSModel)
Compute log prior using distributions in Model.rateprior
called by mhstep() in metropolis_hastings.jl
"""
function logprior(param,model::GRSMmodel)
    r = get_rates(param,model)
    d = model.rateprior
    G = model.G
    R = model.R
    p=0
    j = 1
    #  G rates
    for i in Grange(G)
        p -= logpdf(d[j],r[i])
        j += 1
    end
    # initiation rate
    i = initiation(G)
    p -= logpdf(d[j],r[i])
    j += 1
    # sum of R Steps rates are bounded by length of insertion site to end of gene, i.e sum of reciprocal rates is bounded
    t = sum(1 ./ r[Rrange(G,R)])
    p -= logpdf(d[j],1/t)
    j += 1
    # priors for splice rates
    rs = 0
    for i in Srange(G,R)
        rs += r[i]
        p -= logpdf(d[j],rs)
        j += 1
    end
    return p
end


"""
plot_histogram()

functions to plot data and model predicted histograms

"""
function plot_histogram(data::AbstractRNAData{Array{Array,1}},model)
    h=likelihoodarray(get_param(model),data,model)
    figure(data.gene)
    for i in eachindex(h)
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
    end
    return h
end

function plot_histogram(data::RNALiveCellData,model)
    h=likelihoodtuple(get_param(model),data,model)
    figure(data.gene)
    plot(h[1])
    plot(normalize_histogram(data.OFF))
    plot(h[2])
    plot(normalize_histogram(data.ON))
    figure("FISH")
    plot(h[3])
    plot(normalize_histogram(data.histRNA))
    return h
end

function plot_histogram(data::AbstractRNAData{Array{Float64,1}},model)
    h=likelihoodfn(get_param(model),data,model)
    figure(data.gene)
    plot(h)
    plot(normalize_histogram(data.histRNA))
    return h
end

"""
teststeadystatemodel(model,nhist)

Compare chemical master solution to Gillespie simulation for steadystate mRNA distribution

"""
function teststeadystatemodel(model::AbstractGMmodel,nhist)
    G = model.G
    r = model.rates
    g1 = steady_state(r[1:2*G],G-1,nhist,model.nalleles)
    g2 = telegraph(G-1,r[1:2*G],10000000,1e-5,nhist,model.nalleles)
    return g1,g2
end


"""
residenceprob_G(model)
residenceprob_G(r)

Residence probability of G states
given by exact steady state solution
of master equation

"""
function residenceprob_G(file::String,G)
    r = readdlm(file,',')
    m = size(r)[1]
    p = Array{Any,2}(undef,m,G+1)
    n = G-1
    for i in 1:m
        p[i,1] = r[i,1]
        p[i,2:end] = residenceprob_G(r[i,2:2*n+1],n)
    end
    return p
end

residenceprob_G(model::AbstractGMmodel)=resprob_G(model.rates,model.G-1)
function residenceprob_G(r::Vector,n::Int)
    Gss = Array{Float64,2}(undef,1,n+1)
    Gss[1,1] = 1.
    for k in 1:n
        Gss[1,k+1] = Gss[1,k]*r[2*k-1]/r[2*k]
    end
    Gss ./= sum(Gss)
end
"""
splicesiteusage()

splice site usage probability is the probabilty of ejection times
the probability of not ejection prior to that point
"""
splicesiteusage(model::GRSMmodel) = splicesiteusage(model.rates,model.G-1,model.R)
function splicesiteusage(r::Vector,n::Int,nr::Int)
    nu = get_nu(r,n,nr)
    eta = get_eta(r,n,nr)
	ssf = zeros(nr)
	survival = 1
	for j in 1:nr
		ssf[j] = eta[j]/(nu[j+1]+eta[j])*survival
		survival *= nu[j+1]/(nu[j+1]+eta[j])
	end
	return ssf
end

"""
burstsize(n,nr,r)
Burst size distribution of GRS  model
for total pre-mRNA occupancy and
unspliced (visible) occupancy
obtained by marginalizing over conditional steady state distribution
"""
burstsize(model::GRSMmodel) = burstsize(model.G-1,model.R,model.rates)
function burstsize(n::Int,nr::Int,r::Vector)
    T =  mat_GSR_T(r,n,nr)
    pss = normalized_nullspace(T)
	Rss = zeros(nr)
	Rssvisible = zeros(nr)
	ssf = zeros(nr)
	asum = 0
	for w=1:nr
		for i = 1:n+1, z = 1:3^nr
			zdigits = digits(z-1,base=3,pad=nr)
			a = i + (n+1)*(z-1)
			if sum(zdigits .== 2) == w
				Rssvisible[w] += pss[a]
			end
			if sum(zdigits .> 0) == w
				Rss[w] += pss[a]
			end
		end
	end
	Rss ./= sum(Rss), Rssvisible ./= sum(Rssvisible)
end

# Functions for saving and loading data and models

"""
write_log(file,datafile,data,model)
write all information necessary for rerunning
"""
function save_data(file::String,data::TransientRNAData)
    f = open(file,"w")
    writedlm(f, [typeof(data)])
    writedlm(f,[data.name])
    writedlm(f,[data.gene])
    writedlm(f,[data.nRNA])
    writedlm(f,[data.time])
    writedlm(f,[data.histRNA])
    close(f)
end


function load_data(file::String,model::AbstractGMmodel)


end

function save_model(file::String,model::AbstractGMmodel)
    f = open(file,"w")
    write(f, model.G)
    writedlm(f,model.nalleles)
    writedlm(f,model.ratepriors)
    writedlm(f,model.proposal)
    writedlm(f,model.fittedparam)
    writedlm(f,model.method)
    close(f)

end

function load_model(file::String,model::AbstractGRMmodel)


end
