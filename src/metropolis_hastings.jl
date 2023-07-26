# metropolis_hastings.jl


"""
struct MHOptions <: Options

Structure for MH options
"""
struct MHOptions <: Options
    samplesteps::Int64
    warmupsteps::Int64
    annealsteps::Int64
    maxtime::Float64
    temp::Float64
    tempanneal::Float64
end
"""
struct Fit <: Results

Structure for MH results
"""
struct Fit <: Results
    param::Array
    ll::Array
    parml::Array
    llml::Float64
    ppd::Array
    pwaic::Array
    prior::Float64
    accept::Int
    total::Int
end

"""
struct Stats

Structure for parameter statistics results
"""
struct Stats
    meanparam::Array
    stdparam::Array
    medparam::Array
    madparam::Array
    qparam::Array
    corparam::Array
    covparam::Array
    covlogparam::Array
end

"""
Structure for convergence criteria

"""
struct Measures
    waic::Tuple
    rhat::Vector
end

"""
run_mh(data,model,options)

Run Metropolis-Hastings MCMC algorithm and compute statistics of results
data is an AbstractExperimentalData structure
model is a StochasticGR model with a logprior function
model and data must have a likelihoodfn function
"""
function run_mh(data::AbstractHistogramData,model::AbstractStochasticGRmodel,options::MHOptions)
    fit,waic = metropolis_hastings(data,model,options)
    if options.samplesteps > 0
        stats = compute_stats(fit.param,model)
        rhat = compute_rhat([fit])
        measures = Measures(waic,vec(rhat))
    else
        stats = 0
        measures = 0
    end
    return fit, stats, measures
end
"""
run_mh(data,model,options,nchains)

"""
function run_mh(data::AbstractHistogramData,model::AbstractStochasticGRmodel,options::MHOptions,nchains)
    if nchains == 1
        return run_mh(data,model,options)
    else
        sd = run_chains(data,model,options,nchains)
        chain = extract_chain(sd)
        waic = pooled_waic(chain)
        fit = merge_fit(chain)
        stats = compute_stats(fit.param,model)
        rhat = compute_rhat(chain)
        return fit, stats, Measures(waic,vec(rhat))
    end
end
"""
run_chains(data,model,options,nchains)

"""
function run_chains(data,model,options,nchains)
    sd = Array{Future,1}(undef,nchains)
    for chain in 1:nchains
        sd[chain] = @spawn metropolis_hastings(data,model,options)
    end
    return sd
end

"""
metropolis_hastings(param,data,model,options)

Metropolis-Hastings MCMC algorithm
param = array of parameters to be fit
model, data, and options are structures

stochastic kinetic transcription models
Data can include ON and OFF MS2/PP7 reporter time distributions from live cell recordings,
scRNA or FISH mRNA data, and burst correlations between alleles
"""
function metropolis_hastings(data,model,options)
    param,d = initial_proposal(model)
    ll,predictions = loglikelihood(param,data,model)
    maxtime = options.maxtime
    totalsteps = options.warmupsteps + options.samplesteps + options.annealsteps
    parml = param
    llml = ll
    proposalcv = model.proposal
    if options.annealsteps > 0
        param,parml,ll,llml,predictions,temp = anneal(predictions,param,parml,ll,llml,d,model.proposal,data,model,options.annealsteps,options.temp,options.tempanneal,time(),maxtime*options.annealsteps/totalsteps)
    end
    if options.warmupsteps > 0
        param,parml,ll,llml,d,proposalcv,predictions = warmup(predictions,param,param,ll,ll,d,model.proposal,data,model,options.warmupsteps,options.temp,time(),maxtime*options.warmupsteps/totalsteps)
    end
    fit=sample(predictions,param,parml,ll,llml,d,proposalcv,data,model,options.samplesteps,options.temp,time(),maxtime*options.samplesteps/totalsteps)
    waic = compute_waic(fit.ppd,fit.pwaic,data)
    return fit, waic
end

"""
function anneal(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)

"""
function anneal(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,tempanneal,t1,maxtime)
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    prior = logprior(param,model)
    step = 0
    annealrate = 3/samplesteps
    anneal = 1 - annealrate  # annealing rate with time constant of samplesteps/3
    proposalcv = .02
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        _,predictions,param,ll,prior,d = mhstep(predictions,param,ll,prior,d,proposalcv,model,data,tempanneal)
        if ll < llml
            llml,parml = ll,param
        end
        parout[:,step] = param
        tempanneal = anneal * tempanneal + annealrate * temp
    end
    return param,parml,ll,llml,predictions,temp
end

"""
warmup(predictions,param,rml,ll,llml,d,sigma,data,model,samplesteps,temp,t1,maxtime)

"""
function warmup(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    prior = logprior(param,model)
    step = 0
    accepttotal = 0
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        accept,predictions,param,ll,prior,d = mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
        if ll < llml
            llml,parml = ll,param
        end
        parout[:,step] = param
        accepttotal += accept
    end
    covparam = cov(parout[:,1:step]')*(2.4)^2 / length(param)
    if isposdef(covparam) && step > 1000 && accepttotal/step > .25
        d=proposal_dist(param,covparam,model)
        proposalcv = covparam
        println(proposalcv)
    end
    return param,parml,ll,llml,d,proposalcv,predictions
end

"""
sample(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)

"""
function sample(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)
    llout = Array{Float64,1}(undef,samplesteps)
    pwaic = (0,log.(max.(predictions,eps(Float64))),zeros(length(predictions)))
    ppd = zeros(length(predictions))
    accepttotal = 0
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    prior = logprior(param,model)
    step = 0
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        accept,predictions,param,ll,prior,d = mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
        if ll < llml
            llml,parml = ll,param
        end
        parout[:,step] = param
        llout[step] = ll
        accepttotal += accept
        pwaic = update_waic(ppd,pwaic,predictions)
    end
    pwaic = step > 1 ? pwaic[3]/(step -1) :  pwaic[3]
    ppd /= step
    Fit(parout[:,1:step],llout[1:step],parml,llml,ppd,pwaic,prior,accepttotal,step)
end

"""
mhstep(predictions,param,ll,prior,d,sigma,model,data,temp)

ll is negative log likelihood
"""
function mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
    paramt,dt = proposal(d,proposalcv,model)
    priort = logprior(paramt,model)
    llt,predictionst = loglikelihood(paramt,data,model)
    mhstep(predictions,predictionst,ll,llt,param,paramt,prior,priort,d,dt,temp)
end

function mhstep(predictions,predictionst,ll,llt,param,paramt,prior,priort,d,dt,temp)
    # if rand() < exp((ll + prior - llt - priort + mhfactor(param,d,paramt,dt))/temp)
    if rand() < exp((ll + prior - llt - priort)/temp)
        return 1,predictionst,paramt,llt,priort,dt
    else
        return 0,predictions,param,ll,prior,d
    end
end
"""
computewaic(ppd::Array{T},pwaic::Array{T},data) where {T}

Compute WAIC and SE of WAIC using ppd and pwaic computed in metropolis_hastings()
"""
function compute_waic(ppd::Array{T},pwaic::Array{T},data) where {T}
    hist = datahistogram(data)
    se = sqrt(sum(hist)*(var(ppd,weights(hist),corrected=false)+var(pwaic,weights(hist),corrected=false)))
    lppd = hist'*log.(max.(ppd,eps(T)))
    pwaic = hist'*pwaic
    return -2*(lppd-pwaic), 2*se
end
"""
aic(fit)

AIC
"""
aic(fit::Fit) = 2*length(fit.parml) + 2*fit.llml
aic(nparams::Int,llml) = 2*nparams + 2*llml

"""
initial_proposal(model)
return parameters to be fitted and an initial proposal distribution

"""
function initial_proposal(model)
    param = get_param(model)
    d=proposal_dist(param,model.proposal,model)
    return param,d
end
"""
proposal(d,cv)
return rand(d) and proposal distribution for cv (vector or covariance)

"""
function proposal(d::Distribution,cv,model)
    param = rand(d)
    return param, proposal_dist(param,cv,model)
end

"""
proposal_dist(param,cv)
return proposal distribution specified by location and scale

"""

proposal_scale(cv::Float64,model::AbstractStochasticGRmodel) = sqrt(log(1+cv^2))

proposal_dist(param::Float64,cv::Float64,model) = Normal(param,proposal_scale(cv,model))
function proposal_dist(param::Vector,cv::Float64,model)
    d = Vector{Normal{Float64}}(undef,0)
    for i in eachindex(param)
        push!(d,Normal(param[i],proposal_scale(cv,model)))
    end
    product_distribution(d)
end
function proposal_dist(param::Vector,cv::Vector,model)
    d = Vector{Normal{Float64}}(undef,0)
    for i in eachindex(param)
        push!(d,Normal(param[i],proposal_scale(cv[i],model)))
    end
    product_distribution(d)
end
function proposal_dist(param::Vector,cov::Matrix,model)
    # c = (2.4)^2 / length(param)
    if isposdef(cov)
        return MvNormal(param,cov)
    else
        return MvNormal(param,sqrt.(abs.(diag(cov))))
    end
end


"""
mhfactor(param,d,paramt,dt)
Metropolis-Hastings correction factor for asymmetric proposal distribution
"""
mhfactor(param,d,paramt,dt) = logpdf(dt,param)-logpdf(d,paramt) #pdf(dt,param)/pdf(d,paramt)

"""
update_waic(ppd,pwaic,predictions)

"""
function update_waic(ppd,pwaic,predictions)
    ppd .+= predictions
    var_update(pwaic,log.(max.(predictions,eps(Float64))))
end
# Functions for parallel computing on multiple chains
"""
extract_chain(sdspawn)

"""
function extract_chain(sdspawn)
    chain = Array{Tuple,1}(undef,length(sdspawn))
    for i in eachindex(sdspawn)
        chain[i] = fetch(sdspawn[i])
    end
    return chain
end
"""
collate_fit(chain)
collate fits from multiple metropolis_hastings chains
into and array of Fit structures

"""
function collate_fit(chain)
    fits = Array{Fit,1}(undef,length(chain))
    for i in eachindex(chain)
        fits[i] = chain[i][1]
    end
    return fits
end
"""
collate_waic(chain)
collate waic results from multiple chains
into a 2D array
"""
function collate_waic(chain)
    waic = Array{Float64,2}(undef,length(chain),3)
    for i in eachindex(chain)
        waic[i,1] = chain[i][2][1]
        waic[i,2] = chain[i][2][2]
        waic[i,3] = chain[i][1].total
    end
    return waic
    # mean(waic,dims=1),median(waic,dims=1), mad(waic[:,1],normalize=false),waic
end
function pooled_waic(chain)
    waic = collate_waic(chain)
    return pooled_mean(waic[:,1],waic[:,3]), pooled_std(waic[:,2],waic[:,3])
end

function merge_param(fits::Vector)
    param = fits[1].param
    for i in 2:length(fits)
        param = [param fits[i].param]
    end
    return param
end
function merge_ll(fits::Vector)
    ll = fits[1].ll
    for i in 2:length(fits)
        ll = [ll; fits[i].ll]
    end
    return ll
end

merge_fit(chain::Array{Tuple,1}) = merge_fit(collate_fit(chain))

function merge_fit(fits::Array{Fit,1})
    param = merge_param(fits)
    ll = merge_ll(fits)
    parml,llml = find_ml(fits)
    ppd = fits[1].ppd
    pwaic = fits[1].pwaic
    prior = fits[1].prior
    accept = fits[1].accept
    total = fits[1].total
    for i in 2:length(fits)
        accept += fits[i].accept
        total += fits[i].total
        prior += fits[i].prior
        ppd += fits[i].ppd
        pwaic += fits[i].pwaic
    end
    Fit(param,ll,parml,llml,ppd/length(fits),pwaic/length(fits),prior/length(fits),accept,total)
end

"""
compute_stats(fit::Fit)

Compute mean, std, median, mad, quantiles and correlations, covariances of parameters
"""
function compute_stats(paramin::Array{Float64,2},model)
    param = inverse_transform_rates(paramin,model)
    meanparam = mean(param,dims=2)
    stdparam = std(param,dims=2)
    corparam = cor(param')
    covparam = cov(param')
    # covlogparam = cov(log.(param'))
    covlogparam = cov(paramin')
    medparam = median(param,dims=2)
    np = size(param,1)
    madparam = Array{Float64,1}(undef,np)
    qparam = Matrix{Float64}(undef,3,0)
    for i in 1:np
        madparam[i] = mad(param[i,:],normalize=false)
        qparam = hcat(qparam,quantile(param[i,:],[.025;.5;.975]))
    end
    Stats(meanparam,stdparam,medparam,madparam,qparam,corparam,covparam,covlogparam)
end
"""
compute_rhat(chain::Array{Tuple,1})
compute_rhat(fits::Array{Fit,1})
compute_rhat(params::Vector{Array})


"""
compute_rhat(chain::Array{Tuple,1}) = compute_rhat(collate_fit(chain))

function compute_rhat(fits::Vector{Fit})
    M = length(fits)
    params = Vector{Array}(undef,M)
    for i in 1:M
        params[i] = fits[i].param
    end
    compute_rhat(params)
end


function compute_rhat(params::Vector{Array})
    N = chainlength(params)
    M = length(params)
    m = Matrix{Float64}(undef,size(params[1],1),2*M)
    s = similar(m)
    for i in 1:M
        m[:,2*i-1] = mean(params[i][:,1:N],dims=2)
        m[:,2*i] = mean(params[i][:,N+1:2*N],dims=2)
        s[:,2*i-1] = var(params[i][:,1:N],dims=2)
        s[:,2*i] = var(params[i][:,N+1:2*N],dims=2)
    end
    B = N*var(m,dims=2)
    W = mean(s,dims=2)
    sqrt.((N-1)/N .+ B./W/N)
end


"""
 chainlength(params)

"""
function chainlength(params)
    N = minimum(size.(params,2))
    N = iseven(N) ? N : N-1
    div(N,2)
end

"""
find_ml(fits)
find the maximum likelihood out of multiple chains
"""
function find_ml(fits::Array)
    llml = Array{Float64,1}(undef,length(fits))
    for i in eachindex(fits)
        llml[i] = fits[i].llml
    end
    return fits[argmin(llml)].parml, llml[argmin(llml)]
end

"""
loglikelihood(param,data,model)

Calls likelihoodfn and datahistogram
provided for each data and model type
"""
function loglikelihood(param,data::AbstractHistogramData,model)
    predictions = likelihoodfn(param,data,model)
    hist = datahistogram(data)
    return crossentropy(predictions,hist),predictions
end

"""
crossentropy(predictions::Array{T},data::Array{T}) where {T}

compute crosscentropy between data histogram and likelilhood
"""
function crossentropy(predictions::Array{T1},hist::Array{T2}) where {T1,T2}
        lltot = hist' * log.(max.(predictions,eps(T1)))
        isfinite(lltot) ? -lltot : Inf
end

"""
sumloglikelihood(predictions,data)

sumloglikelihood(predictions,data::RNACountData) = -sum(log.(max.(predictions[data.CountRNA],eps(T1))))
"""

"""
hist_entropy(hist)

Compute entropy of a normalized histogram
"""
function hist_entropy(hist::Array{Float64,1})
    -hist' * (log.(normalize_histogram(max.(hist,eps(Float64)))))
end
function hist_entropy(hist::Array{Array,1})
    l = 0
    for h in hist
        l +=hist_entropy(h)
    end
    return l
end
"""
deviance(predictions,data)
deviance(fit,data,model)

Deviance
"""
# deviance(predictions::Array,data::AbstractHistogramData) = -2*datahistogram(data)'*(log.(max.(predictions,eps(Float64))) - log.(max.(datapdf(data),eps(Float64))))
deviance(predictions::Array,data::AbstractHistogramData) = deviance(predictions,datapdf(data))

# (crossentropy(predictions,datahistogram(data)) - hist_entropy(datahistogram(data)))
function deviance(fit::Fit,data::AbstractHistogramData,model)
    h=likelihoodfn(fit.parml,data,model)
    deviance(h,data)
end

function deviance(data::AbstractHistogramData,model::AbstractStochasticGRmodel)
    h = likelihoodfn(model.rates[model.fittedparam],data,model)
    deviance(h,data)
end

deviance(predictions::Array,hist::Array) = 2*hist'*(log.(max.(hist,eps(Float64)))-log.(max.(predictions,eps(Float64))))
