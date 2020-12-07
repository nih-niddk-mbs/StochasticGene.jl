"""
Structure for MH options
"""
struct MHOptions <: Options
    samplesteps::Int64
    annealsteps::Int64
    warmupsteps::Int64
    maxtime::Float64
    temp::Float64
end
"""
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
run_mh(data,model,options)

Run Metropolis-Hastings MCMC algorithm and compute statistics of results
data is a HistogramData structure
and must have a datahistogram function
model is a StochasticGR model with a logprior function
model and data must have a likelihoodfn function
"""
function run_mh(data::HistogramData,model::StochasticGRmodel,options::MHOptions)
    fit,waic = metropolis_hastings(data,model,options)
    stats = compute_stats(fit)
    # write_results(fit,model)
    return fit, stats, waic, (fit.parml,fit.llml)
end
"""
run_mh(data,model,options,nchains)

"""
function run_mh(data::HistogramData,model::StochasticGRmodel,options::MHOptions,nchains)
    sd = run_chains(data,model,options,nchains)
    chain = extract_chain(sd)
    fits,stats,waics = combine_chains(chain)
    waic = combine_waic(waics)
    parml = find_ml(fits)
    return fits,stats,waic,parml
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
    # if options.annealsteps > 0
    #     anneal()
    # end
    if options.warmupsteps > 0
        param,parml,ll,llml,d,proposalcv,predictions = warmup(predictions,param,param,ll,ll,d,model.proposal,data,model,options.warmupsteps,options.temp,time(),options.maxtime)
    else
        parml = copy(param)
        llml = ll
        proposalcv = model.proposal
    end
    fit=sample(predictions,param,parml,ll,llml,d,proposalcv,data,model,options.samplesteps,options.temp,time(),options.maxtime)
    waic = compute_waic(fit.ppd,fit.pwaic,data)
    return fit, waic
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
initial_proposal(model)
return parameters to be fitted and an initial proposal distribution

"""
function initial_proposal(model)
    param = model.rates[model.fittedparam]
    d=proposal(param,model.proposal)
    return param,d
end
"""
proposal(d,cv)
return rand(d) and proposal distribution for cv (vector or covariance)

"""
function proposal(d::MultivariateDistribution,cv::Matrix)
    param = rand(d)
    return param, proposal(param,cv)
end
function proposal(d::Product,cv::Vector)
    param = rand(d)
    return param, proposal(param,cv)
end
"""
proposal(param,cv)
return proposal distribution specified by location and scale

"""
function proposal(param,cv::Vector)
    d = Array{Truncated{Normal{Float64},Continuous}}(undef,0)
    for i in eachindex(param)
        push!(d,truncated(Normal(param[i],cv[i]*param[i]),0.,1000.))
    end
    product_distribution(d)
end

function proposal(param,cv::Matrix)
    c = (2.4)^2 / length(param)
    return MvLogNormal(log.(param),c*cv)
end

"""
sample(predictions,param,rml,ll,llml,d,sigma,data,model,samplesteps,temp,t1,maxtime)

"""
function sample(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    llout = Array{Float64,1}(undef,samplesteps)
    prior = logprior(param,model)
    pwaic = (0,log.(max.(predictions,eps(Float64))),zeros(length(predictions)))
    ppd = zeros(length(predictions))
    step = 0
    accepttotal = 0
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        accept,predictions,param,ll,prior = mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
        accepttotal += accept
        pwaic = update_waic(ppd,pwaic,predictions)
        if ll < llml
            llml,parml = ll,param
        end
        parout[:,step] = param
        llout[step] = ll
    end
    pwaic = step > 1 ? pwaic[3]/(step -1) :  pwaic[3]
    ppd /= step
    Fit(parout[:,1:step],llout[1:step],parml,llml,ppd,pwaic,prior,accepttotal,step)
end
"""
warmup(predictions,param,rml,ll,llml,d,sigma,data,model,samplesteps,temp,t1,maxtime)

"""
function warmup(predictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    prior = logprior(param,model)
    step = 0
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        _,predictions,param,ll,prior = mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
        if ll < llml
            llml,parml = ll,param
        end
        parout[:,step] = param
    end
    covlogparam = cov(log.(parout[:,1:step]'))
    if isposdef(covlogparam)
        proposalcv = covlogparam
    end
    return param,parml,ll,llml,d,proposalcv,predictions
end

"""
mhstep(predictions,param,ll,prior,d,sigma,model,data,temp)
"""
function mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
    paramt,dt = proposal(d,proposalcv)
    priort = logprior(paramt,model)
    llt,predictionst = loglikelihood(paramt,data,model)

    if rand() < exp((ll + prior - llt - priort + mhfactor(param,d,paramt,dt))/temp)
        return 1,predictionst,paramt,llt,priort,dt
    else
        return 0,predictions,param,ll,prior,d
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
combine_chains(chain)
collate results from multiple metropolis_hastings chains

"""
function combine_chains(chain)
    fits = Array{Fit,1}(undef,length(chain))
    stats = Array{Tuple,1}(undef,length(chain))
    waics = Array{Tuple,1}(undef,length(chain))
    for i in eachindex(chain)
        fits[i] = chain[i][1]
        stats[i] = compute_stats(chain[i][1])
        waics[i] = chain[i][2]
    end
    return fits,stats,waics
end

"""
combine_waic(waics)
transform array of array in to a 2D array for
waic and sd of waic
"""
function combine_waic(waics)
    waic = Array{Float64,2}(undef,length(waics),2)
    for i in eachindex(waics)
        waic[i,1] = waics[i][1]
        waic[i,2] = waics[i][2]
    end
    return mean(waic,dims=1),median(waic,dims=1), mad(waic[:,1],normalize=false),waic
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
aic(fit) = 2*length(fit.parml) + 2*fit.llml

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
compute_stats(fit::Fit)

Compute mean, std, median, mad, quantiles and correlations, covariances of parameters
"""
function compute_stats(fit::Fit)
    meanparam = mean(fit.param,dims=2)
    stdparam = std(fit.param,dims=2)
    corparam = cor(fit.param')
    covparam = cov(fit.param')
    covlogparam = cov(log.(fit.param'))
    medparam = median(fit.param,dims=2)
    np = size(fit.param,1)
    madparam = Array{Float64,1}(undef,np)
    qparam = Array{Array,1}(undef,np)
    for i in 1:np
        madparam[i] = mad(fit.param[i,:],normalize=false)
        qparam[i] = quantile(fit.param[i,:],[.025;.5;.975])
    end
    Stats(meanparam,stdparam,medparam,madparam,qparam,corparam,covparam,covlogparam)
end


"""
loglikelihood(param,data,model)

Calls likelihoodfn and datahistogram
provided for each data and model type
"""
function loglikelihood(param,data,model)
    predictions = likelihoodfn(param,data,model)
    hist = datahistogram(data)
    return crossentropy(predictions,hist),predictions
end
"""
crossentropy(predictions::Array{T},data::Array{T}) where {T}

compute crosscentropy between data histogram and likelilhood
"""
function crossentropy(predictions::Array{T1},data::Array{T2}) where {T1,T2}
    lltot = data' * log.(max.(predictions,eps(T1)))
    isfinite(lltot) ? -lltot : Inf
end
"""
hist_entropy(hist)

Compute entropy of an unnormalized histogram
"""
function hist_entropy(hist::Array{Float64,1})
    -hist' * (log.(max.(hist,eps(Float64))) .- log(sum(hist)))
end
function hist_entropy(hist::Array{Array,1})
    l = 0
    for h in hist
        l +=-h' * (log.(max.(h,eps(Float64))) .- log(sum(h)))
    end
    return l
end
"""
write_results(file::String,x)
"""
function write_all(path::String,fit,stats,waic,data,model)
    if ~isdir(path)
        mkpath(path)
    end
    name = "_$(data.name)" * "_$(data.gene)" * "_$(model.G)" * "_$(model.nalleles)" * txtstr
    write_rates(joinpath(path,"rates" * name ),fit,stats,model)
    write_measures(joinpath(path,"measures" * name),fit,waic)
    write_param_stats(joinpath(path,"param_stats" * name),stats)

end
"""
write_rates(file::String,fit)
"""
function write_rates(file::String,fit,stats,model)
    f = open(file,"w")
    writedlm(f,[get_rates(fit.parml,model)],',')
    writedlm(f,[get_rates(stats.meanparam,model)],',')
    writedlm(f,[get_rates(stats.medparam,model)],',')
    close(f)
end
"""
write_measures(file,fit,waic,model)
"""
function write_measures(file::String,fit,waic)
    f = open(file,"w")
    writedlm(f,[fit.llml mean(fit.ll) std(fit.ll) quantile(fit.ll,[.025;.5;.975])' waic[1] waic[2] aic(fit)],',')
    close(f)

end
"""
write_param_stats(stats,waic,data,model)

"""
function write_param_stats(file,stats)
    # file=joinpath(path,"param_stats" * data.gene * txtstr)
    f = open(file,"w")
    # writedlm(f,[model.G model.nalleles])
    writedlm(f,stats.meanparam')
    writedlm(f,stats.stdparam')
    writedlm(f,stats.medparam')
    writedlm(f,stats.madparam')
    writedlm(f,stats.corparam)
    writedlm(f,stats.covparam)
    writedlm(f,stats.covlogparam)
    close(f)
    # f = open("results" * data.gene * txtstr)
    # y=readdlm(f)
    # # x=read(f,Float64)
    # close(f)
    # return y
end

"""
function anneal()

"""

"""
function rhat()
"""
