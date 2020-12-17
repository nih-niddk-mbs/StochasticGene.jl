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
    stats = compute_stats(fit.param)
    return fit, stats, waic
end
"""
run_mh(data,model,options,nchains)

"""
function run_mh(data::HistogramData,model::StochasticGRmodel,options::MHOptions,nchains)
    sd = run_chains(data,model,options,nchains)
    chain = extract_chain(sd)
    waic = pooled_waic(chain)
    fit = merge_fit(chain)
    stats = compute_stats(fit.param)
    return fit, stats, waic
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
    # if options.annealsteps > 0
    #     anneal()
    # end
    if options.warmupsteps > 0
        param,parml,ll,llml,d,proposalcv,predictions = warmup(predictions,param,param,ll,ll,d,model.proposal,data,model,options.warmupsteps,options.temp,time(),options.maxtime*options.warmupsteps/options.samplesteps)
    else
        parml = param
        llml = ll
        proposalcv = model.proposal
    end
    fit=sample(predictions,param,parml,ll,llml,d,proposalcv,data,model,options.samplesteps,options.temp,time(),options.maxtime)
    waic = compute_waic(fit.ppd,fit.pwaic,data)
    return fit, waic
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
    llout = Array{Float64,1}(undef,samplesteps)
    pwaic = (0,log.(max.(predictions,eps(Float64))),zeros(length(predictions)))
    ppd = zeros(length(predictions))
    accepttotal = 0
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    prior = logprior(param,model)
    step = 0
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        accept,predictions,param,ll,prior = mhstep(predictions,param,ll,prior,d,proposalcv,model,data,temp)
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
function compute_stats(param::Array{Float64,2})
    meanparam = mean(param,dims=2)
    stdparam = std(param,dims=2)
    corparam = cor(param')
    covparam = cov(param')
    covlogparam = cov(log.(param'))
    medparam = median(param,dims=2)
    np = size(param,1)
    madparam = Array{Float64,1}(undef,np)
    qparam = Array{Array,1}(undef,np)
    for i in 1:np
        madparam[i] = mad(param[i,:],normalize=false)
        qparam[i] = quantile(param[i,:],[.025;.5;.975])
    end
    Stats(meanparam,stdparam,medparam,madparam,qparam,corparam,covparam,covlogparam)
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
function loglikelihood(param,data,model)
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
hist_entropy(hist)

Compute entropy of an unnormalized histogram
"""
function hist_entropy(hist::Array{Float64,1})
    -hist' * (log.(max.(hist,eps(Float64))) .- log(sum(hist)))
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
deviance(predictions,data) = 2*(crossentropy(predictions,datahistogram(data)) - hist_entropy(data.histRNA))
function deviance(fit::Fit,data,model)
    h=StochasticGene.likelihoodfn(fit.parml,data,model)
    deviance(h,data)
end

"""
write_results(file::String,x)
"""
function write_all(path::String,fit,stats,waic,data,model::StochasticGRmodel)
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

Write rate parameters, rows in order are
maximum likelihood
mean
median
last accepted
"""
function write_rates(file::String,fit::Fit,stats,model)
    f = open(file,"w")
    writedlm(f,[get_rates(fit.parml,model)],',')
    writedlm(f,[get_rates(stats.meanparam,model)],',')
    writedlm(f,[get_rates(stats.medparam,model)],',')
    writedlm(f,[get_rates(fit.param[:,end],model)],',')
    close(f)
end
"""
write_measures(file,fit,waic,model)
"""
function write_measures(file::String,fit::Fit,waic)
    f = open(file,"w")
    writedlm(f,[fit.llml mean(fit.ll) std(fit.ll) quantile(fit.ll,[.025;.5;.975])' waic[1] waic[2] aic(fit)],',')
    writedlm(f,[fit.accept fit.total])
    close(f)
end
"""
write_param_stats(stats,waic,data,model)

"""
function write_param_stats(file,stats::Stats)
    f = open(file,"w")
    writedlm(f,stats.meanparam')
    writedlm(f,stats.stdparam')
    writedlm(f,stats.medparam')
    writedlm(f,stats.madparam')
    writedlm(f,stats.corparam)
    writedlm(f,stats.covparam)
    writedlm(f,stats.covlogparam)
    close(f)
end

function read_covlogparam(file)
    in = readdlm(file)
    in[end-size(in)[2]+1:end,1:size(in)[2]]
end

"""
read_rates(file::String,type::Int)

type
1       maximum likelihood
2       mean
3       median
4       last value of previous run
"""
function read_rates(file::String,type)
    rates = read_rates(file)
    rates[type,:]
end
read_rates(file::String) = readdlm(file,',')

"""
function anneal()

"""

"""
function rhat()
"""
