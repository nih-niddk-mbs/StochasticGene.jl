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
metropolis_hastings(param,data,model,options)

Metropolis-Hastings MCMC algorithm
r = parameter array
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
write_results(file::String,x)
"""
function write_results(file::String,x)
    f = open(file,"w")
    writedlm(f,x)
    close(f)
end
"""
write_results(stats,waic,data,model)

"""
function write_results(stats,waic,data,model)
    f = open("results" * data.gene * txtstr,"w")
        writedlm(f,model.nalleles)
        writedlm(f,stats[1]')
        writedlm(f,stats[2]')
        writedlm(f,stats[3]')
        writedlm(f,stats[4])
        writedlm(f,stats[5])
        writedlm(f,stats[6])
        writedlm(f,[waic[1] waic[2]])
    close(f)
    f = open("results" * data.gene * txtstr)
        y=readdlm(f)
        # x=read(f,Float64)
    close(f)
    return y
end

"""
run_mh(data,model,options)

"""
function run_mh(data,model,options)
    fit,waic = metropolis_hastings(data,model,options)
    stats = compute_stats(fit)
    # write_results(fit,model)
    return fit, stats, waic, (fit.parml,fit.llml)
end
"""
run_mh(data,model,options,nchains)

"""
function run_mh(data,model,options,nchains)
    sd = run_chains(data,model,options,nchains)
    chain = extract_chain(sd)
    fits,stats,waics = combine_chains(chain)
    waic = combine_waic(waics)
    parml = find_ml(fits)
    return fits,stats,waic,parml
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
            llml = ll
            parml = copy(param)
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
            llml = ll
            parml = copy(param)
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
computewaic(ppd::Array{T},pwaic::Array{T},data) where {T}
"""
function compute_waic(ppd::Array{T},pwaic::Array{T},data) where {T}
    hist = datahistogram(data)
    se = sqrt(sum(hist)*(var(ppd,weights(hist),corrected=false)+var(pwaic,weights(hist),corrected=false)))
    lppd = hist'*log.(max.(ppd,eps(T)))
    pwaic = hist'*pwaic
    return -2*(lppd-pwaic), 2*se
end
"""
update_waic(ppd,pwaic,predictions)

"""
function update_waic(ppd,pwaic,predictions)
    ppd .+= predictions
    var_update(pwaic,log.(max.(predictions,eps(Float64))))
end
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

"""
function combine_chains(chain)
    fits = Array{Fit,1}(undef,length(chain))
    stats = Array{Tuple,1}(undef,length(chain))
    waics = Array{Tuple,1}(undef,length(chain))
    for i in eachindex(chain)
        stats[i] = compute_stats(chain[i][1])
        waics[i] = chain[i][2]
        fits[i] = chain[i][1]
    end
    return fits,stats, waics
end

"""
combine_waic(waics)
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
find_ml(fits)

"""
function find_ml(fits)
    llml = Array{Float64,1}(undef,length(fits))
    for i in eachindex(fits)
        llml[i] = fits[i].llml
    end
    return fits[argmin(llml)].parml, llml[argmin(llml)]
end
"""
compute_stats(fit::Fit)

"""
function compute_stats(fit::Fit)
    meanparam = mean(fit.param,dims=2)
    stdparam = std(fit.param,dims=2)
    corparam = cor(fit.param')
    covparam = cov(fit.param')
    covlogparam = cov(log.(fit.param'))
    medparam = median(fit.param,dims=2)
    madparam = mad(fit.param[1,:],normalize=false)
    qparam = Array{Array,1}(undef,size(fit.param,1))
    for i in eachindex(fit.param[:,1])
        qparam[i] = quantile(fit.param[i,:],[.025;.5;.975])
    end
    return meanparam,stdparam,medparam,madparam,qparam,corparam,covparam,covlogparam
end

"""

"""
function rhat()


end

"""
loglikelihood(param,data,model)

"""
function loglikelihood(param,data,model)
    predictions = likelihoodfn(param,data,model)
    hist = datahistogram(data)
    return crossentropy(predictions,hist),predictions
end
"""
crossentropy(predictions::Array{T},data::Array{T}) where {T}
"""
function crossentropy(predictions::Array{T1},data::Array{T2}) where {T1,T2}
    lltot = data' * log.(max.(predictions,eps(T1)))
    isfinite(lltot) ? -lltot : Inf
end

"""
function anneal()

end
"""
