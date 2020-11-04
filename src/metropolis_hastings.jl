
struct MHOptions <: Options
    samplesteps::Int64
    annealsteps::Int64
    warmupsteps::Int64
    maxtime::Float64
    temp::Float64
end

struct Fit <: Results
    params::Array
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
model, data, and options are struct

stochastic kinetic transcription models
Data can include ON and OFF MS2/PP7 reporter time distributions from live cell recordings
smRNA FISH data, and burst correlations between alleles
"""

"""
run_mhparallel(data,model,options,nchains)
"""

function run_mhparallel(data,model,options,nchains)
    r = model.rates
    sdspawn = Array{Future,1}(undef,nchains)
    for chain in 1:nchains
        sdspawn[chain] = @spawn metropolis_hastings(param,data,model,options)
    end
    fit=combinechains(sdspawn)
    write_results(fit)
end

function run_mh(data,model,options)
    fit=metropolis_hastings(data,model,options)
    write_results(fit,model)
    # compute_stats(fit)
end

function metropolis_hastings(data,model,options)
    param,d = initial_proposal(model)
    ll,predictions = loglikelihood(param,data,model)
    # if options.annealsteps > 0
    #     anneal()
    # end
    # if options.warmupsteps > 0
    #     warmup()
    # end
    fit=sample(predictions,param,param,ll,ll,d,data,model,options.samplesteps,options.temp,time(),options.maxtime)
    waic = computewaic(fit.ppd,fit.pwaic,data)
    return fit, waic
end

function initial_proposal(model)
    param = model.rates[model.fittedparam]
    d=proposal(param,model.proposal)
    return param,product_distribution(d)
end
function proposal(d::Product,cv::Vector)
    param = rand(d)
    dt = proposal(param,cv)
    return param, product_distribution(dt)
end
function proposal(param,cv)
    d = Array{Truncated{Normal{Float64},Continuous}}(undef,0)
    for i in eachindex(param)
        push!(d,truncated(Normal(param[i],cv[i]*param[i]),0.,1000.))
    end
    return d
end

"""
sample(predictions,param,rml,ll,llml,d,sigma,data,model,samplesteps,temp,t1,maxtime)

"""
function sample(predictions,param,parml,ll,llml,d,data,model,samplesteps,temp,t1,maxtime)
    parout = Array{Float64,2}(undef,length(param),samplesteps)
    llout = Array{Float64,1}(undef,samplesteps)
    prior = logprior(param,model)
    priorml = prior
    pwaic = (0,log.(max.(predictions,eps(Float64))),zeros(length(predictions)))
    ppd = zeros(length(predictions))
    step = 0
    accepttotal = 0
    while  step < samplesteps && time() - t1 < maxtime
        step += 1
        accept,predictions,param,ll,prior = mhstep(predictions,param,ll,prior,d,model,data,temp)
        accepttotal += accept
        pwaic = waicupdate(ppd,pwaic,predictions)
        if ll < llml
            llml = ll
            parml = param
            priorml = prior
        end
        parout[:,step] = param
        llout[step] = ll
    end
    pwaic = step > 1 ? pwaic[3]/(step -1) :  pwaic[3]
    ppd /= step
    Fit(parout[:,1:step],llout[1:step],parml,llml,ppd,pwaic,prior,accepttotal,step)
    # return (r=parout[:,1:step],ll=llout[1:step],parml=parml,llml=llml,ppd=ppd,pwaic=pwaic,prior=prior,accept=accepttotal,total=step)
end
"""
mhstep(predictions,param,ll,prior,d,sigma,model,data,temp)
"""
function mhstep(predictions,param,ll,prior,d,model,data,temp)
    paramt,dt = proposal(d,model.proposal)
    priort = logprior(paramt,model)
    llt,predictionst = loglikelihood(paramt,data,model)
    # println(exp((ll + prior - llt - priort)/temp))# * mhfactor(param,d,paramt,dt))
    if rand() < exp((ll + prior - llt - priort)/temp) * mhfactor(param,d,paramt,dt)
        return 1,predictionst,paramt,llt,priort,dt
    else
        return 0,predictions,param,ll,prior,d
    end
end

function mhfactor(param,d,paramt,dt)
    exp(logpdf(dt,param)-logpdf(d,paramt)) #pdf(dt,r)/pdf(d,rt)
end

"""
computewaic(ppd::Array{T},pwaic::Array{T},data) where {T}
"""
function computewaic(ppd::Array{T},pwaic::Array{T},data) where {T}
    hist = datahistogram(data)
    se = sqrt(sum(hist)*(var(ppd,weights(hist),corrected=false)+var(pwaic,weights(hist),corrected=false)))
    lppd = hist'*log.(max.(ppd,eps(T)))
    pwaic = hist'*pwaic
    return -2*(lppd-pwaic), 2*se
end

function waicupdate(ppd,pwaic,predictions)
    ppd .+= predictions
    var_update(pwaic,log.(max.(predictions,eps(Float64))))
end

function combinechains(sdspawn)

    # return rout,llout,rml,llml,ppd,pwaic,accepttotal,step
    return fit,waic
end

function compute_stats(fit::Fit)
    cor(fit.params)
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

function warmup(param,rml,ll,llml,d,cv2,model,data,temp,t1,maxtime,totalsteps,predictions)

    post(param,rml,ll,llml,d,cv2,model,data,temp,predictions)

    waic = WAIC(predictions,log.(predictions),0.)
    sigma = log.(1+model.cv)
    step = 0
    while time() - t1 < maxtime && step < totalsteps
        step += 1
        accept += mhstep!(param,rml,ll,llml,d,sigma,model,data,temp,predictions)
        waicupdate!(waic,predictions)
    end
    return waic

end

function post(param,rml,ll,llml,d,sigma,predictions,model,data,temp)
    paramt,dt = proposal(d,sigma)
    priort = logprior(paramt,model)
    llt,predictionst = loglikelihood(paramt,model,data)
    MHfactor = pdf(dt,param)/pdf(d,paramt)
    exp((ll - llt)/temp - priort + prior) * MHfactor
end

"""
