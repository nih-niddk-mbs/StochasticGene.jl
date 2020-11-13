
struct MHOptions <: Options
    samplesteps::Int64
    annealsteps::Int64
    warmupsteps::Int64
    maxtime::Float64
    temp::Float64
end

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
model, data, and options are struct

stochastic kinetic transcription models
Data can include ON and OFF MS2/PP7 reporter time distributions from live cell recordings
smRNA FISH data, and burst correlations between alleles
"""

"""
run_mhparallel(data,model,options,nchains)
"""
function write_results(file::String,x)
    f = open(file,"w")
    writedlm(f,x)
    close(f)
end

function run_mh(data,model,options,nchains)
    sd = run_chains(data,model,options,nchains)
    chain = extract_chain(sd)
    fits,stats,waics = combine_chains(chain)
    waic = combine_waic(waics)
    parml = find_ml(fits)
    return fits,stats,waic,parml
end

function run_mh(data,model,options)
    fit,waic = metropolis_hastings(data,model,options)
    stats = compute_stats(fit)
    # write_results(fit,model)
    return fit, stats, waic, (fit.parml,fit.llml)
end

function run_chains(data,model,options,nchains)
    sd = Array{Future,1}(undef,nchains)
    for chain in 1:nchains
        sd[chain] = @spawn metropolis_hastings(data,model,options)
    end
    return sd
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
    waic = compute_waic(fit.ppd,fit.pwaic,data)
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
        pwaic = update_waic(ppd,pwaic,predictions)
        if ll < llml
            llml = ll
            parml = copy(param)
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
    #
    # println(ll," ",llt)
    # println(ll + prior - llt - priort)
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
function compute_waic(ppd::Array{T},pwaic::Array{T},data) where {T}
    hist = datahistogram(data)
    se = sqrt(sum(hist)*(var(ppd,weights(hist),corrected=false)+var(pwaic,weights(hist),corrected=false)))
    lppd = hist'*log.(max.(ppd,eps(T)))
    pwaic = hist'*pwaic
    return -2*(lppd-pwaic), 2*se
end

function update_waic(ppd,pwaic,predictions)
    ppd .+= predictions
    var_update(pwaic,log.(max.(predictions,eps(Float64))))
end

function extract_chain(sdspawn)
    chain = Array{Tuple,1}(undef,length(sdspawn))
    for i in eachindex(sdspawn)
        chain[i] = fetch(sdspawn[i])
    end
    return chain
end

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

function combine_waic(waics)
    waic = Array{Float64,2}(undef,length(waics),2)
    for i in eachindex(waics)
        waic[i,1] = waics[i][1]
        waic[i,2] = waics[i][2]
    end
    return mean(waic,dims=1),median(waic,dims=1), mad(waic[:,1],normalize=false),waic
end

function find_ml(fits)
    llml = Array{Float64,1}(undef,length(fits))
    for i in eachindex(fits)
        llml[i] = fits[i].llml
    end
    return fits[argmin(llml)].parml, llml[argmin(llml)]
end



function compute_stats(fit::Fit)
    meanparam = mean(fit.param,dims=2)
    stdparam = std(fit.param,dims=2)
    corparam = cor(fit.param')
    medparam = median(fit.param,dims=2)
    madparam = mad(fit.param[1,:],normalize=false)
    qparam = Array{Array,1}(undef,size(fit.param,1))
    for i in eachindex(fit.param[:,1])
        qparam[i] = quantile(fit.param[i,:],[.025;.5;.975])
    end
    return meanparam,stdparam,medparam,madparam,corparam,qparam
end



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
