# This file is part of StochasticGene.jl  

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
    lppd::Array
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
Structure for computed measures

-`waic`: Watanabe-Akaike information criterion
-`rhat`: r-hat convergence diagnostic
"""
struct Measures
    waic::Tuple
    rhat::Vector
end

struct MeasuresEnhanced
    waic::Tuple
    rhat::Vector
    ess::Vector
    geweke::Vector
    mcse::Vector
    trace_analysis::Dict
end

function compute_measures(fits)
    # Existing code
    waic = compute_waic(fits.lppd, fits.pwaic, data)
    rhat = compute_rhat([fits])

    # New diagnostics
    ess = compute_ess(fits.param)
    geweke = compute_geweke(fits.param)
    mcse = compute_mcse(fits.param)
    trace_analysis = analyze_traces(fits.param)

    return Measures(waic, vec(rhat), ess, geweke, mcse, trace_analysis)
end

"""
run_mh(data,model,options)

returns fits, stats, measures

Run Metropolis-Hastings MCMC algorithm and compute statistics of results

-`data`: AbstractExperimentalData structure
-`model`: AbstractGeneTransitionModel structure with a logprior function
-`options`: MHOptions structure

model and data must have a likelihoodfn function
"""
function run_mh(data::AbstractExperimentalData, model::AbstractGeneTransitionModel, options::MHOptions)
    fits, waic = metropolis_hastings(data, model, options)
    if options.samplesteps > 0
        stats = compute_stats(fits.param, model)
        rhat = compute_rhat([fits])
        measures = Measures(waic, vec(rhat))
    else
        stats = 0
        measures = 0
    end
    return fits, stats, measures
end
"""
run_mh(data,model,options,nchains)

"""
function run_mh(data::AbstractExperimentalData, model::AbstractGeneTransitionModel, options::MHOptions, nchains)
    if CUDA.functional()
        println("CUDA is functional")
        return run_mh_gpu(data, model, options, nchains)
    else
        if nchains == 1
            return run_mh(data, model, options)
        else
            sd = run_chains(data, model, options, nchains)
            chain = extract_chain(sd)
            waic = pooled_waic(chain)
            fits = merge_fit(chain)
            stats = compute_stats(fits.param, model)
            rhat = compute_rhat(chain)
            return fits, stats, Measures(waic, vec(rhat))
        end
    end
end

"""
run_mh_gpu(data, model, options, nchains=1)

Run Metropolis-Hastings MCMC algorithm on GPU(s) and compute statistics of results.

# Arguments
- `data`: AbstractExperimentalData structure
- `model`: AbstractGeneTransitionModel structure with a logprior function
- `options`: MHOptions structure
- `nchains`: Number of chains to run (default: 1)

# Returns
- fits, stats, measures as in run_mh
"""
function run_mh_gpu(data::AbstractExperimentalData, model::AbstractGeneTransitionModel, options::MHOptions, nchains=1)
    # Check if CUDA is available
    if !CUDA.functional()
        error("CUDA is not available. Please use run_mh instead.")
    end
    
    # Get number of available GPUs
    devices = collect(CUDA.devices())
    n_gpus = length(devices)
    if n_gpus == 0
        error("No CUDA devices available")
    end
    
    # # Print GPU information
    # println("Running on $(n_gpus) GPU(s):")
    # for i in 1:n_gpus
    #     device = devices[i]
    #     println("  GPU $i: $(device) with $(device.total_memory / 1024^3) GB memory")
    # end
    
    # If only one chain, run on a single GPU
    if nchains == 1
        # Set the active GPU
        CUDA.device!(0)
        println("Running single chain on GPU 0")
        
        # Run the MCMC chain
        # The GPU-accelerated functions in hmm.jl will be used automatically
        fits, waic = metropolis_hastings(data, model, options)
        
        # Compute statistics if needed
        if options.samplesteps > 0
            stats = compute_stats(fits.param, model)
            rhat = compute_rhat([fits])
            measures = Measures(waic, vec(rhat))
        else
            stats = 0
            measures = 0
        end
        
        return fits, stats, measures
    else
        # For multiple chains, run them sequentially on different GPUs
        chain_results = Vector{Tuple}(undef, nchains)
        
        # Run chains sequentially, each on a different GPU
        for chain in 1:nchains
            # Select GPU for this chain
            gpu_id = (chain - 1) % n_gpus
            
            # Set the active GPU for this chain
            CUDA.device!(gpu_id)
            println("Running chain $chain on GPU $gpu_id")
            
            # Run the MCMC chain
            # The GPU-accelerated functions in hmm.jl will be used automatically
            chain_results[chain] = metropolis_hastings(data, model, options)
        end
        
        # Process results
        waic = pooled_waic(chain_results)
        fits = merge_fit(chain_results)
        stats = compute_stats(fits.param, model)
        rhat = compute_rhat(chain_results)
        return fits, stats, Measures(waic, vec(rhat))
    end
end

"""
run_chains(data,model,options,nchains)

returns an array of Futures

Runs multiple chains of MCMC alorithm
"""
function run_chains(data, model, options, nchains)
    sd = Array{Future,1}(undef, nchains)
    for chain in 1:nchains
        sd[chain] = @spawn metropolis_hastings(data, model, options)
    end
    return sd
end

"""
metropolis_hastings(param,data,model,options)

returns fit structure and waic tuple

Metropolis-Hastings MCMC algorithm
param = array of parameters to be fit
model, data, and options are structures

stochastic kinetic transcription models
Data can include ON and OFF MS2/PP7 reporter time distributions from live cell recordings,
scRNA or FISH mRNA data, and burst correlations between alleles
"""
function metropolis_hastings(data, model, options)
    param, d = initial_proposal(model)
    ll, logpredictions = loglikelihood(param, data, model)
    maxtime = options.maxtime
    totalsteps = options.warmupsteps + options.samplesteps + options.annealsteps
    parml = param
    llml = ll
    proposalcv = model.proposal
    if options.annealsteps > 0
        param, parml, ll, llml, logpredictions, temp = anneal(logpredictions, param, parml, ll, llml, d, model.proposal, data, model, options.annealsteps, options.temp, options.tempanneal, time(), maxtime * options.annealsteps / totalsteps)
    end
    if options.warmupsteps > 0
        param, parml, ll, llml, d, proposalcv, logpredictions = warmup(logpredictions, param, param, ll, ll, d, model.proposal, data, model, options.warmupsteps, options.temp, time(), maxtime * options.warmupsteps / totalsteps)
    end
    fits = sample(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, options.samplesteps, options.temp, time(), maxtime * options.samplesteps / totalsteps)
    waic = compute_waic(fits.lppd, fits.pwaic, data)
    return fits, waic
end

"""
function anneal(logpredictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)

returns variables necessary to continue MCMC (param,parml,ll,llml,logpredictions,temp)

runs MCMC with temperature dropping from tempanneal towards temp

"""
function anneal(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, tempanneal, t1, maxtime)
    parout = Array{Float64,2}(undef, length(param), samplesteps)
    prior = logprior(param, model)
    step = 0
    annealrate = 3 / samplesteps
    anneal = 1 - annealrate  # annealing rate with time constant of samplesteps/3
    proposalcv = 0.02
    while step < samplesteps && time() - t1 < maxtime
        step += 1
        _, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, tempanneal)
        if ll < llml
            llml, parml = ll, param
        end
        parout[:, step] = param
        tempanneal = anneal * tempanneal + annealrate * temp
    end
    return param, parml, ll, llml, logpredictions, temp
end

"""
warmup(logpredictions,param,rml,ll,llml,d,sigma,data,model,samplesteps,temp,t1,maxtime)

returns param,parml,ll,llml,d,proposalcv,logpredictions

runs MCMC for options.warmupstep number of samples and does not collect results

"""
function warmup(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime)
    # parout = Array{Float64,2}(undef, length(param), samplesteps)
    prior = logprior(param, model)
    step = 0
    # accepttotal = 0
    while step < samplesteps && time() - t1 < maxtime
        step += 1
        accept, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)
        if ll < llml
            llml, parml = ll, param
        end
        # parout[:, step] = param
        # accepttotal += accept
    end
    # covparam = cov(parout[:,1:step]')*(2.4)^2 / length(param)
    # if isposdef(covparam) && step > 1000 && accepttotal/step > .25
    #     d=proposal_dist(param,covparam,model)
    #     proposalcv = covparam
    #     println(proposalcv)
    # end
    return param, parml, ll, llml, d, proposalcv, logpredictions
end

"""
sample(logpredictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)

returns Fit structure

ll is negative loglikelihood

"""
function sample(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime)
    llout = Array{Float64,1}(undef, samplesteps)
    pwaic = (0, log.(max.(logpredictions, eps(Float64))), zeros(length(logpredictions)))
    lppd = fill(-Inf, length(logpredictions))
    accepttotal = 0
    parout = Array{Float64,2}(undef, length(param), samplesteps)
    prior = logprior(param, model)
    step = 0
    while step < samplesteps && time() - t1 < maxtime
        step += 1
        accept, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)
        if ll < llml
            llml, parml = ll, param
        end
        parout[:, step] = param
        llout[step] = ll
        accepttotal += accept
        lppd, pwaic = update_waic(lppd, pwaic, logpredictions)
    end
    pwaic = step > 1 ? pwaic[3] / (step - 1) : pwaic[3]
    lppd .-= log(step)
    Fit(parout[:, 1:step], llout[1:step], parml, llml, lppd, pwaic, prior, accepttotal, step)
end

function sample_with_thinning(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime, thin_interval=10)
    # Number of samples to keep after thinning
    kept_samples = div(samplesteps, thin_interval)

    llout = Array{Float64,1}(undef, kept_samples)
    parout = Array{Float64,2}(undef, length(param), kept_samples)

    # WAIC components
    pwaic = (0, log.(max.(logpredictions, eps(Float64))), zeros(length(logpredictions)))
    lppd = fill(-Inf, length(logpredictions))

    accepttotal = 0
    prior = logprior(param, model)
    total_step = 0
    saved_step = 0

    while total_step < samplesteps && time() - t1 < maxtime
        total_step += 1

        # MCMC step
        accept, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)

        # Update maximum likelihood
        if ll < llml
            llml, parml = ll, param
        end

        # Save only every thin_interval steps
        if total_step % thin_interval == 0
            saved_step += 1
            parout[:, saved_step] = param
            llout[saved_step] = ll
            lppd, pwaic = update_waic(lppd, pwaic, logpredictions)
        end

        accepttotal += accept
    end

    # Adjust for actual number of samples saved
    if saved_step < kept_samples
        parout = parout[:, 1:saved_step]
        llout = llout[1:saved_step]
    end

    # Adjust WAIC components
    pwaic = saved_step > 1 ? pwaic[3] / (saved_step - 1) : pwaic[3]
    lppd .-= log(saved_step)

    return Fit(parout, llout, parml, llml, lppd, pwaic, prior, accepttotal, total_step)
end

function compute_ess(params)
    n_params = size(params, 1)
    n_samples = size(params, 2)
    ess = zeros(n_params)

    for i in 1:n_params
        # Calculate autocorrelation
        acf = autocor(params[i, :], 0:min(1000, n_samples - 1))
        # Find where autocorrelation drops below 0.05
        cutoff = findfirst(x -> abs(x) < 0.05, acf)
        cutoff = isnothing(cutoff) ? length(acf) : cutoff

        # Sum of autocorrelations with appropriate cutoff
        tau = 1 + 2 * sum(acf[2:cutoff])
        ess[i] = n_samples / tau
    end

    return ess
end

function compute_geweke(params; first_frac=0.1, last_frac=0.5)
    n_params = size(params, 1)
    n_samples = size(params, 2)

    first_end = Int(floor(first_frac * n_samples))
    last_start = Int(floor((1 - last_frac) * n_samples))

    z_scores = zeros(n_params)

    for i in 1:n_params
        first_mean = mean(params[i, 1:first_end])
        last_mean = mean(params[i, last_start:end])

        # Spectral density estimates for variance
        first_var = spectrum0(params[i, 1:first_end])
        last_var = spectrum0(params[i, last_start:end])

        # Z-score
        z_scores[i] = (first_mean - last_mean) /
                      sqrt(first_var / first_end + last_var / (n_samples - last_start + 1))
    end

    return z_scores
end

# Helper for spectral density estimate
function spectrum0(x)
    n = length(x)
    spec = abs.(fft(x .- mean(x))) .^ 2 / n
    return sum(spec[2:Int(floor(n / 2))]) * 2 / n
end

function compute_mcse(params)
    n_params = size(params, 1)
    n_samples = size(params, 2)
    mcse = zeros(n_params)

    for i in 1:n_params
        # Calculate autocorrelation
        acf = autocor(params[i, :], 0:min(1000, n_samples - 1))
        # Find where autocorrelation drops below 0.05
        cutoff = findfirst(x -> abs(x) < 0.05, acf)
        cutoff = isnothing(cutoff) ? length(acf) : cutoff

        # Sum of autocorrelations
        tau = 1 + 2 * sum(acf[2:cutoff])

        # Standard error
        mcse[i] = sqrt(tau * var(params[i, :]) / n_samples)
    end

    return mcse
end

"""
mhstep(logpredictions,param,ll,prior,d,sigma,model,data,temp)

returns 1,logpredictionst,paramt,llt,priort,dt if accept
        1,logpredictionst,paramt,llt,priort,dt if not

ll is negative log likelihood
"""
function mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)
    paramt, dt = proposal(d, proposalcv, model)
    priort = logprior(paramt, model)
    llt, logpredictionst = loglikelihood(paramt, data, model)
    mhstep(logpredictions, logpredictionst, ll, llt, param, paramt, prior, priort, d, dt, temp)
end

function mhstep(logpredictions, logpredictionst, ll, llt, param, paramt, prior, priort, d, dt, temp)
    # mhfactor not needed for symmetric proposal distribution
    # if rand() < exp((ll + prior - llt - priort + mhfactor(param,d,paramt,dt))/temp)
    # if rand() < exp((ll + prior - llt - priort) / temp)
    # println(mhfactor(param,d,paramt,dt))
    # println(ll + prior - llt - priort)
    if log(rand()) < (ll + prior - llt - priort) / temp
        return 1, logpredictionst, paramt, llt, priort, dt
    else
        return 0, logpredictions, param, ll, prior, d
    end
end
"""
update_waic(lppd,pwaic,logpredictions)

returns lppd and pwaic, which are the running sum and variance of -logpredictions
(- sign needed because logpredictions is negative loglikelihood)
"""
function update_waic(lppd, pwaic, logpredictions)
    lppd = logsumexp(lppd, -logpredictions)
    lppd, var_update(pwaic, -logpredictions)
end
"""
computewaic(lppd::Array{T},pwaic::Array{T},data) where {T}

returns  WAIC and SE of WAIC 

using lppd and pwaic computed in metropolis_hastings()
"""
function compute_waic(lppd::Array{T}, pwaic::Array{T}, data) where {T}
    se = sqrt(length(lppd) * var(lppd - pwaic))
    return -2 * sum(lppd - pwaic), 2 * se
end

"""
    compute_waic(lppd::Array{T},pwaic::Array{T},data::AbstractHistogramData) where {T}

TBW
"""
function compute_waic(lppd::Array{T}, pwaic::Array{T}, data::AbstractHistogramData) where {T}
    hist = datahistogram(data)
    se = sqrt(sum(hist) * var(lppd - pwaic, weights(hist), corrected=false))
    return -2 * hist' * (lppd - pwaic), 2 * se
end

function compute_waic(lppd::Array{T}, pwaic::Array{T}, data::RNACountData) where {T}
    se =  var(lppd - pwaic)
    return -2 * mean(lppd - pwaic), 2 * se
end

"""
    compute_waic(lppd::Array{T},pwaic::Array{T},data::AbstractTraceHistogramData) where {T}

histogram data precedes trace data in vectors
"""
function compute_waic(lppd::Array{T}, pwaic::Array{T}, data::AbstractTraceHistogramData) where {T}
    hist = datahistogram(data)
    N = length(hist)
    elppd = lppd - pwaic
    vh = sum(hist) * var(elppd[1:N], weights(hist), corrected=false)
    vt = length(elppd[N+1:end]) * var(elppd[N+1:end])
    se = sqrt(vh + vt)
    return -2 * hist' * (elppd[1:N]) - 2 * sum(elppd[N+1:end]), 2 * se
end

"""
aic(fits::Fit)
aic(nparams::Int,llml)

-`llml`: negative of maximum loglikelihood

returns AIC
"""
aic(fits::Fit) = 2 * length(fits.parml) + 2 * fits.llml
aic(nparams::Int, llml) = 2 * nparams + 2 * llml

"""
initial_proposal(model)

return parameters to be fitted and an initial proposal distribution

"""
function initial_proposal(model)
    param = get_param(model)
    d = proposal_dist(param, model.proposal, model)
    return param, d
end
"""
proposal(d,cv)

return rand(d) and proposal distribution for cv (vector or covariance)

"""
function proposal(d::Distribution, cv, model)
    param = rand(d)
    return param, proposal_dist(param, cv, model)
end

"""
proposal_dist(param,cv)

return proposal distribution specified by location and scale

"""
proposal_dist(param::Float64, cv::Float64, model) = Normal(param, proposal_scale(cv, model))
function proposal_dist(param::Vector, cv::Float64, model)
    d = Vector{Normal{Float64}}(undef, 0)
    for i in eachindex(param)
        push!(d, Normal(param[i], proposal_scale(cv, model)))
    end
    product_distribution(d)
end
function proposal_dist(param::Vector, cv::Vector, model)
    d = Vector{Normal{Float64}}(undef, 0)
    for i in eachindex(param)
        push!(d, Normal(param[i], proposal_scale(cv[i], model)))
    end
    product_distribution(d)
end
function proposal_dist(param::Vector, cov::Matrix, model)
    # c = (2.4)^2 / length(param)
    if isposdef(cov)
        return MvNormal(param, cov)
    else
        return MvNormal(param, sqrt.(abs.(diag(cov))))
    end
end

"""
proposal_scale(cv::Float64,model::AbstractGeneTransitionModel)

return variance of normal distribution of log(parameter)
"""
proposal_scale(cv::Float64, model::AbstractGeneTransitionModel) = sqrt(log(1 + cv^2))


"""
mhfactor(param,d,paramt,dt)
Metropolis-Hastings correction factor for asymmetric proposal distribution
"""
mhfactor(param, d, paramt, dt) = logpdf(dt, param) - logpdf(d, paramt) #pdf(dt,param)/pdf(d,paramt)


# Functions for parallel computing on multiple chains
"""
extract_chain(sdspawn)

returns array of tuples of fetched futures from multiple chains
"""
function extract_chain(sdspawn)
    chain = Array{Tuple,1}(undef, length(sdspawn))
    for i in eachindex(sdspawn)
        chain[i] = fetch(sdspawn[i])
    end
    return chain
end
"""
collate_fit(chain)

returns array of Fit structures from multiple chains
collate fits from multiple metropolis_hastings chains
into and array of Fit structures

"""
function collate_fit(chain)
    fits = Array{Fit,1}(undef, length(chain))
    for i in eachindex(chain)
        fits[i] = chain[i][1]
    end
    return fits
end
"""
collate_waic(chain)

return 2D array of collated waic from multiple chains

"""
function collate_waic(chain)
    waic = Array{Float64,2}(undef, length(chain), 3)
    for i in eachindex(chain)
        waic[i, 1] = chain[i][2][1]
        waic[i, 2] = chain[i][2][2]
        waic[i, 3] = chain[i][1].total
    end
    return waic
    # mean(waic,dims=1),median(waic,dims=1), mad(waic[:,1],normalize=false),waic
end


"""
    pooled_waic(chain)

returns pooled waic and se from multiple chains
"""
function pooled_waic(chain)
    waic = collate_waic(chain)
    return pooled_mean(waic[:, 1], waic[:, 3]), pooled_std(waic[:, 2], waic[:, 3])
end

"""
    merge_param(fits::Vector)

returns pooled params from multiple chains
"""
function merge_param(fits::Vector)
    param = fits[1].param
    for i in 2:length(fits)
        param = [param fits[i].param]
    end
    return param
end
"""
    merge_ll(fits::Vector)

returns pooled negative loglikelihood from multiple chains
"""
function merge_ll(fits::Vector)
    ll = fits[1].ll
    for i in 2:length(fits)
        ll = [ll; fits[i].ll]
    end
    return ll
end

"""
    merge_fit(chain::Array{Tuple,1})
    merge_fit(fits::Array{Fit,1})

returns Fit structure merged from multiple runs
"""
merge_fit(chain::Array{Tuple,1}) = merge_fit(collate_fit(chain))

function merge_fit(fits::Array{Fit,1})
    param = merge_param(fits)
    ll = merge_ll(fits)
    parml, llml = find_ml(fits)
    lppd = fits[1].lppd
    pwaic = fits[1].pwaic
    prior = fits[1].prior
    accept = fits[1].accept
    total = fits[1].total
    for i in 2:length(fits)
        accept += fits[i].accept
        total += fits[i].total
        prior += fits[i].prior
        lppd = logsumexp(lppd, fits[i].lppd)
        pwaic += fits[i].pwaic
    end
    Fit(param, ll, parml, llml, lppd .- log(length(fits)), pwaic / length(fits), prior / length(fits), accept, total)
end

"""
compute_stats(fits::Fit)

returns Stats structure

Compute mean, std, median, mad, quantiles and correlations, covariances of parameters
"""
function compute_stats(paramin::Array{Float64,2}, model)
    param = inverse_transform_rates(paramin, model)
    meanparam = mean(param, dims=2)
    stdparam = std(param, dims=2)
    medparam = median(param, dims=2)
    np = size(param, 1)
    madparam = Array{Float64,1}(undef, np)
    qparam = Matrix{Float64}(undef, 3, 0)
    if typeof(model) <: AbstractGRSMmodel && hastrait(model, :hierarchical)
        nrates = min(num_all_parameters(model), size(param, 1))
        corparam = cor(param[1:nrates, :]')
        covparam = cov(param[1:nrates, :]')
        covlogparam = cov(paramin[1:nrates, :]')
    else
        corparam = cor(param')
        covparam = cov(param')
        covlogparam = cov(paramin')
    end
    for i in 1:np
        madparam[i] = mad(param[i, :], normalize=false)
        qparam = hcat(qparam, quantile(param[i, :], [0.025; 0.5; 0.975]))
    end
    Stats(meanparam, stdparam, medparam, madparam, qparam, corparam, covparam, covlogparam)
end

# function compute_stats(paramin::Array{Float64,2}, model::AbstractGRSMhierarchicalmodel)
#     param = inverse_transform_rates(paramin, model)
#     meanparam = mean(param, dims=2)
#     stdparam = std(param, dims=2)
#     medparam = median(param, dims=2)
#     np = size(param, 1)
#     madparam = Array{Float64,1}(undef, np)
#     qparam = Matrix{Float64}(undef, 3, 0)
#     nrates = num_all_parameters(model)
#     corparam = cor(param[1:nrates, :]')
#     covparam = cov(param[1:nrates, :]')
#     covlogparam = cov(paramin[1:nrates, :]')
#     for i in 1:np
#         madparam[i] = mad(param[i, :], normalize=false)
#         qparam = hcat(qparam, quantile(param[i, :], [0.025; 0.5; 0.975]))
#     end
#     Stats(meanparam, stdparam, medparam, madparam, qparam, corparam, covparam, covlogparam)
# end
"""
compute_rhat(chain::Array{Tuple,1})
compute_rhat(fits::Array{Fit,1})
compute_rhat(params::Vector{Array})

returns r-hat measure
"""
compute_rhat(chain::Array{Tuple,1}) = compute_rhat(collate_fit(chain))

function compute_rhat(fits::Vector{Fit})
    M = length(fits)
    params = Vector{Array}(undef, M)
    for i in 1:M
        params[i] = fits[i].param
    end
    compute_rhat(params)
end


"""
    compute_rhat(params::Vector{Array})


"""
function compute_rhat(params::Vector{Array})
    N = chainlength(params)
    M = length(params)
    m = Matrix{Float64}(undef, size(params[1], 1), 2 * M)
    s = similar(m)
    for i in 1:M
        m[:, 2*i-1] = mean(params[i][:, 1:N], dims=2)
        m[:, 2*i] = mean(params[i][:, N+1:2*N], dims=2)
        s[:, 2*i-1] = var(params[i][:, 1:N], dims=2)
        s[:, 2*i] = var(params[i][:, N+1:2*N], dims=2)
    end
    B = N * var(m, dims=2)
    W = mean(s, dims=2)
    sqrt.((N - 1) / N .+ B ./ W / N)
end

"""
 chainlength(params)

 returns integer half of the minimum number of MCMC runs
"""
function chainlength(params)
    N = minimum(size.(params, 2))
    div(N, 2)
end

"""
find_ml(fits)

returns the negative maximum likelihood out of all the chains
"""
function find_ml(fits::Array)
    llml = Array{Float64,1}(undef, length(fits))
    for i in eachindex(fits)
        llml[i] = fits[i].llml
    end
    return fits[argmin(llml)].parml, llml[argmin(llml)]
end

"""
crossentropy(logpredictions::Array{T},data::Array{T}) where {T}

compute crosscentropy between data histogram and likelilhood
"""
function crossentropy(logpredictions::Array{T1}, hist::Array{T2}) where {T1,T2}
    lltot = hist' * logpredictions
    isfinite(lltot) ? -lltot : Inf
end
"""
hist_entropy(hist)

returns entropy of a normalized histogram
"""
function hist_entropy(hist::Array{Float64,1})
    -hist' * log.(normalize_histogram(max.(hist, 0)))
end
function hist_entropy(hist::Array{Array,1})
    l = 0
    for h in hist
        l += hist_entropy(h)
    end
    return l
end

# function run_mcmc_parallel(model::HMM, observations::AbstractArray, n_chains::Int, n_steps::Int; gpu_ids::Vector{Int}=Int[])
    # if !CUDA.functional()
    #     error("CUDA is not available. Please install CUDA.jl and ensure you have a compatible GPU.")
    # end
    
    # # Get available GPUs
    # available_gpus = CUDA.devices()
    # if isempty(gpu_ids)
    #     gpu_ids = collect(1:length(available_gpus))
    # end
    
    # # Distribute chains across GPUs
    # chains_per_gpu = ceil(Int, n_chains / length(gpu_ids))
    
    # Run chains in parallel
    # results = Vector{Any}(undef, n_chains)
    # Threads.@threads for i in 1:n_chains
    #     # gpu_id = gpu_ids[mod1(i, length(gpu_ids))]
    #     # device!(available_gpus[gpu_id])
    #     results[i] = run_mcmc_chain(model, observations, n_steps)
    # end
    
    # return results
# end


