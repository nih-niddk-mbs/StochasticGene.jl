# This file is part of StochasticGene.jl  

# metropolis_hastings.jl


"""
struct MHOptions <: Options

Options for Metropolis-Hastings MCMC.

# Fields
- `samplesteps::Int64`: Number of MCMC samples to collect.
- `warmupsteps::Int64`: Number of warmup (burn-in) steps.
- `annealsteps::Int64`: Number of annealing steps.
- `maxtime::Float64`: Maximum allowed runtime (seconds).
- `temp::Float64`: Final temperature for annealing.
- `tempanneal::Float64`: Initial temperature for annealing.
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

Results of a Metropolis-Hastings MCMC run.

# Fields
- `param::Array`: Matrix of sampled parameter values (parameters × samples).
- `ll::Array`: Vector of negative log-likelihoods for each sample.
- `parml::Array`: Parameter vector with maximum likelihood.
- `llml::Float64`: Maximum (minimum value) negative log-likelihood.
- `lppd::Array`: Log pointwise predictive density for WAIC.
- `pwaic::Array`: Pointwise WAIC penalty terms.
- `prior::Float64`: Sum of log prior probabilities for all samples.
- `accept::Int`: Number of accepted proposals.
- `total::Int`: Total number of proposals.
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

Summary statistics for MCMC parameters.

# Fields
- `meanparam::Array`: Mean of each parameter.
- `stdparam::Array`: Standard deviation of each parameter.
- `medparam::Array`: Median of each parameter.
- `madparam::Array`: Median absolute deviation of each parameter.
- `qparam::Array`: Quantiles (2.5%, 50%, 97.5%) for each parameter.
- `corparam::Array`: Correlation matrix of parameters.
- `covparam::Array`: Covariance matrix of parameters.
- `covlogparam::Array`: Covariance matrix of log-parameters.
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
# struct Measures
#     waic::Tuple
#     rhat::Vector
# end

"""
    struct Measures

Computed diagnostic measures for MCMC results.

# Fields
- `waic::Tuple`: Watanabe-Akaike information criterion and its standard error.
- `rhat::Vector`: R-hat convergence diagnostic for each parameter.
- `ess::Vector`: Effective sample size for each parameter.
- `geweke::Vector`: Geweke z-scores for each parameter.
- `mcse::Vector`: Monte Carlo standard error for each parameter.
"""
struct Measures
    waic::Tuple
    rhat::Vector
    ess::Vector
    geweke::Vector
    mcse::Vector
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
        rhat = vec(compute_rhat([fits]))  # ensure rhat is a vector
        ess, geweke, mcse = compute_measures(fits)
        measures = Measures(waic, rhat, ess, geweke, mcse)
    else
        stats = 0
        measures = 0
    end
    return fits, stats, measures
end
"""
run_mh(data,model,options,nchains)

Run Metropolis-Hastings MCMC algorithm with multiple chains
"""
function run_mh(data::AbstractExperimentalData, model::AbstractGeneTransitionModel, options::MHOptions, nchains)
    if false && CUDA.functional()
        println("CUDA is functional")
        return run_mh_gpu(data, model, options, nchains)
    else
        if nchains == 1
            return run_mh(data, model, options)
        else
            Distributed.ENV["JULIA_WORKER_TIMEOUT"] = "60"  # 60 second timeout
            sd = run_chains(data, model, options, nchains)
            chain = extract_chain(sd)
            waic = pooled_waic(chain)
            fits = merge_fit(chain)
            stats = compute_stats(fits.param, model)
            rhat = vec(compute_rhat(chain))
            ess, geweke, mcse = compute_measures(chain)
            return fits, stats, Measures(waic, rhat, ess, geweke, mcse)
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
            ess, geweke, mcse = compute_measures(fits)
            measures = Measures(waic, rhat, ess, geweke, mcse)
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
        ess, geweke, mcse = compute_measures(chain_results)
        return fits, stats, Measures(waic, rhat, ess, geweke, mcse)
    end
end

"""
run_chains(data,model,options,nchains)

returns an array of Futures

Runs multiple chains of MCMC algorithm with retry logic for worker initialization
"""
function run_chains(data, model, options, nchains)
    sd = Array{Future,1}(undef, nchains)
    max_retries = 3
    retry_delay = 5  # seconds
    
    for chain in 1:nchains
        retry_count = 0
        while retry_count < max_retries
            try
                sd[chain] = @spawn metropolis_hastings(data, model, options)
                break  # Success, exit retry loop
            catch e
                retry_count += 1
                if retry_count == max_retries
                    @error "Failed to initialize worker after $max_retries attempts" exception=(e, catch_backtrace())
                    rethrow(e)
                end
                @warn "Worker initialization failed, attempt $retry_count of $max_retries. Retrying in $retry_delay seconds..."
                sleep(retry_delay)
            end
        end
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
        println("Annealing")
        param, parml, ll, llml, logpredictions, temp = anneal(logpredictions, param, parml, ll, llml, d, model.proposal, data, model, options.annealsteps, options.temp, options.tempanneal, time(), maxtime * options.annealsteps / totalsteps)
    end
    if options.warmupsteps > 0
        println("Warmup")
        param, parml, ll, llml, d, proposalcv, logpredictions = warmup(logpredictions, param, param, ll, ll, d, model.proposal, data, model, options.warmupsteps, options.temp, time(), maxtime * options.warmupsteps / totalsteps)
    end
    println("Sampling")
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
    proposalcv = model.proposal
    while step < samplesteps && time() - t1 < maxtime
        step += 1
        _, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, tempanneal)
        if ll > llml
            llml, parml = ll, param
        end
        parout[:, step] = param
        tempanneal = anneal * tempanneal + annealrate * temp
    end
    return param, parml, ll, llml, logpredictions, temp
end

"""
warmup(logpredictions,param,rml,ll,llml,d,sigma,data,model,samplesteps,temp,t1,maxtime)

Run a warmup phase for the Metropolis-Hastings MCMC algorithm to adapt proposal distributions.

# Arguments
- `logpredictions`: Initial log-likelihood predictions.
- `param`: Initial parameter vector.
- `parml`: Initial maximum likelihood parameter vector.
- `ll`: Initial negative log-likelihood.
- `llml`: Initial maximum likelihood value.
- `d`: Initial proposal distribution.
- `proposalcv`: Initial proposal covariance or scale.
- `data`: Experimental data structure.
- `model`: Model structure.
- `samplesteps`: Number of warmup steps to run.
- `temp`: Temperature for MCMC.
- `t1`: Start time.
- `maxtime`: Maximum allowed time for warmup.

# Returns
- `param`: Final parameter vector after warmup.
- `parml`: Maximum likelihood parameter vector found during warmup.
- `ll`: Final negative log-likelihood.
- `llml`: Minimum negative log-likelihood found during warmup.
- `d`: Final proposal distribution.
- `proposalcv`: Final proposal covariance or scale (can be adapted during warmup).
- `logpredictions`: Final log-likelihood predictions.

This function can be used to adapt the proposal distribution (e.g., empirical covariance) before the main MCMC sampling phase.
"""
function warmup(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime)
    parout = Array{Float64,2}(undef, length(param), samplesteps)
    prior = logprior(param, model)
    step = 0
    accepttotal = 0
    while step < samplesteps && time() - t1 < maxtime
        step += 1
        accept, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)
        if ll > llml
            llml, parml = ll, param
        end
        parout[:, step] = param
        accepttotal += accept
    end
    # Compute and scale covariance for proposal adaptation
    d_param = length(param)
    covparam = cov(parout[:,1:step]')
    scaling = (2.38^2) / d_param
    if step > 1000 && accepttotal/step > .15
        if isposdef(covparam)
            proposalcv = covparam * scaling
            d = proposal_dist(param, proposalcv, model)
        else
            # Fallback: use scaled diagonal covariance
            diag_cov = Diagonal(diag(covparam) * scaling)
            proposalcv = diag_cov
            d = proposal_dist(param, diag_cov, model)
        end
        @info "Updated proposal covariance (scaled)" proposalcv
    end
    return param, parml, ll, llml, d, proposalcv, logpredictions
end

"""
sample(logpredictions,param,parml,ll,llml,d,proposalcv,data,model,samplesteps,temp,t1,maxtime)

returns Fit structure

ll is negative loglikelihood

"""
function sample(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime, SLAB=10000)
    llout = Array{Float64,1}(undef, SLAB)
    parout = Array{Float64,2}(undef, length(param), SLAB)
    pwaic = (0, log.(max.(logpredictions, eps(Float64))), zeros(length(logpredictions)))
    lppd = fill(-Inf, length(logpredictions))
    accepttotal = 0
    prior = logprior(param, model)
    step = 0
    while step < samplesteps && time() - t1 < maxtime
        step += 1
        accept, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)

        if step > length(llout)
            parout = hcat(parout, Matrix{Float64}(undef, length(param), SLAB))
            resize!(llout, length(llout) + SLAB)
        end

        if ll > llml
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


# """
#     sample_with_thinning(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime, thin_interval=10)

# Run the MCMC sampling phase, but only keep every `thin_interval`-th sample (thinning).
# Returns a `Fit` structure with thinned samples.
# """
# function sample_with_thinning(logpredictions, param, parml, ll, llml, d, proposalcv, data, model, samplesteps, temp, t1, maxtime, thin_interval=10)
#     # Number of samples to keep after thinning
#     kept_samples = div(samplesteps, thin_interval)

#     llout = Array{Float64,1}(undef, kept_samples)
#     parout = Array{Float64,2}(undef, length(param), kept_samples)

#     # WAIC components
#     pwaic = (0, log.(max.(logpredictions, eps(Float64))), zeros(length(logpredictions)))
#     lppd = fill(-Inf, length(logpredictions))

#     accepttotal = 0
#     prior = logprior(param, model)
#     total_step = 0
#     saved_step = 0

#     while total_step < samplesteps && time() - t1 < maxtime
#         total_step += 1

#         # MCMC step
#         accept, logpredictions, param, ll, prior, d = mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)

#         # Update maximum likelihood
#         if ll > llml
#             llml, parml = ll, param
#         end

#         # Save only every thin_interval steps
#         if total_step % thin_interval == 0
#             saved_step += 1
#             parout[:, saved_step] = param
#             llout[saved_step] = ll
#             lppd, pwaic = update_waic(lppd, pwaic, logpredictions)
#         end

#         accepttotal += accept
#     end

#     # Adjust for actual number of samples saved
#     if saved_step < kept_samples
#         parout = parout[:, 1:saved_step]
#         llout = llout[1:saved_step]
#     end

#     # Adjust WAIC components
#     pwaic = saved_step > 1 ? pwaic[3] / (saved_step - 1) : pwaic[3]
#     lppd .-= log(saved_step)

#     return Fit(parout, llout, parml, llml, lppd, pwaic, prior, accepttotal, total_step)
# end


"""
mhstep(logpredictions,param,ll,prior,d,sigma,model,data,temp)

returns 1,logpredictionst,paramt,llt,priort,dt if accept
        0,logpredictions,param,ll,prior,d if not

ll is log-likelihood (NOT negative log-likelihood)
"""
function mhstep(logpredictions, param, ll, prior, d, proposalcv, model, data, temp)
    paramt, dt = proposal(d, proposalcv, model)
    if instant_reject(paramt, model)
        return 0, logpredictions, param, ll, prior, d
    end
    priort = logprior(paramt, model)
    llt, logpredictionst = loglikelihood(paramt, data, model)
    mhstep(logpredictions, logpredictionst, ll, llt, param, paramt, prior, priort, d, dt, temp)
end

function mhstep(logpredictions, logpredictionst, ll, llt, param, paramt, prior, priort, d, dt, temp)
    # mhfactor not needed for symmetric proposal distribution
    # if rand() < exp((llt + priort - ll - prior + mhfactor(param,d,paramt,dt))/temp)
    
    # println(mhfactor(param,d,paramt,dt))
    # println(llt + priort - ll - prior)

    # Accept if new log-likelihood is higher (or by probability if lower)
    if log(rand()) < (llt + priort - ll - prior) / temp
        return 1, logpredictionst, paramt, llt, priort, dt
    else
        return 0, logpredictions, param, ll, prior, d
    end
end

"""
    instant_reject(paramt, model; minval=1e-12, maxval=1e8)

Return `true` if any parameter in `paramt` is outside [minval, maxval], or if any parameter changes by more than `reltol` (default 0.5 = 50%) relative to `param` (for |param| > abstol). For parameters near zero, only the absolute bound applies.
"""
function instant_reject(paramt, model; minval=-Inf, maxval=Inf)
    n_rates = num_rates(model)
    if n_rates isa Vector
        n_rates = sum(n_rates)
    end
    # Expand paramt to full rate vector (all rates, in model order)
    # Assume get_rates(paramt, model) returns the full rate vector in model order
    all_rates = get_rates(paramt, model)
    rates = all_rates[1:n_rates]
    return any(rates .< minval) || any(rates .> maxval) || any(isnan.(rates)) || any(isinf.(rates))
end

# function instant_reject(paramt, param; reltol=0.5, abstol=1e-6)
#     rel_change = abs.(paramt .- param) ./ max.(abs.(param), abstol)
#     return any(rel_change .> reltol)
# end

"""
update_waic(lppd,pwaic,logpredictions)

returns lppd and pwaic, which are the running sum and variance of logpredictions
(logpredictions is log-likelihood)
"""
function update_waic(lppd, pwaic, logpredictions)
    lppd = logsumexp(lppd, logpredictions)
    lppd, var_update(pwaic, logpredictions)
end
"""
    compute_waic(lppd::Array{T}, pwaic::Array{T}, data) where {T}

Compute the Watanabe-Akaike Information Criterion (WAIC) and its standard error.
This is a unified implementation that works for all data types.

# Arguments
- `lppd`: Log pointwise predictive density (log-likelihoods)
- `pwaic`: Pointwise WAIC penalty terms
- `data`: The data object (can be any type)

# Returns
- Tuple of (WAIC, standard error)

# Notes
- For histogram data: 
  - lppd and pwaic are weighted by histogram counts
  - Each bin's contribution is scaled by its count
- For RNA count data: 
  - lppd and pwaic are raw log-likelihoods
  - Each count is an observation
- For trace data: 
  - lppd and pwaic include both trace and RNA components
  - Each time point is an observation
- The standard error returned is for the **total WAIC** (not per observation), and is scaled by the square root of the number of observations (n_obs).
"""
function compute_waic(lppd::Array{T}, pwaic::Array{T}, data) where {T}
    # Calculate total WAIC
    waic = -2 * sum(lppd - pwaic)
    # Calculate standard error for total WAIC
    n_obs = length(pwaic)  # number of observations
    se = 2 * sqrt(max(sum(pwaic), 0.0)) * sqrt(n_obs)  # scale by sqrt(n_obs)
    return waic, se
end

"""
    aic(fits::Fit)

Compute the Akaike Information Criterion (AIC) for a given fit.
"""
aic(fits::Fit) = 2 * length(fits.parml) - 2 * fits.llml

"""
    aic(nparams::Int, llml)

Compute the Akaike Information Criterion (AIC) given number of parameters and maximum log-likelihood.
"""
aic(nparams::Int, llml) = 2 * nparams - 2 * llml

"""
initial_proposal(model)

return parameters to be fitted and an initial proposal distribution

"""
function initial_proposal(model)
    param = get_param(model)
    d = proposal_dist(param, model.proposal, model)
    return param, d
end


function proposal(d::Distribution, cv, model)
    param = rand(d)
    return param, proposal_dist(param, cv, model)
end

function proposal(d::Vector{T}, cv, model) where {T<:Distribution}
    param = Float64[]
    for i in eachindex(d)
        param = vcat(param..., rand(d[i]))
    end
    return param, proposal_dist(param, cv, model)
end


"""
proposal_dist(param, cv, model)

Construct a proposal distribution centered at `param` with scale/covariance `cv`.
Returns a distribution object.
"""
proposal_dist(param::Float64, cv::Float64, model) = Normal(param, proposal_scale(cv, model))
"""
    proposal_dist(param::Vector, cv::Float64, model)

Construct a product of Normal proposal distributions for each parameter in `param` with shared scale `cv`.
Returns a product distribution object.
"""
function proposal_dist(param::Vector, cv::Float64, model)
    d = Vector{Normal{Float64}}(undef, 0)
    for i in eachindex(param)
        push!(d, Normal(param[i], proposal_scale(cv, model, i, param[i])))
    end
    product_distribution(d)
end
function proposal_dist(param::Vector, cv::Vector, model)
    d = Vector{Normal{Float64}}(undef, 0)
    for i in eachindex(param)
        push!(d, Normal(param[i], proposal_scale(cv[i], model, i, param[i])))
    end
    product_distribution(d)
end
"""
    proposal_dist(param::Vector, cov::Matrix, model)

Construct a multivariate Normal proposal distribution for `param` with covariance matrix `cov`.
Returns an `MvNormal` distribution object.
"""
function proposal_dist(param::Vector, cv::Matrix, model, indiv=0.001)
    n_hyper = size(cv, 1)
    n_param = length(param)
    if hastrait(model, :hierarchical) && n_param > n_hyper
        # Joint proposal for hyperparameters
        hyper_proposal = MvNormal(param[1:n_hyper], cv)
        # Independent proposals for individual-level parameters
        indiv_proposals = proposal_dist(param[n_hyper+1:n_param], indiv, model)
        # Combine into a product distribution
        return [hyper_proposal; indiv_proposals]
    else
        return MvNormal(param, cv)
    end
end

function proposal_dist(param::Vector, cv::Tuple, model)
    proposal_dist(param, cv[1], model, cv[2])
end
"""
proposal_scale(cv::Float64, model::AbstractGeneTransitionModel)

Compute the standard deviation for the proposal distribution given coefficient of variation `cv`.
Returns a Float64.
"""
proposal_scale(cv::Float64, model::AbstractGeneTransitionModel, i=1, param=1.) = sqrt(log(1 + cv^2))

function proposal_scale(cv::Float64, model::GRSMmodel, i, param)
        f = model.transforms.f_cv[i]
        if f == sigmanormal
            return f(param, cv)
        else
            return f(cv)
        end
end
"""
    mhfactor(param, d, paramt, dt)

Compute the Metropolis-Hastings correction factor for asymmetric proposal distributions.
Returns the log ratio of proposal densities.
"""
mhfactor(param, d, paramt, dt) = logpdf(dt, param) - logpdf(d, paramt)


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

Concatenate parameter samples from multiple `Fit` objects into a single array.
Returns a matrix of parameters.
"""
function merge_param_old(fits::Vector)
    param = fits[1].param
    for i in 2:length(fits)
        param = [param fits[i].param]
    end
    return param
end
"""
    merge_param(fits::Vector)

Horizontally concatenate the `param` matrices from a collection of
`Fit` objects (possibly of differing lengths).  Returns a
`d × total_samples` matrix.
"""
function merge_param(fits::Vector)
    d = size(fits[1].param, 1)
    @assert all(size(f.param, 1) == d for f in fits) "parameter dimension differs"

    n_tot = sum(size(f.param, 2) for f in fits)
    out = Matrix{Float64}(undef, d, n_tot)

    col = 1
    for f in fits
        n = size(f.param, 2)
        @inbounds copyto!(out, (col - 1) * d + 1, f.param, 1, n * d)
        col += n
    end
    return out
end

"""
    merge_ll(fits::Vector)

Concatenate log-likelihood samples from multiple `Fit` objects into a single array.
Returns a vector of log-likelihoods.
"""
function merge_ll_old(fits::Vector)
    ll = fits[1].ll
    for i in 2:length(fits)
        ll = [ll; fits[i].ll]
    end
    return ll
end
"""
    merge_ll(fits::Vector)

Concatenate the `ll` vectors from a collection of `Fit` objects into
one long `Vector{Float64}`.
"""
function merge_ll(fits::Vector)
    n_tot = sum(length(f.ll) for f in fits)
    out = Vector{Float64}(undef, n_tot)

    pos = 1
    for f in fits
        n = length(f.ll)
        @inbounds copyto!(out, pos, f.ll, 1, n)
        pos += n
    end
    return out
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
    param = inverse_transform_params(paramin, model)
    meanparam = mean(param, dims=2)
    stdparam = std(param, dims=2)
    medparam = median(param, dims=2)
    np = size(param, 1)
    madparam = Array{Float64,1}(undef, np)
    qparam = Matrix{Float64}(undef, 3, 0)
    for i in 1:np
        madparam[i] = mad(param[i, :], normalize=false)
        qparam = hcat(qparam, quantile(param[i, :], [0.025; 0.5; 0.975]))
    end
    # Thin for covariance/correlation calculations
    param_thin = thin_columns(param)
    if typeof(model) <: AbstractGRSMmodel && hastrait(model, :hierarchical)
        nrates = num_fitted_core_params(model)
        corparam = cor(param_thin[1:nrates, :]')
        covparam = cov(param_thin[1:nrates, :]')
        covlogparam = cov(paramin[1:nrates, :]')
    else
        corparam = cor(param_thin')
        covparam = cov(param_thin')
        covlogparam = cov(paramin')
    end
    Stats(meanparam, stdparam, medparam, madparam, qparam, corparam, covparam, covlogparam)
end

# function compute_stats(paramin::Array{Float64,2}, model::AbstractGRSMhierarchicalmodel)
#     param = inverse_transform_params(paramin, model)
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
function compute_rhat_old(params::Vector{Array})
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


function compute_measures(chain::Array{Tuple,1})
    fits = collate_fit(chain)
    return compute_measures(fits)
end

function compute_measures(fits::Vector{Fit})
    # Each of these is a vector of vectors (one per chain)
    ess_list = [compute_ess(f.param) for f in fits]
    geweke_list = [compute_geweke(f.param) for f in fits]
    mcse_list = [compute_mcse(f.param) for f in fits]

    # Stack into matrices (parameters x chains)
    ess_mat = reduce(hcat, ess_list)
    geweke_mat = reduce(hcat, geweke_list)
    mcse_mat = reduce(hcat, mcse_list)

    # Pooled/summarized diagnostics
    pooled_ess = sum(ess_mat, dims=2)[:, 1]           # sum across chains for each parameter
    mean_geweke = mean(geweke_mat, dims=2)[:, 1]      # mean across chains for each parameter
    mean_mcse = mean(mcse_mat, dims=2)[:, 1]          # mean across chains for each parameter

    return pooled_ess, mean_geweke, mean_mcse
end

function compute_measures(fits::Fit)
    return compute_ess(fits.param), compute_geweke(fits.param), compute_mcse(fits.param)
end

# function Measures(fits::Fit, waic, rhat)
#     ess = compute_ess(fits.param)
#     geweke = compute_geweke(fits.param)
#     mcse = compute_mcse(fits.param)
#     return Measures(waic, rhat, ess, geweke, mcse)
# end
"""
    compute_ess(params)

Compute the effective sample size (ESS) for each parameter in the MCMC chain.
Returns a vector of ESS values.
"""
# function compute_ess(params)
#     n_params = size(params, 1)
#     n_samples = size(params, 2)
#     ess = zeros(n_params)

#     for i in 1:n_params
#         # Calculate autocorrelation
#         acf = autocor(params[i, :], 0:min(1000, n_samples - 1))
#         # Find where autocorrelation drops below 0.05
#         cutoff = findfirst(x -> abs(x) < 0.05, acf)
#         cutoff = isnothing(cutoff) ? length(acf) : cutoff

#         # Sum of autocorrelations with appropriate cutoff
#         tau = 1 + 2 * sum(acf[2:cutoff])
#         ess[i] = n_samples / tau
#     end

#     return ess
# end

"""
    compute_geweke(params; first_frac=0.1, last_frac=0.5)

Compute Geweke z-scores for each parameter to assess MCMC convergence.
Returns a vector of z-scores.
"""
# function compute_geweke(params; first_frac=0.1, last_frac=0.5)
#     n_params = size(params, 1)
#     n_samples = size(params, 2)

#     first_end = Int(floor(first_frac * n_samples))
#     last_start = Int(floor((1 - last_frac) * n_samples))

#     z_scores = zeros(n_params)

#     for i in 1:n_params
#         first_mean = mean(params[i, 1:first_end])
#         last_mean = mean(params[i, last_start:end])

#         # Spectral density estimates for variance
#         first_var = spectrum0(params[i, 1:first_end])
#         last_var = spectrum0(params[i, last_start:end])

#         # Z-score
#         z_scores[i] = (first_mean - last_mean) /
#                       sqrt(first_var / first_end + last_var / (n_samples - last_start + 1))
#     end

#     return z_scores
# end

function compute_geweke(params; first_frac=0.1, last_frac=0.5)
    n_params = size(params, 1)
    n_samples = size(params, 2)

    min_samples_needed = 10 # You can adjust this minimum as needed

    if n_samples < min_samples_needed
        @warn "Insufficient number of samples ($n_samples) to compute Geweke diagnostic. Skipping calculation."
        return fill(NaN, n_params) # Return NaN for all parameters
    end

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
        denominator = sqrt(first_var / first_end + last_var / (n_samples - last_start + 1))
        z_scores[i] = iszero(denominator) ? NaN : (first_mean - last_mean) / denominator
    end

    return z_scores
end

"""
    spectrum0(x)

Estimate the spectral density at frequency zero for a time series `x`.
Used in Geweke diagnostic.
"""
function spectrum0(x)
    n = length(x)
    spec = abs.(fft(x .- mean(x))) .^ 2 / n
    return sum(spec[2:Int(floor(n / 2))]) * 2 / n
end

"""
    compute_mcse(params)

Compute the Monte Carlo standard error (MCSE) for each parameter.
Returns a vector of MCSE values.
"""
# function compute_mcse(params)
#     n_params = size(params, 1)
#     n_samples = size(params, 2)
#     mcse = zeros(n_params)

#     for i in 1:n_params
#         # Calculate autocorrelation
#         acf = autocor(params[i, :], 0:min(1000, n_samples - 1))
#         # Find where autocorrelation drops below 0.05
#         cutoff = findfirst(x -> abs(x) < 0.05, acf)
#         cutoff = isnothing(cutoff) ? length(acf) : cutoff

#         # Sum of autocorrelations
#         tau = 1 + 2 * sum(acf[2:cutoff])

#         # Standard error
#         val = tau * var(params[i, :]) / n_samples
#         mcse[i] = sqrt(max(val, 0.0))
#     end

#     return mcse
# end



function compute_ess(params)
    n_params = size(params, 1)
    n_samples = size(params, 2)
    min_samples_needed = 10 # Adjust as needed
    ess = zeros(n_params)

    if n_samples < min_samples_needed
        @warn "Insufficient number of samples ($n_samples) to compute ESS. Returning NaN."
        return fill(NaN, n_params)
    end

    for i in 1:n_params
        acf = autocor(params[i, :], 0:min(1000, n_samples - 1))
        cutoff = findfirst(x -> abs(x) < 0.05, acf)
        cutoff = isnothing(cutoff) ? length(acf) : cutoff
        tau = 1 + 2 * sum(acf[2:cutoff])
        ess[i] = n_samples / tau
    end
    return ess
end

function compute_mcse(params)
    n_params = size(params, 1)
    n_samples = size(params, 2)
    min_samples_needed = 10 # Adjust as needed
    mcse = zeros(n_params)

    if n_samples < min_samples_needed
        @warn "Insufficient number of samples ($n_samples) to compute MCSE. Returning NaN."
        return fill(NaN, n_params)
    end

    for i in 1:n_params
        acf = autocor(params[i, :], 0:min(1000, n_samples - 1))
        cutoff = findfirst(x -> abs(x) < 0.05, acf)
        cutoff = isnothing(cutoff) ? length(acf) : cutoff
        tau = 1 + 2 * sum(acf[2:cutoff])
        val = tau * var(params[i, :]) / n_samples
        mcse[i] = sqrt(max(val, 0.0))
    end
    return mcse
end

function compute_rhat(params::Vector{Array})
    N = chainlength(params)
    M = length(params)
    n_params = size(params[1], 1)
    min_samples_needed = 10 # Adjust as needed

    if M < 2
        @warn "Insufficient number of chains ($M) to compute R-hat. Returning NaN."
        return fill(NaN, n_params)
    elseif N < min_samples_needed
        @warn "Insufficient number of samples per chain ($N) to compute R-hat. Returning NaN."
        return fill(NaN, n_params)
    end

    m = Matrix{Float64}(undef, n_params, 2 * M)
    s = similar(m)
    for i in 1:M
        m[:, 2*i-1] = mean(params[i][:, 1:N], dims=2)
        m[:, 2*i] = mean(params[i][:, N+1:2*N], dims=2)
        s[:, 2*i-1] = var(params[i][:, 1:N], dims=2)
        s[:, 2*i] = var(params[i][:, N+1:2*N], dims=2)
    end
    B = N * var(m, dims=2)
    W = mean(s, dims=2)
    return sqrt.((N - 1) / N .+ B ./ W / N)
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
    find_ml(fits::Array)

Find the parameter vector and log-likelihood with the maximum likelihood across all fits.
Returns `(parml, llml)`.
"""
function find_ml(fits::Array)
    llml = Array{Float64,1}(undef, length(fits))
    for i in eachindex(fits)
        llml[i] = fits[i].llml
    end
    return fits[argmin(llml)].parml, llml[argmin(llml)]
end

"""
    crossentropy(logpredictions, hist)

Compute the cross-entropy between a data histogram and model log-likelihood predictions.
Returns a Float64.
"""
function crossentropy(logpredictions::Array{T1}, hist::Array{T2}) where {T1,T2}
    lltot = hist' * logpredictions
    isfinite(lltot) ? -lltot : Inf
end

"""
    hist_entropy(hist::Array{Float64,1})

Compute the entropy of a normalized histogram.
Returns a Float64.
"""
function hist_entropy(hist::Array{Float64,1})
    -hist' * log.(normalize_histogram(max.(hist, 0)))
end

"""
    hist_entropy(hist::Array{Array,1})

Compute the entropy of a list of normalized histograms.
Returns a Float64.
"""
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

"""
    thin_columns(mat; maxcols=5000)

Thins the columns of `mat` so that the number of columns does not exceed `maxcols`.
Returns the thinned matrix.
"""
function thin_columns(mat; maxcols=5000)
    ncols = size(mat, 2)
    if ncols > maxcols
        thin_factor = div(ncols, maxcols)
        return mat[:, 1:thin_factor:end]
    else
        return mat
    end
end


function estimate_memory_usage(fits::Vector{Fit})
    # Estimate memory usage in bytes
    # Each Float64 is 8 bytes
    # We need to account for both param and ll arrays
    total_bytes = 0
    for f in fits
        param_bytes = size(f.param, 1) * size(f.param, 2) * 8
        ll_bytes = length(f.ll) * 8
        total_bytes += param_bytes + ll_bytes
    end
    return total_bytes
end

function memory_aware_thin(fits::Vector{Fit}, max_memory::Int=48 * 1024^3)
    # First estimate total memory usage
    total_bytes = estimate_memory_usage(fits)

    # If we're already below the limit, no need to thin
    if total_bytes < max_memory
        return fits
    end

    # Calculate required thinning factor
    required_thinning = ceil(Int, total_bytes / max_memory)

    # Apply thinning to each chain
    thinned_fits = similar(fits)
    for i in 1:length(fits)
        f = fits[i]
        param = f.param[:, 1:size(f.param, 2):required_thinning]
        ll = f.ll[1:length(f.ll):required_thinning]
        thinned_fits[i] = Fit(param, ll, f.parml, f.llml, f.lppd, f.pwaic, f.prior, f.accept, f.total)
    end

    return thinned_fits
end

function merge_fit_memory_aware(fits::Array{Fit,1}, max_memory::Int=64 * 1024^3)
    # First thin chains if necessary
    thinned_fits = memory_aware_thin(fits, max_memory)

    # Then proceed with normal merge
    param = merge_param(thinned_fits)
    ll = merge_ll(thinned_fits)
    parml, llml = find_ml(thinned_fits)
    lppd = thinned_fits[1].lppd
    pwaic = thinned_fits[1].pwaic
    prior = thinned_fits[1].prior
    accept = thinned_fits[1].accept
    total = thinned_fits[1].total

    for i in 2:length(thinned_fits)
        accept += thinned_fits[i].accept
        total += thinned_fits[i].total
        prior += thinned_fits[i].prior
        lppd = logsumexp(lppd, thinned_fits[i].lppd)
        pwaic += thinned_fits[i].pwaic
    end

    return Fit(param, ll, parml, llml, lppd .- log(length(thinned_fits)),
        pwaic / length(thinned_fits), prior / length(thinned_fits),
        accept, total)
end