
# Unified inference entry point for StochasticGene.jl

"""
	run_inference(data, model, options; rng=Random.default_rng())

Unified inference entry point. Dispatches to the appropriate inference method by multiple dispatch on model, data, and option type.
Returns standardized (fits, stats, measures, diagnostics).
"""
function run_inference(data, model, options; rng=Random.default_rng())
	run_inference(data, model, options, rng)
end

# Default fallback: error if no method matches
function run_inference(data, model, options, rng)
	error("No run_inference method for model=$(typeof(model)), data=$(typeof(data)), options=$(typeof(options))")
end

# MH dispatch
function run_inference(data, model, options::MHOptions, rng)
	return _run_mh_inference(data, model, options, rng)
end

# NUTS dispatch
function run_inference(data, model, options::NUTSOptions, rng)
	return _run_nuts_inference(data, model, options, rng)
end

# ADVI dispatch
function run_inference(data, model, options::ADVIOptions, rng)
	return _run_advi_inference(data, model, options, rng)
end

# --- Method-specific runners ---

# Example MH runner (expand as needed)
function _run_mh_inference(data, model, options::MHOptions, rng)
	if options.device == :gpu
		if options.nchains == 1
			return run_mh_gpu(data, model, options, 1)
		else
			return run_mh_gpu(data, model, options, options.nchains)
		end
	elseif options.parallel == :distributed
		if options.nchains == 1
			return run_mh(data, model, options)
		else
			return Distributed.pmap(1:options.nchains) do _
				run_mh(data, model, options)
			end
		end
	elseif options.parallel == :threaded
		if options.nchains == 1
			return run_mh(data, model, options)
		else
			results = Vector{Any}(undef, options.nchains)
			Threads.@threads for i in 1:options.nchains
				results[i] = run_mh(data, model, options)
			end
			return results
		end
	else
		if options.nchains == 1
			return run_mh(data, model, options)
		else
			return run_mh(data, model, options, options.nchains)
		end
	end
end

function _run_nuts_inference(data, model, options::NUTSOptions, rng)
	if options.device == :gpu
		error("NUTS inference on GPU is not implemented.")
	elseif options.parallel == :distributed
		@info "Running NUTS inference with Distributed parallelism."
		nchains = getfield(options, :nchains, 1)
		chains = Distributed.pmap(1:nchains) do _
			run_nuts(data, model, rng, options)
		end
		return chains
	elseif options.parallel == :threaded
		@info "Running NUTS inference with multithreading."
		nchains = getfield(options, :nchains, 1)
		results = Vector{Any}(undef, nchains)
		Threads.@threads for i in 1:nchains
			results[i] = run_nuts(data, model, rng, options)
		end
		return results
	else
		return run_nuts(data, model, rng, options)
	end
end

function _run_advi_inference(data, model, options::ADVIOptions, rng)
	if options.device == :gpu
		error("ADVI inference on GPU is not implemented.")
	elseif options.parallel == :distributed
		@info "Running ADVI inference with Distributed parallelism."
		nchains = getfield(options, :nchains, 1)
		chains = Distributed.pmap(1:nchains) do _
			run_advi(data, model, rng, options)
		end
		return chains
	elseif options.parallel == :threaded
		@info "Running ADVI inference with multithreading."
		nchains = getfield(options, :nchains, 1)
		results = Vector{Any}(undef, nchains)
		Threads.@threads for i in 1:nchains
			results[i] = run_advi(data, model, rng, options)
		end
		return results
	else
		return run_advi(data, model, rng, options)
	end
end

# Add result collation, diagnostics, and option parsing helpers here as needed.
