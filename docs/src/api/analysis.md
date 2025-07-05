# Analysis Functions

Functions for post-fitting analysis, model comparison, and result visualization.

## Model Prediction Functions

### predictedarray

```julia
predictedarray(parameters::Vector{Float64}, model::AbstractModel, data::AbstractExperimentalData) -> Vector{Float64}
```

Generate model predictions for given parameters.

**Arguments:**
- `parameters`: Model parameters
- `model`: Model structure
- `data`: Experimental data

**Returns:**
- Vector of predicted values

**Example:**
```julia
# Generate predictions
fitted_params = [0.1, 0.2, 0.3]
predictions = predictedarray(fitted_params, model, data)
```

### predictedfn

```julia
predictedfn(model::AbstractModel, data::AbstractExperimentalData) -> Function
```

Create a function that generates predictions for given parameters.

**Arguments:**
- `model`: Model structure
- `data`: Experimental data

**Returns:**
- Function that takes parameters and returns predictions

**Example:**
```julia
# Create prediction function
pred_fn = predictedfn(model, data)
predictions = pred_fn(fitted_params)
```

### make_traces

```julia
make_traces(parameters::Vector{Float64}, model::AbstractModel, data::AbstractTraceData) -> Vector{Vector{Float64}}
```

Generate model-predicted intensity traces.

**Arguments:**
- `parameters`: Model parameters
- `model`: Model structure
- `data`: Trace data

**Returns:**
- Vector of predicted traces

**Example:**
```julia
# Generate predicted traces
pred_traces = make_traces(fitted_params, model, trace_data)
```

### make_traces_dataframe

```julia
make_traces_dataframe(traces::Vector{Vector{Float64}}, interval::Float64) -> DataFrame
```

Convert traces to DataFrame format for analysis.

**Arguments:**
- `traces`: Vector of trace vectors
- `interval`: Time interval between points

**Returns:**
- DataFrame with time and intensity columns

**Example:**
```julia
# Convert traces to DataFrame
df = make_traces_dataframe(traces, 1.0)
```

## Model Comparison Functions

### loglikelihood

```julia
loglikelihood(parameters::Vector{Float64}, data::AbstractExperimentalData, model::AbstractModel) -> Float64
```

Calculate log-likelihood for given parameters.

**Arguments:**
- `parameters`: Model parameters
- `data`: Experimental data
- `model`: Model structure

**Returns:**
- Log-likelihood value

**Example:**
```julia
# Calculate log-likelihood
ll = loglikelihood(fitted_params, data, model)
```

### deviance

```julia
deviance(parameters::Vector{Float64}, data::AbstractExperimentalData, model::AbstractModel) -> Float64
```

Calculate deviance (-2 * log-likelihood).

**Arguments:**
- `parameters`: Model parameters
- `data`: Experimental data
- `model`: Model structure

**Returns:**
- Deviance value

**Example:**
```julia
# Calculate deviance
dev = deviance(fitted_params, data, model)
```

### aic

```julia
aic(loglik::Float64, k::Int) -> Float64
```

Calculate Akaike Information Criterion.

**Arguments:**
- `loglik`: Log-likelihood value
- `k`: Number of parameters

**Returns:**
- AIC value

**Example:**
```julia
# Calculate AIC
aic_value = aic(log_likelihood, num_params)
```

### bic

```julia
bic(loglik::Float64, k::Int, n::Int) -> Float64
```

Calculate Bayesian Information Criterion.

**Arguments:**
- `loglik`: Log-likelihood value
- `k`: Number of parameters
- `n`: Number of observations

**Returns:**
- BIC value

**Example:**
```julia
# Calculate BIC
bic_value = bic(log_likelihood, num_params, num_obs)
```

## Burst Analysis Functions

### burstsize

```julia
burstsize(fits::Results, model::AbstractModel) -> Vector{Float64}
```

Calculate burst sizes from fitted parameters.

**Arguments:**
- `fits`: Fit results
- `model`: Model structure

**Returns:**
- Vector of burst sizes

**Example:**
```julia
# Calculate burst sizes
bs = burstsize(fit_results, model)
```

### burstfrequency

```julia
burstfrequency(fits::Results, model::AbstractModel) -> Vector{Float64}
```

Calculate burst frequencies from fitted parameters.

**Arguments:**
- `fits`: Fit results
- `model`: Model structure

**Returns:**
- Vector of burst frequencies

**Example:**
```julia
# Calculate burst frequencies
bf = burstfrequency(fit_results, model)
```

### burst_duration

```julia
burst_duration(fits::Results, model::AbstractModel) -> Vector{Float64}
```

Calculate burst durations from fitted parameters.

**Arguments:**
- `fits`: Fit results
- `model`: Model structure

**Returns:**
- Vector of burst durations

**Example:**
```julia
# Calculate burst durations
bd = burst_duration(fit_results, model)
```

## Statistical Analysis Functions

### posterior_mean

```julia
posterior_mean(samples::Matrix{Float64}) -> Vector{Float64}
```

Calculate posterior means from MCMC samples.

**Arguments:**
- `samples`: MCMC samples matrix

**Returns:**
- Vector of posterior means

**Example:**
```julia
# Calculate posterior means
means = posterior_mean(mcmc_samples)
```

### posterior_std

```julia
posterior_std(samples::Matrix{Float64}) -> Vector{Float64}
```

Calculate posterior standard deviations from MCMC samples.

**Arguments:**
- `samples`: MCMC samples matrix

**Returns:**
- Vector of posterior standard deviations

**Example:**
```julia
# Calculate posterior standard deviations
stds = posterior_std(mcmc_samples)
```

### posterior_quantiles

```julia
posterior_quantiles(samples::Matrix{Float64}, quantiles::Vector{Float64}) -> Matrix{Float64}
```

Calculate posterior quantiles from MCMC samples.

**Arguments:**
- `samples`: MCMC samples matrix
- `quantiles`: Vector of quantile values (0-1)

**Returns:**
- Matrix of quantiles

**Example:**
```julia
# Calculate credible intervals
quantiles = posterior_quantiles(mcmc_samples, [0.025, 0.975])
```

### effective_sample_size

```julia
effective_sample_size(samples::Vector{Float64}) -> Float64
```

Calculate effective sample size for MCMC chain.

**Arguments:**
- `samples`: MCMC samples for single parameter

**Returns:**
- Effective sample size

**Example:**
```julia
# Calculate effective sample size
ess = effective_sample_size(mcmc_samples[:, 1])
```

### rhat

```julia
rhat(chains::Vector{Vector{Float64}}) -> Float64
```

Calculate R-hat convergence diagnostic.

**Arguments:**
- `chains`: Vector of MCMC chains

**Returns:**
- R-hat value

**Example:**
```julia
# Calculate R-hat
rhat_value = rhat([chain1, chain2, chain3])
```

## Residual Analysis Functions

### residuals

```julia
residuals(data::AbstractExperimentalData, predictions::Vector{Float64}) -> Vector{Float64}
```

Calculate residuals between data and predictions.

**Arguments:**
- `data`: Experimental data
- `predictions`: Model predictions

**Returns:**
- Vector of residuals

**Example:**
```julia
# Calculate residuals
resids = residuals(data, predictions)
```

### standardized_residuals

```julia
standardized_residuals(residuals::Vector{Float64}) -> Vector{Float64}
```

Calculate standardized residuals.

**Arguments:**
- `residuals`: Raw residuals

**Returns:**
- Standardized residuals

**Example:**
```julia
# Calculate standardized residuals
std_resids = standardized_residuals(resids)
```

### qq_plot_data

```julia
qq_plot_data(residuals::Vector{Float64}) -> Tuple{Vector{Float64}, Vector{Float64}}
```

Generate data for Q-Q plots.

**Arguments:**
- `residuals`: Residuals for analysis

**Returns:**
- Tuple of (theoretical quantiles, sample quantiles)

**Example:**
```julia
# Generate Q-Q plot data
theoretical, sample = qq_plot_data(residuals)
```

## Model Diagnostic Functions

### convergence_diagnostics

```julia
convergence_diagnostics(chains::Vector{Matrix{Float64}}) -> Dict{String, Any}
```

Calculate comprehensive convergence diagnostics.

**Arguments:**
- `chains`: Vector of MCMC chains

**Returns:**
- Dictionary with diagnostic results

**Example:**
```julia
# Calculate convergence diagnostics
diagnostics = convergence_diagnostics(mcmc_chains)
```

### trace_plots

```julia
trace_plots(chains::Vector{Matrix{Float64}}, parameter_names::Vector{String}) -> Vector{Plot}
```

Generate trace plots for MCMC chains.

**Arguments:**
- `chains`: Vector of MCMC chains
- `parameter_names`: Names of parameters

**Returns:**
- Vector of trace plots

**Example:**
```julia
# Generate trace plots
plots = trace_plots(mcmc_chains, ["rate1", "rate2", "noise"])
```

### autocorrelation

```julia
autocorrelation(samples::Vector{Float64}, max_lag::Int=50) -> Vector{Float64}
```

Calculate autocorrelation function for MCMC samples.

**Arguments:**
- `samples`: MCMC samples
- `max_lag`: Maximum lag for autocorrelation

**Returns:**
- Vector of autocorrelation values

**Example:**
```julia
# Calculate autocorrelation
acf = autocorrelation(mcmc_samples[:, 1], 100)
```

## Sensitivity Analysis Functions

### parameter_sensitivity

```julia
parameter_sensitivity(model::AbstractModel, data::AbstractExperimentalData, 
                     parameters::Vector{Float64}, delta::Float64=0.01) -> Vector{Float64}
```

Calculate parameter sensitivity.

**Arguments:**
- `model`: Model structure
- `data`: Experimental data
- `parameters`: Parameter values
- `delta`: Perturbation size

**Returns:**
- Vector of sensitivity values

**Example:**
```julia
# Calculate parameter sensitivity
sensitivity = parameter_sensitivity(model, data, fitted_params, 0.05)
```

### local_sensitivity_analysis

```julia
local_sensitivity_analysis(model::AbstractModel, data::AbstractExperimentalData,
                          parameters::Vector{Float64}) -> Matrix{Float64}
```

Perform local sensitivity analysis.

**Arguments:**
- `model`: Model structure
- `data`: Experimental data
- `parameters`: Parameter values

**Returns:**
- Sensitivity matrix

**Example:**
```julia
# Perform local sensitivity analysis
sens_matrix = local_sensitivity_analysis(model, data, fitted_params)
```

## Profile Likelihood Functions

### profile_likelihood

```julia
profile_likelihood(parameter_index::Int, values::Vector{Float64}, 
                  model::AbstractModel, data::AbstractExperimentalData) -> Vector{Float64}
```

Calculate profile likelihood for a parameter.

**Arguments:**
- `parameter_index`: Index of parameter to profile
- `values`: Values to evaluate
- `model`: Model structure
- `data`: Experimental data

**Returns:**
- Vector of profile likelihood values

**Example:**
```julia
# Calculate profile likelihood
profile = profile_likelihood(1, [0.05, 0.1, 0.15, 0.2], model, data)
```

### confidence_intervals

```julia
confidence_intervals(profile_ll::Vector{Float64}, values::Vector{Float64}, 
                    confidence_level::Float64=0.95) -> Tuple{Float64, Float64}
```

Calculate confidence intervals from profile likelihood.

**Arguments:**
- `profile_ll`: Profile likelihood values
- `values`: Parameter values
- `confidence_level`: Confidence level (0-1)

**Returns:**
- Tuple of (lower bound, upper bound)

**Example:**
```julia
# Calculate 95% confidence intervals
lower, upper = confidence_intervals(profile_ll, param_values, 0.95)
```

## Model Selection Functions

### model_comparison

```julia
model_comparison(results::Vector{Results}, models::Vector{AbstractModel}) -> DataFrame
```

Compare multiple models using information criteria.

**Arguments:**
- `results`: Vector of fit results
- `models`: Vector of models

**Returns:**
- DataFrame with comparison results

**Example:**
```julia
# Compare models
comparison = model_comparison([results1, results2], [model1, model2])
```

### cross_validation

```julia
cross_validation(model::AbstractModel, data::AbstractExperimentalData, 
                k_folds::Int=5) -> Vector{Float64}
```

Perform k-fold cross-validation.

**Arguments:**
- `model`: Model structure
- `data`: Experimental data
- `k_folds`: Number of folds

**Returns:**
- Vector of cross-validation scores

**Example:**
```julia
# Perform 5-fold cross-validation
cv_scores = cross_validation(model, data, 5)
```

## Visualization Support Functions

### plot_fits

```julia
plot_fits(data::AbstractExperimentalData, predictions::Vector{Float64}, 
         residuals::Vector{Float64}) -> Plot
```

Generate fit quality plots.

**Arguments:**
- `data`: Experimental data
- `predictions`: Model predictions
- `residuals`: Residuals

**Returns:**
- Plot object

**Example:**
```julia
# Generate fit plots
plot = plot_fits(data, predictions, residuals)
```

### plot_traces

```julia
plot_traces(traces::Vector{Vector{Float64}}, interval::Float64) -> Plot
```

Generate trace plots.

**Arguments:**
- `traces`: Vector of traces
- `interval`: Time interval

**Returns:**
- Plot object

**Example:**
```julia
# Plot traces
plot = plot_traces(predicted_traces, 1.0)
```

### plot_distributions

```julia
plot_distributions(data::AbstractHistogramData, predictions::Vector{Float64}) -> Plot
```

Generate distribution comparison plots.

**Arguments:**
- `data`: Histogram data
- `predictions`: Model predictions

**Returns:**
- Plot object

**Example:**
```julia
# Plot distributions
plot = plot_distributions(rna_data, predictions)
```

## Complete Analysis Examples

### RNA Data Analysis

```julia
using StochasticGene

# Load data and fit model
data = load_data("rna", String[], "data/", "exp", "MYC", "ctrl", (), 1)
fits, stats, measures, data, model, options = fit(nchains=4)

# Generate predictions
predictions = predictedarray(stats.medparam, model, data)

# Calculate residuals
resids = residuals(data, predictions)

# Model diagnostics
ll = loglikelihood(stats.medparam, data, model)
dev = deviance(stats.medparam, data, model)
aic_val = aic(ll, length(stats.medparam))

# Burst analysis
bs = burstsize(fits, model)
bf = burstfrequency(fits, model)

# Convergence diagnostics
diagnostics = convergence_diagnostics(fits.chains)
```

### Trace Data Analysis

```julia
# Load trace data
trace_data = load_data("trace", String[], "data/traces/", "exp", "SOX2", "ctrl", 
                      (1.0, 0.0, 100.0, 0.9), 1)

# Fit model
fits, stats, measures, data, model, options = fit(
    datatype="trace", 
    nchains=4,
    noisepriors=[20.0, 10.0]
)

# Generate predicted traces
pred_traces = make_traces(stats.medparam, model, trace_data)

# Convert to DataFrame
trace_df = make_traces_dataframe(pred_traces, 1.0)

# Sensitivity analysis
sensitivity = parameter_sensitivity(model, trace_data, stats.medparam)
```

### Hierarchical Model Analysis

```julia
# Fit hierarchical model
fits, stats, measures, data, model, options = fit(
    datatype="trace",
    hierarchical=(2, [1,2]),  # 2 hyperparameter sets, fit rates 1,2
    nchains=4
)

# Analyze hierarchical results
individual_params = extract_individual_parameters(fits, model)
population_params = extract_population_parameters(fits, model)

# Population-level analysis
pop_means = posterior_mean(population_params)
pop_stds = posterior_std(population_params)
```

### Model Comparison Workflow

```julia
# Compare different models
models = [
    (G=2, R=0, S=0),  # Simple telegraph
    (G=2, R=1, S=0),  # With pre-RNA
    (G=3, R=0, S=0),  # Three states
]

results = []
for (G, R, S) in models
    fits, stats, measures, data, model, options = fit(
        G=G, R=R, S=S, 
        nchains=4
    )
    push!(results, (fits, stats, measures, model))
end

# Compare using information criteria
comparison = model_comparison(
    [r[2] for r in results],  # stats
    [r[4] for r in results]   # models
)
```

## See Also

- [`fit`](@ref): Main fitting function
- [`simulator`](@ref): Model simulation
- [`write_traces`](@ref): Trace generation
- [`utilities`](@ref): Utility functions