# Utility Functions

Collection of utility functions for data processing, model construction, and analysis.

## Data Processing Functions

### normalize_histogram

```julia
normalize_histogram(histogram::Vector{Float64}) -> Vector{Float64}
```

Normalize a histogram to create a probability distribution.

**Arguments:**
- `histogram`: Vector of histogram counts

**Returns:**
- Normalized probability distribution (sums to 1.0)

**Example:**
```julia
# Normalize RNA count histogram
raw_counts = [100, 150, 200, 120, 80, 50]
prob_dist = normalize_histogram(raw_counts)
println("Sum: ", sum(prob_dist))  # Should be 1.0
```

### make_array

```julia
make_array(data) -> Array
```

Convert various data types to arrays for processing.

**Arguments:**
- `data`: Data to convert (Vector, DataFrame, etc.)

**Returns:**
- Array representation of the data

**Example:**
```julia
# Convert DataFrame column to array
df = DataFrame(counts = [1, 2, 3, 4, 5])
arr = make_array(df.counts)
```

### make_mat

```julia
make_mat(data) -> Matrix
```

Convert data to matrix format for analysis.

**Arguments:**
- `data`: Data to convert

**Returns:**
- Matrix representation of the data

**Example:**
```julia
# Convert vector of vectors to matrix
traces = [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]]
mat = make_mat(traces)
```

### digit_vector

```julia
digit_vector(number::Int, base::Int=10) -> Vector{Int}
```

Convert a number to a vector of digits in specified base.

**Arguments:**
- `number`: Number to convert
- `base`: Base for conversion (default: 10)

**Returns:**
- Vector of digits

**Example:**
```julia
# Convert number to digits
digits = digit_vector(12345)  # [1, 2, 3, 4, 5]
binary = digit_vector(10, 2)  # [1, 0, 1, 0]
```

## Model Construction Functions

### prepare_rates

```julia
prepare_rates(rates::Vector{Float64}, model::AbstractModel) -> Vector{Float64}
```

Prepare rate parameters for model fitting by applying transformations.

**Arguments:**
- `rates`: Raw rate parameters
- `model`: Model structure

**Returns:**
- Transformed rates ready for fitting

**Example:**
```julia
# Prepare rates for fitting
raw_rates = [0.1, 0.2, 0.3]
prepared = prepare_rates(raw_rates, model)
```

### get_rates

```julia
get_rates(parameters::Vector{Float64}, model::AbstractModel, log_scale::Bool=true) -> Vector{Float64}
```

Extract rate parameters from fitted parameters.

**Arguments:**
- `parameters`: Fitted parameter vector
- `model`: Model structure
- `log_scale`: Whether parameters are in log scale

**Returns:**
- Rate parameters

**Example:**
```julia
# Extract rates from fitted parameters
fitted_params = [log(0.1), log(0.2), log(0.3)]
rates = get_rates(fitted_params, model, true)
```

### get_param

```julia
get_param(parameters::Vector{Float64}, model::AbstractModel, param_type::String) -> Vector{Float64}
```

Extract specific parameter types from fitted parameters.

**Arguments:**
- `parameters`: Fitted parameter vector
- `model`: Model structure
- `param_type`: Type of parameters to extract ("rates", "noise", "coupling")

**Returns:**
- Requested parameters

**Example:**
```julia
# Extract noise parameters
noise_params = get_param(fitted_params, model, "noise")
```

### num_rates

```julia
num_rates(transitions::Tuple, R::Int, S::Int, insertstep::Int) -> Int
```

Count the number of transition rates in a model.

**Arguments:**
- `transitions`: Model transitions
- `R`: Number of pre-RNA steps
- `S`: Number of splice sites
- `insertstep`: Reporter insertion step

**Returns:**
- Total number of rates

**Example:**
```julia
# Count rates for a GRS model
n_rates = num_rates(([1,2], [2,1]), 3, 2, 1)
```

### num_all_parameters

```julia
num_all_parameters(transitions::Tuple, R::Int, S::Int, insertstep::Int, 
                  reporter, coupling::Tuple, grid) -> Int
```

Count total parameters including rates, noise, coupling, and grid parameters.

**Arguments:**
- `transitions`: Model transitions
- `R`, `S`, `insertstep`: Model structure
- `reporter`: Reporter configuration
- `coupling`: Coupling specification
- `grid`: Grid specification

**Returns:**
- Total parameter count

**Example:**
```julia
# Count all parameters
n_total = num_all_parameters(transitions, 2, 1, 1, reporter, coupling, grid)
```

## Statistical Functions

### prob_Gaussian

```julia
prob_Gaussian(y::Float64, μ::Float64, σ::Float64) -> Float64
```

Calculate Gaussian probability density.

**Arguments:**
- `y`: Observed value
- `μ`: Mean
- `σ`: Standard deviation

**Returns:**
- Probability density

**Example:**
```julia
# Calculate Gaussian probability
prob = prob_Gaussian(2.0, 1.5, 0.5)
```

### prob_Gaussian_grid

```julia
prob_Gaussian_grid(y::Vector{Float64}, μ::Vector{Float64}, σ::Vector{Float64}) -> Vector{Float64}
```

Calculate Gaussian probabilities on a grid.

**Arguments:**
- `y`: Observed values
- `μ`: Mean values
- `σ`: Standard deviations

**Returns:**
- Vector of probabilities

**Example:**
```julia
# Calculate probabilities on grid
y_vals = [1.0, 2.0, 3.0]
mu_vals = [1.1, 2.1, 2.9]
sigma_vals = [0.2, 0.3, 0.4]
probs = prob_Gaussian_grid(y_vals, mu_vals, sigma_vals)
```

### mean_elongationtime

```julia
mean_elongationtime(rates::Vector{Float64}, R::Int) -> Float64
```

Calculate mean elongation time for pre-RNA steps.

**Arguments:**
- `rates`: Rate parameters
- `R`: Number of pre-RNA steps

**Returns:**
- Mean elongation time

**Example:**
```julia
# Calculate elongation time
rates = [0.1, 0.2, 0.3, 0.4]
elong_time = mean_elongationtime(rates, 2)
```

### on_states

```julia
on_states(G::Int, R::Int, S::Int, insertstep::Int) -> Vector{Int}
```

Identify transcriptionally active states.

**Arguments:**
- `G`: Number of gene states
- `R`: Number of pre-RNA steps
- `S`: Number of splice sites
- `insertstep`: Reporter insertion step

**Returns:**
- Vector of active state indices

**Example:**
```julia
# Find active states
active = on_states(2, 3, 1, 1)
```

### source_states

```julia
source_states(transitions::Tuple) -> Vector{Int}
```

Identify source states for transitions.

**Arguments:**
- `transitions`: Model transitions

**Returns:**
- Vector of source state indices

**Example:**
```julia
# Find source states
sources = source_states(([1,2], [2,1]))
```

## File I/O Functions

### folder_path

```julia
folder_path(folder::String, root::String=".", subfolder::String=""; make::Bool=false) -> String
```

Construct folder paths with optional creation.

**Arguments:**
- `folder`: Folder name
- `root`: Root directory
- `subfolder`: Subfolder name
- `make`: Whether to create the folder

**Returns:**
- Full folder path

**Example:**
```julia
# Create results folder path
results_path = folder_path("results", ".", "analysis", make=true)
```

### folder_setup

```julia
folder_setup(base_path::String, folders::Vector{String}) -> Nothing
```

Set up directory structure for analysis.

**Arguments:**
- `base_path`: Base directory
- `folders`: Vector of folder names to create

**Example:**
```julia
# Set up analysis folders
folder_setup("project", ["data", "results", "figures"])
```

### datapdf

```julia
datapdf(data::AbstractExperimentalData, filename::String) -> Nothing
```

Generate PDF summary of data.

**Arguments:**
- `data`: Data structure
- `filename`: Output filename

**Example:**
```julia
# Generate data summary PDF
datapdf(rna_data, "data_summary.pdf")
```

### make_dataframes

```julia
make_dataframes(results, model::AbstractModel) -> DataFrame
```

Create DataFrames from analysis results.

**Arguments:**
- `results`: Analysis results
- `model`: Model structure

**Returns:**
- DataFrame with results

**Example:**
```julia
# Convert results to DataFrame
df = make_dataframes(fit_results, model)
```

### write_dataframes_only

```julia
write_dataframes_only(dataframes::Vector{DataFrame}, folder::String, label::String) -> Nothing
```

Write DataFrames to CSV files.

**Arguments:**
- `dataframes`: Vector of DataFrames
- `folder`: Output folder
- `label`: File label

**Example:**
```julia
# Write results to CSV
write_dataframes_only([df1, df2], "results", "analysis")
```

## Advanced Utility Functions

### zero_median

```julia
zero_median(traces::Vector{Vector{Float64}}, zero::Bool) -> Tuple{Vector{Vector{Float64}}, Vector{Float64}}
```

Zero-center traces by subtracting median values.

**Arguments:**
- `traces`: Vector of trace vectors
- `zero`: Whether to zero-center

**Returns:**
- Tuple of (zero-centered traces, scale factors)

**Example:**
```julia
# Zero-center traces
centered_traces, scales = zero_median(traces, true)
```

### fix

```julia
fix(parameters::Vector{Float64}, fixedeffects::Tuple) -> Vector{Float64}
```

Apply fixed effects to parameters.

**Arguments:**
- `parameters`: Parameter vector
- `fixedeffects`: Fixed effects specification

**Returns:**
- Parameters with fixed effects applied

**Example:**
```julia
# Apply fixed effects
fixed_params = fix(parameters, fixedeffects)
```

### fix_filenames

```julia
fix_filenames(folder::String, pattern::String, replacement::String) -> Nothing
```

Fix filenames in a folder by pattern replacement.

**Arguments:**
- `folder`: Folder path
- `pattern`: Pattern to replace
- `replacement`: Replacement string

**Example:**
```julia
# Fix filenames
fix_filenames("data", "old_prefix", "new_prefix")
```

## Model Analysis Functions

### large_deviance

```julia
large_deviance(chains::Vector, threshold::Float64=1000.0) -> Vector{Int}
```

Identify chains with unusually large deviance.

**Arguments:**
- `chains`: MCMC chains
- `threshold`: Deviance threshold

**Returns:**
- Indices of chains with large deviance

**Example:**
```julia
# Find problematic chains
bad_chains = large_deviance(mcmc_chains, 500.0)
```

### large_rhat

```julia
large_rhat(stats, threshold::Float64=1.1) -> Vector{Int}
```

Identify parameters with large R-hat values.

**Arguments:**
- `stats`: MCMC statistics
- `threshold`: R-hat threshold

**Returns:**
- Indices of parameters with large R-hat

**Example:**
```julia
# Find non-converged parameters
unconverged = large_rhat(mcmc_stats, 1.05)
```

### assemble_measures_model

```julia
assemble_measures_model(measures::Vector, model::AbstractModel) -> Dict
```

Assemble model measures for analysis.

**Arguments:**
- `measures`: Vector of measure objects
- `model`: Model structure

**Returns:**
- Dictionary of assembled measures

**Example:**
```julia
# Assemble measures
assembled = assemble_measures_model(measures, model)
```

### assemble_all

```julia
assemble_all(fits, stats, measures, model::AbstractModel) -> Dict
```

Assemble all analysis results.

**Arguments:**
- `fits`: Fit results
- `stats`: Statistics
- `measures`: Measures
- `model`: Model structure

**Returns:**
- Dictionary of all results

**Example:**
```julia
# Assemble all results
all_results = assemble_all(fits, stats, measures, model)
```

## Performance Utilities

### set_indices

```julia
set_indices(model::AbstractModel) -> Vector{Int}
```

Set parameter indices for efficient computation.

**Arguments:**
- `model`: Model structure

**Returns:**
- Vector of parameter indices

### T_dimension

```julia
T_dimension(model::AbstractModel) -> Int
```

Calculate transition matrix dimension.

**Arguments:**
- `model`: Model structure

**Returns:**
- Matrix dimension

### sparse

```julia
sparse(data::Matrix) -> SparseMatrixCSC
```

Convert matrix to sparse format for efficiency.

**Arguments:**
- `data`: Dense matrix

**Returns:**
- Sparse matrix

## Usage Examples

### Complete Analysis Pipeline

```julia
using StochasticGene

# Load and process data
data = load_data("rna", String[], "data/", "exp1", "MYC", "ctrl", (), 1)
normalized = normalize_histogram(data.histRNA)

# Set up model
rates = prepare_rates([0.1, 0.2], model)
n_params = num_all_parameters(transitions, 0, 0, 1, reporter, (), nothing)

# Analyze results
active_states = on_states(2, 0, 0, 1)
elong_time = mean_elongationtime(rates, 0)

# Output results
folder_setup("analysis", ["results", "figures"])
results_path = folder_path("results", "analysis", make=true)
```

### Parameter Extraction

```julia
# Extract different parameter types
all_rates = get_rates(fitted_params, model, true)
noise_params = get_param(fitted_params, model, "noise")
coupling_params = get_param(fitted_params, model, "coupling")

# Count parameters
n_rates = num_rates(transitions, 2, 1, 1)
n_total = num_all_parameters(transitions, 2, 1, 1, reporter, coupling, grid)
```

### Data Processing Pipeline

```julia
# Process trace data
traces, scales = zero_median(raw_traces, true)
trace_matrix = make_mat(traces)
normalized_traces = normalize_histogram.(traces)

# Statistical analysis
probs = [prob_Gaussian(t, μ, σ) for t in trace_values]
active = on_states(G, R, S, insertstep)
```

## See Also

- [`fit`](@ref): Main fitting function
- [`simulator`](@ref): Model simulation
- [`load_data`](@ref): Data loading
- [`write_traces`](@ref): Trace generation