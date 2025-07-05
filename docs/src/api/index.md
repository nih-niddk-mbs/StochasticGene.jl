# API Reference

## Core Functions

### Model Fitting
- [`fit`](fit.md): Main function for fitting models to data
- [`loglikelihood`](@ref): Calculate log-likelihood for model parameters
- [`run_mh`](@ref): Run Metropolis-Hastings MCMC sampling
- [`run_mcmc_parallel`](@ref): Run parallel MCMC chains

### Model Simulation
- [`simulator`](simulator.md): Simulate stochastic gene expression models
- [`simulate_trace`](@ref): Generate intensity traces from models
- [`simulate_trace_data`](@ref): Generate trace data with metadata
- [`simulate_trace_vector`](@ref): Generate vectors of traces
- [`simulate_trials`](@ref): Simulate multiple model realizations

### Data Loading and Management
- [`load_data`](load_data.md): Load experimental data from files
- [`load_model`](@ref): Load model parameters from files
- [`rna_setup`](@ref): Set up project directory structure
- [`readrates`](@ref): Read rate parameters from files
- [`readfile`](@ref): Read data files with error handling

### Analysis and Visualization
- [`write_traces`](write_traces.md): Generate model-predicted intensity traces
- [`write_ONOFFhistograms`](write_ONOFFhistograms.md): Generate ON/OFF dwell time histograms
- [`write_residency_G_folder`](write_residency_G_folder.md): Generate G state residency probabilities
- [`write_histograms`](@ref): Write RNA histogram predictions
- [`write_dataframes`](@ref): Write results to CSV files

## Comprehensive Function Documentation

### Core Function Libraries
- [`Utilities`](utilities.md): Data processing, model construction, and utility functions
- [`Analysis`](analysis.md): Post-fitting analysis, model comparison, and visualization functions

## Data Types

### RNA Data Structures

#### RNAData
```julia
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    label::String     # Data set label
    gene::String      # Gene name (case sensitive)
    nRNA::nType       # Length of histogram
    histRNA::hType    # RNA histogram data
end
```

A structure for storing steady-state RNA count distributions from techniques like smFISH or scRNA-seq.

**Example:**
```julia
# Create RNA data from histogram
data = RNAData(
    "control",          # label
    "MYC",             # gene
    50,                # nRNA
    [10,20,30,25,15]   # histRNA
)
```

#### RNACountData
```julia
struct RNACountData <: AbstractRNAData{Vector{Int}}
    label::String
    gene::String
    nRNA::Int
    countsRNA::Vector{Int}
    yieldfactor::Vector{Float64}
end
```

A structure for storing individual RNA count measurements with yield factors.

**Example:**
```julia
# Create RNA count data with yield correction
data = RNACountData(
    "single_cell",     # label
    "ACTB",           # gene
    100,              # nRNA
    [1,2,3,4,5],      # countsRNA
    [0.8,0.9,0.85]    # yieldfactor
)
```

#### RNAOnOffData
```julia
struct RNAOnOffData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Vector
    bins::Vector
    ON::Vector        # ON time probability density
    OFF::Vector       # OFF time probability density
end
```

A structure for storing combined RNA count and ON/OFF state duration data.

**Example:**
```julia
# Create combined RNA and ON/OFF data
data = RNAOnOffData(
    "live_cell",      # label
    "SOX2",          # gene
    30,              # nRNA
    [5,10,15,20],    # histRNA
    [0,1,2,3,4],     # bins
    [0.1,0.3,0.4,0.2], # ON
    [0.2,0.4,0.3,0.1]  # OFF
)
```

#### RNADwellTimeData
```julia
struct RNADwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Array
    bins::Vector{Vector}
    DwellTimes::Vector{Vector}
    DTtypes::Vector
end
```

A structure for storing RNA counts with dwell time distributions.

### Trace Data Structures

#### TraceData
```julia
struct TraceData{labelType,geneType,traceType} <: AbstractTraceData
    label::labelType    # Data set label
    gene::geneType      # Gene name
    interval::Float64   # Time interval between trace points
    trace::traceType    # Trace data
end
```

A structure for storing fluorescence intensity time series data.

**Example:**
```julia
# Create trace data
traces = [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]]
data = TraceData(
    "live_imaging",   # label
    "MYC",           # gene
    1.0,             # interval (minutes)
    traces           # trace data
)
```

#### TraceRNAData
```julia
struct TraceRNAData{traceType,hType} <: AbstractTraceHistogramData
    label::String
    gene::String
    interval::Float64
    trace::traceType
    nRNA::Int
    histRNA::hType
end
```

A structure for storing both trace and RNA histogram data.

#### DwellTimeData
```julia
struct DwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    bins::Vector
    DwellTimes::Vector
    DTtypes::Vector
end
```

A structure for storing dwell time distributions only.

## Model Types

### Gene Model (GM)

#### GMmodel
```julia
struct GMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGMmodel
    rates::RateType              # Transition rates
    Gtransitions::Tuple          # G state transitions
    G::Int                       # Number of G states
    nalleles::Int               # Number of alleles
    rateprior::PriorType        # Rate prior distributions
    proposal::ProposalType      # MCMC proposal distribution
    fittedparam::ParamType      # Fitted parameter indices
    fixedeffects::Tuple         # Fixed effects specification
    method::MethodType          # Solution method
    components::ComponentType   # Model components
    reporter::ReporterType      # Reporter configuration
end
```

A structure for Gene (G) models with arbitrary numbers of gene states.

**Example:**
```julia
# Create a simple two-state gene model
model = GMmodel(
    rates = [0.1, 0.2],                    # G1->G2, G2->G1
    Gtransitions = ([1,2], [2,1]),         # State transitions
    G = 2,                                 # Two gene states
    nalleles = 2,                          # Diploid
    # ... other parameters
)
```

### Gene-Reporter-Splice Model (GRSM)

#### GRSMmodel
```julia
struct GRSMmodel{TraitType,RateType,nratesType,GType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{TraitType}
    trait::TraitType            # Model traits (hierarchical, coupling, etc.)
    rates::RateType             # Transition rates
    transforms::Transformation  # Rate transformations
    nrates::nratesType          # Number of rates
    Gtransitions::Tuple         # G state transitions
    G::GType                    # Number of G states
    R::GType                    # Number of R (pre-RNA) steps
    S::GType                    # Number of S (splice) sites
    insertstep::GType           # Reporter insertion step
    nalleles::Int               # Number of alleles
    splicetype::String          # Splicing type
    rateprior::PriorType        # Rate prior distributions
    proposal::ProposalType      # MCMC proposal distribution
    fittedparam::ParamType      # Fitted parameter indices
    fixedeffects::Tuple         # Fixed effects specification
    method::MethodType          # Solution method
    components::ComponentType   # Model components
    reporter::ReporterType      # Reporter configuration
end
```

A comprehensive structure for Gene-Reporter-Splice models with arbitrary numbers of gene states (G), pre-RNA steps (R), and splice sites (S).

**Example:**
```julia
# Create a GRS model with 2 gene states, 3 pre-RNA steps, 2 splice sites
model = GRSMmodel(
    rates = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],  # G, R, S transitions
    Gtransitions = ([1,2], [2,1]),            # G state transitions
    G = 2,                                    # Two gene states
    R = 3,                                    # Three pre-RNA steps
    S = 2,                                    # Two splice sites
    insertstep = 1,                           # Reporter visible from step 1
    nalleles = 2,                             # Diploid
    splicetype = "offeject",                  # Splicing type
    # ... other parameters
)
```

### Supporting Types

#### HMMReporter
```julia
struct HMMReporter
    n::Int                      # Number of noise parameters
    per_state::Vector           # Reporters per state
    probfn::Function           # Noise distribution function
    weightind::Int             # Mixture weight index
    offstates::Vector{Int}     # Off states
    noiseparams::Vector{Int}   # Noise parameter indices
end
```

A structure for configuring reporter properties in Hidden Markov Models.

#### Transformation
```julia
struct Transformation
    f::Vector{Function}        # Forward transformations
    f_inv::Vector{Function}    # Inverse transformations
    f_cv::Vector{Function}     # CV transformations
end
```

A structure for parameter transformations during fitting.

## Utility Functions

### Data Processing
- [`normalize_histogram`](@ref): Normalize probability distributions
- [`make_array`](@ref): Convert data to arrays
- [`make_mat`](@ref): Convert data to matrices
- [`digit_vector`](@ref): Convert numbers to digit vectors

### Model Construction
- [`prepare_rates`](@ref): Prepare rate parameters for fitting
- [`get_rates`](@ref): Extract rates from fitted parameters
- [`get_param`](@ref): Extract specific parameters
- [`num_rates`](@ref): Count number of rates in model
- [`num_all_parameters`](@ref): Count total parameters

### Statistical Functions
- [`prob_Gaussian`](@ref): Gaussian probability density
- [`prob_Gaussian_grid`](@ref): Gaussian probability on grid
- [`mean_elongationtime`](@ref): Calculate mean elongation time
- [`on_states`](@ref): Identify transcriptionally active states
- [`source_states`](@ref): Identify source states for transitions

### File I/O
- [`folder_path`](@ref): Construct folder paths
- [`folder_setup`](@ref): Set up directory structure
- [`datapdf`](@ref): Generate data PDFs
- [`make_dataframes`](@ref): Create DataFrames from results
- [`write_dataframes_only`](@ref): Write DataFrames to files

## Analysis Functions

### Model Diagnostics
- [`large_deviance`](@ref): Identify chains with large deviance
- [`large_rhat`](@ref): Identify parameters with large R-hat
- [`assemble_measures_model`](@ref): Assemble model measures
- [`assemble_all`](@ref): Assemble all results

### Post-fitting Analysis
- [`predictedarray`](@ref): Generate model predictions
- [`predictedfn`](@ref): Generate prediction functions
- [`make_traces`](@ref): Generate trace predictions
- [`make_traces_dataframe`](@ref): Convert traces to DataFrames

### Visualization Support
- [`write_cov`](@ref): Write covariance matrices
- [`write_augmented`](@ref): Write augmented results
- [`write_winners`](@ref): Write best-fit results

## Hierarchical and Coupling Models

### Hierarchical Models
For hierarchical models, use the `hierarchical` parameter in `fit()`:

```julia
# Fit hierarchical model
fits = fit(
    hierarchical = (2, [1,2]),  # 2 hyperparameter sets, fit rates 1,2
    # ... other parameters
)
```

### Coupled Models
For coupled transcriptional units, use tuples for model parameters:

```julia
# Fit coupled model with different G states
fits = fit(
    G = (2, 3),                          # Unit 1: 2 states, Unit 2: 3 states
    transitions = (([1,2], [2,1]), ([1,2], [2,3], [3,1])),
    coupling = (1, 2, [1,2], [3,4], 2),  # Coupling specification
    # ... other parameters
)
```

## Model Components

### Gene States (G)
- Arbitrary number of gene states
- User-specified transitions between states
- One active state for transcription initiation
- Support for multiple alleles

### Pre-RNA Steps (R)
- Irreversible forward transitions
- mRNA ejection from final R step
- Optional reporter insertion at any step
- Support for elongation dynamics

### Splicing (S)
- Up to R splice sites
- PreRNA with or without spliced intron
- Multiple configurations per R step
- Support for different splicing types:
  - `""`: No splicing
  - `"offeject"`: Splice then eject
  - `"offdecay"`: Splice then decay

## Supported Data Types

The package handles multiple experimental data types:

1. **mRNA Count Distributions**
   - Single molecule FISH (smFISH)
   - Single cell RNA sequencing (scRNA-seq)
   - Bulk RNA measurements

2. **Intensity Traces**
   - Live cell fluorescence microscopy
   - Time-lapse imaging
   - Multiple reporter constructs

3. **Dwell Time Distributions**
   - ON/OFF state durations
   - Transcriptional burst analysis
   - Reporter visibility periods

4. **Combined Data Types**
   - RNA + trace data
   - RNA + dwell time data
   - Multi-modal experiments

## Error Handling

The package includes comprehensive error handling:

- **Data validation**: Checks for consistent data formats
- **Parameter validation**: Ensures valid model specifications
- **Convergence monitoring**: Tracks MCMC convergence
- **Memory management**: Handles large datasets efficiently

## Performance Considerations

- **Parallel processing**: Use `nchains > 1` for parallel MCMC
- **Memory usage**: Large datasets may require cluster computing
- **Convergence**: Monitor R-hat values for chain convergence
- **Optimization**: Use appropriate priors for faster convergence
