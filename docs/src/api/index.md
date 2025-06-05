# API Reference

## Types

```@docs
StochasticGene.AbstractModel
StochasticGene.HMMModel
StochasticGene.KalmanModel
StochasticGene.CoupledSystemModel
```

## Macros

```@docs
StochasticGene.@define_model
StochasticGene.@register_model
StochasticGene.@dispatch_by_type
```

## Functions

### Model Creation and Management

```@docs
StochasticGene.create_model
```

### Model Registry

```@docs
StochasticGene.MODEL_REGISTRY
```

## Core Functions

### `fit`

```julia
fit(; <keyword arguments>)
```

Fit steady state or transient GM/GRSM model to RNA data for a single gene, write the result (through function finalize), and return fit results and diagnostics.

#### Arguments

- `annealsteps=0`: Number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `burst=false`: If true then compute burst frequency
- `cell::String=""`: Cell type for halflives and allele numbers
- `coupling=tuple()`: If nonempty, a 4-tuple where elements are:
  1. Tuple of model indices corresponding to each unit
  2. Tuple of vectors indicating source units for each unit
  3. Source states: tuple of vectors of strings
  4. Target transitions: tuple
  5. Int indicating number of coupling parameters
- `datacol=3`: Column of data to use, default is 3 for rna data
- `datatype::String=""`: String that describes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracerna", "tracejoint", "tracegrid"
- `datacond=""`: String or vector of strings describing data
- `datapath=""`: Path to data file or folder or array of files or folders
- `decayrate=1.0`: Decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `ejectnumber=1`: Number of mRNAs produced per burst, default is 1
- `dttype=String[]`: Dwelltime types, choices are "OFF", "ON", for R states and "OFFG", "ONG" for G states
- `elongationtime=6.0`: Average time for elongation, vector of times for coupled model
- `fittedparam::Vector=Int[]`: Vector of rate indices to be fit
- `fixedeffects::String`: If "fixed" is included after a hyphen, then fixedeffects Tuple will be created such that R transitions are fixed to be identical
- `fixedeffects::Tuple=tuple()`: Tuple of vectors of rates that are fixed where first index is fit and others are fixed to first
- `gene::String="MYC"`: Gene name
- `grid=nothing`: Int number of grid points for grid model
- `G=2`: Number of gene states, for coupled models G, R, S, and insertstep are vectors
- `hierarchical=tuple()`: Empty tuple for nonhierarchical model; 3-tuple for hierarchical
- `infolder::String=""`: Result folder used for initial parameters
- `inlabel::String=""`: Label of files used for initial conditions
- `insertstep=1`: R step where reporter is inserted
- `label::String=""`: Label of output files produced
- `maxtime=60`: Maximum wall time for run (in minutes)
- `method=lsoda()`: DifferentialEquations.jl numerical method
- `nalleles=1`: Number of alleles
- `nchains::Int=2`: Number of MCMC chains = number of processors called by Julia
- `noisepriors=[]`: Priors of observation noise
- `onstates=Int[]`: Vector of on or sojourn states
- `optimize=false`: Use optimizer to compute maximum likelihood value
- `priormean=Float64[]`: Mean rates of prior distribution
- `priorcv=10.`: Coefficient of variation(s) for the rate prior distributions
- `probfn=prob_Gaussian`: Probability function for HMM observation probability
- `propcv=0.01`: Coefficient of variation of proposal distribution
- `resultfolder::String=test`: Folder for results of MCMC run
- `R=0`: Number of pre-RNA steps
- `root="."`: Name of root directory for project
- `samplesteps::Int=1000000`: Number of MCMC sampling steps
- `S=0`: Number of splice sites
- `splicetype=""`: RNA pathway for GRS models
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: Annealing temperature
- `temprna=1.`: Reduce RNA counts by temprna compared to dwell times
- `traceinfo=(1.0, 1., -1, 1., [100.,10.])`: 5-tuple of trace parameters
- `TransitionType=""`: String describing G transition type
- `transitions::Tuple=([1,2],[2,1])`: Tuple of vectors that specify state transitions for G states
- `warmupsteps=0`: Number of MCMC warmup steps
- `writesamples=false`: Write out MH samples if true
- `zeromedian=false`: If true, subtract the median of each trace from each trace, then scale by the maximum of the medians

#### Returns

- `fits`: MCMC fit results (posterior samples, log-likelihoods, etc.)
- `stats`: Summary statistics for parameters
- `measures`: Diagnostic measures (including WAIC and its standard error)
- `data`, `model`, `options`: The data, model, and options structures used

### `makeswarm`

```julia
makeswarm(; <keyword arguments>)
```

Write swarm and fit files used on Biowulf.

#### Arguments

- `nthreads=1`: Number of Julia threads per processor
- `swarmfile::String="fit"`: Name of swarmfile to be executed by swarm
- `juliafile::String="fitscript"`: Name of file to be called by julia in swarmfile
- `src=""`: Path to folder containing StochasticGene.jl/src
- All keyword arguments of function `fit`

### `makeswarm_genes`

```julia
makeswarm_genes(genes::Vector{String}; <keyword arguments>)
```

Write a swarmfile and fit files to run each gene in vector genes.

#### Arguments

- `genes`: Vector of genes
- `batchsize=1000`: Number of jobs per swarmfile
- All arguments in `makeswarm`

### `simulator`

```julia
simulator(r, transitions, G, R, S, insertstep; <keyword arguments>)
```

Simulate any GRSM model. Returns steady state mRNA histogram. If bins not a null vector will return a vector of the mRNA histogram and ON and OFF time histograms. If traceinterval > 0, it will return a vector containing the mRNA histogram and the traces.

#### Arguments

- `r`: Vector of rates
- `transitions`: Tuple of vectors that specify state transitions for G states
- `G`: Number of gene states
- `R`: Number of pre-RNA steps
- `S`: Number of splice sites
- `insertstep`: Reporter insertion step

#### Keyword Arguments

- `bins::Vector=Float64[]`: Vector of time bin vectors for each set of ON and OFF histograms
- `coupling=tuple()`: If nonempty, a 4-tuple for coupling configuration
- `nalleles`: Number of alleles
- `nhist::Int`: Size of mRNA histogram
- `onstates::Vector`: Vector of vector of ON states
- `probfn=prob_Gaussian`: Reporter distribution
- `reporterfn=sum`: How individual reporters are combined
- `splicetype::String`: Splice action
- `tol::Float64=1e-6`: Convergence error tolerance
- `totalsteps::Int=10000000`: Maximum number of simulation steps
- `totaltime::Float64=0.0`: Total time of simulation
- `traceinterval`: Interval in minutes between frames for intensity traces
- `verbose::Bool=false`: Flag for printing state information

### `write_dataframes`

```julia
write_dataframes(resultfolder::String, datapath::String; measure::Symbol=:AIC, assemble::Bool=true, fittedparams=Int[])
```

Collates run results into a csv file.

#### Arguments

- `resultfolder`: Name of folder with result files
- `datapath`: Name of folder where data is stored
- `measure`: Measure used to assess winner
- `assemble`: If true then assemble results into summary files

## Utility Functions

### `rna_setup`

```julia
rna_setup(root::String=".")
```

Creates the required directory structure for StochasticGene.

### `simulate_trace_data`

```julia
simulate_trace_data(datafolder::String; <keyword arguments>)
```

Creates simulated trace files in datafolder.

### `write_traces`

```julia
write_traces(resultfolder::String, datapath::String, label::String, interval::Float64)
```

Generates model predicted traces using the Viterbi algorithm.

### `write_ONOFFhistograms`

```julia
write_ONOFFhistograms(resultfolder::String)
```

Generates ON and OFF dwelltime histograms.

### `write_residency_G_folder`

```julia
write_residency_G_folder(resultfolder::String)
```

Generates G state residence time probabilities.

### `write_augmented`

```julia
write_augmented(summaryfile::String, resultfolder::String)
```

Supplements the Summary file with more information.

## Usage Examples

### Creating a Model

```julia
using StochasticGene

# Define a new model type
@define_model MyModel param1::Int, param2::Float64

# Register the model
@register_model MyModel

# Create an instance
model = create_model(:MyModel, 5, 2.0)
```

### Running a Model

```julia
# Create a model instance
model = HMMModel(5)
data = TimeSeriesData()

# Run the model
run_model(model, data)
```

### Data Handling

```@docs
StochasticGene.load_data
StochasticGene.save_data
StochasticGene.rna_setup
```

### Analysis

```@docs
StochasticGene.analyze
StochasticGene.plot
StochasticGene.summarize
```

### Utilities

```@docs
StochasticGene.version
StochasticGene.check_installation
``` 