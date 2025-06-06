# API Reference

## Core Functions

### Model Fitting
- [`fit`](@ref): Main function for fitting models to data
- [`fit_parallel`](@ref): Parallel version of fit for large datasets

### Data Analysis
- [`analyze_results`](@ref): Analyze fitting results
- [`plot_results`](@ref): Plot fitting results
- [`save_results`](@ref): Save results to disk

### Data Types

#### RNAData
```julia
struct RNAData
    histRNA::Vector{Float64}  # Histogram of RNA counts
    gene::String             # Gene name
    condition::String        # Experimental condition
    nalleles::Int           # Number of alleles
end
```

A concrete struct for storing RNA histogram data. This is used for steady-state RNA count distributions from techniques like smFISH or scRNA-seq.

Example:
```julia
data = RNAData(
    histRNA=[1,2,3,4,5],
    gene="GENE1",
    condition="control",
    nalleles=2
)
```

- [`TraceData`](@ref): Live cell imaging data structure
- [`DwellTimeData`](@ref): Dwell time data structure

### Model Types
- [`GMmodel`](@ref): Gene model structure
- [`Transformation`](@ref): Data transformation structure

### Utility Functions
- [`rna_setup`](@ref): Set up project directory structure
- [`write_traces`](@ref): Generate model-predicted intensity traces
- [`write_ONOFFhistograms`](@ref): Generate ON/OFF dwell time histograms
- [`write_residency_G_folder`](@ref): Generate G state residency probabilities

## Model Components

### Gene States (G)
- Arbitrary number of gene states
- User-specified transitions between states
- One active state for transcription initiation

### Pre-RNA Steps (R)
- Irreversible forward transitions
- mRNA ejection from final R step
- Optional reporter insertion step

### Splicing (S)
- Up to R splice sites
- PreRNA with or without spliced intron
- Multiple configurations per R step

## Data Types

The package can handle:
- mRNA count distributions (smFISH, scRNA)
- Image intensity traces (live cell imaging)
- Dwell time distributions
- Combined data types
