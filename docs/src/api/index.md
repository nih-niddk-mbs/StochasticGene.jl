# API Reference

## Core Functions

### Model Fitting
- [`fit`](@ref): Main function for fitting models to data
- [`fit_parallel`](@ref): Parallel version of fit for large datasets
- [`fit_hierarchical`](@ref): Hierarchical model fitting

### Data Analysis
- [`analyze_results`](@ref): Analyze fitting results
- [`plot_results`](@ref): Plot fitting results
- [`save_results`](@ref): Save results to disk

### Data Types
- [`RNAData`](@ref): RNA count data structure
- [`TraceData`](@ref): Live cell imaging data structure
- [`DwellTimeData`](@ref): Dwell time data structure

### Model Types
- [`HMMReporter`](@ref): Hidden Markov model for reporter data
- [`GMmodel`](@ref): Gene model structure
- [`Transformation`](@ref): Data transformation structure

### Traits
- [`HierarchicalTrait`](@ref): Hierarchical model trait
- [`GridTrait`](@ref): Grid search trait
- [`CouplingTrait`](@ref): Coupled model trait

## Utility Functions

### Data Management
- [`rna_setup`](@ref): Set up project directory structure
- [`load_data`](@ref): Load data from files
- [`save_data`](@ref): Save data to files

### Model Management
- [`create_model`](@ref): Create a new model instance
- [`register_model`](@ref): Register a model type
- [`print_model`](@ref): Print model details

### Analysis Tools
- [`calculate_statistics`](@ref): Calculate basic statistics
- [`generate_plots`](@ref): Generate standard plots
- [`export_results`](@ref): Export results to various formats

## Advanced Features

### Parallel Processing
- [`setup_parallel`](@ref): Set up parallel processing
- [`cleanup_parallel`](@ref): Clean up parallel resources

### Bayesian Analysis
- [`bayesian_fit`](@ref): Bayesian model fitting
- [`mcmc_analysis`](@ref): MCMC analysis tools

### Model Validation
- [`validate_model`](@ref): Validate model assumptions
- [`cross_validate`](@ref): Cross-validation tools
