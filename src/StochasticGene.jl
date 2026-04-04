# StochasticGene.jl
# Author: Carson C. Chow

"""
    module StochasticGene

A Julia package for stochastic modeling of gene transcription and Bayesian inference.

**Version:** 1.8.3

# Overview
StochasticGene.jl is a comprehensive Julia package for simulating and fitting stochastic models of gene transcription. It provides tools for:

1. **Model Simulation**
   - Generalized telegraph models (GRSM)
   - Multiple gene states (G)
   - Pre-RNA steps (R)
   - Splice sites (S)
   - Reporter insertion
   - Multiple alleles
   - Coupled gene models

2. **Parameter Inference**
   - Bayesian parameter estimation via MCMC
   - Maximum likelihood estimation
   - Hierarchical modeling
   - Parallel processing support
   - Adaptive proposal distributions

3. **Data Types**
   - mRNA count distributions (smFISH, scRNA-seq)
   - Live cell imaging traces
   - Dwell time distributions
   - Combined data types
   - Time-series data

4. **Analysis Tools**
   - Model fitting and comparison
   - Hidden Markov model analysis
   - Burst size analysis
   - ON/OFF state analysis
   - Model diagnostics
   - Posterior analysis

# Key Features

## Model Flexibility
- **Gene States (G)**: Arbitrary number of states
- **Pre-RNA Steps (R)**: Multiple elongation steps
- **Splice Sites (S)**: Support for splicing dynamics
- **Coupled Models**: Multiple interacting alleles/genes; multi-connection coupling support
- **Hierarchical Models**: Share information across conditions

## Inference Methods
- **MCMC Sampling**: Robust parameter estimation
- **Adaptive Proposals**: Efficient exploration of parameter space
- **Parallel Chains**: Multi-core and distributed computing
- **Convergence Diagnostics**: Ensure reliable results

## Performance
- **Optimized Backend**: For fast simulations
- **Automatic Differentiation**: For gradient-based methods
- **Memory Efficiency**: Handles large datasets

# Quick Start

```julia
using StochasticGene

# Set up a simple two-state model
fits = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/example_data/",
    gene = "MYC"
)
```

# Module Structure

- `common.jl`: Type system and common functions
- `transition_rate_structures.jl`, `transition_rate_elements.jl`, `transition_rate_functions.jl`, `transition_rate_make.jl`: Transition rate matrices and master equations
- `likelihoods.jl`: Likelihood functions for fitting
- `metropolis_hastings.jl`: MCMC parameter estimation
- `io.jl`: Input/output operations (including combined-rate file helpers: `create_combined_file`, `create_combined_files`, …)
- `chemical_master.jl`: Chemical master equation solutions
- `utilities.jl`: Common utility functions
- `simulator_coupled.jl`: Stochastic simulation algorithms
- `fit.jl`: Model fitting functions
- `analysis.jl`: Post-fit analysis tools
- `biowulf.jl`: NIH Biowulf cluster support (`makeswarm`, `makeswarmfiles`, `makeswarm_models`, …)
- `coupled_csv.jl`: `Coupled_models_to_test.csv` → coupling specs (`csv_row_to_connections_simple`, `build_coupled_fit_spec_from_csv_cells`, …)
- `hmm.jl`: Hidden Markov model functions
- `test.jl`: Test functions

Batch workflows (swarms, combined starts, key-based runs) are documented in the manual: **[Cluster and batch workflows](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/cluster_batch_workflows.html)**. Repository layout, `results/` conventions, and nomenclature are in **[Package overview](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/package_overview.html)**.

# Dependencies

- `CSV`, `DataFrames`: Data handling and manipulation
- `Distributed`: Parallel processing
- `Distributions`, `StatsBase`: Statistical distributions and functions
- `LSODA`, `OrdinaryDiffEq`, `SciMLSensitivity`: ODE solvers and adjoints for Zygote through `solve`
- `MultivariateStats`: Statistical analysis
- `Optim`: Optimization algorithms
- `FFTW`: FFT operations
- `JLD2`, `JSON`, `TOML`: Serialization and configuration
- `SparseArrays`, `LinearAlgebra`: Numerical linear algebra

# Documentation

Hosted manual: [stable](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/) · [dev](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/). Stable docs track the released version (v1.8.3).
"""
module StochasticGene
# __precompile__(true)

using AdvancedHMC
using CSV
using ChainRulesCore
using DataFrames
using Dates
using DelimitedFiles
using Distributed
using Distributions
using Downloads
using FFTW
using JLD2
using JSON
using LinearAlgebra
using LogDensityProblems
using LSODA
using MultivariateStats
using Optim
using OrdinaryDiffEq
using SciMLSensitivity  # reverse-mode AD through `solve` (Zygote + trace HMM / `kolmogorov_forward`)
using Printf
using SparseArrays
using Statistics
using StatsBase
using TOML
using Zygote
using Random
# using CUDA

export
    assemble_all,
    assemble_measures_model,
    CSV,
    DataFrame,
    datapdf,
    digit_vector,
    fix,
    fix_filenames,
    find_best_models,
    fit,
    fit_default_spec,
    fit_coupled_default_spec,
    maxtime_seconds,
    INFERENCE_ADVI,
    INFERENCE_CHOICES,
    INFERENCE_MH,
    INFERENCE_NUTS,
    normalize_trace_specs_legacy_t_end!,
    AbstractPriorContext,
    PriorContextCoupled,
    PriorContextTraceSingleUnit,
    PriorContextGenericSingle,
    prior_context,
    set_priormean_empty,
    prior_ratemean_trace,
    trace_prior_variables,
    folder_path,
    folder_setup,
    get_param,
    get_rates,
    get_rates_ad,
    fixed_rates_ad,
    GMmodel,
    GRSMmodel,
    AbstractGRSMmodel,
    AbstractGeneTransitionModel,
    AbstractGMmodel,
    large_deviance,
    large_rhat,
    load_data,
    load_model,
    loglikelihood,
    loglikelihood_ad,
    HMM_STACK_MH,
    HMM_STACK_AD,
    set_hmm_zygote_checkpoint_steps!,
    hmm_zygote_checkpoint_steps,
    lsoda,
    make_array,
    make_dataframes,
    make_mat,
    makeswarm,
    makeswarm_models,
    csv_row_to_connections_simple,
    build_coupled_fit_spec_from_csv_cells,
    parse_coupling_sign_csv,
    replace_csv_cell_legacy_r,
    csv_row_has_legacy_r,
    makeswarm_coupled,
    makeswarm_genes,
    write_fitfile_coupled,
    make_traces,
    make_traces_dataframe,
    mean,
    mean_elongationtime,
    metropolis_hastings,
    GenePosteriorLogDensity,
    NUTSOptions,
    ADVIOptions,
    run_nuts,
    longest_trace_timesteps,
    check_ad_gradient_feasibility,
    run_nuts_fit,
    run_NUTS,
    run_advi,
    logposterior,
    norm,
    normalize_histogram,
    num_all_parameters,
    n_observed_trace_units,
    num_rates,
    on_states,
    predictedarray,
    predictedfn,
    prepare_rates,
    prepare_rates_ad,
    plot_empirical_vs_theory,
    prob_Gaussian,
    prob_Gaussian_grid,
    CorrelationAlgorithm,
    CorrelationTrait,
    StandardCorrelation,
    WindowedCorrelation,
    MultiTauCorrelation,
    IDLCorrelation,
    DEFAULT_CORRELATION_ALGORITHM,
    hastrait,
    readfile,
    RNAData,
    score_models_from_traces,
    RNADwellTimeData,
    RNAOnOffData,
    readrates,
    read_run_spec,
    read_run_spec_for_rates_file,
    info_toml_path_for_rates_file,
    read_rates_folder,
    rate_distribution_summary,
    rna_setup,
    run_mh,
    run_mcmc_parallel,
    set_indices,
    set_elements_TCoupledUnit,
    simulate_trace,
    simulate_trace_data,
    simulate_trace_vector,
    simulate_trials,
    get_3unit_model_params,
    simulator,
    simulator_dwell_specs,
    source_states,
    sparse,
    summarize_model_scores,
    joint_residence_prob_onoff,
    write_joint_residence_prob_onoff,
    test_compare,
    test_compare_3unit,
    test_compare_coupling,
    test_num_reporters_consistency,
    diagnose_sim_vs_cme,
    test_fit_rna,
    test_fit_rnadwelltime,
    test_load_model_keyword_compatibility,
    test_fit_rnaonoff,
    test_fit_simrna,
    test_compare_mh_nuts_posterior,
    test_fit_simrna_mh_nuts,
    test_simulate_trials,
    test_fit_trace,
    test_fit_trace_hierarchical,
    test_fit_tracejoint,
    test_fit_tracejoint_hierarchical,
    test_ad_gradient_smoke,
    test_get_rates_ad_consistency,
    test_run_nuts_fit_smoke,
    test_normalized_nullspace_augmented_pullback_fd,
    test_trace_specs_utilities,
    test_trace_subset_benchmark_keyword_bundle,
    test_benchmark_trace_joint_fit_stacks,
    benchmark_inference_simrna_small,
    benchmark_inference_trace_gr2r2,
    benchmark_inference_trace_coupled_3x3,
    benchmark_scenario_coupled_3x3_10traces_220frames,
    benchmark_inference_trace_coupled_3x3_g3r0,
    benchmark_inference_ensure_workers,
    benchmark_inference_setup_parallel_workers,
    benchmark_inference_run_mh,
    benchmark_inference_run_nuts_parallel,
    benchmark_trace_zygote_gradient,
    benchmark_trace_forwarddiff_gradient,
    benchmark_trace_finitediff_gradient,
    benchmark_trace_compare_forwarddiff_vs_finitediff,
    benchmark_trace_zygote_subset_gradient,
    compare_trace_subset_gradient_benchmarks,
    benchmark_inference_run_advi,
    benchmark_inference_print_summary,
    makeswarmfiles_driver,
    makeswarmfiles,
    makeswarmfiles_coupled,
    makeswarmfiles_h3_latent,
    parse_h3_transition_key,
    COUPLING_MODE_RECIPROCAL_DEFAULT,
    coupling_ranges,
    default_coupling_prior_mean,
    default_coupling_prior_means,
    make_coupling_hidden_latent,
    default_trace_specs_for_coupled,
    coupling_parameter_labels,
    TComponents,
    T_dimension,
    TraceData,
    TraceRNAData,
    Tsit5,
    unbiased_crosscov,
    correlation_function,
    correlation_function_windowed,
    correlation_function_multitau,
    write_augmented,
    write_correlation_functions,
    write_joint_residence_prob_onoff_key,
    write_correlation_functions_file,
    write_correlation_functions_empirical,
    correlation_functions_uncertainty_mc,
    correlation_functions_uncertainty_delta,
    compute_correlation_functions,
    write_dataframes,
    write_dataframes_only,
    write_histograms,
    write_ONOFFhistograms,
    write_residency_G_folder,
    write_RNAhistogram,
    write_traces,
    write_traces_coupling,
    write_traces_coupling_spawn,
    write_winners,
    create_combined_file,
    create_combined_file_mult,
    read_rates_table,
    write_rates_table,
    merge_coupled_two_unit_rates,
    merge_coupled_stacked_units,
    combined_rates_key,
    read_combined_file_specs_csv,
    create_combined_files_driver,
    create_combined_files,
    create_combined_files_h3_latent,
    write_run_spec_preset,
    write_fitscript_tracejoint_key,
    default_model_key,
    zero_median

    

### Source files

# Type system and common functions
include("common.jl")

# Transition rate matrices of stochastic models defining master equations
include("transition_rate_structures.jl")
include("transition_rate_elements.jl")
include("transition_rate_functions.jl")
include("transition_rate_make.jl")

# Likelihood functions for fitting algorithms
include("likelihoods.jl")

# Metropolis Hastings MCMC for computing posterior distributions of model parameters
include("metropolis_hastings.jl")

# NUTS / ADVI (gradient-based inference on transformed parameters)
include("gradient_inference.jl")

# Input output functions
include("io.jl")

# Chemical master equation solutions of stochastic models for likelihood functions in fitting algorithms
include("chemical_master.jl")

# commonly used functions
include("utilities.jl")

# Probability distributions by direct simulation of stochastic models using Gillespie and Gibson-Bruck algorithms
# include("simulator.jl")
include("simulator_coupled.jl")

# functions for fitting models to data
include("fit.jl")

# functions for post fit analysis and plots
include("analysis.jl")

# functions for hidden markov models (before coupled_csv: prob_Gaussian, etc.)
include("hmm.jl")

# Coupled_models_to_test.csv → coupling (used by biowulf batch helpers)
include("coupled_csv.jl")

# functions for use on NIH cluster Biowulf
include("biowulf.jl")

# test functions
include("test.jl")

end # module StochasticGene
