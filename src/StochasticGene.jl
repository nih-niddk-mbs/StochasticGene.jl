# StochasticGene.jl
# Author: Carson C. Chow

"""
    module StochasticGene

A Julia package for stochastic modeling of gene transcription and Bayesian inference.

**Version:** 1.11.0

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
- `stage.jl`: results staging (`stage_label_to_key`, legacy → key-based filenames; `stage_combine_rates` merges rate files)
- `coupled_csv.jl`: `Coupled_models_to_test.csv` → coupling specs (`csv_row_to_connections_simple`, `build_coupled_fit_spec_from_csv_cells`, …)
- `hmm.jl`: Hidden Markov model functions
- `test.jl`, `test_features.jl`, `benchmarks.jl`: full-stack tests, feature/regression checks, and benchmarking helpers

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

Hosted manual: [stable](https://nih-niddk-mbs.github.io/StochasticGene.jl/stable/) · [dev](https://nih-niddk-mbs.github.io/StochasticGene.jl/dev/). Stable docs track the released version (v1.11.0).
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
using FiniteDiff
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
using Profile
using SparseArrays
using Statistics
using StatsBase
using TOML
using Zygote
using Random
# using CUDA

# Dependencies intentionally re-exported for legacy scripts.
export CSV, DataFrame, Tsit5, lsoda, mean, norm, sparse

# Fit and inference entry points.
export fit, fit_default_spec, fit_coupled_default_spec, run_inference, load_options
export metropolis_hastings, run_mh, run_nuts_fit, logposterior
export GenePosteriorLogDensity, NUTSOptions, ADVIOptions
export INFERENCE_ADVI, INFERENCE_CHOICES, INFERENCE_MH, INFERENCE_NUTS
export maxtime_seconds

# Data and model types.
export RNAData, RNADwellTimeData, RNAOnOffData, TraceData, TraceRNAData
export AbstractCombinedData, CombinedData, MultiDatasetData
export AbstractGMmodel, AbstractGRSMmodel, AbstractGeneTransitionModel, GMmodel, GRSMmodel
export DatasetSpec, HierarchyNode, HierarchySpec, HierarchyAssignment
export SharedParameterCachePlan, SharedParameterTrait
export RecursiveHierarchyCachePlan, RecursiveHierarchyTrait
export PriorContextCoupled, PriorContextTraceSingleUnit, PriorContextGenericSingle
export AbstractPriorContext
export TComponents, TransientMasterProblem

# Data loading, model setup, and likelihood bridges.
export load_data, load_model, loglikelihood, loglikelihood_ad
export get_param, get_rates, get_rates_ad, fixed_rates_ad
export prepare_rates, prepare_rates_ad, num_rates, num_all_parameters
export set_indices, set_elements_TCoupledUnit, source_states, on_states
export prior_context, set_priormean_empty, prior_ratemean_trace, trace_prior_variables
export normalize_trace_specs_legacy_t_end!, mean_elongationtime, n_observed_trace_units
export combined_modalities, reconstruct_tracerna
export HMM_STACK_MH, HMM_STACK_AD, set_hmm_zygote_checkpoint_steps!, hmm_zygote_checkpoint_steps

# Shared and recursive hierarchy helpers.
export hierarchy_assignments, hierarchy_node_names, hierarchy_parameter_groups
export hierarchy_path, hierarchy_parameter_levels, n_emission_groups, n_transition_rate_groups
export compile_shared_parameter_trait, shared_parameter_cache_plan
export shared_assignment_rates, shared_transition_group_rates, shared_emission_group_rates
export compile_recursive_hierarchy_trait, recursive_hierarchy_cache_plan
export recursive_assignment_rates, recursive_transition_group_rates, recursive_emission_group_rates

# Probability, simulation, and transient solvers.
export prob_Gaussian, prob_Gaussian_grid, predictedarray, predictedfn
export simulate_trace, simulate_trace_data, simulate_trace_vector, simulate_trials
export get_3unit_model_params, simulator, simulator_dwell_specs
export transient_master_problem, transient_master_initial, transient_master_strang
export transient_master_marginal, transient_A_flow!, transient_B_flow!
export normalize_histogram, datapdf

# Correlation and burst-prediction analysis.
export CorrelationAlgorithm, CorrelationTrait, StandardCorrelation, WindowedCorrelation
export MultiTauCorrelation, IDLCorrelation, DEFAULT_CORRELATION_ALGORITHM, hastrait
export unbiased_crosscov, correlation_function, correlation_function_windowed
export correlation_function_multitau, build_correlation_context, correlation_observable
export correlation_observable_label, correlate_observables, compute_correlation_functions
export posterior_state_probabilities, posterior_burst_probability, auc_binary_score
export burst_prediction_auc, simulate_burst_observability_auc
export simulate_fit_burst_identifiability_auc

# Analysis output writers and dataframe assembly.
export assemble_all, assemble_all_key, assemble_measures_key, assemble_measures_model
export assemble_rates_key, assemble_stats_key
export make_array, make_dataframes, make_dataframes_key, make_mat, make_traces
export make_traces_dataframe, plot_empirical_vs_theory, score_models_from_traces
export summarize_model_scores, joint_residence_prob_onoff
export write_augmented, write_correlation_functions, write_correlation_functions_key
export write_correlation_general_csv, write_correlation_functions_file
export write_correlation_functions_empirical, correlation_functions_uncertainty_mc
export correlation_functions_uncertainty_delta, write_joint_residence_prob_onoff
export write_joint_residence_prob_onoff_key, write_shared_group_measures_from_samples
export write_dataframes, write_dataframes_only, write_dataframes_only_key
export write_histograms, write_ONOFFhistograms, write_ONOFFhistograms_key
export write_residency_G_folder, write_RNAhistogram, write_traces, write_traces_coupling
export write_traces_coupling_spawn, write_winners

# File IO, run specs, and staged batch generation.
export create_label, folder_path, folder_setup, fix, fix_filenames, readfile, readrates
export read_run_spec, read_run_spec_for_rates_file, info_jld2_path
export info_toml_path_for_rates_file, read_rates_folder, rate_distribution_summary
export rna_setup, make_fitscript, write_fitfile_full, build_julia_script_command
export write_julia_command_file, make_commandfile, make_commandfile_from_csv
export make_swarmfile_from_csv, make_fitscripts_from_csv
export make_fitscripts_and_commandfile_from_csv, make_fitscripts_and_swarm_from_csv
export stage_label_to_key, label_to_key, staging_key_segment, stage_write_run_specs
export write_run_spec_preset, write_fitscript_tracejoint_key, default_model_key

# Biowulf and CSV-driven coupled helpers.
export makeswarm, makeswarm_models, makeswarm_coupled, makeswarm_genes
export makeswarmfiles_driver, makeswarmfiles, makeswarmfiles_coupled, makeswarmfiles_h3_latent
export write_fitfile_coupled, csv_row_to_connections_simple
export build_coupled_fit_spec_from_csv_cells, parse_coupling_sign_csv
export replace_csv_cell_legacy_r, csv_row_has_legacy_r
export stage_combine_rates, stage_combine_rates_specs_from_csv, stage_combine_rates_from_csv

# Coupled-model staging and rate-table utilities.
export COUPLING_MODE_RECIPROCAL_DEFAULT, coupling_ranges, coupling_parameter_labels
export default_coupling_prior_mean, default_coupling_prior_means, default_trace_specs_for_coupled
export make_coupling_hidden_latent, parse_h3_transition_key
export create_combined_file, create_combined_file_mult, create_combined_files_driver
export create_combined_files, create_combined_files_h3_latent
export read_rates_table, write_rates_table, merge_coupled_two_unit_rates
export merge_coupled_stacked_units, combined_rates_key, read_combined_file_specs_csv

# Miscellaneous helpers retained for compatibility.
export digit_vector, find_best_models, large_deviance, large_rhat
export longest_trace_timesteps, check_ad_gradient_feasibility, T_dimension, zero_median

# Test and benchmark helpers exported for compatibility with existing scripts.
export test_compare, test_compare_3unit, test_compare_coupling, test_num_reporters_consistency
export diagnose_sim_vs_cme, test_fit_rna, test_fit_rnadwelltime
export test_load_model_keyword_compatibility, test_fit_rnaonoff, test_fit_simrna
export test_compare_mh_nuts_posterior, test_fit_simrna_mh_nuts, test_simulate_trials
export test_fit_trace, test_fit_trace_hierarchical, test_fit_tracejoint
export test_fit_tracejoint_hierarchical, test_ad_gradient_smoke, test_get_rates_ad_consistency
export test_run_nuts_fit_smoke, test_run_advi_fit_smoke, test_run_advi_trace_smoke
export test_nuts_vs_advi_consistency, test_normalized_nullspace_augmented_pullback_fd
export test_trace_specs_utilities, test_trace_subset_benchmark_keyword_bundle
export test_benchmark_trace_joint_fit_stacks
export benchmark_inference_simrna_small, benchmark_inference_trace_gr2r2
export benchmark_inference_trace_coupled_3x3, benchmark_scenario_coupled_3x3_10traces_220frames
export benchmark_inference_trace_coupled_3x3_g3r0, benchmark_inference_ensure_workers
export benchmark_inference_setup_parallel_workers, benchmark_inference_run_mh
export benchmark_combined_likelihood_stack, benchmark_inference_run_nuts_parallel
export benchmark_trace_zygote_gradient, benchmark_trace_forwarddiff_gradient
export benchmark_trace_finitediff_gradient, benchmark_trace_compare_forwarddiff_vs_finitediff
export benchmark_trace_zygote_subset_gradient, compare_trace_subset_gradient_benchmarks
export benchmark_inference_run_advi, benchmark_inference_compare_mh_nuts_advi
export benchmark_inference_print_summary

    

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

# Unified inference dispatcher (run_inference)
include("inference_common.jl")

# Input output functions
include("io.jl")

# Chemical master equation solutions of stochastic models for likelihood functions in fitting algorithms
include("chemical_master.jl")

# Transient RNA closure / operator-splitting solvers. Kept separate from the
# steady-state CME stack in chemical_master.jl.
include("transient_master.jl")

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

# staged run directories and companion file writers (batch / reproducibility)
include("stage.jl")

# test functions
include("test.jl")

end # module StochasticGene
