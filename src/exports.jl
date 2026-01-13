# Core Model Types
export GMmodel, GRSMmodel, AbstractGRSMmodel, AbstractGeneTransitionModel, AbstractGMmodel

# Model Creation and Analysis
export fit, fit_parallel, loglikelihood, predictedarray, predictedfn

# Model Simulation
export simulator, simulate_trace, simulate_trace_data, simulate_trace_vector, simulate_trials

# Data Types
export RNAData, RNADwellTimeData, RNAOnOffData, TraceData, TraceRNAData

# Data Handling
export load_data, load_model, write_histograms, write_traces, write_dataframes

# Analysis Functions
export mean, norm, normalize_histogram, mean_elongationtime, on_states, source_states

# Utility Functions
export make_array, make_mat, prepare_rates, prob_Gaussian, prob_GaussianMixture, prob_Gaussian_grid

# Correlation Algorithms
export CorrelationAlgorithm, StandardCorrelation, WindowedCorrelation, MultiTauCorrelation, IDLCorrelation, DEFAULT_CORRELATION_ALGORITHM

# Note: Test functions are not exported as they are for internal use only 