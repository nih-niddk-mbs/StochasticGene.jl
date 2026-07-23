# Examples

This section contains examples of using StochasticGene.jl for various analysis
scenarios.

The most current, runnable workflow examples are:

- [RNA Histogram Analysis](rna_histogram.md)
- [Multi-State Model](multi_state_rna.md)
- [Shared Parameter Trace Fits](shared_parameter_fits.md)
- [Cluster and batch workflows](../cluster_batch_workflows.md)

Several older pages in this section are retained as conceptual sketches while
the v1.11/2.0 API is being finalized. Those pages are marked with a legacy
notice and should not be copied verbatim until refreshed.

## Basic Examples

### RNA Count Analysis
- [Basic Telegraph Model](basic_telegraph.md): Legacy conceptual sketch for a simple two-state model
- [Multi-State Model](multi_state_rna.md): Current workflow for fitting a model with multiple gene states
- [RNA Histogram Analysis](rna_histogram.md): Current workflow for RNA histogram/scRNA sweeps

### Live Cell Imaging
- [MS2 Reporter Analysis](ms2_analysis.md): Legacy conceptual sketch for MS2 reporter data
- [Trace Analysis](trace_analysis.md): Legacy conceptual sketch for time-series data

### Dwell Time Analysis
- [Basic Dwell Times](dwell_time_analysis.md): Legacy conceptual sketch for dwell-time distributions
- [RNA Dwell Time Analysis](rna_dwell_time.md): Legacy conceptual sketch for RNA+dwell-time analysis

## Advanced Examples

Cluster / Biowulf batch generation and combined rate files are documented on the top-level page [Cluster and batch workflows](../cluster_batch_workflows.md) (`makeswarmfiles`, `create_combined_file`, etc.), not in this examples tree.

### Parallel Processing
- [Parallel Processing](parallel_processing.md): Legacy conceptual sketch; use [Cluster and batch workflows](../cluster_batch_workflows.md) for current batch generation
- [Shared Parameter Trace Fits](shared_parameter_fits.md): Current workflow for exact parameter sharing across grouped datasets

## Contributing Examples

We welcome contributions of new examples! If you have an interesting use case or analysis workflow, please consider sharing it with the community. See the Contributing guide for more information.
