# Getting Started

## Folder Structure

StochasticGene expects a specific directory structure for data and results. For **naming conventions**, **rate ordering**, and why **`results/`** is typically git-ignored, see [Package overview](package_overview.md).

```bash
project_root/
├── data/
│   ├── alleles/      # Contains allele count information
│   └── halflives/    # Contains mRNA half-life information
└── results/          # Output directory for analysis results (add to .gitignore locally)
```

## Setting Up

Use the built-in helper to set up a new project folder with the required data structure and reference files:

```julia
using StochasticGene
rna_setup()           # Sets up folders in the current directory
# or, to set up in a subfolder:
rna_setup("scRNA")    # Sets up folders in ./scRNA
cd("scRNA")           # Move into the project directory (if you want to work there)
```

This creates the necessary folders and seeds bundled reference CSV files for alleles, halflives, and the small `HCT116_testdata` mock RNA histogram set. If a reference file is not bundled in the local checkout, setup falls back to downloading it.

To generate Biowulf scripts for the mock RNA histogram data in the folder created by `rna_setup`, scan the data folder with `makeswarm_genes`:

```julia
using StochasticGene

root = "scRNA"
rna_setup(root)

makeswarm_genes(
    root=root,
    datapath="HCT116_testdata",
    datacond="MOCK",
    resultfolder="HCT116_test",
    filedir=joinpath(root, "run-HCT116-testdata-rna"),
    cell="HCT116",
    datatype="rna",
    G=2,
    R=0,
    S=0,
    insertstep=1,
    transitions=([1, 2], [2, 1]),
    nchains=8,
)
```

This writes one shared fitscript with `gene=ARGS[1]` and one command-file line
per gene found in `root/data/HCT116_testdata`. For large Biowulf RNA sweeps,
keep the default single `fit.swarm` file and use `swarm -b <bundle>` when
submitting to bundle several gene commands into each Slurm array task.

After the jobs finish, assemble the fit outputs into analysis CSVs:

```julia
make_dataframes(
    resultfolder = joinpath(root, "results", "HCT116_test"),
    datapath = joinpath(root, "data", "HCT116_testdata"),
    datatype = "rna",
)

write_dataframes_only(
    joinpath(root, "results", "HCT116_test"),
    joinpath(root, "data", "HCT116_testdata");
    datatype = "rna",
)
```

This creates files such as `Summary_rna-HCT116_2.csv`, with fitted rates,
parameter summaries, condition/model metadata, and observed RNA moments.

For key-based result folders with `rates_<key>.txt`, `param-stats_<key>.txt`,
and `info_<key>.jld2`, use:

```julia
make_dataframes_key(joinpath(root, "results", "HCT116_test"); datatype = "rna")
write_dataframes_only_key(joinpath(root, "results", "HCT116_test"); datatype = "rna")
```

The key version writes `Summary_key.csv` and keeps the original fit key in a
`Key` column.

## Data Preparation

### RNA Count Data
RNA count data should be provided as a text file with the following format:
```text
# gene_condition.txt
count1
count2
count3
...
```

The data will be loaded into an `RNAData` structure with the following fields:
- `histRNA`: Vector of RNA counts
- `gene`: Gene name
- `condition`: Experimental condition
- `nalleles`: Number of alleles

### Live Cell Imaging Data
For live cell imaging data, provide a text file with time series data:
```text
# gene_condition.txt
time1 intensity1
time2 intensity2
time3 intensity3
...
```

## Fitting a Model

### Basic Two-State Model
Fit a simple two-state model to RNA count data:

```julia
fits, stats, measures, data, model, options = fit(
    G = 2,                    # Two gene states
    R = 0,                    # No pre-RNA steps
    transitions = ([1,2], [2,1]),  # Transitions between states
    datatype = "rna",         # RNA count data
    datapath = "data/HCT116_testdata/",
    gene = "MYC",            # Gene name
    datacond = "MOCK"        # Experimental condition
    # Optional: inference_method = :mh (default), :nuts, :advi; see docs/src/api/fit.md
)
```

### Model with Pre-RNA Steps
Fit a model with pre-RNA processing steps:

```julia
fits, stats, measures, data, model, options = fit(
    G = 2,                    # Two gene states
    R = 3,                    # Three pre-RNA steps
    S = 2,                    # Two splice sites
    insertstep = 1,           # Reporter insertion at step 1
    transitions = ([1,2], [2,1]),  # Gene state transitions
    datatype = "trace",       # Live cell imaging data
    datapath = "data/testtraces",
    gene = "MS2",            # Gene name
    datacond = "testtrace",   # Experimental condition
    trace_specs = [(unit = 1, interval = 1.0, start = 1.0, t_end = -1.0, zeromedian = true, active_fraction = 1.0, background = 0.0)],
    noisepriors = [40., 20., 200., 10.],  # Noise parameters
    nchains = 4              # Number of parallel chains (pooling for MH; see `run_inference` for NUTS/ADVI)
)
```

## Model Components

### Gene States (G)
- Arbitrary number of gene states
- User-specified transitions between states
- One active state for transcription initiation
- Example transitions:
  - Two-state: `([1,2], [2,1])` (telegraph model)
  - Three-state: `([1,2], [2,1], [2,3], [3,1])` (cyclic model)

### Pre-RNA Steps (R)
- Irreversible forward transitions
- mRNA ejection from final R step
- Optional reporter insertion step
- Example: R=3 with reporter at step 1:
  ```julia
  R = 3
  insertstep = 1
  ```

### Splicing (S)
- Up to R splice sites
- PreRNA with or without spliced intron
- Multiple configurations per R step
- Example: S=2 with R=3:
  ```julia
  R = 3
  S = 2
  ```

## Model Types

StochasticGene supports several model types:

### G Models (Telegraph Models)
- Basic two-state model
- Multiple gene states
- Example:
  ```julia
  G = 2
  R = 0
  transitions = ([1,2], [2,1])
  ```

### GR Models (with Pre-RNA Steps)
- Gene states plus pre-RNA processing
- Example:
  ```julia
  G = 2
  R = 3
  transitions = ([1,2], [2,1])
  ```

### GRS Models (with Splicing)
- Gene states, pre-RNA steps, and splicing
- Example:
  ```julia
  G = 2
  R = 3
  S = 2
  insertstep = 1
  ```

### Coupled Models (Multiple Alleles)
- Multiple alleles with coupling
- Example:
  ```julia
  G = (2,2)  # Two alleles, each with 2 states
  R = (3,3)  # Three pre-RNA steps for each
  coupling = (
      (1,1),  # Model indices
      ([2], [1]),  # Source units
      (["G2"], ["G1"]),  # Source states
      ([2], [1])  # Target transitions
  )
  ```

## Data Types

The package can handle:

### RNA Count Data
- Single molecule FISH (smFISH)
- Single cell RNA sequencing (scRNA-seq)
- Example:
  ```julia
  datatype = "rna"
  datapath = "data/rna_counts/"
  ```

### Live Cell Imaging
- MS2 reporter data
- PP7 reporter data
- Example:
  ```julia
  datatype = "trace"
  datapath = "data/ms2_traces/"
  traceinfo = (1.0, 1., -1, 1.)  # Frame interval, start, end, active fraction
  ```

### Dwell Time Analysis
- ON/OFF state durations
- Example:
  ```julia
  datatype = "dwelltime"
  datapath = "data/dwell_times/"
  ```

### Combined Data
- RNA counts with dwell times
- Example:
  ```julia
  datatype = (:rna, :dwelltime)
  datapath = (
      rna = "data/rna_counts/",
      dwelltime = ["data/dwell_times/on.csv", "data/dwell_times/off.csv"],
  )
  dwell_specs = [(unit = 1, onstates = [[2], Int[]], dttype = ["ON", "OFF"])]
  ```

Use tuple/vector `datatype` values for new multimodal fits. StochasticGene builds a `CombinedData` object with independent modality legs, evaluates likelihoods separately, and combines the results for total likelihood and WAIC. Legacy strings such as `"rnadwelltime"` and `"tracerna"` remain available for older scripts.

Older `infolder` and `inlabel` inputs have been retired. Use `datapath` for inputs, `label` for output/data labels, `resultfolder` for outputs, and `root` for the project root.

## Next Steps

- **Cluster / Biowulf / batch fits:** read [Cluster and batch workflows](cluster_batch_workflows.md) for `makeswarm`, `makeswarmfiles`, and the **coupled-model workflow** (fit units separately, merge rates with `create_combined_file`, then run the coupled fit).
- Check the Examples section for more complex usage scenarios
- Read the API Reference for detailed function documentation
- Join the [GitHub discussions](https://github.com/nih-niddk-mbs/StochasticGene.jl/discussions) for community support
