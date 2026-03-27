# Package overview

This page complements the [home page](index.md) and [Getting started](getting_started.md) with **repository conventions**, **nomenclature**, and **I/O** details that apply across workflows.

## Project root, `data/`, and `results/`

StochasticGene expects a **project root** (`root` keyword, often `"."`). Under it:

- **`data/`** holds experimental inputs: histograms, trace folders, condition labels, and optional tables such as `CellType_alleles.csv` and `CellType_halflife.csv`.
- **`results/<resultfolder>/`** holds outputs from [fit](api/fit.md): rate files, diagnostics, and optional run-spec files.

The canonical resolution uses `folder_path` (see [Utilities](api/utilities.md)): if `joinpath(root, resultfolder)` exists it is used; otherwise `joinpath(root, "results", resultfolder)` is used (and created if needed). So fits usually land in **`root/results/<name>/`**.

**Version control:** the repository’s `.gitignore` excludes **`results/`** so large MCMC outputs stay local. Archive what you need for papers or collaboration separately (e.g. exported CSV summaries, small `info_*.toml` markers without huge binaries).

Use `rna_setup("dirname")` (see the [API reference](api/index.md)) to create a minimal tree with example data for learning the layout.

## Run specification files (`info_<key>`)

For key-based workflows, each run can be described by:

- **`info_<key>.toml`** — small marker pointing at the JLD2 companion.
- **`info_<key>.jld2`** — full keyword dict (including types like `method`, `probfn`) read by `read_run_spec` (see [Run specification (info TOML)](run_spec_toml.md)).

[fit](api/fit.md) with keyword `key` loads this spec when present and merges with explicit keywords (keywords win). Batch helpers such as `write_run_spec_preset`, `makeswarm`, and `makeswarmfiles` generate these files (see [Cluster and batch workflows](cluster_batch_workflows.md)).

## States vs steps (nomenclature)

- **G states** are **mutually exclusive** configurations of the promoter / gene regulatory model (the gene is in exactly one of **G** states at a time).
- **R steps** are positions along the transcription unit; many **occupancy patterns** are possible at once (which steps are loaded), so the “state space” of the elongation ladder is combinatorial. For **R** steps, the full GR state is a **tensor product** of G and the R occupancy pattern (and mRNA counts when tracked).
- With **splicing**, each R step can be represented in a higher alphabet (e.g. unoccupied / spliced / unspliced); state-space size grows with **G**, **R**, and **S**.

This distinction matters when setting **`onstates`**, interpreting **`simulator`** output, and reading coupling specs.

## Rate vector ordering (typical)

For telegraph-style models, the rate vector stacks contributions in a fixed order (see docstrings and papers): **G transitions**, then **R** transitions, **S** if present, then **decay** (and noise parameters for traces). Not every index is necessarily **fitted**; use **`fittedparam`** to select indices. Coupled models use **tuple** layouts for per-unit parameters and a **coupling** block.

## Data types (`datatype`)

Common values for [fit](api/fit.md):

| `datatype` | Meaning (short) |
|------------|-----------------|
| `"rna"` | Stationary RNA count histogram |
| `"rnaonoff"` | RNA + ON/OFF dwell histograms |
| `"rnadwelltime"` | RNA + multiple dwell-time types |
| `"trace"` | Intensity traces (HMM likelihood) |
| `"tracerna"` | Traces + RNA histogram |
| `"tracejoint"` | Coupled / joint traces between units |
| `"tracegrid"` | Grid-based trace likelihood |

`datapath` may be a file, folder, or vector of paths depending on type. **`traceinfo`** encodes frame interval, time window, active fraction, and background handling for traces.

## Output file families

Under each **`resultfolder`**, names encode cell, condition, gene, model string (**G**, **R**, **S**, **insertstep**), and alleles. Typical prefixes:

- **`rates_*.txt`** — posterior summaries and ML row (see file header / docs).
- **`measures_*.txt`**, **`param-stats_*.txt`** — diagnostics and parameter summaries.
- **`burst_*.txt`** — optional burst statistics when requested.
- **`optimized_*.txt`** — optional optimizer output.

Underscore **`_`** separates fields in filenames; avoid **`_`** inside user labels where the naming convention would become ambiguous.

## Cluster workflows (pointer)

For NIH Biowulf **swarm** generation, **`makeswarm`** / **`makeswarm_genes`** / **`makeswarmfiles`**, and the **coupled** pipeline (single-unit fits → **`create_combined_file`** → coupled **`fit`**), read [Cluster and batch workflows](cluster_batch_workflows.md).

## Coupled CSV (`Coupled_models_to_test`)

Batch coupled jobs can read a seven-column (or wider) CSV processed in **`coupled_csv.jl`**: model name plus six cells describing enhancer→gene and gene→enhancer connections and signs. See docstrings for `csv_row_to_connections_simple`, `build_coupled_fit_spec_from_csv_cells`, and `makeswarmfiles_coupled` in the Julia package (use `?` in the REPL after `using StochasticGene`), and the *Coupled CSV* section in the long docstring of `makeswarmfiles` (source: `biowulf.jl`; overview: [Cluster and batch workflows](cluster_batch_workflows.md)).
