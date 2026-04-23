# Package overview

This page complements the [home page](index.md) and [Getting started](getting_started.md) with **repository conventions**, **nomenclature**, and **I/O** details that apply across workflows.

## Project root, `data/`, and `results/`

StochasticGene expects a **project root** (`root` keyword, often `"."`). Under it:

- **`data/`** holds experimental inputs: histograms, trace folders, condition labels, and optional tables such as `CellType_alleles.csv` and `CellType_halflife.csv`.
- **`results/<resultfolder>/`** holds outputs from [fit](api/fit.md): rate files, diagnostics, and optional run-spec files.

The canonical resolution uses `folder_path` (see [Utilities](api/utilities.md)): if `joinpath(root, resultfolder)` exists it is used; otherwise `joinpath(root, "results", resultfolder)` is used (and created if needed). So fits usually land in **`root/results/<name>/`**.

**Version control:** the repository’s `.gitignore` excludes **`results/`** so large fit outputs (MCMC, NUTS, ADVI) stay local. Archive what you need for papers or collaboration separately (e.g. exported CSV summaries, small `info_*.toml` markers without huge binaries).

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
- **`proposal-cov_*.jld2`** — proposal covariance matrix and metadata (see MCMC proposal & warmup below).
- **`burst_*.txt`** — optional burst statistics when requested.
- **`optimized_*.txt`** — optional optimizer output.

Underscore **`_`** separates fields in filenames; avoid **`_`** inside user labels where the naming convention would become ambiguous.

## MCMC proposal covariance and warmup

When fitting expensive models (e.g., with ODE-based likelihood evaluation that takes minutes per step), proposal covariance reuse can significantly speed up workflows:

### Proposal Covariance Reuse

The `propcv` keyword controls the proposal distribution:

- **`propcv=0.01` (positive)**: Use fixed coefficient of variation. MCMC will compute empirical covariance during warmup (if `warmupsteps > 0`) and save it to `proposal-cov_*.jld2`.
- **`propcv=-0.01` (negative)**: Attempt to load covariance from `proposal-cov_*.jld2` if it exists and model parameters match exactly (G, R, S, transitions, fittedparam, nalleles all equal). 
  - If **loading succeeds**: Warmup is **automatically skipped** (even if `warmupsteps > 0`), and sampling proceeds immediately with the loaded proposal.
  - If **loading fails**: Falls back to `abs(propcv)` and warmup proceeds normally.

**Workflow example:**
```julia
# First run: compute and save covariance
fits1 = fit(; G=2, R=0, transitions=([1,2],[2,1]), 
            propcv=0.01, warmupsteps=10000, ...)

# Subsequent run: reuse covariance (warmup skipped automatically)
fits2 = fit(; G=2, R=0, transitions=([1,2],[2,1]), 
            propcv=-0.01, warmupsteps=10000, ...)  # warmup still specified but skipped
```

No need to remember to set `warmupsteps=0` on the second run — the presence of a loaded covariance automatically prevents warmup from running.

### Adaptive Warmup

When `warmupsteps > 0` and no covariance is loaded, the warmup phase adapts the proposal covariance to improve MCMC efficiency:

- **Periodic Adaptation**: Adapts every `max(1000, samplesteps ÷ 3)` steps (typically 2–3 times per warmup).
- **Acceptance Rate Targeting**: 
  - Acceptance target scales with problem dimension: 44% (d=1) → 30% (d=5–20) → 23.4% (d>>1).
  - If current rate < 15%, shrinks proposals; if > 40%, expands them.
- **Time Allocation**: Warmup time is proportional to step count: `warmup_time = maxtime × (warmupsteps / total_steps)`. For expensive steps, increase `maxtime` or reduce `warmupsteps` if warmup times out before adaptation triggers.

The saved covariance is stored as `proposal-cov_<name>.jld2` with metadata validation, ensuring proposals are only reused when model structure matches.

## Cluster workflows (pointer)

For NIH Biowulf **swarm** generation, **`makeswarm`** / **`makeswarm_genes`** / **`makeswarmfiles`**, and the **coupled** pipeline (single-unit fits → **`create_combined_file`** → coupled **`fit`**), read [Cluster and batch workflows](cluster_batch_workflows.md).

## Coupled CSV (`Coupled_models_to_test`)

Batch coupled jobs can read CSVs processed in **`coupled_csv.jl`**. Each row defines a coupled model by specifying connections between units.

**Column naming** (by pattern; order-independent):
- `Model_name` — required, key for the model
- `enhancer_to_gene_1`, `enhancer_sign_1` — optional, enhancer→gene connections for unit 1
- `enhancer_to_gene_2`, `enhancer_sign_2` — optional, enhancer→gene connections for unit 2  
- `gene_to_enhancer` or `gene_to_enhancer_sign` — optional, gene→enhancer connections
- `background_gene`, `background_gene_sign` —optional, genetic background connections
- Sign columns (`*_sign`) accept: `">0"` (activate), anything else except empty/`0`/`"free"` (inhibit), or empty/`0`/`"free"` (**:free** mode)

Minimal example: `Model_name`, `enhancer_to_gene_1`, `enhancer_sign_1` (other couplings default to `:free` mode).

See docstrings for `csv_row_to_connections_simple`, `build_coupled_fit_spec_from_csv_cells`, and `makeswarmfiles_coupled` (use `?` in REPL after `using StochasticGene`), and the *Coupled CSV* section in the long docstring of `makeswarmfiles` (source: `biowulf.jl`; overview: [Cluster and batch workflows](cluster_batch_workflows.md)).
