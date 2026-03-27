# Run specification (info TOML)

Every fit run can write an **info TOML** file next to the rates and measures files (same stem, `.toml` extension). For example, `rates_myrun.txt` is accompanied by `info_myrun.toml`. This file records everything needed to reproduce or vary the run.

## Key-based naming

When you pass `key = "identifier"` to `fit`, all output files use that stem:

- `rates_identifier.txt`
- `info_identifier.toml`
- `measures_identifier.txt`, `param-stats_identifier.txt`, etc.

Without `key`, the existing long naming (e.g. `rates_label_gene_cond_model_nalleles.txt`) is used unchanged. Both conventions are supported; no migration is required.

## Contents of the info TOML

The file has four top-level sections:

| Section       | Description |
|---------------|-------------|
| `[run]`       | Full set of fit arguments (datapath, G, R, S, transitions, coupling, nchains, samplesteps, noisepriors, etc.). Keys match `fit(; kwargs...)` names. |
| `[output]`    | Key results from the run (e.g. `llml`, `accept`, `total`, `median_param`). |
| `[model_info]`| Model metadata (e.g. `rate_labels`, `interval`). |
| `[environment]`| `julia_version`, `threads` (for replication). |

Unset or optional arguments are stored as the string `"nothing"`; when the file is read back, `"nothing"` is interpreted as Julia `nothing`.

Coupling is stored as `coupling_unit_model` and `coupling_connections` in `[run]` (format `(unit_model, connections)` with each connection `(β, s, α, t)`).

## Starting a fit from a TOML file

Pass `key = "identifier"` to load the run specification from `info_<key>.toml` in the results folder. The results folder is `folder_path(resultfolder, root, "results")`; set `resultfolder` and `root` in the call if needed. All run options are read from the TOML's `[run]` section; any keyword arguments you pass to `fit` override those values.

- If the rates file `rates_<key>.txt` exists in that results folder, it is used as initial rates (warm start).

```julia
# Run using spec from info_myrun.toml (if it exists in resultfolder/root/results)
fit(key = "myrun", resultfolder = "myfolder", root = ".")

# Override options (same key, same output stem)
fit(key = "myrun", resultfolder = "myfolder", nchains = 4, samplesteps = 50_000)
```

## Programmatic use

```julia
using StochasticGene

# Load the run spec as a dictionary
spec = read_run_spec("results/myfolder/info_myrun.toml")
fit(; spec..., nchains = 4)

# Find the info TOML path for a given rates file
toml_path = info_toml_path_for_rates_file("results/myfolder/rates_myrun.txt")
# returns "results/myfolder/info_myrun.toml"

# Load spec from the rates file's companion TOML (if it exists)
spec = read_run_spec_for_rates_file("results/myfolder/rates_myrun.txt")
```

## See also

- [Model Fitting (fit)](@ref): full list of `fit` arguments.
- [Key-based naming](@ref): `fit(; key = "id", ...)`.
- [Cluster and batch workflows](cluster_batch_workflows.md): generating `info_<key>` presets with [`makeswarmfiles`](@ref) / [`write_run_spec_preset`](@ref) and combining single-unit rates for coupled starts.
