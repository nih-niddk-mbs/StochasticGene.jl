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

Use `spec_file` to load the run specification from an info TOML. All run options are read from the TOML's `[run]` section; any keyword arguments you pass to `fit` override those values.

- **Warm start**: If the TOML has `infolder` and (for key-based runs) `key`, and the corresponding rates file exists, it is used as initial rates.
- **Cold start**: Pass `cold = true` to ignore `infolder`/`inlabel` and start from the prior.

```julia
# Reproduce a run exactly
fit(spec_file = "results/myfolder/info_myrun.toml")

# Reproduce but override options and write to a new stem
fit(
    spec_file = "results/myfolder/info_myrun.toml",
    nchains = 4,
    samplesteps = 50_000,
    resultfolder = "myfolder_long_run",
    key = "long_run"
)

# Cold start from TOML (ignore any existing rates file)
fit(spec_file = "results/myfolder/info_myrun.toml", cold = true, key = "cold_run")
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
