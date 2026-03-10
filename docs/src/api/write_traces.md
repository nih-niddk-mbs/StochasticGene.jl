# write_traces Function

Generate **model-predicted** intensity traces from fitted rates: for each rates file in a folder, load trace data, run the model with those rates, and write predicted traces (and optional state columns) to CSV.

## Syntax

```julia
write_traces(folder; datapath=nothing, interval=1.0, kwargs...)
write_traces(folder, datapath; interval=1.0, kwargs...)   # legacy
```

## Arguments

- `folder::String`: Folder containing `rates_*.txt` (and optionally `info_*.toml`) files.
- `datapath`: (optional) Path to trace data. Required for **legacy** runs when no info TOML exists for a rates file; can be omitted when all rates files have an info TOML (key-based).
- `interval::Float64 = 1.0`: Default frame interval (minutes). Ignored for key-based runs (interval comes from info TOML).
- `ratetype::String = "median"`: Which rate row to use (`"median"`, `"ml"`, etc.).
- `start`, `stop`: Trace frame range (defaults `1`, `-1`).
- `probfn`, `noiseparams`, `splicetype`: Model options (passed to prediction).
- `state::Bool = true`: Include G state / R state / reporter columns in output.
- `grid`, `zeromedian`, `datacol`: Data and output options.
- `hlabel::String = "-h"`: (legacy only) Label substring used to detect hierarchical in filename.

## Behavior

1. **Discovery**: Finds all `rates_*.txt` in `folder`.
2. **Per file**:
   - If an **info TOML** exists for that rates file (key-based): reads `read_run_spec_and_interval_for_rates_file`; uses `write_trace_dataframe_from_info` so **all** parameters (datapath, datacond, interval, hierarchical, coupling, G, R, S, etc.) come from the TOML. Hierarchical is taken from the TOML `hierarchical` field (nonempty → hierarchical stack).
   - **Legacy** (no info TOML): requires `datapath` (and uses `interval` kwarg). Uses `write_trace_dataframe(file, datapath, interval, ...)`, which extracts model parameters from the **filename** via `parse_filename` and `make_coupling` in io.jl.
3. **Output**: Writes `predictedtraces_<key>.csv` (or legacy stem) next to each rates file. Returns list of output CSV paths.

## Returns

- `Vector{String}`: Paths to the written predicted-trace CSV files.

## Examples

```julia
# Key-based: folder has rates_*.txt and info_*.toml; no datapath needed
write_traces("results/3Prime-coupled-2026-01-13/")

# Legacy: same folder, single datapath and interval for all rate files
write_traces("results/trace-test/", "data/testtraces"; interval=1.0, zeromedian=true)
```

## See also

- [Run specification (info TOML)](../run_spec_toml.md): how `info_*.toml` is read and what fields (datapath, datacond, interval, hierarchical, etc.) it supplies.
- `read_run_spec_and_interval_for_rates_file`, `write_trace_dataframe_from_info`, `write_trace_dataframe` in the source (analysis.jl, io.jl).
