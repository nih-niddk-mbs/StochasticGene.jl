# Pre-2026 trace-reading code (reference)

Trace reading is implemented in **io.jl** and **fit.jl**. The following is the established API as of the last pre-2026 commit (`3396160`, 2025-08-19, "keep nonunique traces"). Use these entry points; do not reimplement path or trace-load logic in analysis.

---

## io.jl – low-level trace file reading

- **`read_tracefile(path, start, stop, col=3)`**  
  Read one trace file into a vector. `path` = file path; `start`/`stop` = frame indices (-1 = end); `col` = column.

- **`read_tracefiles(path, label, start, stop, col=3, uniquetrace=false)`**  
  Single condition: `path` = folder, `label` = string to match in filename. Returns `Vector{Vector}`.

- **`read_tracefiles(path, label::Vector{String}, start, stop, col=3)`**  
  Joint (tracejoint): `label` = vector of condition names; files are paired by name. Returns `Vector{Matrix}`.

- **`read_tracefiles(path, label, traceinfo::Tuple, col=3)`**  
  Same as above but `traceinfo = (dt, start_time, stop_time, ...)`; converts to frame indices via `start = max(round(Int, traceinfo[2]/traceinfo[1]), 1)` and `stop` from `traceinfo[3]`.

- **`read_tracefiles_unbalanced(path, label::Vector{String}, start, stop, col=3; backup_path="")`**  
  Joint traces when file counts per condition differ; can use `backup_path` for missing partners.

- **`read_tracefiles_grid(path, label, traceinfo)`** / **`read_tracefiles_grid(path, label, start, stop)`**  
  Grid traces; returns unique traces by sum.

- **`read_tracefile_grid(path, start, stop, header=true)`**  
  Single grid file via `readdlm`.

Path convention: **`path`** is the directory to search; **`label`** (or **`datacond`** in fit.jl) is the condition name(s) used to match filenames. No extra path rewriting is done in these functions.

---

## fit.jl – trace data loading (pre-2026 API)

- **`load_data_trace(datapath, label, gene, datacond, traceinfo, datatype::Symbol, col=3, zeromedian=false)`**  
  - Calls **`read_tracefiles(datapath, datacond, traceinfo, col)`** (or `datapath[1]` when `datapath` is a vector).  
  - `datacond`: single condition string or vector of strings for joint.  
  - Applies `zero_median`, builds `TraceData` (or `TraceRNAData` for `:tracerna`).  
  - **Does not** change or reinterpret `datapath`; it is passed straight through to `read_tracefiles`.

- **`load_data_tracegrid(datapath, label, gene, datacond, traceinfo)`**  
  Uses **`read_tracefiles_grid(datapath, datacond, traceinfo)`**.

- **`load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, ...)`**  
  Dispatches to `load_data_trace` or `load_data_tracegrid` for trace types.

So: **datapath** and **datacond** come from the run spec (or caller). Trace reading is entirely **io.jl** (`read_tracefiles` / `read_tracefiles_grid`); **fit.jl** only wires datatype, zeromedian, and traceinfo.

---

## analysis.jl – pre-2026 trace writing pattern

- **`write_trace_dataframe(outfile, datapath, datacond, interval, r, transitions, G, R, S, insertstep, ...)`**  
  Core: builds `traceinfo = (interval, start_time, stop_time)` and calls  
  **`data = load_data_trace(datapath, "", "", datacond, traceinfo, :trace, datacol, zeromedian)`**  
  then uses `data` and `r` to write predictions. **datapath** and **datacond** are taken as given (e.g. from run spec).

- **`write_trace_dataframe(file, datapath, interval, ratetype, ...)`**  
  Legacy file-based: uses **`parse_filename(file)`** and **`make_coupling`** from io.jl to get `datacond`, transitions, G, R, S, etc.; then one call that rewrote `datapath` for coupled runs:  
  `if G isa Int && !isempty(coupling) datapath = joinpath(datapath, string(coupling_field[1]))`.

- **`write_traces(folder, datapath, interval, ...)`**  
  Scans folder for rate files and calls the file-based `write_trace_dataframe(..., datapath, datacond, ...)` per file.

So the pre-2026 pattern is: **get (datapath, datacond, interval) from run spec or file parsing → call `load_data_trace(datapath, "", "", datacond, traceinfo, :trace, ...)` → no extra path logic in analysis.**

---

## Summary

| Layer    | Pre-2026 entry points |
|----------|------------------------|
| **io.jl**  | `read_tracefile`, `read_tracefiles` (String; Vector{String}; traceinfo), `read_tracefiles_unbalanced`, `read_tracefiles_grid`, `read_tracefile_grid` |
| **fit.jl** | `load_data_trace`, `load_data_tracegrid`, `load_data` (for trace types) |
| **analysis** | `write_trace_dataframe(outfile, datapath, datacond, interval, r, ...)` using `load_data_trace`; file-based overload uses `parse_filename` / `make_coupling` |

For key-based (folder + key) flows, the intended approach is: **resolve (datapath, datacond, interval) in io.jl** (e.g. `trace_load_spec_for_rates_file` or equivalent that reads the info TOML and optional root), then call the **same** `load_data_trace` / `write_trace_dataframe(out, datapath, datacond, interval, r, ...)` with no new path logic in analysis.
