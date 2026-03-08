# Notes

## 2026-03 – write_traces / write_trace_dataframe

- **Verified working:** The current `write_traces` / `write_trace_dataframe` flow works with the correct data path.
- **Design:** Single `write_traces(folder; datapath=nothing, interval=1.0, ...)` discovers all `rates_*.txt`; for each file, if an info TOML exists (key-based) it uses `write_trace_dataframe_from_info` (all parameters from TOML, including hierarchical); otherwise (legacy) it uses `write_trace_dataframe(file, datapath, interval, ...)` with parameters from the filename. Legacy call `write_traces(folder, datapath; ...)` is supported via keyword forwarding.
- **Hierarchical:** Key-based path reads `hierarchical` from the TOML; if the field is present and nonempty, the computation is sent through the hierarchical stack.
- **Data path:** Using the wrong `datapath` (or wrong root in TOML) can cause trace length mismatches (e.g. BoundsError when model output length differs from data); ensure the path points at the intended trace data.
