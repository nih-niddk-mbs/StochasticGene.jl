# Notes

## Doc updates (since last doc update)

- **docs/src/api/fit.md**: Added "Trace and dwell specs (coupled models, hidden units)" with `trace_specs` and `dwell_specs`; noted storage in info TOML and use for 3-unit / hidden-unit runs.
- **docs/src/api/write_traces.md**: Rewritten to match current API: `write_traces(folder; datapath=nothing, interval=1.0, ...)` and legacy `write_traces(folder, datapath; ...)`; key-based vs legacy behavior; hierarchical from TOML; returns list of output paths.
- **docs/src/run_spec_toml.md**: Documented `trace_specs`, `dwell_specs`, and `interval` in TOML; added `read_run_spec_and_interval_for_rates_file` example.
- **README.md**: Updated write_traces example to current signature and key-based vs legacy.
- **fit.jl**: Docstrings for `fit(; ...)` and `make_reporter_components(::AbstractTraceData, ...)` document `trace_specs`, `dwell_specs`, and 3-unit hidden-unit pattern.
- **docs/src/trace_and_dwell_specs.md**: New guide for using trace and dwell specs (TraceSpec/DwellSpec fields, building from legacy or by hand, passing to fit, hidden units, helper functions, examples). Linked from index, api/index, fit.md, and run_spec_toml.md.

## 2026-03 – write_traces / write_trace_dataframe

- **Verified working:** The current `write_traces` / `write_trace_dataframe` flow works with the correct data path.
- **Design:** Single `write_traces(folder; datapath=nothing, interval=1.0, ...)` discovers all `rates_*.txt`; for each file, if an info TOML exists (key-based) it uses `write_trace_dataframe_from_info` (all parameters from TOML, including hierarchical); otherwise (legacy) it uses `write_trace_dataframe(file, datapath, interval, ...)` with parameters from the filename. Legacy call `write_traces(folder, datapath; ...)` is supported via keyword forwarding.
- **Hierarchical:** Key-based path reads `hierarchical` from the TOML; if the field is present and nonempty, the computation is sent through the hierarchical stack.
- **Data path:** Using the wrong `datapath` (or wrong root in TOML) can cause trace length mismatches (e.g. BoundsError when model output length differs from data); ensure the path points at the intended trace data.

## fit.jl – coupled models, hidden units, trace_specs

- **Docstrings updated** (fit.jl): Main `fit(; ...)` docstring now documents `trace_specs` and `dwell_specs`. `make_reporter_components(::AbstractTraceData, ...)` docstring documents `trace_specs` and the 3-unit-with-hidden pattern.
- **Coupled 3 units, 1 hidden:** Use `trace_specs` with one spec per *observed* unit (e.g. two specs for two observed units); the third unit is hidden (no observed trace). Pass coupling and tuple G for all three units; `trace_specs` length can be less than the number of units.
