# Trace and dwell specs

**TraceSpec** and **DwellSpec** are the canonical way to define which units are observed and how for trace (intensity time-series) and dwell-time data. They support single-unit and coupled multi-unit models, including **hidden units** (units with no observed trace).

## When to use specs

- **trace_specs**: A vector of **`TraceSpec`** structures. Use with `datatype = "trace"` or `"tracejoint"`. When not nothing, `trace_specs` are used and legacy `onstates`/`traceinfo` are ignored. Supports per-unit control and **hidden units** (e.g. 3 coupled units but only 2 observed traces).
- **dwell_specs**: A vector of **`DwellSpec`** structures. Use with `datatype = "dwelltime"` or `"rnadwelltime"`. When not nothing, `dwell_specs` are used and legacy `onstates`/`bins`/`dttype` are ignored.

If you do not pass specs (or pass nothing), `fit` uses the legacy arguments. When you pass a non-nothing vector of specs, they are used and the corresponding legacy quantities are ignored. Specs are stored in `info_*.toml` when using key-based runs.

---

## TraceSpec

A **TraceSpec** describes one trace (intensity) channel: which unit it observes, what emission it reports, and the time window.

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `unit` | `Int` | Which unit this trace observes (1-based). |
| `emission` | `Symbol` or `NamedTuple` | `:reporter_sum` (default) = sum over R steps from insertstep; or `(kind=:G, state=k)` for a single G state k. |
| `frame_interval` | `Float64` | Time between frames (minutes). |
| `start_frame` | `Float64` | Start time/frame (minutes). |
| `end_frame` | `Float64` or `Int` | End time/frame; use `-1` for last. |
| `active_fraction` | `Float64` | Fraction of observed active traces (default `1.0`). |
| `background_mean` | optional | Background mean (scalar or vector per channel). |
| `output_index` | `Int` or `nothing` | Optional ordering (default `nothing`). |

One trace channel per TraceSpec; order of specs = order of columns in trace data.

### Construction

**By hand (single unit, reporter sum):**

```julia
using StochasticGene

# Unit 1, reporter sum, 1 min interval, frames 1 to end, active fraction 1
spec = TraceSpec(1, :reporter_sum, 1.0, 1.0, -1)
trace_specs = [spec]
```

**From legacy (single unit):**

```julia
# Same as legacy onstates=Int[], traceinfo=(1.0, 1.0, -1)
trace_specs = trace_specs_from_legacy(Int[], (1.0, 1.0, -1), G, R, S, insertstep; unit=1)
```

**From legacy (coupled, two units, same traceinfo):**

```julia
# onstates per unit: empty = reporter_sum for that unit
onstates = [Int[], Int[]]  # both units: reporter sum
trace_specs = trace_specs_from_legacy(onstates, (1.0, 1.0, -1), G, R, S, insertstep, coupling)
```

**Coupled with hidden unit (3 units, only units 1 and 2 observed):**

```julia
# G, R, S, insertstep, coupling for 3 units
# Only build specs for units 1 and 2; unit 3 has no spec = hidden
onstates = [Int[], Int[], Int[]]  # reporter sum for units 1 and 2; unit 3 unused
trace_specs = trace_specs_from_legacy(onstates, (1.0, 1.0, -1), G, R, S, insertstep, coupling; units=[1, 2])
```

Then pass to fit:

```julia
fit(;
    datatype = "tracejoint",
    G = (3, 3, 3),
    R = (3, 3, 3),
    # ... transitions, coupling, datapath, datacond, etc.
    trace_specs = trace_specs,
)
```

---

## DwellSpec

A **DwellSpec** describes one dwell-time (sojourn) histogram: which unit, gene vs R space, in-state vs out-of-state, on-states, and bins.

### Fields

| Field | Type | Description |
|-------|------|-------------|
| `unit` | `Int` | Which unit (1-based). |
| `kind` | `Symbol` | `:G` (gene-state space) or `:R` (R-step / full state space). |
| `sojourn_kind` | `Symbol` | `:in_state` (time in the "on" set) or `:out_of_state` (time out). |
| `onstates` | `Vector{Int}` | The set of states considered "on" for this histogram (empty for R = any step). |
| `bins` | `Vector{Float64}` | Histogram bin edges. |
| `output_index` | `Int` or `nothing` | Optional ordering (default `nothing`). |

**Mapping to legacy dttype:**

- `kind=:R`, `sojourn_kind=:out_of_state` → `"OFF"`
- `kind=:R`, `sojourn_kind=:in_state` → `"ON"`
- `kind=:G`, `sojourn_kind=:out_of_state` → `"OFFG"`
- `kind=:G`, `sojourn_kind=:in_state` → `"ONG"`

### Construction

**By hand:**

```julia
# Unit 1, R-space, OFF (out of state), onstates e.g. [2,3], bins
spec_off = DwellSpec(1, :R, :out_of_state, [2, 3], [0.0, 1.0, 2.0, 10.0])
spec_on  = DwellSpec(1, :R, :in_state,  [2, 3], [0.0, 1.0, 2.0, 10.0])
dwell_specs = [spec_off, spec_on]
```

**From legacy (single unit):**

```julia
onstates = [[2, 3], [2, 3]]   # ON set for each dwell type
bins = [[0.0, 1.0, 2.0, 10.0], [0.0, 1.0, 2.0, 10.0]]
dttype = ["OFF", "ON"]
dwell_specs = dwell_specs_from_legacy(onstates, bins, dttype; unit=1)
```

Then pass to fit:

```julia
fit(;
    datatype = "dwelltime",
    G = 2, R = 3, S = 0, insertstep = 1,
    # ... datapath, etc.
    dwell_specs = dwell_specs,
)
```

---

## Passing specs to fit

- **trace_specs**: `fit(; trace_specs = [TraceSpec(...), ...], ...)`. Vector of `TraceSpec` structures. When not nothing, used in place of legacy `onstates`/`traceinfo`. For coupled models, one HMMReporter is built per spec; fewer specs than units means the missing units are **hidden**.
- **dwell_specs**: `fit(; dwell_specs = [DwellSpec(...), ...], ...)`. Vector of `DwellSpec` structures. When not nothing, used in place of legacy `onstates`/`bins`/`dttype`.

When using key-based runs (`fit(; key = "myrun", ...)`), `trace_specs` and `dwell_specs` are read from and written to `info_myrun.toml` in the `[run]` section (as arrays of tables).

---

## Hidden units (trace)

For coupled trace/tracejoint models with **3 units and 1 hidden** (only 2 observed traces):

1. Define `G`, `R`, `S`, `insertstep`, and `coupling` for **all 3 units**.
2. Build **trace_specs** with **2** specs (one per observed unit). Omit a spec for the hidden unit.
3. Pass `trace_specs` to `fit`. The model still has 3 units; the likelihood uses only the 2 observed channels.

Example (conceptually):

```julia
# 3 units; only units 1 and 2 have trace data
trace_specs = [
    TraceSpec(1, :reporter_sum, 1.0, 1.0, -1),
    TraceSpec(2, :reporter_sum, 1.0, 1.0, -1),
]
fit(;
    datatype = "tracejoint",
    G = (3, 3, 3),
    R = (3, 3, 3),
    S = (0, 0, 0),
    insertstep = (1, 1, 1),
    transitions = (...),
    coupling = ((1, 2), [...]),  # connections between units
    datapath = "...",
    datacond = ["label1", "label2"],
    trace_specs = trace_specs,
)
```

`datacond` should have one label per **observed** channel (2 here). See `trace_specs_from_legacy(..., units=[1, 2])` for building from legacy onstates/traceinfo.

---

## Helper functions (specs.jl)

| Function | Purpose |
|----------|---------|
| `trace_specs_from_legacy(onstates, traceinfo, G, R, S, insertstep; unit=1)` | Single-unit trace specs from legacy. |
| `trace_specs_from_legacy(onstates, traceinfo, G, R, S, insertstep, coupling; units=nothing)` | Coupled; `units` = which unit indices get specs (omit for hidden). |
| `dwell_specs_from_legacy(onstates, bins, dttype; unit=1)` | Single-unit dwell specs from legacy. |
| `legacy_onstates_traceinfo(trace_specs; nunits=1)` | Convert trace_specs back to (onstates, traceinfo). |
| `legacy_dwell(dwell_specs)` | Convert dwell_specs back to (onstates, bins, dttype). |

---

## See also

- [Model Fitting (fit)](api/fit.md): full `fit` arguments including `trace_specs` and `dwell_specs`.
- [Run specification (info TOML)](run_spec_toml.md): how specs are stored in `info_*.toml`.
- `src/specs.jl`: struct definitions and conversion helpers.
