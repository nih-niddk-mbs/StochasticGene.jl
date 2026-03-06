# Plan: Measurement Spec redesign and multi-unit multi-allele RNA

Naming: **TraceSpec** and **DwellSpec** (replacing "parcel" in the API).

## 1. Current scope (single-entity trace and dwelltime)

### 1.1 TraceSpec (single entity) — implemented
- **unit**: single unit index
- **emission**: `:reporter_sum` or `(kind=:G, state=k)` for a single G state
- **frame_interval**, **start_frame**, **end_frame**
- Optional: **active_fraction**, **background_mean**, **output_index**

### 1.2 DwellSpec (single entity) — implemented
- **unit**: single unit index
- **kind**: `:G` or `:R` (gene-state vs R-step space)
- **sojourn_kind**: `:in_state` | `:out_of_state` (replaces ON/OFF/ONG/OFFG)
- **onstates**: the set of "on" states for this histogram
- **bins**: histogram bins for this dwell type
- Optional: **output_index**

### 1.3 Implementation status
- **Done**: `TraceSpec` and `DwellSpec` structs in `src/specs.jl`; conversion helpers `trace_specs_from_legacy`, `dwell_specs_from_legacy`, `legacy_onstates_traceinfo`, `legacy_dwell`. `fit(; trace_specs=..., dwell_specs=...)` and `make_structures` accept specs and derive legacy onstates/traceinfo/bins/dttype for the current pipeline.
- **Done**: Refactor make_reporter_components to build from trace_specs/dwell_specs when provided: trace path builds one `HMMReporter` per `TraceSpec` for coupled models (enables hidden units); dwell path uses `legacy_dwell(dwell_specs)` for single-unit. Helper `hmm_reporter_from_trace_spec` in `specs.jl`.
- **Done**: I/O: `run_spec_to_toml` and `read_run_spec` in `io.jl` serialize/deserialize `trace_specs` and `dwell_specs` to/from info_*.toml (array-of-tables for specs).
- **Done**: Tests: `test_spec_conversion`, `test_spec_io_roundtrip`, `test_spec_trace_in_fit`, `test_spec_dwell_in_fit`.
- **Done (entry-point)**: Simulator accepts optional `trace_specs` and `dwell_specs`; derives onstates, bins, traceinterval via `legacy_onstates_traceinfo` / `legacy_dwell` and runs with existing logic. Full internal "loop over specs" (set_onoff, set_before, make_trace from spec.unit/emission) remains optional for later.
- **Done**: Likelihood/prediction is spec-driven via make_structures: when trace_specs/dwell_specs are passed to fit, reporter and components are built from specs (or legacy derived), so ll_hmm_trace and predictedfn use the spec-derived model.
- **Hidden units**: supported in trace coupled path (fewer reporter channels than units when fewer trace_specs than units).
- Leave RNA and nalleles unchanged for now.

---

## 2. Deferred: multi-unit multi-allele RNA

### 2.1 RNA histogram parcel (entity)
- **units**: which units contribute, e.g. `[1]` or `[1, 2]`
- **nalleles**: per contributing unit, e.g. `2` or `[2, 1]`
- **yield**: per contributing unit, e.g. `1.0` or `[1.0, 0.8]`
- Other RNA-specific fields (bins, nRNA, etc.) as needed.

### 2.2 Prediction pipeline (chemical_master.jl)
1. **Joint steady state**: solve combined master equation for coupled units → **P(n₁, n₂, …)**.
2. **Marginals**: **Pᵢ(nᵢ)** from joint for each unit.
3. **Per-unit allele convolution**: for unit i, convolve in that coordinate **(kᵢ − 1)** times with **Pᵢ** (1D convolutions per dimension).
4. **Aggregation**: **Y = ∑ᵢ yᵢ Sᵢ**; compute **P(Y = n)** from the convolved joint.
5. API: accept RNA parcel(s); from each parcel use **units**, **nalleles**, **yield** and run the above.

### 2.3 Simulator
- For each RNA parcel: use **units**, **nalleles**, **yield** to define what is summed into each histogram (already compatible with current simulator idea).

### 2.4 Optional later: multi-unit trace/dwell
- **Source** for trace/dwell: `SingleUnit(unit, nalleles=1)` vs `AggregatedUnits(units, aggregation)` if we ever need combined trace or dwell from several units.

---

## 3. Order of work (later)

1. **RNA parcel type**: define struct with `units`, `nalleles`, `yield` (and existing RNA fields).
2. **chemical_master.jl**:
   - Joint steady state for coupled units.
   - Marginals from joint.
   - Convolution loop: per unit, (kᵢ − 1) 1D convolutions with Pᵢ.
   - Aggregation: form Y and P(Y = n) from convolved joint.
3. **Wire RNA prediction**: call the new pipeline from RNA parcel(s); deprecate or replace old nalleles/units wiring.
4. **Simulator and I/O**: use RNA parcel(s) so multi-unit multi-allele runs match the new prediction.
