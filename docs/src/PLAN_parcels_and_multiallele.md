# Plan: Parcel redesign and multi-unit multi-allele RNA

## 1. Current scope (single-entity trace and dwelltime)

### 1.1 Trace parcel (single entity)
- **unit**: single unit index
- **emission**: reporter sum, specific R step, G state, or small struct
- **frame_interval**, **start_frame**, **end_frame**
- Optional: **active_fraction**, **background_mean**, **reporter** / output index

### 1.2 Dwell-time parcel (single entity)
- **unit**: single unit index
- **state_spec**: `(kind = :G, state = k)` or `(kind = :R, step = k)` (same idea as trace)
- **sojourn_kind**: `:in_state` | `:out_of_state` (replaces ON/OFF/ONG/OFFG)
- **bins**: histogram bins for this dwell type
- Optional: output index

### 1.3 Implementation tasks (now)
- Define structs/named tuples for trace and dwell parcels.
- Replace or wrap current onstates / bins / dttype with lists of these parcels for trace and dwell.
- Simulator: loop over parcels; use `parcel.unit` to index state; use parcel fields for emission, state_spec, sojourn_kind, bins, frames.
- Likelihood / prediction: same parcel-driven logic for trace and dwell.
- Hidden units: no parcel for that unit ⇒ no indexing in trace/dwell.
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
