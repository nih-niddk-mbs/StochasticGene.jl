# Multi-Connection Coupling Design

Design for extending the coupling model to support multiple (source_state, target_transition) pairs per (source_unit, target_unit) direction.

**Current scope:** 4-character reciprocal coupling (models 1–8, format "s1t1s2t2") is supported. Models 9–12 (extended multi-connection format) will be implemented later.

---

## Source vs target (don’t get it backwards)

- **Target unit α:** The unit whose dynamics we are updating. Its transition rates are modulated by other units. When we are “focused on unit α,” we are building/updating the rates for **α**.
- **Source unit β:** The unit whose **state** we read. That state (e.g. “is β in state s?”) multiplies a coupling constant and is added into the **target**’s transition rate.

So for a connection **(β, α, s, t)**:
- **β = source** (influencer; we read β’s state).
- **α = target** (influenced; we modify α’s transition **t** using β’s state **s**).

When focused on **unit α**, we sum over **all sources to α**: i.e. over all **β** in `sources[α]`. The 5-tuple is indexed by **target**: `sources[α]` = list of units β that couple **into** α.

---

## 0. Design principles for models 9–12

1. **(α, β, s, t) everywhere:** Implement multiple source states **s** and target transitions **t** for a given direction (source unit β → target unit α). Each connection is (β, α, s, t) with its own coupling constant γ.

2. **Pairing by index:** For each target unit α, list all connections into α as parallel tuples: `sources[α]`, `source_state[α]`, `target_transition[α]` have the same length and are aligned by index. If one source state couples to multiple target transitions, repeat that source (and its s) so each (s, t) is its own entry. One coupling constant per connection; coupling constants are ordered to align with this list.

3. **Fit function input:** The `coupling` argument to `fit` must be able to represent this extended format (see §2). We need a clear way to construct and pass it for models 9–12.

4. **File coupling model identifier:** We need a way to encode the coupling topology in the filename (e.g. the condition/label field that currently holds "31" or "3131") so that rate files and downstream code (e.g. `write_correlation_functions`, `write_traces`) can parse it. Example: extended string for 9–12 such as "24,33|33" or similar.

5. **Field 1 – unit_model:** The first element of the coupling tuple identifies which **model** (transition structure, G, R, S, etc.) corresponds to which **unit**. Example: `(1, 2)` = unit 1 uses model 1, unit 2 uses model 2. This allows "repeated" models: e.g. `(1, 1)` means both units use the same model type, and parameters can overlap or be shared as needed.

6. **Source state numbering:** Use a single, consistent numbering scheme for source states (e.g. G states 1..G, then R states G+1..G+R per unit, or a global scheme) so that **s** and **t** are unambiguous when building U(s) and V(t) and when specifying coupling in the fit and in file identifiers.

---

## 1. Models 9–12 (require multiple connections; planned)

| Model | Connections |
|-------|-------------|
| 9 | (1→2): state 2→rate 4, state 3→rate 3; (2→1): state 3→rate 3 |
| 10 | (1→2): state 2→rate 4, state 3→rate 5; (2→1): state 3→rate 5 |
| 11 | (1→2): state 2→rate 3, state 3→rate 3; (2→1): state 3→rate 3 |
| 12 | (1→2): state 2→rate 3, state 3→rate 5; (2→1): state 3→rate 5 |

Each connection has its own coupling strength γ in the rate vector.

---

## 2. Extended 5-Tuple Format (backward compatible)

**Legacy** (one connection per (β,α)):
```
(unit_model, sources, source_state, target_transition, ncoupling)
sources[α] = (β,)           # single source
source_state[α] = s         # scalar or (s,)
target_transition[α] = t    # scalar or (t,)
```

**Extended** (multiple connections per target):
```
sources[α] = (β1, β2, ...)           # can repeat β
source_state[α] = (s1, s2, ...)      # same length as sources[α]
target_transition[α] = (t1, t2, ...) # same length as sources[α]
```

**Invariant:** For each target α, `sources[α]`, `source_state[α]`, and `target_transition[α]` must have the same length and be pairwise aligned. If one source state couples to multiple target transitions, list it once per target (i.e. repeat the source/source_state so each connection is its own entry). Example: unit 1 state 2→trans 4 and state 3→trans 3: `sources[2]=(1,1)`, `source_state[2]=(2,3)`, `target_transition[2]=(4,3)`.

Canonical order of coupling params:
```
[(β, α, s, t) for α in 1:n for k in 1:length(sources[α])
    with β=sources[α][k], s=source_state[α][k], t=target_transition[α][k]]
```

---

## 3. Connections List (internal representation)

```
connections = [(β, α, source_state, target_transition), ...]
ncoupling = length(connections)
```

- `source_state` can be `Int` (single G state) or `Vector{Int}` (e.g. R states)
- `target_transition` is `Int`
- k-th coupling param = k-th element of `connections`

---

## 4. Code Changes Required

| Component | Change |
|-----------|--------|
| **io.jl** | `to_connections(coupling, G, R)` normalizes to flat connection list; `make_coupling_extended` for models 9–12; **file identifier**: extended format (e.g. "24,33\|33") parsed in `extract_source_target` (analysis.jl) and `make_coupling` / `make_coupling_extended` (io.jl); `parse_filename` already passes through full condition suffix as `coupling_field` |
| **analysis.jl** | **write_traces** / **write_correlation_functions**: `extract_source_target` returns extended coupling string when filename segment contains `,` or `\|`; `make_coupling` builds extended 5-tuple via `make_coupling_extended`; `write_trace_dataframe` uses full `coupling_field` as datapath subpath when extended (else first character for legacy) |
| **fit.jl** | `coupling` argument accepts extended 5-tuple; **input**: document or helper to build (unit_model, sources, source_state, target_transition, ncoupling) for 9–12; `coupling_indices` uses connection order so γ aligns 1:1 with connections |
| **transition_rate_make.jl** | `make_mat_TC` loops over all connection indices; for each (β, α, s, t) build U_β(s), V_α(t) and add γ[k]*term |
| **transition_rate_structures.jl** | Components stay in matrix coords; optionally store what’s needed to build U(s), V(t) per connection |
| **simulator_coupled.jl** | `prepare_rates_sim`, `update_coupling!` iterate over connection order |
| **likelihoods.jl** | `prepare_rates_coupled` uses connection order for γ |
| **Source state numbering** | Use one consistent scheme (e.g. G states 1..G, R states G+1..G+R per unit) everywhere: coupling spec, file identifier, and when building U/V |

---

## 5. Building U(β,s) and V(α,t) per Connection

For each connection (β, α, s, t):
1. Get unit β transitions, G, R, S, insertstep
2. `elementsSource = set_elements_Source(s, G, R, S)`; `U = make_mat_S(elementsSource, nT)`
3. Get unit α transitions, G, R, S, insertstep  
4. `elementsTarget = set_elements_Target(t, transitions, G, R, S, insertstep, indices, nT, splicetype)`; `V = make_mat(elementsTarget, rates, nT)`
5. Build coupling term with U, V in correct Kronecker order
6. Add γ[k] * term to Tc

---

## 6. Backward Compatibility

- Legacy 5-tuple (one connection per (β,α)) unchanged
- `to_connections()` converts both legacy and extended → flat list
- All consumers use `connections`; legacy format is converted at boundaries
