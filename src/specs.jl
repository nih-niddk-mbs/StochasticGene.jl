# This file is part of StochasticGene.jl
# Measurement specs: TraceSpec and DwellSpec for trace and dwell-time observables.
# Specs are the canonical way to define which units are observed and how (supports hidden units).

# --------------- Region selectors (sigma-algebra elements) ---------------

abstract type RegionSelector end

"""
    AllStates()

Region selector meaning "all states" in the given space (:R, :G, :T, ...).

For `space = :R`, this is interpreted consistently with existing code as
"all R states starting at insertstep" when passed through the current
`num_reporters_per_state` / dwell-time helpers.
"""
struct AllStates <: RegionSelector
end

"""
    ExplicitStates(states)

Region selector for an explicit set of state indices in the given space.
"""
struct ExplicitStates <: RegionSelector
    states::Vector{Int}
end

"""
    ReporterNonzero()

Region selector meaning "all states where the reporter count is nonzero"
in the given space (typically :R or :T). Interpretation is delegated to the
existing reporter helper functions so semantics stay consistent.
"""
struct ReporterNonzero <: RegionSelector
end

function _infer_trace_region(emission)
    if emission === :reporter_sum
        return :R, AllStates()
    elseif emission isa NamedTuple && hasproperty(emission, :kind) && hasproperty(emission, :state) && emission.kind === :G
        return :G, ExplicitStates([Int(emission.state)])
    else
        throw(ArgumentError("Unsupported TraceSpec emission: $(emission)"))
    end
end

function _infer_dwell_space_side(dttype::Symbol)
    dttype === :OFF  && return (:R, :out_of_region)
    dttype === :ON   && return (:R, :in_region)
    dttype === :OFFG && return (:G, :out_of_region)
    dttype === :ONG  && return (:G, :in_region)
    throw(ArgumentError("Unknown dwell dttype: $(dttype); expected :OFF, :ON, :OFFG, or :ONG"))
end

"""
    TraceSpec

Specification for a single trace (intensity time-series) measurement.

# Fields
- `unit::Int`: Which unit this trace observes (1-based index).
- `emission`: What is reported. `:reporter_sum` (default) = sum over R steps from insertstep;
  or `(kind=:G, state=k)` for a single G state k.
- `frame_interval::Float64`: Time between frames (minutes).
- `start_frame::Float64`: Start time/frame (minutes).
- `end_frame`: End time/frame (minutes); use -1 for last index).
- `active_fraction::Float64`: Fraction of observed active traces (default 1.0).
- `background_mean`: Optional background mean (scalar or vector per channel).
- `output_index::Union{Int,Nothing}`: Optional index for ordering (default nothing).

One trace channel per `TraceSpec`; order of specs = order of columns in trace data.
"""
struct TraceSpec
    unit::Int
    space::Symbol
    region::RegionSelector
    emission::Union{Symbol,NamedTuple}
    frame_interval::Float64
    start_frame::Float64
    end_frame::Union{Float64,Int}
    active_fraction::Float64
    background_mean::Any  # Nothing, scalar, or vector
    output_index::Union{Int,Nothing}
    function TraceSpec(unit::Int, emission, frame_interval::Real, start_frame::Real, end_frame;
                       active_fraction::Real=1.0, background_mean=nothing, output_index=nothing)
        space, region = _infer_trace_region(emission)
        new(unit, space, region, emission, Float64(frame_interval), Float64(start_frame), end_frame,
            Float64(active_fraction), background_mean, output_index)
    end
end

"""
    DwellSpec

Specification for a single dwell-time (sojourn) histogram measurement.

# Fields
- `unit::Int`: Which unit this dwell type refers to (1-based index).
- `dttype::Symbol`: Dwell type: `:OFF`, `:ON` (R space), `:OFFG`, `:ONG` (G space).
- `onstates::Vector{Int}`: The set of states considered "on" for this histogram (empty for R = any step).
- `bins::Vector{Float64}`: Histogram bin edges for this dwell type.
- `output_index::Union{Int,Nothing}`: Optional index for ordering (default nothing).
"""
struct DwellSpec
    unit::Int
    space::Symbol
    region::RegionSelector
    side::Symbol
    dttype::Symbol
    onstates::Vector{Int}
    bins::Vector{Float64}
    output_index::Union{Int,Nothing}
    function DwellSpec(unit::Int, dttype::Symbol, onstates::Vector{Int}, bins::Vector{<:Real};
                       output_index=nothing)
        (dttype === :OFF || dttype === :ON || dttype === :OFFG || dttype === :ONG) ||
            throw(ArgumentError("dttype must be :OFF, :ON, :OFFG, or :ONG"))
        space, side = _infer_dwell_space_side(dttype)
        region = ExplicitStates(Int.(onstates))
        new(unit, space, region, side, dttype, Int.(onstates), Float64.(bins), output_index)
    end
end

# --------------- Legacy → Spec conversion -----------------

"""
    trace_specs_from_legacy(onstates, traceinfo, G, R, S, insertstep; unit=1)

Build a vector of `TraceSpec` from legacy `onstates` and `traceinfo`.

- If `onstates` is empty (or Int[]): one spec with `emission = :reporter_sum`.
- If `onstates` is a vector of G states (e.g. [2, 3]): one spec per state with `emission = (kind=:G, state=k)`.
- `traceinfo`: (frame_interval, start_frame, end_frame) or (..., active_fraction, background_mean).
"""
function trace_specs_from_legacy(onstates, traceinfo, G, R, S, insertstep; unit::Int=1)
    fi = Float64(traceinfo[1])
    sf = Float64(traceinfo[2])
    ef = traceinfo[3]
    af = length(traceinfo) >= 4 ? Float64(traceinfo[4]) : 1.0
    bg = length(traceinfo) >= 5 ? traceinfo[5] : nothing
    specs = TraceSpec[]
    if isempty(onstates)
        push!(specs, TraceSpec(unit, :reporter_sum, fi, sf, ef; active_fraction=af, background_mean=bg))
    else
        for (i, k) in enumerate(onstates)
            push!(specs, TraceSpec(unit, (kind=:G, state=Int(k)), fi, sf, ef;
                                  active_fraction=af, background_mean=bg, output_index=i))
        end
    end
    specs
end

"""
    trace_specs_from_legacy(onstates, traceinfo, G::Tuple, R, S, insertstep, coupling; units=nothing)

Coupled case: build trace specs for each unit. `onstates` is Vector{Vector{Int}} (per unit).
If `units` is given, only those unit indices get specs (hidden units omitted).
"""
function trace_specs_from_legacy(onstates::Vector{Vector{Int}}, traceinfo, G::Tuple, R, S, insertstep, coupling;
                                 units=nothing)
    nunits = length(onstates)
    units_use = units === nothing ? (1:nunits) : units
    specs = TraceSpec[]
    for u in units_use
        i = u
        fi = Float64(traceinfo[1])
        sf = Float64(traceinfo[2])
        ef = traceinfo[3]
        af = length(traceinfo) >= 4 ? (traceinfo[4] isa AbstractVector ? traceinfo[4][i] : traceinfo[4]) : 1.0
        bg = length(traceinfo) >= 5 ? (traceinfo[5] isa AbstractVector ? traceinfo[5][i] : traceinfo[5]) : nothing
        ons = i <= length(onstates) ? onstates[i] : Int[]
        if isempty(ons)
            push!(specs, TraceSpec(i, :reporter_sum, fi, sf, ef; active_fraction=af, background_mean=bg))
        else
            for (j, k) in enumerate(ons)
                push!(specs, TraceSpec(i, (kind=:G, state=Int(k)), fi, sf, ef;
                                      active_fraction=af, background_mean=bg, output_index=j))
            end
        end
    end
    specs
end

"""
    dwell_specs_from_legacy(onstates, bins, dttype; unit=1)

Build a vector of `DwellSpec` from legacy `onstates`, `bins`, and `dttype`.

- `onstates`: Vector of "on" state sets (one per dwell type); empty means all R steps / all G.
- `bins`: Vector of bin vectors (one per dwell type).
- `dttype`: e.g. ["OFF", "ON", "OFFG", "ONG"].
"""
function dwell_specs_from_legacy(onstates, bins, dttype; unit::Int=1)
    length(onstates) == length(bins) == length(dttype) ||
        throw(ArgumentError("onstates, bins, dttype must have same length"))
    specs = DwellSpec[]
    for i in eachindex(dttype)
        dt = dttype[i]
        ons = eltype(onstates) <: Vector ? onstates[i] : copy(onstates)
        b = eltype(bins) <: Vector && eltype(bins) != Float64 ? bins[i] : bins
        if dt == "OFF"
            push!(specs, DwellSpec(unit, :OFF, ons, b; output_index=i))
        elseif dt == "ON"
            push!(specs, DwellSpec(unit, :ON, ons, b; output_index=i))
        elseif dt == "OFFG"
            push!(specs, DwellSpec(unit, :OFFG, ons, b; output_index=i))
        elseif dt == "ONG"
            push!(specs, DwellSpec(unit, :ONG, ons, b; output_index=i))
        else
            throw(ArgumentError("Unknown dttype: $dt"))
        end
    end
    specs
end

"""
    dwell_specs_from_legacy(onstates::Vector{Vector{Vector{Int}}}, bins, dttype, G::Tuple, R, S, insertstep)

Coupled case: onstates is Vector{Vector{Vector{Int}}} (per unit, per dwell type). Returns flat list of
DwellSpec with correct unit index.
"""
function dwell_specs_from_legacy(onstates::Vector{Vector{Vector{Int}}}, bins, dttype, G::Tuple, R, S, insertstep)
    specs = DwellSpec[]
    for i in eachindex(onstates)
        # bins and dttype for coupled: bins[i], dttype[i] are the vectors for unit i
        bi = bins[i]
        dti = dttype[i]
        length(onstates[i]) == length(bi) == length(dti) ||
            throw(ArgumentError("Per-unit onstates, bins, dttype length mismatch"))
        for j in eachindex(dti)
            dt = dti[j]
            ons = onstates[i][j]
            b = bi[j]
            if dt == "OFF"
                push!(specs, DwellSpec(i, :OFF, ons, b; output_index=j))
            elseif dt == "ON"
                push!(specs, DwellSpec(i, :ON, ons, b; output_index=j))
            elseif dt == "OFFG"
                push!(specs, DwellSpec(i, :OFFG, ons, b; output_index=j))
            elseif dt == "ONG"
                push!(specs, DwellSpec(i, :ONG, ons, b; output_index=j))
            else
                throw(ArgumentError("Unknown dttype: $dt"))
            end
        end
    end
    specs
end

# --------------- Spec → legacy (for incremental refactor) -----------------

"""
    legacy_onstates_traceinfo(trace_specs::Vector{TraceSpec}; nunits=1)

Return (onstates, traceinfo) equivalent to the given trace specs.
- Single unit (nunits=1): onstates = Int[] or Vector{Int} of G states; traceinfo = (frame_interval, start, end, af, bg).
- Coupled (nunits>1): onstates = Vector{Vector{Int}} with one entry per unit (empty = no trace for that unit);
  traceinfo from first spec. Units with no specs get empty onstates.
"""
function legacy_onstates_traceinfo(trace_specs::Vector{TraceSpec}; nunits::Int=1)
    isempty(trace_specs) && return nunits == 1 ? (Int[], (1.0, 0.0, -1)) : (fill(Int[], nunits), (1.0, 0.0, -1))
    s = trace_specs[1]
    traceinfo = (s.frame_interval, s.start_frame, s.end_frame, s.active_fraction, s.background_mean)
    if nunits == 1
        onstates = Int[]
        for spec in trace_specs
            spec.unit != 1 && continue
            if spec.emission === :reporter_sum
                return Int[], traceinfo
            elseif spec.emission isa NamedTuple && spec.emission.kind === :G
                push!(onstates, spec.emission.state)
            end
        end
        return onstates, traceinfo
    else
        onstates_per_unit = [Int[] for _ in 1:nunits]
        for spec in trace_specs
            u = spec.unit
            (u < 1 || u > nunits) && continue
            if spec.emission === :reporter_sum
                # single reporter-sum for this unit; leave others as-is
                onstates_per_unit[u] = Int[]
                break
            elseif spec.emission isa NamedTuple && spec.emission.kind === :G
                push!(onstates_per_unit[u], spec.emission.state)
            end
        end
        return onstates_per_unit, traceinfo
    end
end

"""
    legacy_dwell(dwell_specs::Vector{DwellSpec})

Return (onstates, bins, dttype) equivalent to the given dwell specs (single unit).
"""
function legacy_dwell(dwell_specs::Vector{DwellSpec})
    onstates = Vector{Int}[]
    bins = Vector{Float64}[]
    dttype = String[]
    for s in dwell_specs
        push!(onstates, s.onstates)
        push!(bins, s.bins)
        push!(dttype, string(s.dttype))
    end
    onstates, bins, dttype
end

"""
    legacy_dwell_coupled(dwell_specs::Vector{DwellSpec}, nunits::Int)

Return (onstates, bins, dttype) in the coupled-format expected by `set_onoff`:

- `onstates::Vector{Vector{Vector{Int}}}`: per-unit, per-dwell-type ON-state sets
- `bins::Vector{Vector{Vector{Float64}}}`: per-unit, per-dwell-type bin edges
- `dttype::Vector{Vector{String}}`: per-unit, per-dwell-type dwell type strings
"""
function legacy_dwell_coupled(dwell_specs::Vector{DwellSpec}, nunits::Int)
    onstates = [Vector{Vector{Int}}() for _ in 1:nunits]
    bins = [Vector{Vector{Float64}}() for _ in 1:nunits]
    dttype = [String[] for _ in 1:nunits]
    for s in dwell_specs
        u = s.unit
        (u < 1 || u > nunits) && throw(ArgumentError("DwellSpec unit $(u) out of range for $(nunits) units"))
        push!(onstates[u], s.onstates)
        push!(bins[u], s.bins)
        push!(dttype[u], string(s.dttype))
    end
    return onstates, bins, dttype
end

# --------------- Build HMMReporter from TraceSpec (coupled; one channel per spec) ---------------

"""
    _onstates_for_trace_spec(spec::TraceSpec) -> Vector{Int}

Return the legacy onstates vector for this spec's emission (:reporter_sum => Int[], (kind=:G, state=k) => [k]).
"""
function _onstates_for_trace_spec(spec::TraceSpec)
    if spec.emission === :reporter_sum
        return Int[]
    elseif spec.emission isa NamedTuple && spec.emission.kind === :G
        return [spec.emission.state]
    else
        throw(ArgumentError("Unsupported TraceSpec emission: $(spec.emission)"))
    end
end

"""
    hmm_reporter_from_trace_spec(spec::TraceSpec, transitions::Tuple, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, probfn, noisepriors)

Build one HMMReporter for the given TraceSpec in a coupled model (G,R,S,etc. are per-unit tuples).
Used when building reporter from trace_specs to support hidden units (fewer reporter channels than units).
"""
function hmm_reporter_from_trace_spec(spec::TraceSpec, transitions::Tuple, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, probfn, noisepriors)
    i = spec.unit
    (i < 1 || i > length(G)) && throw(ArgumentError("TraceSpec unit $(i) out of range for $(length(G)) units"))
    onstates = _onstates_for_trace_spec(spec)
    Gu, Ru, Su, insertstepu = G[i], R[i], S[i], insertstep[i]
    trans = transitions[i]
    n = num_rates(trans, Ru, Su, insertstepu)
    nnoise = length(noisepriors[i])
    per_state = num_reporters_per_state(Gu, Ru, Su, insertstepu, onstates)
    offs = off_states(Gu, Ru, Su, insertstepu, onstates)
    pfn = probfn isa Union{Tuple,Vector} ? probfn[i] : probfn
    weightind = occursin("Mixture", "$(pfn)") ? n + nnoise : 0
    noiseparams = collect(n+1:n+nnoise)
    HMMReporter(nnoise, per_state, pfn, weightind, offs, noiseparams)
end

# --------------- Spec ↔ Dict (for TOML I/O) ---------------

function _trace_spec_to_dict(spec::TraceSpec)
    em = spec.emission === :reporter_sum ? "reporter_sum" : Dict("kind" => string(spec.emission.kind), "state" => spec.emission.state)
    Dict{String,Any}(
        "unit" => spec.unit,
        "emission" => em,
        "frame_interval" => spec.frame_interval,
        "start_frame" => spec.start_frame,
        "end_frame" => spec.end_frame,
        "active_fraction" => spec.active_fraction,
        "background_mean" => spec.background_mean,
        "output_index" => spec.output_index,
    )
end

function _dict_to_trace_spec(d::Dict)
    em = d["emission"]
    emission = (em isa String && em == "reporter_sum") ? :reporter_sum : (kind=Symbol(em["kind"]), state=Int(em["state"]))
    bg = get(d, "background_mean", nothing)
    (bg isa String && bg == "nothing") && (bg = nothing)
    oi = get(d, "output_index", nothing)
    (oi isa String && oi == "nothing") && (oi = nothing)
    oi !== nothing && (oi = Int(oi))
    TraceSpec(
        Int(d["unit"]),
        emission,
        Float64(d["frame_interval"]),
        Float64(d["start_frame"]),
        d["end_frame"] isa Int ? d["end_frame"] : Int(d["end_frame"]);
        active_fraction=Float64(d["active_fraction"]),
        background_mean=bg,
        output_index=oi,
    )
end

function _dwell_spec_to_dict(spec::DwellSpec)
    Dict{String,Any}(
        "unit" => spec.unit,
        "dttype" => string(spec.dttype),
        "onstates" => spec.onstates,
        "bins" => spec.bins,
        "output_index" => spec.output_index,
    )
end

function _dict_to_dwell_spec(d::Dict)
    oi = get(d, "output_index", nothing)
    (oi isa String && oi == "nothing") && (oi = nothing)
    oi !== nothing && (oi = Int(oi))
    dt = haskey(d, "dttype") ? Symbol(d["dttype"]) : _legacy_dttype_from_kind(d)
    DwellSpec(
        Int(d["unit"]),
        dt,
        Int.(d["onstates"]),
        Float64.(d["bins"]);
        output_index=oi,
    )
end

function _legacy_dttype_from_kind(d::Dict)
    k = Symbol(get(d, "kind", "R"))
    s = Symbol(get(d, "sojourn_kind", "out_of_state"))
    (k === :R && s === :out_of_state) && return :OFF
    (k === :R && s === :in_state) && return :ON
    (k === :G && s === :out_of_state) && return :OFFG
    (k === :G && s === :in_state) && return :ONG
    error("Invalid legacy kind/sojourn_kind in dwell_spec")
end
