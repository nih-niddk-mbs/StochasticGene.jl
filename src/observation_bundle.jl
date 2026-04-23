# observation_bundle.jl — multi-modality experimental data container
#
# Legacy combined structs (e.g. TraceRNAData) are refactored toward separate
# per-modality payloads held in an ObservationBundle. Likelihood code paths
# that require joint dispatch (trace + RNA FISH) reassemble a TraceRNAData
# representative until HMM entry points accept decoupled trace views.

"""
    AbstractObservationBundle <: AbstractExperimentalData

Supertype for **joint experimental data**: several co-registered modalities, each
stored as its own [`AbstractExperimentalData`](@ref) leg (e.g. live-cell trace +
RNA FISH histogram).

See also [`ObservationBundle`](@ref).
"""
abstract type AbstractObservationBundle <: AbstractExperimentalData end

"""
    ObservationBundle(legs::NamedTuple)

Container holding named modalities. Typical keys: `:rna`, `:trace`, `:dwell`.
Use `nothing` for absent legs, e.g. `(rna=RNAData(...), trace=TraceData(...), dwell=nothing)`.

# Construction

- [`ObservationBundle(::TraceRNAData)`](@ref) splits trace + histogram into
  [`TraceData`](@ref) and [`RNAData`](@ref) legs (same on-disk semantics as `tracerna`).

# Likelihood

[`loglikelihood`](@ref) / [`loglikelihood_ad`](@ref) on `ObservationBundle` currently
supports the **trace + RNA** joint shape by delegating to the existing
[`TraceRNAData`](@ref) implementation (reconstructed from legs). Other shapes
throw until independent-leg joint likelihoods are implemented.

# Note

[`RNADwellTimeData`](@ref) is **not** split here: its histogram likelihood uses a
single joint `predictedarray` / MT stack; splitting would change the math until
that path is refactored.
"""
struct ObservationBundle{N<:NamedTuple} <: AbstractObservationBundle
    legs::N
end

ObservationBundle(pairs::Vararg{Pair{Symbol,<:Any}}) = ObservationBundle(NamedTuple(pairs))

"""Return modality keys with non-`nothing` legs (declaration order)."""
function observation_modalities(b::ObservationBundle)
    Tuple(k for (k, v) in pairs(b.legs) if v !== nothing)
end

function ObservationBundle(d::TraceRNAData)
    td = TraceData(d.label, d.gene, d.interval, d.trace, d.units)
    rna = RNAData(d.label, d.gene, d.nRNA, d.histRNA, d.yield, d.units)
    ObservationBundle((rna=rna, trace=td, dwell=nothing))
end

function _is_tracerna_shaped_bundle(b::ObservationBundle)
    legs = b.legs
    hasproperty(legs, :rna) || return false
    hasproperty(legs, :trace) || return false
    legs.rna isa AbstractRNAData || return false
    legs.trace isa AbstractTraceData || return false
    if hasproperty(legs, :dwell) && legs.dwell !== nothing
        return false
    end
    return true
end

"""
    reconstruct_tracerna(b::ObservationBundle) -> TraceRNAData

Rebuild a [`TraceRNAData`](@ref) from `:rna` + `:trace` legs (inverse of
[`ObservationBundle(::TraceRNAData)`](@ref)). Used internally for likelihood dispatch.
"""
function reconstruct_tracerna(b::ObservationBundle)
    _is_tracerna_shaped_bundle(b) || throw(ArgumentError("ObservationBundle is not a trace+RNA (tracerna) shape; got modalities $(observation_modalities(b))"))
    r = b.legs.rna
    t = b.legs.trace
    return TraceRNAData(t.label, t.gene, t.interval, t.trace, r.nRNA, r.histRNA, r.yield, t.units)
end
