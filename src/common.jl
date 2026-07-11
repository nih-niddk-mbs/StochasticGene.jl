# This file is part of StochasticGene.jl   

###  common.jl

# Data types
"""
    AbstractExperimentalData

Abstract type for experimental data
"""
abstract type AbstractExperimentalData end
"""
    AbstractSampleData

Abstract type for data in the form of samples
"""
abstract type AbstractSampleData <: AbstractExperimentalData end
"""
    AbstractHistogramData

Abstract type for data in the form of a histogram (probability distribution)
"""
abstract type AbstractHistogramData <: AbstractExperimentalData end
"""
    AbstractTraceData
    
Abstract type for intensity time series data
"""
abstract type AbstractTraceData <: AbstractExperimentalData end

"""
    AbstractRNAData{hType}

Abstract type for steady state RNA histogram data.
"""
abstract type AbstractRNAData{hType} <: AbstractHistogramData end


"""
    AbstractTraceHistogramData
    
Abstract type for intensity time series data with RNA histogram data
"""

abstract type AbstractTraceHistogramData <: AbstractTraceData end

"""
    HMM_STACK_MH

Likelihood **MH / MCMC track**: in-place `forward`, `kolmogorov_forward` (`fkf!`). This is the default for
[`loglikelihood`](@ref) and [`ll_hmm_trace`](@ref) when you do not pass `hmm_stack`.

See also [`HMM_STACK_AD`](@ref).
"""
const HMM_STACK_MH = :fast

"""
    HMM_STACK_AD

Likelihood **AD track** (NUTS, ADVI, Zygote): [`forward_ad`](@ref), [`kolmogorov_forward_ad`](@ref), etc.
Default for [`loglikelihood_ad`](@ref) on trace data. Use the same tag with `ll_hmm_trace(...; hmm_stack=HMM_STACK_AD)`.

See also [`HMM_STACK_MH`](@ref).
"""
const HMM_STACK_AD = :zygote

"""
    set_hmm_zygote_checkpoint_steps!(n::Union{Nothing,Integer})
    hmm_zygote_checkpoint_steps()

Configure **gradient checkpointing** for reverse-mode AD (Zygote) through the trace HMM forward pass
([`forward_ad`](@ref) / [`forward_grid_ad`](@ref)): when `n` is a positive integer, time steps are processed in
chunks of `n` using [`Zygote.checkpointed`](@ref), trading extra forward work for lower memory during the
backward pass. Set to `nothing` to disable (default).

For **multiple discrete traces**, [`loglikelihood_ad`](@ref) / [`ll_hmm_trace`](@ref) with [`HMM_STACK_AD`](@ref)
also applies [`Zygote.checkpointed`](@ref) **per trace** so reverse-mode does not retain one tape over all traces
(only the per-trace forward pass is recomputed on the backward pass).

The same frame setting is preferred on inference **options** as
[`NUTSOptions`](@ref) / [`ADVIOptions`](@ref) field `gradient_checkpoint_length` (run-spec keys
`gradient_checkpoint_length` or legacy `hmm_checkpoint_steps`). Per-call `hmm_checkpoint_steps` on
[`run_nuts`](@ref), [`ll_hmm_trace`](@ref), or [`loglikelihood_ad`](@ref) still overrides for that call.

Does not affect [`ForwardDiff`](@ref) or finite-difference gradients.
"""
const HMM_ZYGOTE_CHECKPOINT_STEPS = Ref{Union{Nothing,Int}}(nothing)

function set_hmm_zygote_checkpoint_steps!(n::Union{Nothing,Integer})
    HMM_ZYGOTE_CHECKPOINT_STEPS[] = n === nothing ? nothing : Int(n)
end

hmm_zygote_checkpoint_steps() = HMM_ZYGOTE_CHECKPOINT_STEPS[]

function with_hmm_zygote_checkpoint(f::Function, ck::Union{Nothing,Integer})
    ck === nothing && return f()
    prev = HMM_ZYGOTE_CHECKPOINT_STEPS[]
    HMM_ZYGOTE_CHECKPOINT_STEPS[] = Int(ck)
    try
        return f()
    finally
        HMM_ZYGOTE_CHECKPOINT_STEPS[] = prev
    end
end

# Data structures 
#
# Do not use underscore "_" in label

"""
    RNAData{nType,hType}

Structure for storing RNA histogram data.

# Fields
- `label`: Label for the data set.
- `gene`: Gene name (case sensitive).
- `nRNA`: Length of the histogram (type varies).
- `histRNA`: RNA histograms (type varies).
- `yield`: Detection efficiency (Float64 when = 1.0) or (yield, nRNA_true) tuple when < 1.0.
"""
struct RNAData{nType,hType} <: AbstractRNAData{hType}
    label::String
    gene::String
    nRNA::nType
    histRNA::hType
    yield::Union{Float64, Tuple{Float64, Int}}
    units::Vector{Int}
end
RNAData(label, gene, nRNA, histRNA, yield) = RNAData(label, gene, nRNA, histRNA, yield, [1])

struct RNACountData <: AbstractRNAData{Vector{Int}}
    label::String
    gene::String
    nRNA::Int
    countsRNA::Vector{Int}
    yieldfactor::Vector{Float64}
    units::Vector{Int}
end
RNACountData(label, gene, nRNA, countsRNA, yieldfactor) = RNACountData(label, gene, nRNA, countsRNA, yieldfactor, [1])

struct DwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    bins::Vector
    DwellTimes::Vector
    DTtypes::Vector
    units::Vector{Int}
end
DwellTimeData(label, gene, bins, DwellTimes, DTtypes) = DwellTimeData(label, gene, bins, DwellTimes, DTtypes, [1])

"""
    RNAOnOffData

Structure for storing RNA ON/OFF time data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::Int`: Length of the histogram.
- `histRNA::Vector`: RNA histograms.
- `bins::Vector`: Number of live cell recording time bins.
- `ON::Vector`: ON time probability density.
- `OFF::Vector`: OFF time probability density.
"""
struct RNAOnOffData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Vector
    bins::Vector
    ON::Vector
    OFF::Vector
    yield::Union{Float64, Tuple{Float64, Int}}
    units::Vector{Int}
end
RNAOnOffData(label, gene, nRNA, histRNA, bins, ON, OFF, yield) = RNAOnOffData(label, gene, nRNA, histRNA, bins, ON, OFF, yield, [1])
"""
    RNADwellTimeData

Structure for storing RNA dwell time data.

# Fields
- `label::String`: Label for the data set.
- `gene::String`: Gene name (case sensitive).
- `nRNA::Int`: Length of the histogram.
- `histRNA::Array`: RNA histograms.
- `bins::Vector{Vector}`: Number of live cell recording time bins.
- `DwellTimes::Vector{Vector}`: Dwell times.
- `DTtypes::Vector`: Types of dwell times.
"""
struct RNADwellTimeData <: AbstractHistogramData
    label::String
    gene::String
    nRNA::Int
    histRNA::Array
    bins::Vector{Vector}
    DwellTimes::Vector{Vector}
    DTtypes::Vector
    yield::Union{Float64, Tuple{Float64, Int}}
    units::Vector{Int}
end
RNADwellTimeData(label, gene, nRNA, histRNA, bins, DwellTimes, DTtypes, yield) = RNADwellTimeData(label, gene, nRNA, histRNA, bins, DwellTimes, DTtypes, yield, [1])
"""
    TraceData{labelType,geneType,traceType}

Structure for storing trace data.

# Fields
- `label::labelType`: Label for the data set.
- `gene::geneType`: Gene name (case sensitive).
- `interval::Float64`: Time interval between trace points.
- `trace::traceType`: Trace data.
- `units::Vector{Int}`: Unit index per observation (trace column); empty means legacy 1:1 (observation i = unit i).
"""
struct TraceData{labelType,geneType,traceType} <: AbstractTraceData
    label::labelType
    gene::geneType
    interval::Float64
    trace::traceType
    units::Vector{Int}
end
# Legacy constructor: omit units (default empty = 1:1 mapping).
TraceData(label, gene, interval, trace) = TraceData(label, gene, interval, trace, Int[])
"""
    TraceRNAData{traceType,hType}

Joint **trace + RNA FISH histogram** payload (legacy single-struct layout).

# Fields
- `label`: Label for the data set.
- `gene`: Gene name.
- `interval`: Time between trace points.
- `trace`: Trace data (type varies).
- `nRNA`: Histogram length.
- `histRNA`: RNA histogram (type varies).
- `yield`: Detection efficiency (Float64 or (yield, nRNA_true) tuple).
- `units::Vector{Int}`: Unit index per observation; empty means legacy 1:1.

For a **split** representation (separate [`TraceData`](@ref) and [`RNAData`](@ref) legs
in one container), use [`CombinedData`](@ref) and [`CombinedData(::TraceRNAData)`](@ref).
"""
struct TraceRNAData{traceType,hType} <: AbstractTraceHistogramData
    label::String
    gene::String
    interval::Float64
    trace::traceType
    nRNA::Int
    histRNA::hType
    yield::Union{Float64, Tuple{Float64, Int}}
    units::Vector{Int}
end
# Legacy constructor: omit units (default empty = 1:1 mapping).
TraceRNAData(label, gene, interval, trace, nRNA, histRNA, yield) = TraceRNAData(label, gene, interval, trace, nRNA, histRNA, yield, Int[])

# --- Combined experimental data (elementary legs only; one type per modality set) ---
#
# Elementary legs: e.g. TraceData, RNAData, RNACountData, DwellTimeData. Legacy joint products
# (RNAOnOffData, RNADwellTimeData, TraceRNAData, …) are not valid legs; use CombinedData once split.
# CombinedData{NT} dispatches on NamedTuple key type; constructor sorts names lexicographically.

"""
    AbstractCombinedData <: AbstractExperimentalData

Supertype for experimental data formed by **combining** elementary modalities in one object
(see [`CombinedData`](@ref)).
"""
abstract type AbstractCombinedData <: AbstractExperimentalData end

"""
    CombinedData{NT<:NamedTuple} <: AbstractCombinedData

Single container type for **any** combination of elementary modalities. `NT` is a
`NamedTuple` type whose **names** are the combination signature used for Julia dispatch
(names are stored in **lexicographic** order by the constructor).

# Dispatch pattern

```julia
foo(d::CombinedData{NT}) where {NT<:NamedTuple{(:rna, :trace)}} = ...
```

Use `combined_modalities(d)` for a runtime `Tuple` of symbols (same order as `d.legs`).

# Construction

- `CombinedData((; trace=td, rna=rna))` or `CombinedData(trace=td, rna=rna)` — keys are canonicalized.
- [`CombinedData(::TraceRNAData)`](@ref) splits into `:rna` + `:trace` legs (tracerna on-disk shape).

Legs must not be nested [`CombinedData`](@ref), [`TraceRNAData`](@ref), [`AbstractTraceHistogramData`](@ref),
[`RNAOnOffData`](@ref), or [`RNADwellTimeData`](@ref) (use a composition of simpler histogram + trace legs when those loaders exist).
"""
struct CombinedData{NT<:NamedTuple} <: AbstractCombinedData
    legs::NT
    function CombinedData(legs::NamedTuple{names}) where {names}
        c = _canonical_named_tuple(legs)
        isempty(keys(c)) && throw(ArgumentError("CombinedData requires at least one non-empty leg"))
        _assert_elementary_legs(c)
        return new{typeof(c)}(c)
    end
end

CombinedData(pairs::Vararg{Pair{Symbol,<:Any}}) = CombinedData(NamedTuple(pairs))

@generated function _canonical_named_tuple(legs::NamedTuple{names}) where {names}
    syms = sort(collect(names))
    isempty(syms) && return :(NamedTuple())
    perm = Tuple(syms)
    getters = [:(getfield(legs, $(QuoteNode(s)))) for s in perm]
    return :(NamedTuple{$(QuoteNode(perm))}(($(getters...),)))
end

function _assert_elementary_legs(legs::NamedTuple)
    for (k, v) in pairs(legs)
        _assert_elementary_leg(v, k)
    end
    return nothing
end

function _assert_elementary_leg(v, k::Symbol)
    v isa AbstractCombinedData &&
        throw(ArgumentError("CombinedData leg :$k cannot be nested AbstractCombinedData"))
    v isa TraceRNAData &&
        throw(ArgumentError("CombinedData leg :$k cannot be TraceRNAData; split into :rna + :trace elementary legs"))
    v isa AbstractTraceHistogramData &&
        throw(ArgumentError("CombinedData leg :$k cannot be AbstractTraceHistogramData (e.g. TraceRNAData)"))
    v isa RNAOnOffData &&
        throw(ArgumentError("CombinedData leg :$k cannot be RNAOnOffData; represent ON/OFF as separate elementary legs when that layout is supported"))
    v isa RNADwellTimeData &&
        throw(ArgumentError("CombinedData leg :$k cannot be RNADwellTimeData; use a composition of RNA + dwell elementary legs when that layout is supported"))
    if !(v isa AbstractTraceData || v isa AbstractHistogramData)
        throw(ArgumentError("CombinedData leg :$k must be AbstractTraceData or AbstractHistogramData, got $(typeof(v))"))
    end
    return nothing
end

"""Return modality keys of `d.legs` (canonical order)."""
combined_modalities(d::CombinedData) = Tuple(keys(d.legs))

function CombinedData(d::TraceRNAData)
    td = TraceData(d.label, d.gene, d.interval, d.trace, d.units)
    rna = RNAData(d.label, d.gene, d.nRNA, d.histRNA, d.yield, d.units)
    return CombinedData((rna=rna, trace=td))
end

"""
    reconstruct_tracerna(b::CombinedData{NT}) where {NT<:NamedTuple{(:rna, :trace)}}

Rebuild [`TraceRNAData`](@ref) from `:rna` + `:trace` legs.
"""
function reconstruct_tracerna(b::CombinedData{NT}) where {NT<:NamedTuple{(:rna, :trace)}}
    r = b.legs.rna
    t = b.legs.trace
    return TraceRNAData(t.label, t.gene, t.interval, t.trace, r.nRNA, r.histRNA, r.yield, t.units)
end

# --- Recursive hierarchy/data-collection vocabulary ---
#
# These types are intentionally passive. They describe how separately loaded datasets
# should be grouped for future hierarchical parameter sharing without changing the
# existing `hierarchical=(nhypersets, fittedindividual, fixedindividual)` behavior.

"""
    DatasetSpec

Description of one input dataset before it is loaded.

`DatasetSpec` is a small data-wrangling record: it keeps the dataset's own
`datatype`, `datapath`, `datacond`, and optional per-dataset loader settings
together with a metadata table used by recursive hierarchies. For example,
metadata can record `sample = :ThreePrime`, `replicate = :rep1`, or any other
category that a [`HierarchyNode`](@ref) uses as its `key`.

The type is deliberately permissive so existing `fit` keyword values can be
stored unchanged. It does not load data and does not alter existing fit paths.
"""
struct DatasetSpec
    name::Symbol
    datatype::Any
    datapath::Any
    datacond::Any
    metadata::NamedTuple
    trace_specs::Any
    dwell_specs::Any
end

DatasetSpec(name::Symbol, datatype, datapath, datacond;
            metadata::NamedTuple=NamedTuple(),
            trace_specs=[],
            dwell_specs=[]) =
    DatasetSpec(name, datatype, datapath, datacond, metadata, trace_specs, dwell_specs)

DatasetSpec(name::AbstractString, datatype, datapath, datacond; kwargs...) =
    DatasetSpec(Symbol(name), datatype, datapath, datacond; kwargs...)

"""
    HierarchyNode

One node in a recursive hierarchy used to assign datasets/observations to
parameter-sharing levels.

Fields:
- `name`: stable level name, e.g. `:top`, `:group`, `:individual`.
- `key`: metadata key used to assign children at this level. The root can use
  `nothing`.
- `levels`: allowed level values. An empty vector means levels can be discovered
  from data metadata later.
- `children`: nested [`HierarchyNode`](@ref)s.

This type only describes the hierarchy tree. Parameter ownership is stored on
[`HierarchySpec`](@ref), not on the node, so the same tree can be reused with
different sharing rules.
"""
struct HierarchyNode
    name::Symbol
    key::Union{Nothing,Symbol}
    levels::Vector{Any}
    children::Vector{HierarchyNode}
end

HierarchyNode(name::Symbol; key=nothing, levels=[], children=HierarchyNode[]) =
    HierarchyNode(name, key === nothing ? nothing : Symbol(key), collect(levels), collect(children))

HierarchyNode(name::AbstractString; kwargs...) = HierarchyNode(Symbol(name); kwargs...)

"""
    HierarchySpec

Recursive hierarchy plus a parameter-scope map.

`parameter_scope` maps parameter families or concrete parameter identifiers to
the hierarchy node that owns them. Examples:

```julia
HierarchySpec(
    HierarchyNode(:top, children=[
        HierarchyNode(:group, key=:sample, levels=[:ThreePrime, :FivePrime],
            children=[HierarchyNode(:individual, key=:trace_id)])
    ]);
    parameter_scope=Dict(:G => :top, :R => :group, :noise => :individual),
)
```

The current implementation stores the specification only. Later fitting code can
compile this declarative form into explicit parameter indices while preserving
the existing tuple-based hierarchical path.
"""
struct HierarchySpec
    root::HierarchyNode
    parameter_scope::Dict{Any,Symbol}
end

HierarchySpec(root::HierarchyNode; parameter_scope=Dict{Any,Symbol}()) =
    HierarchySpec(root, Dict{Any,Symbol}(parameter_scope))

"""
    HierarchyAssignment

Resolved hierarchy membership for one dataset.

Fields:
- `dataset`: dataset name.
- `path`: named tuple mapping hierarchy node names to the dataset's level at
  that node, e.g. `(top=:top, group=:ThreePrime, individual=:cell_001)`.
- `parameter_groups`: map from each `HierarchySpec.parameter_scope` key to the
  hierarchy prefix that owns that parameter. For example, if `:R => :group`,
  the group key is `(top=:top, group=:ThreePrime)`.
"""
struct HierarchyAssignment
    dataset::Symbol
    path::NamedTuple
    parameter_groups::Dict{Any,NamedTuple}
end

"""
    RecursiveHierarchyCachePlan

Compiled, data-independent cache plan for recursive hierarchical likelihoods.

The plan is intentionally general: callers choose which parameter families alter
the transition matrix (`transition_families`) and which alter the emission model
(`emission_families`). The compiler then derives the minimal group IDs from the
resolved hierarchy assignments. No 3Prime/5Prime assumptions are embedded here.

Fields:
- `assignments`: one [`HierarchyAssignment`](@ref) per dataset/trace-like item.
- `transition_families`: families whose parameters change transition matrices.
- `emission_families`: families whose parameters change emission distributions.
- `transition_group_keys`: unique group keys for transition-cache entries.
- `emission_group_keys`: unique group keys for emission-cache entries.
- `assignment_transition_group`: transition group index per assignment.
- `assignment_emission_group`: emission group index per assignment.
"""
struct RecursiveHierarchyCachePlan
    assignments::Vector{HierarchyAssignment}
    transition_families::Vector{Any}
    emission_families::Vector{Any}
    transition_group_keys::Vector{NamedTuple}
    emission_group_keys::Vector{NamedTuple}
    assignment_transition_group::Vector{Int}
    assignment_emission_group::Vector{Int}
end

"""
    RecursiveHierarchyTrait

Model trait for recursive hierarchical likelihoods.

This trait is the bridge between `fit.jl` and `likelihoods.jl`: `fit.jl`
attaches a compiled [`RecursiveHierarchyCachePlan`](@ref), and the likelihood
uses it to cache transition matrices by transition-rate group and dispatch
emissions by emission group.
"""
struct RecursiveHierarchyTrait
    cache_plan::RecursiveHierarchyCachePlan
    rate_families::Dict{Any,Vector{Int}}
    initial_values::Dict{Any,Any}
    parameter_base_indices::Vector{Int}
    assignment_parameter_indices::Matrix{Int}
    transition_group_parameter_indices::Vector{Vector{Int}}
    emission_group_parameter_indices::Vector{Vector{Int}}
end

RecursiveHierarchyTrait(cache_plan::RecursiveHierarchyCachePlan, rate_families::Dict{Any,Vector{Int}}) =
    RecursiveHierarchyTrait(cache_plan, rate_families, Dict{Any,Any}(), Int[], Matrix{Int}(undef, 0, 0), Vector{Int}[], Vector{Int}[])

RecursiveHierarchyTrait(cache_plan::RecursiveHierarchyCachePlan, rate_families::Dict{Any,Vector{Int}}, initial_values::AbstractDict) =
    RecursiveHierarchyTrait(cache_plan, rate_families, Dict{Any,Any}(pairs(initial_values)), Int[], Matrix{Int}(undef, 0, 0), Vector{Int}[], Vector{Int}[])

function _normalize_rate_families_for_trait(rate_families)
    d = Dict{Any,Vector{Int}}()
    for (k, v) in pairs(rate_families)
        d[k] = collect(Int, v)
    end
    return d
end

function _base_parameter_families(rate_families::AbstractDict, nbase::Integer)
    owners = Vector{Any}(undef, nbase)
    assigned = falses(nbase)
    for (family, inds) in pairs(rate_families)
        for i in inds
            1 <= i <= nbase || throw(ArgumentError("rate family $(repr(family)) contains index $i outside 1:$nbase"))
            assigned[i] && throw(ArgumentError("base parameter index $i belongs to more than one rate family"))
            owners[i] = family
            assigned[i] = true
        end
    end
    missing = findall(!, assigned)
    isempty(missing) || throw(ArgumentError("recursive_hierarchy rate_families must cover every base parameter index; missing $(repr(missing))"))
    return owners
end

"""
    compile_recursive_hierarchy_trait(cache_plan, rate_families, nbase; transition_indices, emission_indices)

Create the parameter-index layout for a recursive hierarchy. The resulting trait
maps each base parameter slot and assignment to a global parameter index, and
precomputes the indices needed for transition-matrix and emission caches.
"""
function compile_recursive_hierarchy_trait(
    cache_plan::RecursiveHierarchyCachePlan,
    rate_families::AbstractDict,
    nbase::Integer;
    transition_indices=1:nbase,
    emission_indices=Int[],
)
    rate_families_out = _normalize_rate_families_for_trait(rate_families)
    family_for_base = _base_parameter_families(rate_families_out, nbase)
    nassign = length(cache_plan.assignments)
    assignment_parameter_indices = Matrix{Int}(undef, nbase, nassign)
    parameter_base_indices = Int[]
    parameter_lookup = Dict{Tuple{Any,NamedTuple,Int},Int}()
    for aidx in 1:nassign
        groups = cache_plan.assignments[aidx].parameter_groups
        for base_idx in 1:nbase
            family = family_for_base[base_idx]
            haskey(groups, family) ||
                throw(ArgumentError("assignment $(cache_plan.assignments[aidx].dataset) has no parameter group for family $(repr(family))"))
            key = (family, groups[family], base_idx)
            pidx = get(parameter_lookup, key, 0)
            if pidx == 0
                push!(parameter_base_indices, base_idx)
                pidx = length(parameter_base_indices)
                parameter_lookup[key] = pidx
            end
            assignment_parameter_indices[base_idx, aidx] = pidx
        end
    end

    transition_group_parameter_indices = Vector{Int}[]
    for gid in 1:n_transition_rate_groups(cache_plan)
        aidx = findfirst(==(gid), cache_plan.assignment_transition_group)
        aidx === nothing && throw(ArgumentError("transition group $gid has no assignments"))
        push!(transition_group_parameter_indices, assignment_parameter_indices[collect(transition_indices), aidx])
    end

    emission_group_parameter_indices = Vector{Int}[]
    for gid in 1:n_emission_groups(cache_plan)
        aidx = findfirst(==(gid), cache_plan.assignment_emission_group)
        aidx === nothing && throw(ArgumentError("emission group $gid has no assignments"))
        push!(emission_group_parameter_indices, assignment_parameter_indices[collect(emission_indices), aidx])
    end

    return RecursiveHierarchyTrait(
        cache_plan,
        rate_families_out,
        Dict{Any,Any}(),
        parameter_base_indices,
        assignment_parameter_indices,
        transition_group_parameter_indices,
        emission_group_parameter_indices,
    )
end

"""
    MultiDatasetData

Loaded datasets plus their hierarchy metadata.

`MultiDatasetData` is the loaded counterpart to a vector of [`DatasetSpec`](@ref)
objects. Each entry in `datasets` is an existing StochasticGene data object
(`RNAData`, `TraceData`, `CombinedData`, etc.), and `metadata[i]` stores the
hierarchy categories for `datasets[i]`.

No likelihood behavior is attached to this type yet; it is an additive container
for the recursive hierarchy work.
"""
struct MultiDatasetData{D<:AbstractVector,M<:AbstractVector} <: AbstractExperimentalData
    datasets::D
    metadata::M
    names::Vector{Symbol}
end

function MultiDatasetData(datasets::AbstractVector, metadata::AbstractVector; names=nothing)
    length(datasets) == length(metadata) ||
        throw(ArgumentError("MultiDatasetData requires one metadata entry per dataset"))
    names_out = names === nothing ? [Symbol(:dataset_, i) for i in eachindex(datasets)] : Symbol.(names)
    length(names_out) == length(datasets) ||
        throw(ArgumentError("MultiDatasetData requires one name per dataset"))
    return MultiDatasetData(datasets, metadata, names_out)
end

"""
    hierarchy_node_names(node)

Return node names in pre-order. Useful for validation and tests.
"""
function hierarchy_node_names(node::HierarchyNode)
    names = Symbol[node.name]
    for child in node.children
        append!(names, hierarchy_node_names(child))
    end
    return names
end

"""
    hierarchy_parameter_levels(spec)

Return the unique hierarchy levels referenced by `spec.parameter_scope`.
"""
hierarchy_parameter_levels(spec::HierarchySpec) = unique(collect(values(spec.parameter_scope)))

function _metadata_value(metadata, key::Symbol)
    if metadata isa NamedTuple
        haskey(metadata, key) && return metadata[key]
    elseif metadata isa AbstractDict
        haskey(metadata, key) && return metadata[key]
    end
    hasproperty(metadata, key) && return getproperty(metadata, key)
    throw(ArgumentError("hierarchy metadata is missing key $(repr(key))"))
end

function _push_hierarchy_path!(names::Vector{Symbol}, values::Vector{Any}, node::HierarchyNode, metadata)
    node.name in names &&
        throw(ArgumentError("hierarchy node names must be unique; found duplicate $(repr(node.name))"))
    value = node.key === nothing ? node.name : _metadata_value(metadata, node.key)
    if !isempty(node.levels) && value ∉ node.levels
        throw(ArgumentError(
            "metadata value $(repr(value)) for hierarchy node $(repr(node.name)) is not in allowed levels $(repr(node.levels))",
        ))
    end
    push!(names, node.name)
    push!(values, value)
    for child in node.children
        _push_hierarchy_path!(names, values, child, metadata)
    end
    return nothing
end

"""
    hierarchy_path(spec, metadata)

Resolve one metadata row against a recursive hierarchy, returning a named tuple
from hierarchy node names to level values.
"""
function hierarchy_path(spec::HierarchySpec, metadata)
    names = Symbol[]
    values = Any[]
    _push_hierarchy_path!(names, values, spec.root, metadata)
    return NamedTuple{Tuple(names)}(Tuple(values))
end

function _hierarchy_prefix(path::NamedTuple, level::Symbol)
    path_keys = collect(keys(path))
    i = findfirst(==(level), path_keys)
    i === nothing && throw(ArgumentError("hierarchy level $(repr(level)) is not present in path $(repr(path_keys))"))
    prefix_keys = Tuple(path_keys[1:i])
    prefix_values = Tuple(values(path))[1:i]
    return NamedTuple{prefix_keys}(prefix_values)
end

"""
    hierarchy_parameter_groups(spec, path)

For one resolved hierarchy `path`, return the owning group key for each
parameter scope in `spec.parameter_scope`.
"""
function hierarchy_parameter_groups(spec::HierarchySpec, path::NamedTuple)
    groups = Dict{Any,NamedTuple}()
    for (param, level) in spec.parameter_scope
        groups[param] = _hierarchy_prefix(path, level)
    end
    return groups
end

"""
    hierarchy_assignments(spec, multidata)
    hierarchy_assignments(spec, dataset_specs)

Resolve every dataset into hierarchy paths and parameter-sharing groups.

This is still passive data wrangling: it validates metadata and returns
[`HierarchyAssignment`](@ref) records, but does not change likelihoods or the
legacy `hierarchical=(...)` tuple path.
"""
function hierarchy_assignments(spec::HierarchySpec, multidata::MultiDatasetData)
    assignments = HierarchyAssignment[]
    for i in eachindex(multidata.datasets)
        path = hierarchy_path(spec, multidata.metadata[i])
        push!(assignments, HierarchyAssignment(
            multidata.names[i],
            path,
            hierarchy_parameter_groups(spec, path),
        ))
    end
    return assignments
end

function hierarchy_assignments(spec::HierarchySpec, dataset_specs::AbstractVector{<:DatasetSpec})
    assignments = HierarchyAssignment[]
    for ds in dataset_specs
        path = hierarchy_path(spec, ds.metadata)
        push!(assignments, HierarchyAssignment(
            ds.name,
            path,
            hierarchy_parameter_groups(spec, path),
        ))
    end
    return assignments
end

function _family_group_key(groups::AbstractDict, families)
    names = Symbol[]
    values = Any[]
    seen = Set{Symbol}()
    for family in families
        haskey(groups, family) ||
            throw(ArgumentError("parameter family $(repr(family)) is not present in hierarchy assignment groups"))
        group = groups[family]
        for name in keys(group)
            value = group[name]
            if name in seen
                first_i = findfirst(==(name), names)
                values[first_i] == value ||
                    throw(ArgumentError("incompatible hierarchy values for key $(repr(name)): $(repr(values[first_i])) and $(repr(value))"))
            else
                push!(names, name)
                push!(values, value)
                push!(seen, name)
            end
        end
    end
    return NamedTuple{Tuple(names)}(Tuple(values))
end

function _group_indices(keys_in_order)
    unique_keys = NamedTuple[]
    indices = Int[]
    lookup = Dict{NamedTuple,Int}()
    for key in keys_in_order
        idx = get(lookup, key, 0)
        if idx == 0
            push!(unique_keys, key)
            idx = length(unique_keys)
            lookup[key] = idx
        end
        push!(indices, idx)
    end
    return unique_keys, indices
end

"""
    recursive_hierarchy_cache_plan(assignments; transition_families, emission_families)
    recursive_hierarchy_cache_plan(spec, dataset_specs; ...)

Compile hierarchy assignments into transition/emission cache groups.

Families in `transition_families` determine when transition matrices must differ.
Families in `emission_families` determine when emission distributions must differ.
The function returns stable integer group IDs in first-seen order.
"""
function recursive_hierarchy_cache_plan(
    assignments::AbstractVector{<:HierarchyAssignment};
    transition_families,
    emission_families,
)
    assignments_out = collect(assignments)
    transition = collect(transition_families)
    emission = collect(emission_families)
    transition_keys = [_family_group_key(a.parameter_groups, transition) for a in assignments_out]
    emission_keys = [_family_group_key(a.parameter_groups, emission) for a in assignments_out]
    transition_group_keys, assignment_transition_group = _group_indices(transition_keys)
    emission_group_keys, assignment_emission_group = _group_indices(emission_keys)
    return RecursiveHierarchyCachePlan(
        assignments_out,
        transition,
        emission,
        transition_group_keys,
        emission_group_keys,
        assignment_transition_group,
        assignment_emission_group,
    )
end

function recursive_hierarchy_cache_plan(
    spec::HierarchySpec,
    dataset_specs::AbstractVector{<:DatasetSpec};
    transition_families,
    emission_families,
)
    recursive_hierarchy_cache_plan(
        hierarchy_assignments(spec, dataset_specs);
        transition_families=transition_families,
        emission_families=emission_families,
    )
end

const SharedParameterCachePlan = RecursiveHierarchyCachePlan
const SharedParameterTrait = RecursiveHierarchyTrait

shared_parameter_cache_plan(args...; kwargs...) =
    recursive_hierarchy_cache_plan(args...; kwargs...)

compile_shared_parameter_trait(args...; kwargs...) =
    compile_recursive_hierarchy_trait(args...; kwargs...)

"""
    n_transition_rate_groups(plan)

Number of unique transition-rate groups in a recursive hierarchy cache plan.
This is the number of transition probability/steady-state objects `(a, p0)` that
must be computed and cached per likelihood evaluation.
"""
n_transition_rate_groups(plan::RecursiveHierarchyCachePlan) = length(plan.transition_group_keys)

"""
    n_emission_groups(plan)

Number of unique emission/noise groups in a recursive hierarchy cache plan.
"""
n_emission_groups(plan::RecursiveHierarchyCachePlan) = length(plan.emission_group_keys)

# Helper functions for yield Union type
"""
    get_yield_value(yield::Union{Float64, Tuple{Float64, Int}})

Extract the yield factor value from yield (either Float64 or tuple).
"""
get_yield_value(yield::Float64) = yield
get_yield_value(yield::Tuple{Float64, Int}) = yield[1]

"""
    get_nRNA_true(yield::Union{Float64, Tuple{Float64, Int}}, nRNA_observed::Int)

Get the true nRNA size. Returns nRNA_true from tuple if available, otherwise nRNA_observed.
"""
get_nRNA_true(yield::Float64, nRNA_observed::Int) = nRNA_observed
get_nRNA_true(yield::Tuple{Float64, Int}, nRNA_observed::Int) = yield[2]

# --- trace length (for AD feasibility / benchmarks) ---

function _trace_trials_longest(tr::AbstractVector)::Union{Nothing,Int}
    isempty(tr) && return 0
    L = 0
    for t in tr
        t isa AbstractArray || return nothing
        L = max(L, Int(size(t, 1)))
    end
    return L
end

function _trace_trials_longest(raw::Tuple)::Union{Nothing,Int}
    isempty(raw) && return nothing
    t1 = first(raw)
    if t1 isa AbstractVector
        return _trace_trials_longest(t1)
    elseif t1 isa AbstractMatrix
        return _trace_trials_longest([t1])
    else
        return nothing
    end
end

function _trace_trials_longest(raw::AbstractMatrix)::Union{Nothing,Int}
    _trace_trials_longest([raw])
end

function _trace_nframes_hint(raw)::Union{Nothing,Int}
    raw isa Tuple && length(raw) >= 4 || return nothing
    n = raw[4]
    (n isa Integer && n > 0) ? Int(n) : nothing
end

"""
    longest_trace_timesteps(data::AbstractTraceData) -> Union{Nothing,Int}

Largest per-trial frame count inferred from `data.trace` (and optional `nframes` hint in
coupled `tracejoint` tuples). Used for diagnostics and for warning when Zygote reverse-mode
is used on long HMM likelihoods.
"""
function longest_trace_timesteps(data::AbstractTraceData)
    hasproperty(data, :trace) || return nothing
    raw = data.trace
    L_trial = _trace_trials_longest(raw)
    L_nf = _trace_nframes_hint(raw)
    if L_trial !== nothing && L_nf !== nothing
        return max(L_trial, L_nf)
    elseif L_trial !== nothing
        return L_trial
    else
        return L_nf
    end
end

# Model structures



"""
    HMMReporter

Structure for reporters.

# Fields
- `n::Int`: Number of noise parameters.
- `per_state::Vector`: Number of reporters per state.
- `probfn::Function`: Noise distribution function, e.g., `prob_GaussianMixture`.
- `weightind::Int`: Index for mixture model bias parameter (restricted to range [0,1]).
- `offstates::Vector{Int}`: Vector of off states.
"""
struct HMMReporter
    n::Int
    per_state::Vector
    probfn::Function
    weightind::Int
    offstates::Vector{Int}
    noiseparams::Vector{Int}
end



"""
    Abstract model types
"""
abstract type AbstractModel end
abstract type AbstractGeneTransitionModel <: AbstractModel end
abstract type AbstractGMmodel <: AbstractGeneTransitionModel end
abstract type AbstractGRSMmodel{TraitType} <: AbstractGeneTransitionModel end



"""
    Model structures

fields:
- `rates`: transition rates
- `Gtransitions`:  tuple of vectors of G state transitions
- `G`: number of G steps
- `R`: number of R steps
- `S`: indicator for splicing, 0 no splicing, > 1 splicing
- `insertstep`: R step where reporter is inserted (first step where reporter is visible)
- `nalleles`: number of alleles producing RNA
- `splicetype`: choices are "", "offeject", "offdecay"
- `rateprior`: prior distribution for rates
- `proposal`: MCMC proposal distribution
- `fittedparam`: indices of rates to be fitted
- `fixedeffects`: indices of rates that are fixed to each other, in the form of a 2 tuple of vectors
    with index 1 the tied index vector and 2 the corresponding fitted index vector
- `fixedeffects`: tuple of vectors of rates that are locked together
- `method`: method option, for nonhierarchical models 1 indicates solving Master equation directly, otherwise by eigendecomposition, 
            for hierarchical models, 2-tuple, where 1st component is same as above and 2nd is Bool where true means rates are fixed for all individuals
-` reporter`: vector of reporters or sojorn states (onstates) or vectors of vectors depending on model and data

"""

"""
    GMmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GM models.

# Fields
- `rates::RateType`: Transition rates.
- `Gtransitions::Tuple`: Tuple of vectors of G state transitions.
- `G::Int`: Number of G steps.
- `nalleles::Int`: Number of alleles producing RNA.
- `rateprior::PriorType`: Prior distribution for rates.
- `proposal::ProposalType`: MCMC proposal distribution.
- `fittedparam::ParamType`: Indices of rates to be fitted.
"""
struct GMmodel{RateType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGMmodel
    rates::RateType
    Gtransitions::Tuple
    G::Int
    nalleles::Int
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end

"""
    Transformation

Structure for transformation of rates.

# Fields
- `f::Vector{Function}`: Vector of functions.
- `f_inv::Vector{Function}`: Vector of inverse functions.
"""
struct Transformation
    f::Vector{Function}
    f_inv::Vector{Function}
    f_cv::Vector{Function}
end

"""
Trait system
each trait keeps track of its indices
model keeps track of number of model rates
reporter keeps track of noise parameter indices


nallparams it total number of parameters per individual of all types including transition rates, noise parameters, coupling parameters, and grid parameters
nrates is number of transition rates per individual
nindividualparams is number of fitted parameters per individual

"""



"""
    HierarchicalTrait

Structure for hierarchical model traits.

# Fields
- `nhypersets::Int`: Number of hyperparameter sets
- `nrates::Int`: Number of transition rates per individual
- `nindividualparams::Int`: Number of fitted parameters per individual
- `nindividuals::Int`: Number of individuals (traces)
- `individualstart::Int`: Starting index for individual parameters
- `paramstart::Int`: Starting index for hyperparameters
- `hyperindices::Vector{Vector}`: Indices for hyperparameters
- `fittedshared::Vector{Int}`: Indices of shared parameters that are fitted
"""
struct HierarchicalTrait
    nhypersets::Int
    nrates::Int
    nindividualparams::Int
    nindividuals::Int
    individualstart::Int
    paramstart::Int
    hyperindices::Vector{Vector}
    fittedshared::Vector{Int}
    fittedpriors::Vector{Int}
end

"""
    GridTrait

Structure for grid-based parameter space exploration.

# Fields
- `ngrid::Int`: Number of grid points
- `gridindices::Vector`: Indices for grid parameters
"""
struct GridTrait
    ngrid::Int
    gridindices::Vector
end

"""
    CouplingTrait

Structure for coupled model traits.

# Fields
- `ncoupling::Int`: Number of coupling parameters
- `couplingindices::Vector`: Indices for coupling parameters
- `labels::Union{Vector{String}, Nothing}`: Optional canonical labels (e.g. from `coupling_parameter_labels`) for rate file headers
"""
struct CouplingTrait
    ncoupling::Int
    couplingindices::Vector
    labels::Union{Vector{String}, Nothing}
end
CouplingTrait(ncoupling::Int, couplingindices::Vector) = CouplingTrait(ncoupling, couplingindices, nothing)

"""
    hastrait(model, trait)

Check if a trait is present in a model.

# Arguments
- `model::AbstractModel`: Model to check.
- `trait::Symbol`: Trait to check.

"""
function hastrait(model, trait)
    if !(model isa AbstractGRSMmodel)
        return false
    end
    !isnothing(model.trait) && haskey(model.trait, trait)
end

"""
    GRSMmodel{RateType, PriorType, ProposalType, ParamType, MethodType, ComponentType, ReporterType}

Structure for GRSM models.

# Fields
- `trait::TraitType`: Trait type.
- `rates::RateType`: Transition rates.
- `transforms::Transformation`: Transformation of rates.
- `nrates::nratesType`: Number of transition rates.
- `Gtransitions::Tuple`: Tuple of vectors of G state transitions.
- `G::Int`: Number of G steps.
- `R::Int`: Number of R steps.
- `S::Int`: Indicator for splicing, 0 means no splicing, > 1 means splicing.
- `insertstep::Int`: R step where reporter is inserted (first step where reporter is visible).
- `nalleles::Int`: Number of alleles producing RNA.
- `splicetype::String`: Choices are "", "offeject", "offdecay".
- `rateprior::PriorType`: Prior distribution for rates.
- `proposal::ProposalType`: MCMC proposal distribution.
- `fittedparam::ParamType`: Indices of rates to be fitted.
- `fixedeffects::Tuple`: Indices of rates that are fixed to each other, in the form of a 2-tuple of vectors with index 1 being the tied index vector and 2 being the corresponding fitted index vector.
- `method::MethodType`: Method option, for non-hierarchical models 1 indicates solving Master equation directly, otherwise by eigendecomposition.
- `components::ComponentType`: Components of the model.
- `reporter::ReporterType`: Vector of reporters or sojourn states (onstates) or vectors of vectors depending on model and data.
"""

struct GRSMmodel{TraitType,RateType,nratesType,GType,PriorType,ProposalType,ParamType,MethodType,ComponentType,ReporterType} <: AbstractGRSMmodel{TraitType}
    trait::TraitType
    rates::RateType
    transforms::Transformation
    nrates::nratesType
    Gtransitions::Tuple
    G::GType
    R::GType
    S::GType
    insertstep::GType
    nalleles::Int
    splicetype::String
    rateprior::PriorType
    proposal::ProposalType
    fittedparam::ParamType
    fixedeffects::Tuple
    method::MethodType
    components::ComponentType
    reporter::ReporterType
end


"""
    print_model(model::AbstractModel)

Print all fields of model

# Arguments
- `model::AbstractModel`: Model to print.

# Returns
- `nothing`
"""
function print_model(model::AbstractModel)
    for fname in fieldnames(model)
        println("$fname =", getfield(model, fname))
    end
end

"""
Abstract Option types for fitting methods
"""
abstract type Options end
"""
    MHOptions <: Options

Options for Metropolis-Hastings MCMC.

# Fields
- `samplesteps::Int64`: Number of MCMC samples to collect.
- `warmupsteps::Int64`: Number of warmup (burn-in) steps.
- `sample_stride::Int64`: Number of MH proposal steps between stored posterior samples.
- `maxtime::Float64`: Maximum allowed runtime (seconds).
- `temp::Float64`: Temperature for MCMC.
- `merge_max_memory::Int64`: Approximate byte budget for multi-chain merge/statistics arrays.
- `init_jitter::Float64`: Initial sampling-space jitter for each MH chain. `0.0` disables jitter.
- `init_jitter_individual::Float64`: Initial jitter for hierarchical individual-level parameters. `0.0` uses `init_jitter`.
- `init_jitter_noise::Float64`: Initial jitter for noise parameters. `0.0` uses `init_jitter`.
- `proposal_cv_rates::Float64`: Optional scalar proposal CV for rate parameters. `0.0` uses `propcv`.
- `proposal_cv_noise::Float64`: Optional scalar proposal CV for noise parameters. `0.0` uses `proposal_cv_rates` / `propcv`.
- `device::Symbol`: :cpu, :gpu
- `parallel::Symbol`: :single, :threaded, :distributed
- `gradient::Symbol`: Gradient type (:none, :finite, :ForwardDiff, :Zygote)
- `likelihood_executor::Symbol`: Which **likelihood implementation** to use for modalities that support a primal
  vs reverse-diff-friendly path (same values as [`HMM_STACK_MH`](@ref) / [`HMM_STACK_AD`](@ref); see `hmm.jl`).
  Default [`HMM_STACK_MH`](@ref) (`:fast`, in-place / MCMC-oriented). Other differentiable likelihood subgraphs can
  read this field in the future; trace HMMs use it today.
- `gradient_checkpoint_length::Union{Nothing,Int}`: When an integer `> 0`, reverse-mode through **long time-series**
  likelihoods may process time in chunks of this many frames (Zygote checkpointing; see [`with_hmm_zygote_checkpoint`](@ref)).
  `nothing` leaves the process-global setting unchanged except where inference wrappers set it explicitly.
"""
struct MHOptions <: Options
    samplesteps::Int64
    warmupsteps::Int64
    sample_stride::Int64
    maxtime::Float64
    temp::Float64
    merge_max_memory::Int64
    init_jitter::Float64
    init_jitter_individual::Float64
    init_jitter_noise::Float64
    proposal_cv_rates::Float64
    proposal_cv_noise::Float64
    device::Symbol  # :cpu, :gpu
    parallel::Symbol  # :single, :threaded, :distributed
    gradient::Symbol
    likelihood_executor::Symbol
    gradient_checkpoint_length::Union{Nothing,Int}
end

function MHOptions(
    samplesteps::Integer, warmupsteps::Integer, maxtime::Real, temp::Real;
    sample_stride::Integer=1,
    merge_max_memory::Integer=48 * 1024^3,
    init_jitter::Real=0.0,
    init_jitter_individual::Real=0.0,
    init_jitter_noise::Real=0.0,
    proposal_cv_rates::Real=0.0,
    proposal_cv_noise::Real=0.0,
    device::Symbol=:cpu, parallel::Symbol=:single, gradient::Symbol=:none,
    likelihood_executor::Symbol=HMM_STACK_MH,
    gradient_checkpoint_length::Union{Nothing,Integer}=nothing,
)
    gck = gradient_checkpoint_length === nothing ? nothing : Int(gradient_checkpoint_length)
    stride = Int64(sample_stride)
    stride > 0 || throw(ArgumentError("sample_stride must be positive, got $(repr(sample_stride))"))
    merge_mem = Int64(merge_max_memory)
    merge_mem > 0 || throw(ArgumentError("merge_max_memory must be positive, got $(repr(merge_max_memory))"))
    jitter = Float64(init_jitter)
    jitter >= 0.0 || throw(ArgumentError("init_jitter must be nonnegative, got $(repr(init_jitter))"))
    jitter_individual = Float64(init_jitter_individual)
    jitter_individual >= 0.0 || throw(ArgumentError("init_jitter_individual must be nonnegative, got $(repr(init_jitter_individual))"))
    jitter_noise = Float64(init_jitter_noise)
    jitter_noise >= 0.0 || throw(ArgumentError("init_jitter_noise must be nonnegative, got $(repr(init_jitter_noise))"))
    pcv_rates = Float64(proposal_cv_rates)
    pcv_rates >= 0.0 || throw(ArgumentError("proposal_cv_rates must be nonnegative, got $(repr(proposal_cv_rates))"))
    pcv_noise = Float64(proposal_cv_noise)
    pcv_noise >= 0.0 || throw(ArgumentError("proposal_cv_noise must be nonnegative, got $(repr(proposal_cv_noise))"))
    MHOptions(
        Int64(samplesteps), Int64(warmupsteps), stride, Float64(maxtime), Float64(temp),
        merge_mem, jitter, jitter_individual, jitter_noise, pcv_rates, pcv_noise,
        device, parallel, gradient, likelihood_executor, gck,
    )
end

"""
    NUTSOptions

Options for `run_nuts` (NUTS/AdvancedHMC).

# Fields
- `n_samples`, `n_adapts`: post-warmup samples and adaptation steps
- `δ`: target acceptance (NUTS dual averaging)
- `gradient`: :finite (default), :ForwardDiff, :Zygote
- `fd_ε`: finite-difference step when `gradient === :finite`
- `verbose`, `progress`: passed to AdvancedHMC.sample (`progress` defaults to `true`, enabling AdvancedHMC’s ProgressMeter sampling bar; set `nuts_progress=false` in [`fit`](@ref) / run-spec when running many parallel chains)
- `max_depth`: optional NUTS tree-depth cap; `nothing` uses AdvancedHMC's default
- `device::Symbol`: :cpu, :gpu
- `parallel::Symbol`: :single, :threaded, :distributed
- `likelihood_executor::Symbol`: likelihood track for modalities with dual implementations (`:fast` for finite-difference NUTS, [`HMM_STACK_AD`](@ref) for AD gradients).
- `gradient_checkpoint_length::Union{Nothing,Int}`: optional frame chunk size for reverse-mode through long trace likelihoods.
"""
struct NUTSOptions <: Options
    n_samples::Int
    n_adapts::Int
    δ::Float64
    gradient::Symbol
    fd_ε::Float64
    verbose::Bool
    progress::Bool
    max_depth::Union{Nothing,Int}
    device::Symbol  # :cpu, :gpu
    parallel::Symbol  # :single, :threaded, :distributed
    likelihood_executor::Symbol
    gradient_checkpoint_length::Union{Nothing,Int}
end

function NUTSOptions(;
    n_samples=100, n_adapts=100, δ=0.8, gradient=:finite, fd_ε=1e-5, verbose=true, progress=true,
    device::Symbol=:cpu, parallel::Symbol=:single,
    likelihood_executor::Symbol=(gradient === :finite ? HMM_STACK_MH : HMM_STACK_AD),
    gradient_checkpoint_length::Union{Nothing,Integer}=nothing,
    max_depth=nothing,
)
    gck = gradient_checkpoint_length === nothing ? nothing : Int(gradient_checkpoint_length)
    md = max_depth === nothing ? nothing : Int(max_depth)
    NUTSOptions(n_samples, n_adapts, Float64(δ), gradient, Float64(fd_ε), verbose, progress, md, device, parallel, likelihood_executor, gck)
end

"""
    ADVIOptions

Options for mean-field ADVI.

# Fields
- `maxiter`: Optim iterations
- `n_mc`: Monte Carlo draws for reparameterization gradient
- `σ_floor`: lower bound on σ
- `init_s_raw`: initial value for log-scale state s
- `verbose`: Optim show trace
- `gradient`: :Zygote (default), :finite, :ForwardDiff
- `time_limit`: wall-clock limit (seconds)
- `device::Symbol`: :cpu, :gpu
- `parallel::Symbol`: :single, :threaded, :distributed
- `likelihood_executor::Symbol`: likelihood track for AD-friendly modalities (default [`HMM_STACK_AD`](@ref)).
- `gradient_checkpoint_length::Union{Nothing,Int}`: optional frame chunk size for reverse-mode through long trace likelihoods.
"""
struct ADVIOptions <: Options
    maxiter::Int
    n_mc::Int
    σ_floor::Float64
    init_s_raw::Float64
    verbose::Bool
    gradient::Symbol
    time_limit::Union{Nothing,Float64}
    zygote_trace::Bool
    device::Symbol  # :cpu, :gpu
    parallel::Symbol  # :single, :threaded, :distributed
    likelihood_executor::Symbol
    gradient_checkpoint_length::Union{Nothing,Int}
end

function ADVIOptions(;
    maxiter=500, n_mc=8, σ_floor=1e-4, init_s_raw=-4.0, verbose=false, gradient=:Zygote, time_limit=nothing,
    zygote_trace::Bool=false,
    device::Symbol=:cpu, parallel::Symbol=:single,
    likelihood_executor::Symbol=HMM_STACK_AD,
    gradient_checkpoint_length::Union{Nothing,Integer}=nothing,
)
    gck = gradient_checkpoint_length === nothing ? nothing : Int(gradient_checkpoint_length)
    ADVIOptions(
        Int(maxiter), Int(n_mc), Float64(σ_floor), Float64(init_s_raw), verbose, gradient,
        time_limit === nothing ? nothing : Float64(time_limit), zygote_trace, device, parallel,
        likelihood_executor, gck,
    )
end

"""
    INFERENCE_MH, INFERENCE_NUTS, INFERENCE_ADVI, INFERENCE_CHOICES

Symbolic tags for `fit` / run-spec `inference_method` (also accept plain symbols like `:mh`).
"""
const INFERENCE_MH = :mh
const INFERENCE_NUTS = :nuts
const INFERENCE_ADVI = :advi
const INFERENCE_CHOICES = (INFERENCE_MH, INFERENCE_NUTS, INFERENCE_ADVI)

abstract type Results end



# === Behavioral trait functions ===
#
# These functions provide soft behavioral tags for data types
# without relying on type inheritance. They help determine which
# methods are valid for a given data type without enforcing hard
# subtyping relationships.
#
# Usage example:
#     if is_deviance_supported(data)
#         println("Deviance: ", deviance(fits, data, model))
#     end

"""
    is_histogram_compatible(data) -> Bool

Return `true` if the data can be treated as a histogram (e.g. for
plotting, normalization, or PDF comparisons). 

Used to prevent histogram methods from being accidentally applied to
non-histogram types like `RNACountData`.
"""
is_histogram_compatible(::AbstractExperimentalData) = false
is_histogram_compatible(::AbstractHistogramData) = true
is_histogram_compatible(::RNACountData) = false




# ============================================================================
# Correlation Function Abstraction Layer
# ============================================================================

"""
    CorrelationTrait

Structure for correlation algorithm traits (features).

# Fields
- `centering::Symbol`: Centering method:
  - `:none`: No centering (uncentered correlation)
  - `:global_mean`: Center by global mean (constant across all lags)
  - `:windowed_mean`: Center by lag-dependent windowed mean (IDL-style)
- `multitau::Symbol`: Multi-tau binning:
  - `:none`: No multi-tau binning (uniform lag spacing)
  - `:multitau`: Use multi-tau progressive binning (IDL-style)
- `normalization::Symbol`: Normalization method:
  - `:none`: No normalization (return raw covariance/correlation)
  - `:global_mean`: Normalize by global means: (E[XY] - E[X]E[Y])/(E[X]E[Y])
  - `:windowed_mean`: Normalize by windowed means (lag-dependent)
  - `:variance`: Normalize by variance (autocorrelation normalization)
- `m::Int`: Number of points per level for multi-tau (default: 16)

# Examples
Standard correlation (uncentered, no multi-tau, no normalization):
    alg = CorrelationTrait()  # All defaults
    # or explicitly:
    alg = CorrelationTrait(centering=:none, multitau=:none, normalization=:none)

IDL correlation (windowed means, multi-tau, global mean normalization):
    alg = CorrelationTrait(centering=:windowed_mean, multitau=:multitau, normalization=:global_mean)

Windowed correlation (windowed means, no multi-tau, no normalization):
    alg = CorrelationTrait(centering=:windowed_mean)
    # or explicitly:
    alg = CorrelationTrait(centering=:windowed_mean, multitau=:none, normalization=:none)

Windowed centering + global mean normalization:
    alg = CorrelationTrait(centering=:windowed_mean, normalization=:global_mean)
"""
struct CorrelationTrait
    centering::Symbol
    multitau::Symbol
    normalization::Symbol
    m::Int  # Points per level for multi-tau
    biased::Bool  # If true, use fixed divisor (N). If false, use unbiased divisor (N-τ).
    
    function CorrelationTrait(; centering::Symbol=:none, multitau::Symbol=:none, normalization::Symbol=:none, m::Int=16, biased::Bool=false)
        centering in (:none, :global_mean, :windowed_mean, :per_trace_mean) || 
            error("centering must be :none, :global_mean, :windowed_mean, or :per_trace_mean")
        multitau in (:none, :multitau) || 
            error("multitau must be :none or :multitau")
        normalization in (:none, :global_mean, :windowed_mean, :per_trace_mean, :variance) || 
            error("normalization must be :none, :global_mean, :windowed_mean, :per_trace_mean, or :variance")
        new(centering, multitau, normalization, m, biased)
    end
end

"""
    hastrait(alg::CorrelationTrait, feature::Symbol)

Check if a correlation algorithm has a specific feature.

# Arguments
- `alg::CorrelationTrait`: Correlation algorithm to check
- `feature::Symbol`: Feature to check (`:centering`, `:multitau`, `:normalization`, or specific values like `:windowed_mean`)

# Returns
- `Bool`: Whether the algorithm has the specified feature

# Examples
```julia
alg = CorrelationTrait(:windowed_mean, :multitau, :global_mean)
hastrait(alg, :multitau)  # true
hastrait(alg, :windowed_mean)  # true (checks if centering is :windowed_mean)
hastrait(alg, :none)  # false
```
"""
function hastrait(alg::CorrelationTrait, feature::Symbol)
    if feature == :centering
        return alg.centering != :none
    elseif feature == :multitau
        return alg.multitau == :multitau
    elseif feature == :normalization
        return alg.normalization != :none
    elseif feature in (:none, :global_mean, :windowed_mean)
        return alg.centering == feature || alg.normalization == feature
    elseif feature == :variance
        return alg.normalization == :variance
    else
        return false
    end
end

# Convenience constructors for common algorithms
"""
    StandardCorrelation()

Standard correlation algorithm (uncentered, no multi-tau, no normalization).
"""
StandardCorrelation() = CorrelationTrait(centering=:none, multitau=:none, normalization=:none)

"""
    WindowedCorrelation()

Windowed means correlation algorithm (windowed centering, no multi-tau, no normalization).
"""
WindowedCorrelation() = CorrelationTrait(centering=:windowed_mean, multitau=:none, normalization=:none)

"""
    MultiTauCorrelation(; m=16)

Multi-tau correlation algorithm (uncentered, multi-tau binning, no normalization).
"""
MultiTauCorrelation(; m=16) = CorrelationTrait(centering=:none, multitau=:multitau, normalization=:none, m=m)

"""
    IDLCorrelation(; m=16)

Full IDL algorithm (windowed centering, multi-tau binning, global mean normalization).
This exactly matches the IDL Xcor algorithm implementation.
"""
IDLCorrelation(; m=16) = CorrelationTrait(centering=:windowed_mean, multitau=:multitau, normalization=:global_mean, m=m)

# Default correlation algorithm
const DEFAULT_CORRELATION_ALGORITHM = StandardCorrelation()

# Type alias for backward compatibility
const CorrelationAlgorithm = CorrelationTrait
