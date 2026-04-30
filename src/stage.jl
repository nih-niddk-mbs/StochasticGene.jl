# This file is part of StochasticGene.jl
#
# stage.jl — restructure results on disk (legacy names → key-based names for `fit(; key=...)`).
#
# Legacy basenames are parsed with [`fields`](@ref) (`io.jl`): for
# `rates_scRNA_CANX_3331_2.txt` that yields [`Result_Fields`](@ref) with
# `name=\"rates\"` (artifact type prefix), **`label=\"scRNA\"`**, `cond=\"\"`, `gene=\"CANX\"`,
# `model=\"3331\"`, `nalleles=\"2\"`. For `key_mode=:fields`, default **`(:label, :gene, :model, :nalleles)`**
# applies to [`Result_Fields`](@ref) names; four-part **`*.csv`** aggregate summaries use raw tokens (no per-gene — see docstring). Prefer `:label` / `:gene` / `:model` / … — not `:name`,
# which changes between `rates_*`, `info_*`, etc.

"""Strip one leading `results/` (or `results\\`) so `joinpath(root, \"results\", tail)` does not become `.../results/results/...`."""
function _strip_optional_results_prefix(p::AbstractString)::String
    x = String(strip(p))
    if startswith(x, "results/")
        x = String(strip(chopprefix(x, "results/")))
    elseif startswith(x, "results" * "\\")
        x = String(strip(chopprefix(x, "results\\")))
    end
    isempty(x) && throw(ArgumentError("path became empty after stripping a leading results/; pass a name like \"scRNAkey\" instead of \"results/scRNAkey\""))
    return x
end

function _stage_resolve_src(path::AbstractString, root::AbstractString)::String
    s = String(strip(path))
    isempty(s) && throw(ArgumentError("src must be non-empty"))
    isdir(s) && return abspath(s)
    p1 = joinpath(root, s)
    isdir(p1) && return abspath(p1)
    p2 = joinpath(root, "results", _strip_optional_results_prefix(s))
    isdir(p2) && return abspath(p2)
    throw(ArgumentError("results directory does not exist: $(repr(s)) under root $(repr(root))"))
end

function _stage_resolve_dst(dst::AbstractString, root_dst::AbstractString; dry_run::Bool, make::Bool)::String
    d = String(strip(dst))
    isempty(d) && throw(ArgumentError("dst must be non-empty"))
    isdir(d) && return abspath(d)
    if isabspath(d)
        make && !dry_run && mkpath(d)
        return abspath(d)
    end
    p1 = joinpath(root_dst, d)
    isdir(p1) && return abspath(p1)
    tail = _strip_optional_results_prefix(d)
    p2 = joinpath(root_dst, "results", tail)
    isdir(p2) && return abspath(p2)
    make && !dry_run && mkpath(p2)
    return abspath(p2)
end

"""
    folder_setup(root=".")

Create `data/` and `results/` under `root` when missing.
"""
function folder_setup(root::AbstractString=".")
    data = joinpath(root, "data")
    results = joinpath(root, "results")
    !ispath(data) && mkpath(data)
    !ispath(results) && mkpath(results)
    return (; data, results)
end

"""
    rna_setup(root=".")

Set up `data/` and `results/` under `root`, then seed bundled reference RNA test data
into `root/data/`.
"""
function rna_setup(root::AbstractString=".")
    folder_setup(root)
    _stage_seed_reference_data(root)
end

# Longest first so `sampled_rates` does not match before `rates`.
const _STAGE_LEGACY_RESULT_PREFIXES = (
    "ll_sampled_rates", "sampled_rates", "param-stats", "proposal-cov", "advi-measures",
    "measures", "rates", "shared", "optimized", "burst", "info",
)

function _legacy_parse_result_basename(basename_fn::AbstractString)
    for p in _STAGE_LEGACY_RESULT_PREFIXES
        pfx = string(p, "_")
        startswith(basename_fn, pfx) || continue
        sfn = String(basename_fn)
        rest = chopprefix(sfn, pfx)
        stem, ext = splitext(rest)
        isempty(stem) && continue
        return (prefix=String(p), stem=String(stem), ext=String(ext))
    end
    return nothing
end

function _default_legacy_stem_to_keybase(stem::AbstractString)
    s = String(stem)
    s = replace(s, " " => "-")
    s = replace(s, "_" => "-")
    return s
end

"""Map `key_fields` symbols onto `Result_Fields` with explicit slots (kept local so staging does not depend on external helper dispatch)."""
function _stage_result_slot(parts::Result_Fields, slot::Symbol)::Union{Nothing,String}
    v = if slot === :name
        parts.name
    elseif slot === :label
        parts.label
    elseif slot === :cond
        parts.cond
    elseif slot === :gene
        parts.gene
    elseif slot === :model
        parts.model
    elseif slot === :nalleles
        parts.nalleles
    else
        throw(ArgumentError(
            "unknown key_fields slot $(repr(slot)) for Result_Fields artifacts; use one of :name, :label, :cond, :gene, :model, :nalleles",
        ))
    end
    x = String(v)
    isempty(x) && return nothing
    return _staging_norm_seg(x)
end

"""Normalize one segment for aggregate summary CSV keys (four underscore-separated stem tokens before `.csv`)."""
function _staging_agg_csv_segment(tokens::NTuple{4,String}, slot::Symbol)::Union{Nothing,String}
    # Aggregate summaries do not encode a gene; `:gene` only reuses the third token (`cond`) for tuple compatibility with Result_Fields key lists.
    x = if slot === :name
        tokens[1]
    elseif slot === :label
        tokens[2]
    elseif slot === :cond || slot === :gene
        tokens[3]
    elseif slot === :model
        tokens[4]
    elseif slot === :nalleles
        return nothing
    else
        throw(ArgumentError(
            "unknown staging slot $(repr(slot)) for four-part aggregate *.csv; use one of :name, :label, :cond, :model, or :gene/:nalleles (see docstring)",
        ))
    end
    isempty(x) && return nothing
    return _staging_norm_seg(x)
end

function _stage_key_from_agg_csv(tokens::NTuple{4,String}, key_fields, hierarchical::Bool)::String
    isempty(key_fields) && throw(ArgumentError("key_fields must be non-empty"))
    segs = String[]
    for f in key_fields
        seg = _staging_agg_csv_segment(tokens, f)
        seg === nothing || push!(segs, seg)
    end
    isempty(segs) && throw(ArgumentError("empty key after joining aggregate CSV tokens"))
    combined_rates_key(join(segs, "-"); hierarchical=hierarchical)
end

"""Build key from [`fields`](@ref)-style slots: [`Result_Fields`](@ref) for non-CSV; four-token aggregate `*.csv` without using [`Summary_Fields`](@ref)."""
function _stage_key_from_fields(
    fn::AbstractString,
    key_fields,
    hierarchical::Bool,
)
    stem, suffix = get_suffix(String(fn))
    v = split(stem, "_")
    if suffix == "csv"
        length(v) == 4 ||
            throw(ArgumentError(
                "key_mode=:fields expects four underscore-separated stem segments before `.csv` for aggregate summaries (name_label_cond_model.csv), got $(length(v)) segments in $(repr(fn))",
            ))
        tokens = (String(v[1]), String(v[2]), String(v[3]), String(v[4]))
        return _stage_key_from_agg_csv(tokens, key_fields, hierarchical)
    end
    parts = fields(fn)
    parts isa Result_Fields ||
        throw(ArgumentError(
            "key_mode=:fields needs Result_Fields-style names for non-CSV files (e.g. *.{txt,toml,jld2}), not $(typeof(parts)) for $(repr(fn))",
        ))
    isempty(key_fields) && throw(ArgumentError("key_fields must be non-empty"))
    segs = String[]
    for f in key_fields
        seg = _stage_result_slot(parts, f)
        seg === nothing || push!(segs, seg)
    end
    isempty(segs) && throw(ArgumentError("empty key after joining fields for $(repr(fn))"))
    combined_rates_key(join(segs, "-"); hierarchical=hierarchical)
end

function _stage_legacy_file_key(
    fn::AbstractString,
    legacy_stem::AbstractString;
    key_mode::Symbol,
    key_fields,
    hierarchical::Bool,
    key_from_stem::Union{Nothing,Function},
    random_cache::Dict{String,String},
    random_key_length::Int,
)::String
    if key_from_stem !== nothing
        return String(key_from_stem(String(legacy_stem)))
    end
    if key_mode === :stem
        return combined_rates_key(_default_legacy_stem_to_keybase(legacy_stem); hierarchical=hierarchical)
    elseif key_mode === :random
        st = String(legacy_stem)
        if !haskey(random_cache, st)
            random_cache[st] = Random.randstring(random_key_length)
        end
        return combined_rates_key(random_cache[st]; hierarchical=hierarchical)
    elseif key_mode === :fields
        return _stage_key_from_fields(fn, key_fields, hierarchical)
    else
        throw(ArgumentError("key_mode must be :stem, :random, or :fields, got $(repr(key_mode))"))
    end
end

"""
    stage_label_to_key(src, dst; root=\".\", root_dst=root,
        key_mode=:stem, key_fields=(:label, :gene, :model, :nalleles), random_key_length=12,
        hierarchical=false, key_from_stem=nothing, dry_run=false, force=true, copy_unparsed=false) -> (; n, rows)

Copy files from a **legacy** flat results directory into `dst`, renaming known artifacts to the
**key-based** layout (`rates_<key>.txt`, `info_<key>.toml`, …) for `fit(; key=...)` and [`combined_rates_key`](@ref).

`src` / `dst` resolution: existing path, or `joinpath(root, path)`, or `joinpath(root, \"results\", path)` with a **single** `results` segment (a leading `results/` on `path` is stripped before that last join so `dst=\"results/scRNAkey\"` does not become `.../results/results/scRNAkey`).

### How legacy names map to [`Result_Fields`](@ref)

[`fields`](@ref) splits the full basename on `_`. For **`rates_scRNA_CANX_3331_2.txt`** (five tokens after `rates_`):

- **`name`** = `rates` (artifact prefix — **not** the experimental “run name”; it differs for `info_*`, `measures_*`, …)
- **`label`** = `scRNA`
- **`cond`** = `""`
- **`gene`** = `CANX`
- **`model`** = `3331`
- **`nalleles`** = `2`

So in your Biowulf folder, **`scRNA` is the `label`**, not `name`.

### Key selection (`key_from_stem` wins if set)

- **`key_mode=:stem`** (default): hyphenate the legacy stem (`scRNA_CANX_3331_2` → `scRNA-CANX-3331-2`), then [`combined_rates_key`](@ref) (and `-h` if `hierarchical`).
- **`key_mode=:random`**: one random alphanumeric string per legacy stem (stable across `rates`/`info`/…); length `random_key_length`. Still passed through `combined_rates_key` for the hierarchical suffix.
- **`key_mode=:fields`**: **Per-gene / per-run artifacts** (`rates_*`, `info_*`, …) use [`Result_Fields`](@ref) via [`fields`](@ref). **Default `key_fields`** **`(:label, :gene, :model, :nalleles)`** matches those stems (minus the artifact-type prefix). **Aggregate summary `*.csv`** files (four underscore-separated stem segments: name/label/cond/model) summarize many runs and **do not encode a gene**; staging parses them by raw tokens only. For aggregate CSV, `:gene` reuses the third token (same position as `cond`) only so the same `key_fields` tuple can list both result files and batch summaries — it is **not** a gene name. `:nalleles` is ignored for CSV. Example: `measures_scRNA__3331.csv` with the default tuple → key `scRNA-3331`.

# Arguments
- `key_fields`: only for `key_mode=:fields`; see above for per-run vs aggregate CSV behaviour.
- `copy_unparsed`: if `true`, copy non-matching files under their original names.
- `dry_run`, `force`: see implementation.

Returns `(n=n_files, rows)` with rows `(; src, dst, stem, key)`.
"""
function stage_label_to_key(
    src::AbstractString,
    dst::AbstractString;
    root::AbstractString=".",
    root_dst::AbstractString=root,
    key_mode::Symbol=:stem,
    key_fields::Tuple{Vararg{Symbol}} = (:label, :gene, :model, :nalleles),
    random_key_length::Int=12,
    hierarchical::Bool=false,
    key_from_stem::Union{Nothing,Function}=nothing,
    dry_run::Bool=false,
    force::Bool=true,
    copy_unparsed::Bool=false,
)
    key_mode in (:stem, :random, :fields) || throw(ArgumentError("key_mode must be :stem, :random, or :fields"))
    random_key_length >= 4 || throw(ArgumentError("random_key_length should be at least 4"))

    src_abs = _stage_resolve_src(src, String(root))
    dst_abs = _stage_resolve_dst(dst, String(root_dst); dry_run=dry_run, make=true)

    random_cache = Dict{String,String}()
    rows = NamedTuple{(:src, :dst, :stem, :key), Tuple{String,String,String,String}}[]
    seen_dest_first_src = Dict{String,String}()

    for fn in sort(readdir(src_abs))
        startswith(fn, ".") && continue
        full = joinpath(src_abs, fn)
        isfile(full) || continue
        parsed = _legacy_parse_result_basename(fn)
        if parsed === nothing
            if copy_unparsed
                dest_file = joinpath(dst_abs, fn)
                bdest = basename(dest_file)
                if haskey(seen_dest_first_src, bdest)
                    throw(ArgumentError(
                        "duplicate destination $(repr(bdest)): $(repr(seen_dest_first_src[bdest])) and $(repr(fn)); choose different key_fields or key_mode",
                    ))
                end
                seen_dest_first_src[bdest] = fn
                push!(rows, (; src=full, dst=dest_file, stem="", key=""))
            end
            continue
        end
        key = _stage_legacy_file_key(
            fn,
            parsed.stem;
            key_mode=key_mode,
            key_fields=key_fields,
            hierarchical=hierarchical,
            key_from_stem=key_from_stem,
            random_cache=random_cache,
            random_key_length=random_key_length,
        )
        dest_bn = string(parsed.prefix, "_", key, parsed.ext)
        dest_file = joinpath(dst_abs, dest_bn)
        bdest = basename(dest_file)
        if haskey(seen_dest_first_src, bdest)
            prev = seen_dest_first_src[bdest]
            throw(ArgumentError(
                "duplicate destination $(repr(bdest)): $(repr(prev)) and $(repr(fn)) both map to key $(repr(key)) with key_fields=$(repr(key_fields)). Per-gene legacy names usually need :gene (and often :model / :nalleles) in key_fields, or use key_mode=:stem.",
            ))
        end
        seen_dest_first_src[bdest] = fn
        push!(rows, (; src=full, dst=dest_file, stem=parsed.stem, key=key))
    end

    if !dry_run
        for r in rows
            mkpath(dirname(r.dst))
            cp(r.src, r.dst; force=force)
        end
    end
    return (; n=length(rows), rows=rows)
end

"""
    label_to_key(src, dst; kwargs...) -> (; n, rows)

Short alias for [`stage_label_to_key`](@ref).
All keyword arguments are forwarded unchanged.
"""
function label_to_key(src::AbstractString, dst::AbstractString; kwargs...)
    stage_label_to_key(src, dst; kwargs...)
end

function _stage_resolve_rate_file(path::AbstractString, root::AbstractString)::String
    s = String(strip(path))
    isempty(s) && throw(ArgumentError("rate file path must be non-empty"))
    isfile(s) && return abspath(s)
    p1 = joinpath(root, s)
    isfile(p1) && return abspath(p1)
    tail = _strip_optional_results_prefix(s)
    p2 = joinpath(root, "results", tail)
    isfile(p2) && return abspath(p2)
    throw(ArgumentError("rates file does not exist: $(repr(s)) under root $(repr(root))"))
end

function _stage_resolve_output_file(path::AbstractString, root_out::AbstractString; write_out::Bool)::String
    d = String(strip(path))
    isempty(d) && throw(ArgumentError("output file path must be non-empty"))
    out_abs = isabspath(d) ? abspath(d) : abspath(joinpath(root_out, d))
    write_out && mkpath(dirname(out_abs))
    return out_abs
end

"""
Merge rate matrices (one per unit) set-wise.

**Stride:** each unit’s table is sliced in blocks of width `set_widths[i]` (= `number_of_parameters[i]` from
the caller), i.e. the same “columns per set” role as `Nenh` / `Ngene` in [`merge_coupled_two_unit_rates`](@ref)
(`io.jl`). That walk `(set-1)*sw+1:set*sw` is the stride logic; there is no separate hidden stride parameter.

**`hierarchical` here (insert layout only):** this flag is **not** the same as `fit(; hierarchical=...)` or the
`hierarchical` argument in external swarm builders (e.g. `makescriptcoupled.jl` → `makeswarm_coupled`): those
choose priors / hierarchy tuples / labels for **fitting**. Here it only chooses how `new_params` are laid
out after strided unit columns:

- `false`: after **each** set’s unit blocks, append `new_params` (labels suffixed with `_<set>` when
  `nset > 1`). This matches the **per-set** coupling tail pattern of [`merge_coupled_two_unit_rates`](@ref).
- `true`: append `new_params` **once** after **all** sets (shared tail); labels unchanged.

For combined starts that should mirror `create_combined_file` / `merge_coupled_two_unit_rates`, use
`hierarchical=false` and set `number_of_parameters` to the same per-unit widths you would pass as `Nenh`,
`Ngene` (often **larger** when each set already bundles hierarchical hyper + individual columns — see the
docstring of [`create_combined_file`](@ref)).
"""
function _stage_hcat_rate_units(
    data::Vector{Matrix{Float64}},
    hdrs::Vector{Vector{String}},
    set_widths::Vector{Int},
    new_params::Vector{Float64},
    new_labels::Vector{String},
    hierarchical::Bool,
    key::AbstractString,
)::Tuple{Matrix{Float64},Vector{String},Int}
    nrows = size(data[1], 1)
    nsets = Int[]
    for i in eachindex(data)
        size(data[i], 1) == nrows || throw(ArgumentError("row count mismatch for $(repr(key)) across units"))
        w = size(data[i], 2)
        sw = set_widths[i]
        w % sw == 0 || throw(ArgumentError(
            "$(repr(key)) unit $(i): width $w not divisible by set width $sw (number_of_parameters for this unit)",
        ))
        length(hdrs[i]) == w || throw(ArgumentError("$(repr(key)) unit $(i): header length mismatch"))
        push!(nsets, div(w, sw))
    end
    all(==(nsets[1]), nsets) || throw(ArgumentError("$(repr(key)): inconsistent set counts across units: $(nsets)"))
    nset = nsets[1]

    out = zeros(Float64, nrows, 0)
    outh = String[]
    if hierarchical
        for s in 1:nset
            for i in eachindex(data)
                sw = set_widths[i]
                cols = (s - 1) * sw + 1:s * sw
                out = hcat(out, data[i][:, cols])
                append!(outh, hdrs[i][cols])
            end
        end
        if !isempty(new_params)
            out = hcat(out, repeat(reshape(new_params, 1, :), nrows, 1))
            append!(outh, new_labels)
        end
    else
        for s in 1:nset
            for i in eachindex(data)
                sw = set_widths[i]
                cols = (s - 1) * sw + 1:s * sw
                out = hcat(out, data[i][:, cols])
                append!(outh, hdrs[i][cols])
            end
            if !isempty(new_params)
                out = hcat(out, repeat(reshape(new_params, 1, :), nrows, 1))
                for j in eachindex(new_params)
                    base = new_labels[j]
                    push!(outh, nset == 1 ? base : string(base, "_", s))
                end
            end
        end
    end
    return (out, outh, nset)
end

"""Columns per set for each unit; `number_of_parameters` is always a vector of length `n_units`."""
function _stage_set_widths_per_unit(
    number_of_parameters::AbstractVector{<:Integer},
    n_units::Int,
)::Vector{Int}
    n_units > 0 || throw(ArgumentError("n_units must be positive"))
    v = collect(Int, number_of_parameters)
    length(v) == n_units || throw(ArgumentError(
        "number_of_parameters must have length $(n_units) (one per rate file), got $(length(v))",
    ))
    all(>(0), v) || throw(ArgumentError("each number_of_parameters entry must be > 0"))
    return v
end

"""
    stage_combine_rates(rate_files, output_file, number_of_parameters, new_params=Float64[],
                        new_param_labels=String[], hierarchical=false, write_out=true, force=true, root=".", root_out=root) -> (; srcs, dst, nsets)

Merge the tables at `rate_files` (one path per unit) into a single `output_file`. **All arguments are
positional** (optional trailing arguments use defaults as shown).

`number_of_parameters` is a vector of length `length(rate_files)`: entry `i` is the number of rate columns
**per set** in `rate_files[i]` (e.g. `[13, 13]`). This is the **stride** through that file’s columns (same idea
as `Nenh` / `Ngene` in [`merge_coupled_two_unit_rates`](@ref)). Compute from your model / transitions / `G`,
`R`, `S`, noise counts so widths divide evenly.

`new_params` / `new_param_labels`: extra columns (same value on every row). The boolean `hierarchical` only
controls **per-set vs single tail** placement (see the internal `_stage_hcat_rate_units` docstring); it does
not change stride. External scripts such as `makescriptcoupled.jl` use a different notion of “hierarchical”
for **swarms and fits**; for rate **merging**, align `number_of_parameters` with your on-disk template.

Examples:

```julia
# Match coupled two-unit per-set layout (like merge_coupled_two_unit_rates + coupling columns):
stage_combine_rates(rate_files, out, [13, 13], [0.1, 0.2], ["coupling_a", "coupling_b"], false)

# Shared tail once after all strided sets:
stage_combine_rates(rate_files, out, [13, 13], [0.1, 0.2], ["coupling_a", "coupling_b"], true)
```

If `write_out` is false, nothing is written. If `force` is false, an existing output file causes an error.
Relative paths resolve under `root` / `root_out`.
"""
function stage_combine_rates(
    rate_files::AbstractVector{<:AbstractString},
    output_file::AbstractString,
    number_of_parameters::AbstractVector{<:Integer},
    new_params::AbstractVector{<:Real}=Float64[],
    new_param_labels::AbstractVector{<:AbstractString}=String[],
    hierarchical::Bool=false,
    write_out::Bool=true,
    force::Bool=true,
    root::AbstractString=".",
    root_out::AbstractString=root,
)
    n_units = length(rate_files)
    n_units > 0 || throw(ArgumentError("rate_files must be non-empty"))

    set_widths = _stage_set_widths_per_unit(number_of_parameters, n_units)

    np = isempty(new_params) ? Float64[] : collect(Float64, new_params)
    nl = isempty(new_param_labels) ? String[] : [String(strip(x)) for x in new_param_labels]
    if isempty(np)
        isempty(nl) || throw(ArgumentError("new_param_labels must be empty when new_params is empty"))
    else
        length(nl) == length(np) || throw(ArgumentError(
            "new_param_labels must have the same length as new_params ($(length(np))), got $(length(nl))",
        ))
        any(isempty, nl) && throw(ArgumentError("new_param_labels must be non-empty strings"))
    end

    srcs = [_stage_resolve_rate_file(f, String(root)) for f in rate_files]
    dst = _stage_resolve_output_file(output_file, String(root_out); write_out=write_out)
    key_label = basename(dst)

    data = Matrix{Float64}[]
    hdrs = Vector{String}[]
    for fp in srcs
        d, h = read_rates_table(fp)
        push!(data, d)
        push!(hdrs, h)
    end
    out, outh, nset = _stage_hcat_rate_units(data, hdrs, set_widths, np, nl, hierarchical, key_label)
    if write_out
        !force && isfile(dst) && throw(ArgumentError("destination exists: $(repr(dst))"))
        write_rates_table(dst, out, outh)
    end
    return (; srcs=srcs, dst=dst, nsets=nset)
end

"""Normalize coupling mode text (`:activate`, `>0`, etc.) to `:activate`/`:inhibit`/`:free`."""
function _stage_parse_coupling_mode(x)::Symbol
    x === missing && return :free
    if x isa Symbol
        s = lowercase(String(x))
        s == "activate" && return :activate
        s == "inhibit" && return :inhibit
        s == "free" && return :free
    elseif x isa Number
        v = Float64(x)
        v > 0 && return :activate
        v < 0 && return :inhibit
        return :free
    end
    s = lowercase(strip(String(x)))
    isempty(s) && return :free
    s in ("activate", ":activate", "+", "+1", "positive", "pos", ">0") && return :activate
    s in ("inhibit", ":inhibit", "-", "-1", "negative", "neg", "<0") && return :inhibit
    s in ("free", ":free", "0", "none", "neutral") && return :free
    return parse_coupling_sign_csv(s)
end

function _stage_apply_key_flag(key::AbstractString, key_flag::AbstractString)::String
    f = strip(String(key_flag))
    isempty(f) && return String(key)
    # Keep filename/key-safe behavior aligned with existing key normalization style.
    f = replace(f, " " => "-")
    return string(f, "-", key)
end

"""
    stage_combine_rates_specs_from_csv(
        csv_path, coupling_cols, coupling_labels;
        key_col="Model_name", skip_empty=true, key_normalizer=identity,
        base_new_params=Float64[], base_new_labels=String[],
        key_flag="",
        coupling_mode_values=(free=default_coupling_gamma_csv(:free), activate=default_coupling_gamma_csv(:activate), inhibit=default_coupling_gamma_csv(:inhibit)),
    ) -> Vector{NamedTuple}

Read a coupling CSV and build per-key append specs for [`stage_combine_rates`](@ref).

- `coupling_cols`: CSV columns containing coupling sign/mode text per coupling parameter.
- `coupling_labels`: output labels for those coupling parameters (same length as `coupling_cols`).
- `coupling_mode_values`: Dict/NamedTuple providing numeric values for `:free`, `:activate`, `:inhibit`.
- `base_new_params`/`base_new_labels`: optional non-coupling appended parameters (e.g. hidden units), copied
  into every output spec before coupling parameters.
- `key_flag`: if non-empty, the effective key is `"{flag}-{key}"` (spaces in the flag become `-`).
"""
function stage_combine_rates_specs_from_csv(
    csv_path::AbstractString,
    coupling_cols::AbstractVector{<:AbstractString},
    coupling_labels::AbstractVector{<:AbstractString};
    key_col::AbstractString="Model_name",
    skip_empty::Bool=true,
    key_normalizer=identity,
    base_new_params::AbstractVector{<:Real}=Float64[],
    base_new_labels::AbstractVector{<:AbstractString}=String[],
    key_flag::AbstractString="",
    coupling_mode_values=(free=default_coupling_gamma_csv(:free), activate=default_coupling_gamma_csv(:activate), inhibit=default_coupling_gamma_csv(:inhibit)),
)
    length(coupling_cols) == length(coupling_labels) ||
        throw(ArgumentError("coupling_cols and coupling_labels must have the same length"))
    length(base_new_params) == length(base_new_labels) ||
        throw(ArgumentError("base_new_params and base_new_labels must have the same length"))

    mode_vals = if coupling_mode_values isa NamedTuple
        Dict{Symbol,Float64}(k => Float64(getfield(coupling_mode_values, k)) for k in keys(coupling_mode_values))
    else
        Dict{Symbol,Float64}(Symbol(k) => Float64(v) for (k, v) in pairs(coupling_mode_values))
    end
    for k in (:free, :activate, :inhibit)
        haskey(mode_vals, k) || throw(ArgumentError("coupling_mode_values must define $(repr(k))"))
    end

    df = DataFrame(CSV.File(csv_path))
    name_syms = Symbol.(names(df))
    kc = Symbol(key_col)
    kc in name_syms || throw(ArgumentError("missing key column $(repr(key_col)) in $(repr(csv_path))"))
    csyms = Symbol.(coupling_cols)
    for c in csyms
        c in name_syms || throw(ArgumentError("missing coupling column $(repr(String(c))) in $(repr(csv_path))"))
    end

    rows = NamedTuple{(:key, :new_params, :new_labels), Tuple{String,Vector{Float64},Vector{String}}}[]
    base_vals = collect(Float64, base_new_params)
    base_lbls = [String(strip(x)) for x in base_new_labels]
    any(isempty, base_lbls) && throw(ArgumentError("base_new_labels must be non-empty strings"))
    coup_lbls = [String(strip(x)) for x in coupling_labels]
    any(isempty, coup_lbls) && throw(ArgumentError("coupling_labels must be non-empty strings"))

    for row in eachrow(df)
        key_raw = strip(string(row[kc]))
        (skip_empty && isempty(key_raw)) && continue
        key = String(key_normalizer(key_raw))
        key = _stage_apply_key_flag(key, key_flag)
        isempty(strip(key)) && throw(ArgumentError("key became empty after key_normalizer for input $(repr(key_raw))"))

        cvals = Float64[]
        for c in csyms
            mode = _stage_parse_coupling_mode(row[c])
            push!(cvals, mode_vals[mode])
        end
        push!(rows, (;
            key=key,
            new_params=vcat(base_vals, cvals),
            new_labels=vcat(base_lbls, coup_lbls),
        ))
    end
    return rows
end

function _stage_is_coupling_cell(x)::Bool
    x === missing && return true
    if x isa Number
        return true
    end
    s = lowercase(strip(String(x)))
    isempty(s) && return true
    s in ("activate", ":activate", "+", "+1", "positive", "pos", ">0") && return true
    s in ("inhibit", ":inhibit", "-", "-1", "negative", "neg", "<0") && return true
    s in ("free", ":free", "0", "none", "neutral") && return true
    try
        parse(Float64, s)
        return true
    catch
        return false
    end
end

"""Infer coupling sign/mode columns from CSV by excluding `key_col` and keeping columns with sign-like values."""
function _stage_detect_coupling_cols(csv_path::AbstractString, key_col::AbstractString)::Vector{String}
    df = DataFrame(CSV.File(csv_path))
    name_strs = String.(names(df))
    key = String(key_col)
    out = String[]
    for nm in name_strs
        nm == key && continue
        col = df[!, Symbol(nm)]
        has_nonempty = false
        ok = true
        for v in col
            v === missing && continue
            s = strip(string(v))
            isempty(s) || (has_nonempty = true)
            if !_stage_is_coupling_cell(v)
                ok = false
                break
            end
        end
        if ok && has_nonempty
            push!(out, nm)
        end
    end
    isempty(out) && throw(ArgumentError(
        "could not infer coupling columns from $(repr(csv_path)); pass explicit coupling columns or use sign-like values (>0/<0/free)",
    ))
    return out
end

"""
    stage_combine_rates_from_csv(
        csv_path, rate_files, number_of_parameters, output_folder;
        key_col="Model_name", skip_empty=true, key_normalizer=identity,
        base_new_params=Float64[], base_new_labels=String[],
        key_flag="",
        coupling_mode_values=(free=default_coupling_gamma_csv(:free), activate=default_coupling_gamma_csv(:activate), inhibit=default_coupling_gamma_csv(:inhibit)),
        coupling_cols=nothing, coupling_labels=nothing,
        hierarchical=false, write_out=true, force=true, root=".", root_out=root,
    ) -> (; n, rows)

Batch driver over [`stage_combine_rates`](@ref): for each CSV row, generate `new_params/new_labels` from
coupling sign/mode columns and write `rates_<key>.txt` into `output_folder`.
"""
function stage_combine_rates_from_csv(
    csv_path::AbstractString,
    rate_files::AbstractVector{<:AbstractString},
    number_of_parameters::AbstractVector{<:Integer},
    output_folder::AbstractString,
    ;
    key_col::AbstractString="Model_name",
    skip_empty::Bool=true,
    key_normalizer=identity,
    base_new_params::AbstractVector{<:Real}=Float64[],
    base_new_labels::AbstractVector{<:AbstractString}=String[],
    key_flag::AbstractString="",
    coupling_mode_values=(free=default_coupling_gamma_csv(:free), activate=default_coupling_gamma_csv(:activate), inhibit=default_coupling_gamma_csv(:inhibit)),
    coupling_cols=nothing,
    coupling_labels=nothing,
    hierarchical::Bool=false,
    write_out::Bool=true,
    force::Bool=true,
    root::AbstractString=".",
    root_out::AbstractString=root,
)
    cols = coupling_cols === nothing ? _stage_detect_coupling_cols(csv_path, key_col) : String.(coupling_cols)
    labels = if coupling_labels === nothing
        ["Coupling_$(i)" for i in eachindex(cols)]
    else
        String.(coupling_labels)
    end
    specs = stage_combine_rates_specs_from_csv(
        csv_path,
        cols,
        labels;
        key_col=key_col,
        skip_empty=skip_empty,
        key_normalizer=key_normalizer,
        base_new_params=base_new_params,
        base_new_labels=base_new_labels,
        key_flag=key_flag,
        coupling_mode_values=coupling_mode_values,
    )
    out_abs = _stage_resolve_output_file(output_folder, String(root_out); write_out=write_out)
    rows = NamedTuple{(:key, :dst, :result), Tuple{String,String,NamedTuple}}[]
    for s in specs
        outfile = joinpath(out_abs, "rates_" * s.key * ".txt")
        r = stage_combine_rates(
            rate_files,
            outfile,
            number_of_parameters,
            s.new_params,
            s.new_labels,
            hierarchical,
            write_out,
            force,
            root,
            ".",
        )
        push!(rows, (; key=s.key, dst=r.dst, result=r))
    end
    return (; n=length(rows), rows=rows)
end

function _stage_read_keys_from_csv(csv_path::AbstractString; key_col::AbstractString="Model_name", skip_empty::Bool=true)
    df = DataFrame(CSV.File(csv_path))
    kc = Symbol(key_col)
    kc in Symbol.(names(df)) || throw(ArgumentError("missing key column $(repr(key_col)) in $(repr(csv_path))"))
    keys = String[]
    for row in eachrow(df)
        key = strip(string(row[kc]))
        (skip_empty && isempty(key)) && continue
        isempty(key) || push!(keys, key)
    end
    isempty(keys) && throw(ArgumentError("no keys found in $(repr(csv_path)) using key_col=$(repr(key_col))"))
    return unique(keys)
end

function _stage_seed_reference_data(root::AbstractString=".")
    alleles = joinpath(root, "data", "alleles")
    halflives = joinpath(root, "data", "halflives")
    testdata = joinpath(root, "data", "HCT116_testdata")
    rnatest = joinpath(root, "data", "rnatest")
    for d in (alleles, halflives, testdata)
        !ispath(d) && mkpath(d)
    end
    !ispath(rnatest) && mkpath(rnatest)
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/alleles/CH12_alleles.csv", joinpath(alleles, "CH12_alleles.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/alleles/HCT116_alleles.csv", joinpath(alleles, "HCT116_alleles.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/ESC_halflife.csv", joinpath(halflives, "ESC_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/CH12_halflife.csv", joinpath(halflives, "CH12_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/HCT116_halflife.csv", joinpath(halflives, "HCT116_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/OcaB_halflife.csv", joinpath(halflives, "OcaB_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "aB_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "CAST_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "FIBS_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "MAST_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "NK_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "TEC_halflife.csv"))
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", joinpath(halflives, "SKIN_halflife.csv"))
    repo_root = normpath(joinpath(@__DIR__, ".."))
    candidates = (
        joinpath(repo_root, "data", "HCT116_testdata"),
        joinpath(repo_root, "test", "data", "HCT116_testdata"),
    )
    src_testdata = ""
    for c in candidates
        if isdir(c)
            src_testdata = c
            break
        end
    end
    isempty(src_testdata) && throw(ArgumentError("could not locate bundled HCT116_testdata in repository"))
    for f in readdir(src_testdata)
        cp(joinpath(src_testdata, f), joinpath(testdata, f); force=true)
    end
    # Anonymized bundled test set used by default `fit()`.
    src_gene1 = if isfile(joinpath(src_testdata, "MYC_MOCK.txt"))
        joinpath(src_testdata, "MYC_MOCK.txt")
    else
        joinpath(src_testdata, readdir(src_testdata)[1])
    end
    src_gene2 = if isfile(joinpath(src_testdata, "CENPL_MOCK.txt"))
        joinpath(src_testdata, "CENPL_MOCK.txt")
    else
        src_gene1
    end
    cp(src_gene1, joinpath(rnatest, "gene1.txt"); force=true)
    cp(src_gene2, joinpath(rnatest, "gene2.txt"); force=true)
    return (; alleles, halflives, testdata, rnatest)
end

"""
    create_label(label, datatype, datacond, cell; transition_type="")

Return `label` when non-empty; otherwise build default
`datatype * "-" * cell * ("-" * transition_type when set) * "_" * datacond`.
"""
function create_label(label::AbstractString, datatype::AbstractString, datacond::AbstractString, cell::AbstractString; transition_type::AbstractString="")
    if isempty(strip(label))
        out = string(datatype) * "-" * string(cell)
        if !isempty(strip(transition_type))
            out = out * "-" * string(transition_type)
        end
        return out * "_" * string(datacond)
    end
    return string(label)
end

create_label(label::AbstractString, datatype::AbstractString, datacond::Symbol, cell::AbstractString; transition_type::AbstractString="") =
    create_label(label, datatype, string(datacond), cell; transition_type=transition_type)

function create_label(label::AbstractString, datatype::AbstractString, datacond::AbstractVector{<:AbstractString}, cell::AbstractString; transition_type::AbstractString="")
    create_label(label, datatype, join(string.(datacond), "-"), cell; transition_type=transition_type)
end

function create_label(label::AbstractString, datatype::AbstractString, datacond::AbstractVector, cell::AbstractString; transition_type::AbstractString="")
    create_label(label, datatype, join((string(x) for x in datacond), "-"), cell; transition_type=transition_type)
end

"""
    folder_path(folder::String, root::String, folderatetype::String=""; make=false)

Resolve folder path relative to `root` and optional folder type (e.g. `"results"`), optionally creating it.
"""
function folder_path(folder::String, root::String, folderatetype::String=""; make=false)
    f = folder
    if ~ispath(folder) && ~isempty(folder)
        f = joinpath(root, folder)
        if ~ispath(f)
            f = joinpath(root, folderatetype, folder)
            if ~ispath(f) && !make
                println("$folder not found")
            elseif ~ispath(f) && make
                mkpath(f)
            end
        end
    end
    f
end

function folder_path(folder::Vector, root, foldertype)
    fv = folder
    for i in eachindex(fv)
        fv[i] = folder_path(fv[i], root, foldertype)
    end
    fv
end

"""
    make_fitscript(key; juliafile="fitscript", filedir=".", src="", kwargs...) -> String

Write one key-based fit script (`juliafile_<key>.jl`) via [`write_fitfile_key`](@ref) and
return the script path.
"""
function make_fitscript(
    key::AbstractString;
    juliafile::String="fitscript",
    filedir::AbstractString=".",
    src::AbstractString="",
    kwargs...,
)
    !isempty(filedir) && !isdir(filedir) && mkpath(filedir)
    scriptpath = joinpath(String(filedir), juliafile * "_" * sanitize_for_filename(String(key)) * ".jl")
    write_fitfile_key(scriptpath, String(key); src=src, kwargs...)
    return scriptpath
end

"""
    make_fitscripts_from_csv(csv_path; key_col="Model_name", skip_empty=true, juliafile="fitscript", filedir=".", src="", kwargs...) -> Vector{String}

Stage-native script emitter: read keys from a CSV and write one key-based fit script per key
(`juliafile_<key>.jl`) via [`make_fitscript`](@ref). This writes scripts only (no command file).
"""
function make_fitscripts_from_csv(
    csv_path::AbstractString;
    key_col::AbstractString="Model_name",
    skip_empty::Bool=true,
    juliafile::String="fitscript",
    filedir::AbstractString=".",
    src::AbstractString="",
    kwargs...,
)
    keys = _stage_read_keys_from_csv(csv_path; key_col=key_col, skip_empty=skip_empty)
    return [
        make_fitscript(
            k;
            juliafile=juliafile,
            filedir=filedir,
            src=src,
            kwargs...,
        ) for k in keys
    ]
end

"""
    build_julia_script_command(script; julia_bin="julia", nthreads=1, nprocs=2, project="", sysimage="", extra_flags=String[], script_args=String[]) -> String

Build one command line to execute a fit script. This is scheduler-agnostic and can be used for
swarm files, plain command lists, or other launchers.
"""
function build_julia_script_command(
    script::AbstractString;
    julia_bin::AbstractString="julia",
    nthreads::Integer=1,
    nprocs::Integer=2,
    project::AbstractString="",
    sysimage::AbstractString="",
    extra_flags::AbstractVector{<:AbstractString}=String[],
    script_args::AbstractVector{<:AbstractString}=String[],
)::String
    parts = String[String(julia_bin)]
    !isempty(project) && push!(parts, "--project=$(String(project))")
    !isempty(sysimage) && push!(parts, "--sysimage=$(String(sysimage))")
    nthreads > 0 && push!(parts, "-t $(Int(nthreads))")
    nprocs > 0 && push!(parts, "-p $(Int(nprocs))")
    append!(parts, String.(extra_flags))
    push!(parts, String(script))
    append!(parts, String.(script_args))
    return join(parts, " ")
end

"""
    write_julia_command_file(command_path, scripts; julia_bin="julia", nthreads=1, nprocs=2, project="", sysimage="", extra_flags=String[]) -> String

Write one command per script into `command_path`.
"""
function write_julia_command_file(
    command_path::AbstractString,
    scripts::AbstractVector{<:AbstractString};
    julia_bin::AbstractString="julia",
    nthreads::Integer=1,
    nprocs::Integer=2,
    project::AbstractString="",
    sysimage::AbstractString="",
    extra_flags::AbstractVector{<:AbstractString}=String[],
)::String
    f = open(command_path, "w")
    for s in scripts
        cmd = build_julia_script_command(
            s;
            julia_bin=julia_bin,
            nthreads=nthreads,
            nprocs=nprocs,
            project=project,
            sysimage=sysimage,
            extra_flags=extra_flags,
        )
        write(f, cmd * "\n")
    end
    close(f)
    return String(command_path)
end

"""
    make_commandfile(scripts; commandfile="fit.commands", filedir=".", julia_bin="julia", nthreads=1, nprocs=2, project="", sysimage="", extra_flags=String[]) -> String

Write one Julia launch command per script path in `scripts`.
"""
function make_commandfile(
    scripts::AbstractVector{<:AbstractString};
    commandfile::String="fit.commands",
    filedir::AbstractString=".",
    julia_bin::AbstractString="julia",
    nthreads::Integer=1,
    nprocs::Integer=2,
    project::AbstractString="",
    sysimage::AbstractString="",
    extra_flags::AbstractVector{<:AbstractString}=String[],
)::String
    !isempty(filedir) && !isdir(filedir) && mkpath(filedir)
    command_path = joinpath(String(filedir), commandfile)
    return write_julia_command_file(
        command_path,
        scripts;
        julia_bin=julia_bin,
        nthreads=nthreads,
        nprocs=nprocs,
        project=project,
        sysimage=sysimage,
        extra_flags=extra_flags,
    )
end

"""
    make_commandfile_from_csv(csv_path; key_col="Model_name", skip_empty=true, commandfile="fit.commands", juliafile="fitscript", filedir=".", julia_bin="julia", nthreads=1, nprocs=2, project="", sysimage="", extra_flags=String[]) -> String

Read keys from CSV and write a scheduler-agnostic Julia command file with one command per
`juliafile_<key>.jl` script.
"""
function make_commandfile_from_csv(
    csv_path::AbstractString;
    key_col::AbstractString="Model_name",
    skip_empty::Bool=true,
    commandfile::String="fit.commands",
    juliafile::String="fitscript",
    filedir::AbstractString=".",
    julia_bin::AbstractString="julia",
    nthreads::Integer=1,
    nprocs::Integer=2,
    project::AbstractString="",
    sysimage::AbstractString="",
    extra_flags::AbstractVector{<:AbstractString}=String[],
)
    keys = _stage_read_keys_from_csv(csv_path; key_col=key_col, skip_empty=skip_empty)
    scripts = [juliafile * "_" * sanitize_for_filename(k) * ".jl" for k in keys]
    return make_commandfile(
        scripts;
        commandfile=commandfile,
        filedir=filedir,
        julia_bin=julia_bin,
        nthreads=nthreads,
        nprocs=nprocs,
        project=project,
        sysimage=sysimage,
        extra_flags=extra_flags,
    )
end

"""
    make_swarmfile_from_csv(csv_path; key_col="Model_name", skip_empty=true, swarmfile="fit", juliafile="fitscript", filedir=".", nchains=2, nthreads=1, project="", sysimage="") -> String

Biowulf-compatible wrapper: read keys from a CSV and write `<swarmfile>.swarm` with one line
per key script (`juliafile_<key>.jl`). Internally this is a naming/keyword compatibility wrapper
over [`make_commandfile_from_csv`](@ref).
"""
function make_swarmfile_from_csv(
    csv_path::AbstractString;
    key_col::AbstractString="Model_name",
    skip_empty::Bool=true,
    swarmfile::String="fit",
    juliafile::String="fitscript",
    filedir::AbstractString=".",
    nchains::Int=2,
    nthreads=1,
    project::AbstractString="",
    sysimage::AbstractString="",
)
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    return make_commandfile_from_csv(
        csv_path;
        key_col=key_col,
        skip_empty=skip_empty,
        commandfile=sfile,
        juliafile=juliafile,
        filedir=filedir,
        julia_bin="julia",
        nthreads=nthreads,
        nprocs=nchains,
        project=project,
        sysimage=sysimage,
    )
end

"""
    make_fitscripts_and_commandfile_from_csv(csv_path; kwargs...) -> NamedTuple

Convenience wrapper that first writes scripts via [`make_fitscripts_from_csv`](@ref), then writes the
command file via [`make_commandfile_from_csv`](@ref). Returns `(; commandfile, scripts)`.
"""
function make_fitscripts_and_commandfile_from_csv(
    csv_path::AbstractString;
    key_col::AbstractString="Model_name",
    skip_empty::Bool=true,
    commandfile::String="fit.commands",
    juliafile::String="fitscript",
    filedir::AbstractString=".",
    nprocs::Int=2,
    nthreads=1,
    project::AbstractString="",
    sysimage::AbstractString="",
    julia_bin::AbstractString="julia",
    extra_flags::AbstractVector{<:AbstractString}=String[],
    src::AbstractString="",
    kwargs...,
)
    scripts = make_fitscripts_from_csv(
        csv_path;
        key_col=key_col,
        skip_empty=skip_empty,
        juliafile=juliafile,
        filedir=filedir,
        src=src,
        kwargs...,
    )
    command_path = make_commandfile_from_csv(
        csv_path;
        key_col=key_col,
        skip_empty=skip_empty,
        commandfile=commandfile,
        juliafile=juliafile,
        filedir=filedir,
        julia_bin=julia_bin,
        nthreads=nthreads,
        nprocs=nprocs,
        project=project,
        sysimage=sysimage,
        extra_flags=extra_flags,
    )
    return (; commandfile=command_path, scripts)
end

"""
    make_fitscripts_and_swarm_from_csv(csv_path; kwargs...) -> NamedTuple

Compatibility wrapper over [`make_fitscripts_and_commandfile_from_csv`](@ref) that keeps historical
`swarm` naming and returns `(; swarm, scripts)`.
"""
function make_fitscripts_and_swarm_from_csv(
    csv_path::AbstractString;
    swarmfile::String="fit",
    nchains::Int=2,
    kwargs...,
)
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    out = make_fitscripts_and_commandfile_from_csv(
        csv_path;
        commandfile=sfile,
        nprocs=nchains,
        kwargs...,
    )
    return (; swarm=out.commandfile, scripts=out.scripts)
end

