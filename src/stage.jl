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

_stage_rates_key(fn::AbstractString) =
    startswith(fn, "rates_") && endswith(fn, ".txt") ? String(chopprefix(chopsuffix(fn, ".txt"), "rates_")) : nothing

function _stage_rates_map(folder::AbstractString)::Dict{String,String}
    out = Dict{String,String}()
    for fn in readdir(folder)
        k = _stage_rates_key(fn)
        k === nothing && continue
        full = joinpath(folder, fn)
        isfile(full) || continue
        haskey(out, k) && throw(ArgumentError("duplicate rates key $(repr(k)) in folder $(repr(folder))"))
        out[k] = full
    end
    out
end

"""
    stage_combine_rates(folders, outfolder; grsm, noise, hidden_per_set=0, hidden_values=nothing,
                        ncoupling=1, coupling_values=nothing, key_select=:intersection,
                        root=".", root_out=root, dry_run=false, force=true) -> (; n, rows)

Combine key-matched `rates_<key>.txt` files across multiple input folders and write only combined
`rates_*.txt` files into `outfolder`.

`folders` is a list of rate folders (one per unit). For each folder `i`, provide:
- `grsm[i]`: number of GRSM rate columns per set
- `noise[i]`: number of noise columns per set

Per unit, one set width is `grsm[i] + noise[i]`. Each source rates table width must be divisible by
its set width; the inferred number of sets must match across all units for a given key (this is what
makes hierarchical stride insertion deterministic).

For each set `s`, output column order is:
1. unit 1 block for set `s`
2. unit 2 block for set `s`
3. ...
4. hidden block (`hidden_per_set` columns, optional)
5. coupling block (`ncoupling` columns)

So hidden/coupling are inserted per set (not only once at file end).
"""
function stage_combine_rates(
    folders::AbstractVector{<:AbstractString},
    outfolder::AbstractString;
    grsm::AbstractVector{<:Integer},
    noise::AbstractVector{<:Integer},
    hidden_per_set::Integer=0,
    hidden_values=nothing,
    ncoupling::Integer=1,
    coupling_values=nothing,
    key_select::Symbol=:intersection,
    root::AbstractString=".",
    root_out::AbstractString=root,
    dry_run::Bool=false,
    force::Bool=true,
)
    n_units = length(folders)
    n_units > 0 || throw(ArgumentError("folders must be non-empty"))
    length(grsm) == n_units || throw(ArgumentError("grsm length must match folders length"))
    length(noise) == n_units || throw(ArgumentError("noise length must match folders length"))
    key_select in (:intersection, :first) || throw(ArgumentError("key_select must be :intersection or :first"))
    hidden_per_set >= 0 || throw(ArgumentError("hidden_per_set must be >= 0"))
    ncoupling >= 0 || throw(ArgumentError("ncoupling must be >= 0"))

    src_abs = [_stage_resolve_src(f, String(root)) for f in folders]
    dst_abs = _stage_resolve_dst(outfolder, String(root_out); dry_run=dry_run, make=true)
    set_widths = Int.(grsm) .+ Int.(noise)
    any(set_widths .<= 0) && throw(ArgumentError("each grsm[i] + noise[i] must be > 0"))

    maps = [_stage_rates_map(f) for f in src_abs]
    matched_keys = if key_select === :intersection
        ks = isempty(maps) ? String[] : collect(Base.keys(maps[1]))
        for m in maps[2:end]
            ks = [k for k in ks if haskey(m, k)]
        end
        sort(ks)
    else
        sort(collect(Base.keys(maps[1])))
    end
    isempty(matched_keys) && throw(ArgumentError("no matching rates_<key>.txt files found for requested key_select=$(repr(key_select))"))

    hidden_vals = hidden_values === nothing ? fill(0.01, Int(hidden_per_set)) : collect(Float64, hidden_values)
    length(hidden_vals) == Int(hidden_per_set) || throw(ArgumentError("hidden_values length must equal hidden_per_set"))
    coupling_vals = coupling_values === nothing ? fill(0.0, Int(ncoupling)) : collect(Float64, coupling_values)
    length(coupling_vals) == Int(ncoupling) || throw(ArgumentError("coupling_values length must equal ncoupling"))

    rows = NamedTuple{(:key, :srcs, :dst, :nsets), Tuple{String,Vector{String},String,Int}}[]
    for key in matched_keys
        srcs = [haskey(m, key) ? m[key] : throw(ArgumentError("missing key $(repr(key)) in folder $(repr(src_abs[i]))")) for (i, m) in enumerate(maps)]
        data = Matrix{Float64}[]
        hdrs = Vector{String}[]
        for fp in srcs
            d, h = read_rates_table(fp)
            push!(data, d)
            push!(hdrs, h)
        end
        nrows = size(data[1], 1)
        nsets = Int[]
        for i in eachindex(data)
            size(data[i], 1) == nrows || throw(ArgumentError("row count mismatch for key $(repr(key)) across units"))
            w = size(data[i], 2)
            sw = set_widths[i]
            w % sw == 0 || throw(ArgumentError("key $(repr(key)) unit $(i): width $w not divisible by set width $sw (grsm=$(grsm[i]), noise=$(noise[i]))"))
            length(hdrs[i]) == w || throw(ArgumentError("key $(repr(key)) unit $(i): header length mismatch"))
            push!(nsets, div(w, sw))
        end
        all(==(nsets[1]), nsets) || throw(ArgumentError("key $(repr(key)) has inconsistent set counts across units: $(nsets)"))
        nset = nsets[1]

        out = zeros(Float64, nrows, 0)
        outh = String[]
        for s in 1:nset
            for i in eachindex(data)
                sw = set_widths[i]
                cols = (s - 1) * sw + 1:s * sw
                out = hcat(out, data[i][:, cols])
                append!(outh, hdrs[i][cols])
            end
            if hidden_per_set > 0
                out = hcat(out, repeat(reshape(hidden_vals, 1, :), nrows, 1))
                append!(outh, ["Hidden_$(j)_$s" for j in 1:Int(hidden_per_set)])
            end
            if ncoupling > 0
                out = hcat(out, repeat(reshape(coupling_vals, 1, :), nrows, 1))
                append!(outh, Int(ncoupling) == 1 ? ["Coupling_$s"] : ["Coupling_$(j)_$s" for j in 1:Int(ncoupling)])
            end
        end
        dst = joinpath(dst_abs, "rates_" * key * ".txt")
        if !dry_run
            !force && isfile(dst) && throw(ArgumentError("destination exists: $(repr(dst))"))
            write_rates_table(dst, out, outh)
        end
        push!(rows, (; key=key, srcs=srcs, dst=dst, nsets=nset))
    end
    (; n=length(rows), rows=rows)
end

function _stage_unit_from_spec(u)::NamedTuple{(:folder, :grsm, :noise), Tuple{String, Int, Int}}
    (u isa NamedTuple && haskey(u, :folder) && haskey(u, :grsm) && haskey(u, :noise)) ||
        throw(ArgumentError("each unit spec must provide :folder, :grsm, :noise"))
    (folder=String(u.folder), grsm=Int(u.grsm), noise=Int(u.noise))
end

"""
    stage_combine_rates_spec(spec; outfolder, root=".", root_out=root) -> NamedTuple

Resolve an abstract combine spec into the concrete arguments consumed by
[`stage_combine_rates`](@ref).

`spec` must include `units`, where each unit is a named tuple with:
- `folder` (String)
- `grsm` (Int)
- `noise` (Int)

Optional spec fields are forwarded when present:
`hidden_per_set`, `hidden_values`, `ncoupling`, `coupling_values`, `key_select`, `subfolder`.
"""
function stage_combine_rates_spec(
    spec::NamedTuple;
    outfolder::AbstractString,
    root::AbstractString=".",
    root_out::AbstractString=root,
)
    haskey(spec, :units) || throw(ArgumentError("combine spec must include :units"))
    units_raw = collect(spec.units)
    isempty(units_raw) && throw(ArgumentError("spec.units must be non-empty"))
    units = [_stage_unit_from_spec(u) for u in units_raw]
    folders = [u.folder for u in units]
    grsm = [u.grsm for u in units]
    noise = [u.noise for u in units]
    dst = haskey(spec, :subfolder) ? joinpath(String(outfolder), String(spec.subfolder)) : String(outfolder)
    (
        folders=folders,
        outfolder=dst,
        grsm=grsm,
        noise=noise,
        hidden_per_set=haskey(spec, :hidden_per_set) ? Int(spec.hidden_per_set) : 0,
        hidden_values=haskey(spec, :hidden_values) ? spec.hidden_values : nothing,
        ncoupling=haskey(spec, :ncoupling) ? Int(spec.ncoupling) : 1,
        coupling_values=haskey(spec, :coupling_values) ? spec.coupling_values : nothing,
        key_select=haskey(spec, :key_select) ? Symbol(spec.key_select) : :intersection,
        root=String(root),
        root_out=String(root_out),
    )
end

"""
    stage_combine_rates_batch(specs, outfolder; root=".", root_out=root, dry_run=false, force=true) -> (; n, specs)

Loop over abstract combine specs, resolve each with [`stage_combine_rates_spec`](@ref), and invoke
[`stage_combine_rates`](@ref). Returns one result block per spec.
"""
function stage_combine_rates_batch(
    specs::AbstractVector{<:NamedTuple},
    outfolder::AbstractString;
    root::AbstractString=".",
    root_out::AbstractString=root,
    dry_run::Bool=false,
    force::Bool=true,
)
    rows = NamedTuple[]
    for (i, s) in enumerate(specs)
        resolved = stage_combine_rates_spec(s; outfolder=outfolder, root=root, root_out=root_out)
        res = stage_combine_rates(
            resolved.folders,
            resolved.outfolder;
            grsm=resolved.grsm,
            noise=resolved.noise,
            hidden_per_set=resolved.hidden_per_set,
            hidden_values=resolved.hidden_values,
            ncoupling=resolved.ncoupling,
            coupling_values=resolved.coupling_values,
            key_select=resolved.key_select,
            root=resolved.root,
            root_out=resolved.root_out,
            dry_run=dry_run,
            force=force,
        )
        push!(rows, (; index=i, result=res))
    end
    (; n=length(rows), specs=rows)
end
