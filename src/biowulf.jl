# This file is part of StochasticGene.jl  

# biowulf.jl
# NIH Biowulf helpers: swarm files, fit scripts, run-spec presets (`makeswarm`, `makeswarm_models`,
# `makeswarmfiles`, `makeswarmfiles_h3_latent`, `write_run_spec_preset`, …). Run specs and emitted scripts
# align with [`fit`](@ref) / [`make_structures`](@ref) / [`load_options`](@ref) (MH, NUTS, ADVI): include
# `inference_method`, `device`, `parallel`, `gradient`, and shared budgets like `samplesteps` / `warmupsteps`.


"""
    sanitize_for_filename(s::AbstractString)

Replace characters that are unsafe in shell filenames (e.g. `,` and `|`) with `-`.
Underscore is reserved for preset filename fields (see io.jl).
Use when building swarm/.jl filenames from label or coupling spec (e.g. "24,35|35").
"""
sanitize_for_filename(s::AbstractString) = replace(replace(string(s), "," => "-"), "|" => "-")

"""
    default_model_key(G, R, S, insertstep)

Filename-safe key from [`create_modelstring`](@ref) (same default naming as model grids). Override with
`key=...` in each model row when you want arbitrary unique labels.
"""
default_model_key(G, R, S, insertstep) = sanitize_for_filename(create_modelstring(G, R, S, insertstep))

"""
    COUPLING_MODE_RECIPROCAL_DEFAULT

Default coupling sign modes (`:activate` / `:inhibit` / `:free`, or tuples per connection) keyed by
reciprocal coupling spec strings (e.g. `"3333"`, `"24,33|33"`). Used with `make_coupling` /
`makeswarm_coupled_reciprocal`-style workflows.
"""
const COUPLING_MODE_RECIPROCAL_DEFAULT = Dict{String,Union{Symbol,Tuple{Vararg{Symbol}}}}(
    "3333" => :inhibit,
    "3434" => :activate,
    "3535" => :inhibit,
    "R5R5" => :inhibit,
    "R4R4" => :activate,
    "R3R3" => :inhibit,
    "2424" => :activate,
    "2323" => :inhibit,
    "24,33|33" => (:inhibit, :activate, :inhibit),
    "24,35|35" => (:inhibit, :activate, :inhibit),
    "23,33|33" => :inhibit,
    "23,35|35" => :inhibit,
    "24R4" => (:activate, :activate),
    "2434" => (:activate, :activate),
    "24R5" => (:inhibit, :activate),
    "2435" => (:inhibit, :activate),
    "24R3" => (:inhibit, :activate),
    "2433" => (:inhibit, :activate),
    "23R4" => (:activate, :inhibit),
    "2334" => (:activate, :inhibit),
    "23R5" => (:inhibit, :inhibit),
    "2335" => (:inhibit, :inhibit),
    "23R3" => (:inhibit, :inhibit),
    "2333" => (:inhibit, :inhibit),
)

const _SWARM_ONLY_KEYS = (:nthreads, :nchains, :swarmfile, :juliafile, :filedir, :project, :sysimage, :src)

function _split_swarm_fit_kwargs(kwargs::Dict{Symbol,Any})
    swarm = Dict{Symbol,Any}()
    fit = Dict{Symbol,Any}()
    for (k, v) in kwargs
        if k in _SWARM_ONLY_KEYS
            swarm[k] = v
        else
            fit[k] = v
        end
    end
    return swarm, fit
end

"""
    write_run_spec_preset(resultfolder, key, run_spec; root=".")

Write `results/.../info_<key>.jld2` via [`write_run_spec_jld2`](@ref) (required for [`fit`](@ref)`(;
key=...)`) and a short `info_<key>.toml` marker. `run_spec` should contain the same keyword fields you
would pass to `fit`, including `resultfolder` and `root`.

# Notes
- [`read_run_spec`](@ref) loads the **JLD2** companion; keep complex types (e.g. `method`, `probfn`) in
  the dict so restarts match an interactive `fit` call.
- Before writing, [`normalize_trace_specs_legacy_t_end!`](@ref) rewrites legacy huge `t_end` (old `1e30`
  sentinel) to `-1.0` so merged specs from older `info_<key>.jld2` files do not keep invalid values.
"""
function write_run_spec_preset(resultfolder::AbstractString, key::AbstractString, run_spec::Dict{Symbol,Any}; root::AbstractString=".")
    normalize_trace_specs_legacy_t_end!(run_spec)
    rr = folder_path(resultfolder, root, "results")
    mkpath(rr)
    path_toml = joinpath(rr, "info_" * string(key) * ".toml")
    write_run_spec_jld2(path_toml, run_spec)
    j2 = basename(info_jld2_path(path_toml))
    open(path_toml, "w") do io
        println(io, "# Pre-generated run specification. Machine-readable state is in `", j2, "` (read by `fit(; key=...)`).")
        println(io, "key = ", repr(string(key)))
    end
    return path_toml
end

"""
    makeswarm_models(models::Vector{<:NamedTuple}; transitions=nothing, kwargs...)
    makeswarm_models(Gvec, Rvec, Svec, insertvec; combine=:product, kwargs...)

Batch workflow for **many models, few genes**: merge `fit_default_spec` with shared kwargs (leave
`priormean` empty to use `fit` defaults: for trace single-unit models, `make_structures` applies
`prior_ratemean_trace` / structured `priorcv` automatically), write one `write_run_spec_preset` per
model, then write swarm + one `fitscript_<key>.jl` per key. Each script is a multi-line `fit(; …)` with
the **full merged run spec** (same content as the companion `info_<key>.jld2`), not a minimal `fit(;
key=…)` plus a few overrides.

Each element of `models` should include `G`, `R`, `S`, `insertstep` and usually `transitions` (or pass
`transitions=...` once for all rows). Optional `key` in a row overrides [`default_model_key`](@ref).

The `combine=:product` form iterates `Base.product(Gvec, Rvec, Svec, insertvec)`; `combine=:zip` requires
four vectors of equal length.

Use **`transitionsvec`** (same length as **`Gvec`**) so each gene-state count **`Gvec[i]`** gets its own
**`transitionsvec[i]`** (required when sweeping multiple **`G`** with different topologies). **`Gvec`**
must have **unique** entries. Omit **`transitionsvec`** to use a single shared **`transitions=`** for
every row (positional / kw to [`makeswarm_models`](@ref)).

# Keyword split
Keywords `nthreads`, `nchains`, `swarmfile`, `juliafile`, `filedir`, `project`, `sysimage`, `src` are
swarm-only; all others are merged into each run spec (and into `fit_default_spec()`), including inference
options consumed by [`load_options`](@ref) (`inference_method`, `device`, `parallel`, `gradient`, `n_mc`, …).

# Coupled follow-up
For slow coupled fits, fit single units first, combine rates with [`create_combined_file`](@ref), then
run coupled models using those starts (see README Workflow 2).
"""
function makeswarm_models(models::AbstractVector{<:NamedTuple}; transitions=nothing, unique_model_keys::Bool=false, kwargs...)
    swarm_kw, fit_kw = _split_swarm_fit_kwargs(Dict{Symbol,Any}(kwargs))
    base = merge(fit_default_spec(), fit_kw)
    if transitions !== nothing && !haskey(base, :transitions)
        base[:transitions] = transitions
    end
    keys_out = String[]
    specs_by_key = Dict{String, Dict{Symbol,Any}}()
    for m in models
        spec = merge(copy(base), Dict{Symbol,Any}(pairs(m)))
        if !haskey(spec, :transitions) || spec[:transitions] === nothing
            throw(ArgumentError("each model needs transitions= or pass transitions= for all"))
        end
        key = if get(spec, :key, nothing) !== nothing
            string(spec[:key])
        else
            default_model_key(spec[:G], spec[:R], spec[:S], spec[:insertstep])
        end
        if unique_model_keys
            suf = string(hash((key, time_ns(), objectid(models))), base=16)
            key = sanitize_for_filename(string(key, "-", suf[1:min(8, length(suf))]))
        end
        spec[:key] = key
        push!(keys_out, key)
        write_run_spec_preset(spec[:resultfolder], key, spec; root=spec[:root])
        specs_by_key[key] = spec
    end
    # Omit default `nchains` so `_makeswarm_with_run_specs` can match `-p` to each spec (see `_nchains_for_swarm_line`).
    swarm_out = merge(Dict{Symbol,Any}(
        :nthreads => get(swarm_kw, :nthreads, 1),
        :swarmfile => get(swarm_kw, :swarmfile, "fit"),
        :juliafile => get(swarm_kw, :juliafile, "fitscript"),
        :filedir => get(swarm_kw, :filedir, "."),
        :project => get(swarm_kw, :project, ""),
        :sysimage => get(swarm_kw, :sysimage, ""),
        :src => get(swarm_kw, :src, ""),
    ), swarm_kw)
    _makeswarm_with_run_specs(keys_out, swarm_out, specs_by_key)
end

function makeswarm_models(
    Gvec::AbstractVector{<:Integer},
    Rvec::AbstractVector{<:Integer},
    Svec::AbstractVector{<:Integer},
    insertvec::AbstractVector{<:Integer};
    combine::Symbol=:product,
    unique_model_keys::Bool=false,
    transitionsvec=nothing,
    kwargs...,
)
    g_tr = nothing
    if transitionsvec !== nothing
        length(transitionsvec) == length(Gvec) ||
            throw(ArgumentError("transitionsvec length ($(length(transitionsvec))) must match Gvec length ($(length(Gvec)))"))
        length(unique(Gvec)) == length(Gvec) ||
            throw(ArgumentError("transitionsvec requires unique entries in Gvec (duplicate G would be ambiguous)"))
        g_tr = Dict{Int,Any}(Gvec[i] => transitionsvec[i] for i in eachindex(Gvec))
    end
    if combine == :product
        # Array comprehensions over ProductIterator allocate an ndims(A) array matching the
        # product axes; vec(...) yields a Vector{NamedTuple} for makeswarm_models(models; ...).
        if g_tr === nothing
            models = vec([(G=gs[1], R=gs[2], S=gs[3], insertstep=gs[4]) for gs in Base.product(Gvec, Rvec, Svec, insertvec)])
        else
            models = vec([(G=gs[1], R=gs[2], S=gs[3], insertstep=gs[4], transitions=g_tr[gs[1]]) for gs in Base.product(Gvec, Rvec, Svec, insertvec)])
        end
    elseif combine == :zip
        n = length(Gvec)
        length(Rvec) == length(Svec) == length(insertvec) == n || throw(ArgumentError("combine=:zip requires equal-length Gvec, Rvec, Svec, insertvec"))
        if g_tr === nothing
            models = [(G=Gvec[i], R=Rvec[i], S=Svec[i], insertstep=insertvec[i]) for i in 1:n]
        else
            models = [(G=Gvec[i], R=Rvec[i], S=Svec[i], insertstep=insertvec[i], transitions=transitionsvec[i]) for i in 1:n]
        end
    else
        throw(ArgumentError("combine must be :product or :zip"))
    end
    makeswarm_models(models; unique_model_keys=unique_model_keys, kwargs...)
end

"""
    makeswarm(keys::Vector{String}; <keyword arguments>)
    makeswarm(; key::String, <keyword arguments>)

Write a swarm file and one fit script per key for Biowulf. Each swarm line runs one script:
`julia -t nthreads -p nchains fitscript_<key>.jl`. Each script calls `fit(; key=key, ...)` so the
run is defined by `info_<key>.toml` (if present) plus any overrides.

# Arguments
- `keys` or `key`: run key(s); each key gets one swarm line and one script `fitscript_<key>.jl`.
- `nthreads=1`: Julia threads per process.
- `nchains=2`: number of parallel chains (passed to `-p`).
- `swarmfile="fit"`: base name for the swarm file (writes `swarmfile.swarm`).
- `juliafile="fitscript"`: base name for scripts (writes `juliafile_<key>.jl`).
- `filedir="."`: directory for swarm and script files.
- `project=""`, `sysimage=""`: optional `--project` and `--sysimage` for the julia command.
- `src=""`: path to StochasticGene src (for script prolog; empty if package is installed).
- Overrides (optional): `resultfolder`, `root`, `maxtime`, `samplesteps`, `warmupsteps`, `inference_method`,
  `device`, `parallel` / `parallelism`, `gradient`, etc. are written into each script so
  `fit(; key=key, resultfolder=..., ...)` uses them. Serializable types: `String`, `Number`, `Bool`, and `Symbol`
  (symbols are emitted as `repr`, e.g. `inference_method=:nuts`).

# Examples
```julia
makeswarm(["33il", "44il"]; filedir="swarm", resultfolder="HCT116_test", root=".")
makeswarm(; key="33il", resultfolder="HCT116_test", maxtime="120m")
```
"""
function makeswarm(keys::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript", filedir=".", project="", sysimage="", src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, warmupsteps=nothing, kwargs...)
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    isempty(keys) && return
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    write_swarmfile_keys(joinpath(filedir, sfile), nchains, nthreads, juliafile, keys, project, sysimage)
    for k in keys
        scriptpath = joinpath(filedir, juliafile * "_" * sanitize_for_filename(k) * ".jl")
        write_fitfile_key(scriptpath, k; src=src, resultfolder=resultfolder, root=root, maxtime=maxtime, samplesteps=samplesteps, warmupsteps=warmupsteps, kwargs...)
    end
end

function makeswarm(; key::String, kwargs...)
    isempty(key) && error("makeswarm(; key=...) requires a non-empty key")
    makeswarm([key]; kwargs...)
end

"""Like [`folder_path`](@ref) but never prints to stdout (used when probing many candidate paths)."""
function _folder_path_quiet(folder::AbstractString, root::AbstractString, folderatetype::AbstractString="")
    isempty(strip(folder)) && return string(folder)
    f = string(folder)
    ispath(f) && return f
    f = joinpath(root, folder)
    ispath(f) && return f
    parts = splitpath(folder)
    nested_under_type = !isempty(folderatetype) && length(parts) >= 1 && parts[1] == folderatetype
    nested_under_type ? joinpath(root, folder) : joinpath(root, folderatetype, folder)
end

"""
    _find_existing_info_toml(key, resultfolder, root) -> String or nothing

Same resolution as `fit(; key=...)`: `joinpath(folder_path(resultfolder, root, \"results\"), \"info_<key>.toml\")`.
"""
function _find_existing_info_toml(key::AbstractString, resultfolder, root)
    stem = "info_" * string(key) * ".toml"
    isempty(strip(string(resultfolder))) && return nothing
    rr = _folder_path_quiet(resultfolder, root, "results")
    p = abspath(joinpath(rr, stem))
    isfile(p) ? p : nothing
end

"""Parse `H3-t1-t2` or `H3-t1-t2-h` batch keys into `(t1, t2)` for hidden-latent coupling; else `nothing`."""
function parse_h3_transition_key(key::AbstractString)
    s = String(strip(key))
    m = match(r"^H3-(\d+)-(\d+)(?:-h)?$", s)
    m === nothing && return nothing
    return (parse(Int, m[1]), parse(Int, m[2]))
end

"""Shared three-unit transitions for H3 latent (units 1–2 observed GRSM, unit 3 hidden G-only all-pairs)."""
const _H3_LATENT_TRANSITIONS = (
    ([1, 2], [2, 1], [2, 3], [3, 2]),
    ([1, 2], [2, 1], [2, 3], [3, 2]),
    ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]),
)

function _apply_h3_latent_key_overlay!(spec::Dict{Symbol,Any}, t1::Int, t2::Int; coupling_ranges::Tuple=(:activate, :activate))
    spec[:G] = (3, 3, 3)
    spec[:R] = (3, 3, 0)
    spec[:S] = (0, 0, 0)
    spec[:insertstep] = (1, 1, 0)
    spec[:transitions] = _H3_LATENT_TRANSITIONS
    spec[:coupling] = make_coupling_hidden_latent(t1, t2; coupling_ranges=coupling_ranges)
    spec[:noisepriors] = ([0.0, 0.1, 0.5, 0.15], [0.0, 0.1, 0.9, 0.2], Float64[])
    spec[:probfn] = (prob_Gaussian, prob_Gaussian)
    el = get(spec, :elongationtime, (20.0, 5.0))
    if el isa Real
        spec[:elongationtime] = (Float64(el), Float64(el))
    elseif el isa Tuple && length(el) >= 2
        spec[:elongationtime] = (Float64(el[1]), Float64(el[2]))
    else
        spec[:elongationtime] = (20.0, 5.0)
    end
    spec[:zeromedian] = Bool[true, true]
    spec[:trace_specs] = default_trace_specs_for_coupled((5.0 / 3.0, 1.0, -1.0), spec[:zeromedian], 2)
    spec[:decayrate] = 1.0
    spec[:splicetype] = "full"
    spec[:prerun] = get(spec, :prerun, 0.0)
    return spec
end

function _driver_write_specs_and_makeswarm(run_keys::Vector{String}, fit_base::Dict{Symbol,Any}, swarm_kw::Dict{Symbol,Any};
        merge_existing_info::Bool=false,
        fit_kw::Dict{Symbol,Any}=Dict{Symbol,Any}(),
        warn_missing_info::Bool=true,
        apply_h3_latent_overlay::Bool=false,
        h3_coupling_ranges::Tuple=(:activate, :activate),
    )
    n_missing = 0
    specs_by_key = Dict{String, Dict{Symbol,Any}}()
    for key in run_keys
        spec = copy(fit_base)
        if merge_existing_info
            p = _find_existing_info_toml(key, get(spec, :resultfolder, ""), get(spec, :root, "."))
            if p !== nothing
                loaded = read_run_spec(p)
                spec = merge(spec, loaded)
                for k in (
                    :datapath, :resultfolder, :maxtime, :root,
                    :inference_method, :device, :parallel, :parallelism, :gradient,
                    :samplesteps, :warmupsteps, :n_mc, :maxiter, :n_samples, :n_adapts,
                    :nuts_delta, :δ, :fd_ε, :nuts_fd_ε, :verbose, :progress, :nuts_verbose, :nuts_progress,
                    :time_limit, :advi_time_limit, :σ_floor, :advi_σ_floor, :init_s_raw, :advi_init_s_raw,
                    :advi_verbose, :advi_n_mc, :zygote_trace,
                )
                    haskey(fit_base, k) && (spec[k] = fit_base[k])
                end
                # Explicit `makeswarmfiles(...; kw...)` fit keys win over loaded info (merge loaded first for priors).
                # Use `Base.keys` — a parameter named `keys` would shadow `keys(...)`.
                for k in Base.keys(fit_kw)
                    haskey(fit_base, k) && (spec[k] = fit_base[k])
                end
            else
                n_missing += 1
            end
        end
        if apply_h3_latent_overlay
            h3p = parse_h3_transition_key(key)
            if h3p !== nothing
                _apply_h3_latent_key_overlay!(spec, h3p[1], h3p[2]; coupling_ranges=h3_coupling_ranges)
                for k in Base.keys(fit_kw)
                    haskey(fit_base, k) && (spec[k] = fit_base[k])
                end
            end
        end
        spec[:key] = key
        write_run_spec_preset(spec[:resultfolder], key, spec; root=spec[:root])
        specs_by_key[key] = spec
    end
    if merge_existing_info && n_missing > 0 && warn_missing_info
        @warn "makeswarmfiles: no info_<key>.toml (with JLD2) found for $n_missing of $(length(run_keys)) keys under resultfolder; run specs use merged defaults only — add prior `info_<key>` files next to those fits or pass full model kwargs. Silence when expected: `warn_missing_info=false`." n_missing=n_missing nkeys=length(run_keys)
    end
    _makeswarm_with_run_specs(run_keys, swarm_kw, specs_by_key)
    return run_keys
end

"""
    _extract_coupling_columns_from_row(row, col_names::Vector{String})

Extract coupling columns (e1, e1s, e2, e2s, ge, ges) from a DataFrame row,
supporting flexible column naming and missing columns.

Column name patterns searched (in order):
- enhancer_to_gene_{1,2}, enhancer_sign_{1,2}
- gene_to_enhancer, gene_to_enhancer_sign
- e1, e1s, e2, e2s, ge, ges
- 1,2,3,4,5,6 (positional fallback)

Missing columns default to empty string.
"""
function _extract_coupling_columns_from_row(row, col_names::Vector{String})
    function get_col(patterns::Vector{String}; default="")
        for pat in patterns
            pat in col_names && return row[Symbol(pat)]
        end
        return default
    end
    
    e1 = get_col(["enhancer_to_gene_1", "e1"])
    e1s = get_col(["enhancer_sign_1", "e1s"])
    e2 = get_col(["enhancer_to_gene_2", "gene_to_enhancer", "e2"])
    e2s = get_col(["enhancer_sign_2", "gene_to_enhancer_sign", "e2s"])
    ge = get_col(["background_gene", "ge"])
    ges = get_col(["background_gene_sign", "ges"])
    
    return e1, e1s, e2, e2s, ge, ges
end

function _makeswarmfiles_coupled_models_csv(
    csv_path::AbstractString;
    filedir::AbstractString,
    key_col::String,
    skip_empty::Bool,
    hierarchical_modes,
    coupled_csv_cols::NTuple{6,Int},
    coupled_emit_legacy_r_variants::Bool,
    coupled_tie_rsum::Bool,
    G::Tuple,
    R::Tuple,
    S::Tuple,
    insertstep::Tuple,
    transitions::Tuple,
    merged::Dict{Symbol,Any},
    swarm_kw::Dict{Symbol,Any},
)
    df = DataFrame(CSV.File(csv_path))
    key_col in names(df) || throw(ArgumentError("missing key column '$key_col' in $csv_path"))
    
    col_names = String.(names(df))
    
    keys_ordered = String[]
    specs_by_key = Dict{String, Dict{Symbol,Any}}()
    for hierarchical in hierarchical_modes
        for row in eachrow(df)
            kn = strip(string(row[Symbol(key_col)]))
            (skip_empty && isempty(kn)) && continue
            kn == "Model_name" && continue
            key_base = replace(kn, " " => "-")
            
            # Extract coupling columns by name (flexible order/naming)
            e1, e1s, e2, e2s, ge, ges = _extract_coupling_columns_from_row(row, col_names)
            has_lr = coupled_emit_legacy_r_variants && csv_row_has_legacy_r(e1, e2, ge)
            variant_kinds = has_lr ? ((:rsum, "-Rsum"), (:rany, "-Rany")) : ((:none, ""),)
            for (vkind, vsuffix) in variant_kinds
                e1u = vkind === :none ? e1 : replace_csv_cell_legacy_r(e1, vkind)
                e2u = vkind === :none ? e2 : replace_csv_cell_legacy_r(e2, vkind)
                geu = vkind === :none ? ge : replace_csv_cell_legacy_r(ge, vkind)
                key = key_base * (hierarchical ? "-h" : "") * vsuffix
                ts0 = get(merged, :trace_specs, [])
                interval, tracetime = if !isempty(ts0)
                    s1 = ts0[1]
                    Float64(get(s1, :interval, 5.0 / 3.0)), Float64(get(s1, :t_end, -1.0))
                else
                    ti = get(merged, :traceinfo, (5.0 / 3.0, 1.0, -1.0))
                    Float64(ti[1]), (length(ti) >= 3 ? Float64(ti[3]) : -1.0)
                end
                dc = get(merged, :datacond, ["gene", "enhancer"])
                dc_vec = dc isa AbstractVector ? String[string(x) for x in dc] : String[string(dc)]
                pfn0 = get(merged, :probfn, nothing)
                pfn = if pfn0 === nothing || pfn0 isa Function
                    (prob_Gaussian, prob_Gaussian)
                elseif pfn0 isa Tuple && length(pfn0) >= 2
                    pfn0
                else
                    (prob_Gaussian, prob_Gaussian)
                end
                np = get(merged, :noisepriors, nothing)
                if np === nothing || np == [] || (np isa AbstractVector && isempty(np))
                    np = ([0.0, 0.1, 0.5, 0.15], [0.0, 0.1, 0.9, 0.2])
                end
                el = get(merged, :elongationtime, (20.0, 5.0))
                elongationtime2 = if el isa Tuple && length(el) >= 2
                    (Float64(el[1]), Float64(el[2]))
                elseif el isa Real
                    (Float64(el), Float64(el))
                else
                    (20.0, 5.0)
                end
                spec0 = build_coupled_fit_spec_from_csv_cells(
                    e1u, e1s, e2u, e2s, geu, ges, key;
                    G=G,
                    R=R,
                    S=S,
                    insertstep=insertstep,
                    transitions=transitions,
                    datapath=string(merged[:datapath]),
                    resultfolder=string(merged[:resultfolder]),
                    root=string(get(merged, :root, ".")),
                    gene=string(get(merged, :gene, "MYC")),
                    cell=string(get(merged, :cell, "HBEC")),
                    datacond=dc_vec,
                    noisepriors=np,
                    elongationtime=elongationtime2,
                    hierarchical=hierarchical,
                    interval=interval,
                    tracetime=tracetime,
                    zeromedian=get(merged, :zeromedian, Bool[true, true]),
                    probfn=pfn,
                    trace_specs=get(merged, :trace_specs, []),
                    tie_rsum=coupled_tie_rsum,
                    maxtime=Float64(get(merged, :maxtime, 60.0)),
                    nchains=Int(get(merged, :nchains, 16)),
                    samplesteps=Int(get(merged, :samplesteps, 100_000)),
                    warmupsteps=Int(get(merged, :warmupsteps, 0)),
                    propcv=Float64(get(merged, :propcv, 0.05)),
                    ratetype=string(get(merged, :ratetype, "median")),
                    datacol=Int(get(merged, :datacol, 3)),
                    writesamples=Bool(get(merged, :writesamples, false)),
                    prerun=Float64(get(merged, :prerun, 0.0)),
                    splicetype=string(get(merged, :splicetype, "")),
                    decayrate=Float64(get(merged, :decayrate, 1.0)),
                )
                spec0 === nothing && continue
                spec = merge(merged, spec0)
                spec[:key] = key
                write_run_spec_preset(spec[:resultfolder], key, spec; root=spec[:root])
                specs_by_key[key] = spec
                push!(keys_ordered, key)
            end
        end
    end
    isempty(keys_ordered) && throw(ArgumentError("coupled_models_csv: no runnable rows (empty coupling for all rows?)"))
    _makeswarm_with_run_specs(keys_ordered, swarm_kw, specs_by_key)
    return keys_ordered
end

"""
    makeswarmfiles_driver(; transitions=([1, 2], [2, 1]), transitionsvec=nothing, Gvec=[3], Rvec=[3, 4, 5], Svec=[0], insertvec=[1], combine=:product, models=nothing, kwargs...)

Developer driver for **single-unit** GRSM batches using [`makeswarm_models`](@ref): define model
layouts, write `info_<key>.jld2/.toml`, and generate swarm + fit scripts. For **coupled key-first**
workflows (CSV / explicit keys / H3 grid), use [`makeswarmfiles`](@ref) instead.

Use `models` (vector of NamedTuples with `G,R,S,insertstep` and optional `key`) for explicit
control, or leave `models=nothing` to generate from `Gvec/Rvec/Svec/insertvec`.
All `kwargs...` are forwarded to [`makeswarm_models`](@ref) (e.g. data paths, fit options,
`swarmfile`, `juliafile`, `filedir`, `project`, `sysimage`, etc.).
"""
function makeswarmfiles_driver(; transitions=([1, 2], [2, 1]), transitionsvec=nothing, Gvec=[3], Rvec=[3, 4, 5], Svec=[0], insertvec=[1], combine::Symbol=:product, models=nothing, unique_model_keys::Bool=false, kwargs...)
    if models === nothing
        return makeswarm_models(Gvec, Rvec, Svec, insertvec; combine=combine, transitions=transitions, transitionsvec=transitionsvec, unique_model_keys=unique_model_keys, kwargs...)
    else
        transitionsvec !== nothing &&
            throw(ArgumentError("makeswarmfiles_driver: transitionsvec only applies to Gvec/Rvec/Svec/insertvec grids, not models="))
        return makeswarm_models(models; transitions=transitions, unique_model_keys=unique_model_keys, kwargs...)
    end
end

"""
    makeswarmfiles(; filedir=".", csv_file="", datapath, folder, maxtime=30.0, hierarchical_modes=(true,false), key_col="Model_name", skip_empty=true, base_keys=nothing, h3_transition_range=nothing, coupled=nothing, merge_existing_info=true, warn_missing_info=true, h3_latent=false, h3_coupling_ranges=(:activate, :activate), warn_coupled_csv_shape=true, models=nothing, Gvec=nothing, Rvec=nothing, Svec=nothing, insertvec=nothing, combine=:product, transitions=([1,2],[2,1]), transitionsvec=nothing, kwargs...)

Unified entry for [`write_run_spec_preset`](@ref) plus swarm and `fitscript_<key>.jl` under `filedir`
(one `info_<key>.jld2` and marker TOML per key). Single-unit `Gvec`…`insertvec` / `models` batches use
[`write_fitscript_tracejoint_key`](@ref) so each script contains the **full** default-merged `fit(;
…)`; string-key / coupled-CSV paths do the same. Direct [`makeswarm`](@ref) with a vector of keys
still writes a compact one-line `fit(; key=…)` when used alone (no presets).

Pick **exactly one** workflow below (mutually exclusive sources).

# Coupled / key-first models (string keys + [`combined_rates_key`](@ref))

Use when each run is identified by a **predefined key** (e.g. `H3-3-3`, combined-rate filenames).

Provide keys in **one** of these ways:

1. **`csv_file`** (alias **`csv`**): path to a UTF-8 CSV with a **key column** (default `Model_name`; set `key_col`).
   **Only that column is used for batch keys.** This is **not** the same as the
   `makeswarmfiles(csv_path, filedir)` / `makeswarmfiles_from_csv` workflow in **`makescriptcoupled.jl`**:
   there, each row’s **coupling** is built from **seven columns** (model name + six coupling cells:
   enhancer→gene ×2 + signs, gene→enhancer + sign). That logic lives in the script (calls
   `StochasticGene.makeswarm_coupled` / `_csv_row_to_connections_*` per row and writes **full** fit scripts).
   Here, coupling and priors for each key must come from **`kwargs...`**, from merged
   **`info_<key>.jld2`** when **`merge_existing_info=true`**, or from your own precomputed specs—not from
   extra CSV columns.
2. **`base_keys`**: `Vector` of strings, e.g. `["H3-3-3", "H3-5-5"]`.
3. **`h3_transition_range`**: e.g. `1:5` builds all `H3-t1-t2` base keys (same idea as
   [`makeswarmfiles_h3_latent`](@ref)). This automatically enables the **H3 hidden-latent**
   three-unit skeleton (`G=(3,3,3)`, `make_coupling_hidden_latent(t1,t2)`, …) on each key; it is **not**
   the two-unit `fit_coupled_default_spec` layout. For CSV / `base_keys` only, pass **`h3_latent=true`**
   (and optionally **`h3_coupling_ranges=(:activate, :activate)`**) so the same overlay applies.

Set **`coupled`** for this workflow (defaults to **`true`** when omitted):

- **`coupled=true`** or **`coupled=nothing`**: run specs start from [`fit_coupled_default_spec`](@ref) (trace-joint / coupled batch baseline) merged with `kwargs...`.
- **`coupled=false`**: run specs start from [`fit_default_spec`](@ref) (single-unit GRSM baseline). Use when keys refer to **uncoupled** models (e.g. CSV of `default_model_key` strings) but you still want string-key batching.

If **`merge_existing_info=true`** (default), each key also merges an existing
`info_<key>.toml` / JLD2 found under `folder` (`resultfolder`; same resolution as `fit(; key=...)`),
so layouts, priors, and trace settings from a prior fit are preserved; explicit
`folder` (`resultfolder`), `maxtime`, and `root` from this call still win, as do
any fit keywords you pass (including **`datapath`** when you set it here).

Optional **`hierarchical_modes`**: tuple of `Bool` (default `(true, false)`) applies
[`combined_rates_key`](@ref) for each mode (e.g. with/without `-h`). Use `(false,)` for non-hierarchical only.

# Single-unit GRSM sweeps (keys from `G,R,S,insertstep`)

Use when each run is a **different uncoupled** layout and keys come from [`default_model_key`](@ref)
or per-row `key`:

- **`models`**: `Vector` of named tuples with `G`, `R`, `S`, `insertstep`, optional `key`, plus
  `transitions` (globally or per row as supported by [`makeswarm_models`](@ref)).
- **Or** all four **`Gvec`, `Rvec`, `Svec`, `insertvec`**: grid with `combine=:product` or `:zip`.
  With **`transitionsvec`** (same length as **`Gvec`**; **`Gvec`** must have unique entries), use
  **`transitionsvec[i]`** for **`Gvec[i]`**; otherwise one shared **`transitions`** applies to every row.
  Not used with **`models=`** (put **`transitions`** in each named tuple instead).

Delegates to [`makeswarmfiles_driver`](@ref) → [`makeswarm_models`](@ref), which always uses
[`fit_default_spec`](@ref) as the batch baseline (**single-unit**). **Do not** combine with
`csv_file`, `base_keys`, or `h3_transition_range`. Do not pass **`coupled=true`** here.

# Common keywords

- **`filedir`**: directory for `*.swarm` and `fitscript_<key>.jl`.
- **`datapath`, `folder`**: merged into the run spec (`folder` is written as `resultfolder`).
  There is **no** repository default for **`datapath`**—set it to your machine’s data root. Default
  **`folder="coupled"`**; if that path is missing, [`folder_path`](@ref) logs a warning naming the folder (e.g. `coupled`)—not a Julia “variable not found” error.
- **`maxtime`**: wall-time hint for scripts.
- **`coupled`**: only for the **key-first** branch (`csv_file` / `base_keys` / `h3_transition_range`); see above.
- **`merge_existing_info`**: if `true`, merge each key’s existing `info_<key>` run spec when found (see above).
  When some keys have no prior file, a **`@warn`** is emitted unless **`warn_missing_info=false`** (e.g. you supply full model kwargs for every key).
- **`h3_latent`**, **`h3_coupling_ranges`**: for CSV / `base_keys` whose entries look like `H3-t1-t2`, set **`h3_latent=true`** to apply the same three-unit hidden-latent skeleton as **`h3_transition_range`** (see above). **`h3_coupling_ranges`** defaults to **`(:activate, :activate)`** (passed to [`make_coupling_hidden_latent`](@ref)).
- **`warn_coupled_csv_shape`**: if `true` (default), and the CSV has many columns (≥7), emit a **one-time**
  hint that this API does not parse `Coupled_models_to_test.csv`-style coupling columns; set to `false` to
  silence when your file is wide for other reasons.
- **`kwargs...`**: forwarded; swarm-only keys (`nthreads`, `nchains`, `project`, `juliafile`, …) go to
  [`makeswarm`](@ref), the rest merge into each run spec (see [`_split_swarm_fit_kwargs`](@ref)).
  If you omit **`nchains`** in **`kwargs`**, the swarm file’s **`-p`** is set from each run spec’s **`nchains`**
  (e.g. **16** from [`fit_coupled_default_spec`](@ref)), matching the generated **`fit(; …)`** scripts.

# Removed API

Positional `makeswarmfiles(csv_path, filedir; ...)` is **not** supported — use **only** the keyword form below.

# Coupled CSV (per-row coupling, `Coupled_models_to_test` layout)

Prefer **[`makeswarmfiles_coupled`](@ref)**`(; csv=path, filedir=..., ...)` or **`makeswarmfiles(; csv_file=path, coupled_csv=true, ...)`** — **`csv`** and **`csv_file`** are the same argument (path to the CSV).
Legacy **`coupled_models_csv=path`** is still accepted (same as **`coupled_csv=true`** with that path). See `coupled_csv.jl` for columns: **`key_col`** plus six coupling columns (**`coupled_csv_cols`**, default `(2,…,7)`). Shared tuple kwargs **`coupled_G`**, **`coupled_R`**, **`coupled_S`**, **`coupled_insertstep`**, **`coupled_transitions`** set the two-unit skeleton. Mutually exclusive with keys-only **`csv`/`csv_file`** without **`coupled_csv`**, **`base_keys`**, and **`h3_transition_range`**.

# Uncoupled grids and unique keys (`makeswarm_models`)

For **`models=`** or **`Gvec`…`insertvec`**, pass **`unique_model_keys=true`** on [`makeswarm_models`](@ref) / [`makeswarmfiles_driver`](@ref) to append a random suffix so keys stay unique when layouts repeat. For **`Gvec`** grids only, **`transitionsvec`** (same length as **`Gvec`**, **`Gvec`** unique) sets per-**`G`** transitions; otherwise use one shared **`transitions`**.
"""
function _resolve_csv_kw(csv_file::AbstractString, csv::Union{Nothing,AbstractString})
    f = strip(csv_file)
    c = csv === nothing ? "" : strip(string(csv))
    if !isempty(f) && !isempty(c) && f != c
        throw(ArgumentError("csv and csv_file disagree; pass one path (csv_file=$(repr(f)), csv=$(repr(c)))"))
    end
    return !isempty(f) ? f : c
end

function makeswarmfiles(;
        filedir::AbstractString=".",
        csv_file::AbstractString="",
        csv::Union{Nothing,AbstractString}=nothing,
        csv_path::Union{Nothing,AbstractString}=nothing,
        coupled_csv::Bool=false,
        coupled_models_csv::Union{Nothing,AbstractString}=nothing,
        datapath::AbstractString="",
        folder="coupled",
        maxtime=30.0,
        hierarchical_modes=(true, false),
        key_col::String="Model_name",
        skip_empty::Bool=true,
        base_keys=nothing,
        h3_transition_range=nothing,
        models=nothing,
        Gvec=nothing,
        Rvec=nothing,
        Svec=nothing,
        insertvec=nothing,
        combine::Symbol=:product,
        transitions=([1, 2], [2, 1]),
        transitionsvec=nothing,
        coupled::Union{Nothing,Bool}=nothing,
        merge_existing_info::Bool=true,
        warn_missing_info::Bool=true,
        h3_latent::Bool=false,
        h3_coupling_ranges::Tuple=(:activate, :activate),
        warn_coupled_csv_shape::Bool=true,
        coupled_csv_cols::NTuple{6,Int}=(2, 3, 4, 5, 6, 7),
        coupled_emit_legacy_r_variants::Bool=true,
        coupled_tie_rsum::Bool=true,
        coupled_G::Tuple=(3, 3),
        coupled_R::Tuple=(3, 3),
        coupled_S::Tuple=(0, 0),
        coupled_insertstep::Tuple=(1, 1),
        coupled_transitions::Tuple=(([1, 2], [2, 1], [2, 3], [3, 2]), ([1, 2], [2, 1], [2, 3], [3, 2])),
        kwargs...,
    )
    path_csv = _resolve_csv_kw(csv_file, csv)
    if csv_path !== nothing
        !isempty(path_csv) &&
            throw(ArgumentError("makeswarmfiles: pass only one of csv_file/csv or csv_path"))
        path_csv = strip(string(csv_path))
    end
    cms = coupled_models_csv === nothing ? "" : strip(string(coupled_models_csv))
    if !isempty(cms)
        if !isempty(path_csv) && path_csv != cms
            throw(ArgumentError("makeswarmfiles: csv/csv_file/csv_path and coupled_models_csv disagree; use one path"))
        end
        path_csv = isempty(path_csv) ? cms : path_csv
    end
    use_coupled_rows = coupled_csv || !isempty(cms)

    if use_coupled_rows
        models !== nothing && throw(ArgumentError("makeswarmfiles: coupled CSV cannot be used with models="))
        (Gvec !== nothing || Rvec !== nothing || Svec !== nothing || insertvec !== nothing) &&
            throw(ArgumentError("makeswarmfiles: coupled CSV cannot be used with Gvec/Rvec/Svec/insertvec"))
        base_keys !== nothing &&
            throw(ArgumentError("makeswarmfiles: coupled CSV cannot be used with base_keys"))
        h3_transition_range !== nothing &&
            throw(ArgumentError("makeswarmfiles: coupled CSV cannot be used with h3_transition_range"))
        isempty(path_csv) &&
            throw(ArgumentError("makeswarmfiles: coupled CSV requires csv/csv_file, csv_path, or coupled_models_csv"))
        swarm_kw, fit_kw = _split_swarm_fit_kwargs(Dict{Symbol,Any}(kwargs))
        merged = merge(fit_coupled_default_spec(), fit_kw)
        merged[:datapath] = datapath
        merged[:resultfolder] = folder
        merged[:maxtime] = maxtime
        swarm_kw[:filedir] = filedir
        return _makeswarmfiles_coupled_models_csv(
            path_csv;
            filedir=filedir,
            key_col=key_col,
            skip_empty=skip_empty,
            hierarchical_modes=hierarchical_modes,
            coupled_csv_cols=coupled_csv_cols,
            coupled_emit_legacy_r_variants=coupled_emit_legacy_r_variants,
            coupled_tie_rsum=coupled_tie_rsum,
            G=coupled_G,
            R=coupled_R,
            S=coupled_S,
            insertstep=coupled_insertstep,
            transitions=coupled_transitions,
            merged=merged,
            swarm_kw=swarm_kw,
        )
    end

    single_unit = models !== nothing ||
        (Gvec !== nothing && Rvec !== nothing && Svec !== nothing && insertvec !== nothing)
    if single_unit
        coupled === true &&
            throw(ArgumentError("makeswarmfiles: coupled=true is invalid for single-unit workflows (models= or Gvec/Rvec/Svec/insertvec); omit coupled or use coupled=false, or use csv_file/base_keys/h3 for coupled key lists"))
        !isempty(path_csv) &&
            throw(ArgumentError("makeswarmfiles: csv_file/csv_path cannot be used together with single-unit models= or Gvec/Rvec/Svec/insertvec"))
        h3_transition_range !== nothing &&
            throw(ArgumentError("makeswarmfiles: h3_transition_range cannot be used together with single-unit model sweep"))
        base_keys !== nothing &&
            throw(ArgumentError("makeswarmfiles: base_keys cannot be used together with single-unit model sweep"))
        return makeswarmfiles_driver(; transitions=transitions, transitionsvec=transitionsvec, Gvec=Gvec, Rvec=Rvec, Svec=Svec, insertvec=insertvec,
            combine=combine, models=models, datapath=datapath, resultfolder=folder, maxtime=maxtime,
            filedir=filedir, kwargs...)
    end

    path = path_csv
    resolved = if !isempty(path)
        h3_transition_range !== nothing &&
            throw(ArgumentError("makeswarmfiles: cannot use csv_file together with h3_transition_range"))
        !isnothing(base_keys) &&
            @warn "makeswarmfiles: csv_file is set; ignoring base_keys" csv_file=path
        df = DataFrame(CSV.File(path))
        if warn_coupled_csv_shape && ncol(df) >= 7
            @warn """Wide CSV ($(ncol(df)) columns): keys-only mode uses `key_col` ($(repr(key_col))) only. For per-row coupling, use `makeswarmfiles_coupled` or `makeswarmfiles(; coupled_csv=true, csv_file=..., ...)` (or legacy `coupled_models_csv`). Or `merge_existing_info=true` with existing `info_<key>.jld2`. Silence: `warn_coupled_csv_shape=false`.""" path = path
        end
        key_col in names(df) || throw(ArgumentError("missing key column '$key_col' in $path"))
        out = String[]
        for row in eachrow(df)
            k = strip(string(row[Symbol(key_col)]))
            (skip_empty && isempty(k)) && continue
            push!(out, k)
        end
        isempty(out) &&
            throw(ArgumentError("makeswarmfiles: CSV at '$path' produced no model keys (empty column or all rows skipped)."))
        out
    elseif h3_transition_range !== nothing
        !isnothing(base_keys) &&
            throw(ArgumentError("makeswarmfiles: use either h3_transition_range or base_keys, not both"))
        [string("H3-", t1, "-", t2) for t1 in h3_transition_range for t2 in h3_transition_range]
    elseif base_keys !== nothing
        out = String[]
        for k in base_keys
            k2 = strip(string(k))
            (skip_empty && isempty(k2)) && continue
            push!(out, k2)
        end
        isempty(out) &&
            throw(ArgumentError("makeswarmfiles: base_keys produced no non-empty keys (check skip_empty and entries)."))
        out
    else
        throw(ArgumentError("makeswarmfiles: specify one source — csv/csv_file / csv_path, base_keys, h3_transition_range, or single-unit models / Gvec+Rvec+Svec+insertvec"))
    end

    keys = String[]
    for h in hierarchical_modes
        append!(keys, [combined_rates_key(k; hierarchical=h) for k in resolved])
    end
    swarm_kw, fit_kw = _split_swarm_fit_kwargs(Dict{Symbol,Any}(kwargs))
    is_coupled = coupled === nothing ? true : coupled
    fit_base = merge((is_coupled ? fit_coupled_default_spec : fit_default_spec)(), fit_kw)
    fit_base[:datapath] = datapath
    fit_base[:resultfolder] = folder
    fit_base[:maxtime] = maxtime
    swarm_kw[:filedir] = filedir
    apply_h3 = h3_latent || (h3_transition_range !== nothing)
    _driver_write_specs_and_makeswarm(keys, fit_base, swarm_kw; merge_existing_info=merge_existing_info, fit_kw=fit_kw, warn_missing_info=warn_missing_info,
        apply_h3_latent_overlay=apply_h3, h3_coupling_ranges=h3_coupling_ranges)
end

function makeswarmfiles(csv_path::AbstractString, filedir::AbstractString; kwargs...)
    throw(ArgumentError(
        "positional `makeswarmfiles(csv_path, filedir; kw...)` was removed. Use e.g. " *
        "`makeswarmfiles(; csv_file=" * repr(csv_path) * ", filedir=" * repr(filedir) * ", kw...)` (keys-only CSV), or " *
        "`makeswarmfiles_coupled(; csv=" * repr(csv_path) * ", filedir=" * repr(filedir) * ", kw...)` / " *
        "`makeswarmfiles(; csv=" * repr(csv_path) * ", coupled_csv=true, filedir=" * repr(filedir) * ", kw...)` for Coupled_models_to_test (`csv` and `csv_file` are aliases).",
    ))
end

"""
    makeswarmfiles_coupled(; csv | csv_file, filedir, kwargs...)

Convenience entry for **Coupled_models_to_test**-style CSVs: calls [`makeswarmfiles`](@ref) with
**`coupled_csv=true`**, merging [`fit_coupled_default_spec`](@ref) with **`kwargs`**.

Pass the CSV path as **`csv`** or **`csv_file`** (aliases; at least one required, same semantics as [`makeswarmfiles`](@ref)).
"""
function makeswarmfiles_coupled(; csv=nothing, csv_file::AbstractString="", filedir::AbstractString, kwargs...)
    p = _resolve_csv_kw(csv_file, csv)
    isempty(p) && throw(ArgumentError("makeswarmfiles_coupled: pass csv=... or csv_file=... (path to Coupled_models_to_test-style CSV)"))
    makeswarmfiles(; csv_file=p, coupled_csv=true, filedir=filedir, kwargs...)
end

"""
    makeswarmfiles_h3_latent(filedir::String; transition_range=1:5, datapath, folder="coupled", maxtime=43000.0, hierarchical_modes=(true,false), merge_existing_info=true, warn_missing_info=true, kwargs...)

Convenience wrapper: calls [`makeswarmfiles`](@ref) with `h3_transition_range=transition_range`.
**`datapath`** is required (no default)—set your machine’s trace data root. Other defaults
match [`fit_coupled_default_spec`](@ref)-style coupled trace jobs: **`maxtime`** in **seconds**
(~`43000` ≈ 12 h wall for the MCMC engine, same units as [`fit`](@ref) → `MHOptions`),
and **`folder`** (`resultfolder`) aligned with [`makeswarmfiles`](@ref) (`"coupled"`). Override
**`folder`** if you use a separate results directory (e.g. `"coupled-h3"`).
Forwarded kwargs include **`warn_missing_info`** (see [`makeswarmfiles`](@ref)).
"""
function makeswarmfiles_h3_latent(filedir::AbstractString;
        transition_range=1:5,
        datapath::AbstractString,
        folder="coupled",
        maxtime=43000.0,
        hierarchical_modes=(true, false),
        kwargs...,
    )
    makeswarmfiles(; filedir=filedir, h3_transition_range=transition_range, datapath=datapath,
        folder=folder, maxtime=maxtime, hierarchical_modes=hierarchical_modes, kwargs...)
end

"""
    makeswarm_genes(genes::Vector{String}; <keyword arguments>)

Write a **swarm file** and **one shared fit script** so each gene runs as a separate job (genome-scale
scRNA-style: one model specification, many genes). Each swarm line runs:

`julia -t nthreads -p nchains fitscript_<label>_<model>.jl <gene>`

The generated script passes the gene name as `ARGS[1]` into `fit(...)`.

# Arguments
- `genes`: vector of gene names to include (list explicitly; there is no “all genes in folder” mode here).
- `batchsize=1000`: split into multiple `.swarm` files when `length(genes) > batchsize`.
- `filedir="."`: directory for `*.swarm` and the `fitscript_*.jl` file.
- **`method`**: like other Biowulf helpers, this must be a **`String`** in the emitted script (e.g.
  `\"Tsit5()\"`), not a raw `Tsit5()` value, so the script parses when submitted.
- Plus the same keywords as [`fit`](@ref) for the model and data (`datatype`, `datapath`, `cell`,
  `datacond`, `transitions`, `G`, `R`, `S`, `insertstep`, …), inference options forwarded as trailing
  `fit(...; ...)` kwargs (e.g. `inference_method=:nuts`, `parallel=:distributed`, `gradient=:finite`), and
  swarm-only options (`nchains`, `nthreads`, `swarmfile`, `juliafile`, `project`, `sysimage`, `src`).

# Cluster use

Submit the generated swarm with Biowulf `swarm` (or run the swarm file under Bash for sequential local
execution). Allocate wall time and CPUs so `maxtime` and `nchains` are feasible.

# Example
```julia
makeswarm_genes(["MYC", "SOX9"]; cell="HBEC", datatype="rna", datapath="data/", resultfolder="out")
```
"""
function makeswarm_genes(genes::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="rna", datapath="", cell::String="HCT116", datacond="MOCK", resultfolder::String="HCT116_test", label::String="",
    fittedparam::Vector=Int[], fixedeffects=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, coupling=tuple(), grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime=60.0, samplesteps::Int=1000000, warmupsteps=0, temp=1.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian=true, datacol=3, ejectnumber=1, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[], filedir=".", project="", sysimage="", kwargs...)
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    modelstring = create_modelstring(G, R, S, insertstep)
    label = create_label(label, datatype, datacond, cell)
    label_safe = sanitize_for_filename(label)
    model_safe = sanitize_for_filename(modelstring)
    ngenes = length(genes)
    juliafile_full = juliafile * "_" * label_safe * "_" * model_safe * ".jl"
    if ngenes > batchsize
        batches = getbatches(genes, ngenes, batchsize)
        for batch in eachindex(batches)
            sfile = (endswith(swarmfile, ".swarm") ? swarmfile[1:end-6] : swarmfile) * "_" * label_safe * "_" * model_safe * "_" * "$batch" * ".swarm"
            write_swarmfile(joinpath(filedir, sfile), nchains, nthreads, juliafile_full, batches[batch], project, sysimage)
        end
    else
        sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
        write_swarmfile(joinpath(filedir, sfile), nchains, nthreads, juliafile_full, genes, project, sysimage)
    end
    write_fitfile_genes(joinpath(filedir, juliafile_full), nchains, datatype, datapath, cell, datacond, resultfolder, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs; kwargs...)
    end

"""
    write_swarmfile(sfile, nchains, nthreads, juliafile::String, project="")

Writes a swarmfile for a single Julia file.

# Arguments
- `sfile`: The path to the swarmfile to be written.
- `nchains`: The number of chains to run.
- `nthreads`: The number of threads to use.
- `juliafile::String`: The path to the Julia file to be executed.
- `project`: The Julia project to use (optional, default is an empty string).

# Description
This function writes a swarmfile that specifies how to run a single Julia file with the given number of chains and threads. If a project is specified, it includes the `--project` flag in the command.
"""
function write_swarmfile(sfile, nchains, nthreads, juliafile::String, project="", sysimage = "")
    f = open(sfile, "w")
    if isempty(project)
        writedlm(f, ["julia -t $nthreads -p" nchains juliafile])
    elseif isempty(sysimage)
        writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile])
    else
        writedlm(f, ["julia --project=$project --sysimage=$sysimage -t $nthreads -p" nchains juliafile])
    end
    close(f)
end

"""
    write_swarmfile(sfile, nchains, nthreads, juliafile::String, genes::Vector{String}, project="", sysimage="")

Writes a swarm file with one line per gene. Each line runs: `julia ... juliafile gene` (gene passed as argument).
"""
function write_swarmfile(sfile, nchains, nthreads, juliafile::String, genes::Vector{String}, project="", sysimage="")
    f = open(sfile, "w")
    for gene in genes
        gene_safe = check_genename(gene, "(")
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains juliafile gene_safe])
        elseif isempty(sysimage)
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile gene_safe])
        else
            writedlm(f, ["julia --project=$project --sysimage=$sysimage -t $nthreads -p" nchains juliafile gene_safe])
        end
    end
    close(f)
end

"""
    write_swarmfile_keys(sfile, nchains, nthreads, juliafile_base, keys::Vector{String}, project="", sysimage="")

Writes a swarm file with one line per key. Each line runs: `julia -t nthreads -p nchains juliafile_base_<key>.jl`.
"""
function write_swarmfile_keys(sfile, nchains, nthreads, juliafile_base::String, keys::Vector{String}, project="", sysimage="")
    f = open(sfile, "w")
    base = isempty(juliafile_base) ? "fitscript" : juliafile_base
    for k in keys
        scriptname = base * "_" * sanitize_for_filename(k) * ".jl"
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains scriptname])
        elseif isempty(sysimage)
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains scriptname])
        else
            writedlm(f, ["julia --project=$project --sysimage=$sysimage -t $nthreads -p" nchains scriptname])
        end
    end
    close(f)
end

function _format_fit_override(k::Symbol, v)::Union{String,Nothing}
    v === nothing && return nothing
    v isa AbstractString && return "$k=\"$(escape_string(String(v)))\""
    v isa Real && return "$k=$v"
    v isa Bool && return "$k=$v"
    v isa Symbol && return "$k=$(repr(v))"
    return nothing
end

"""Build a `; kw1=..., kw2=...` suffix for emitted `fit(...)` scripts from keyword arguments."""
function _fit_positional_script_kw_suffix(; kwargs...)
    parts = String[]
    for (k, v) in pairs(kwargs)
        sym = k isa Symbol ? k : Symbol(k)
        line = _format_fit_override(sym, v)
        line === nothing || push!(parts, line)
    end
    return isempty(parts) ? "" : ("; " * join(parts, ", "))
end

"""
    write_fitscript_tracejoint_key(scriptpath, key, spec; src="")

Write a standalone fit script: `@everywhere using StochasticGene` then a multi-line
`fit(; key=..., kw=..., ...)` with one keyword per line (same style as `makescriptcoupled.jl`’s
`_write_fitscript_key`). Used for coupled/key-list batches and for **single-unit** [`makeswarm_models`](@ref)
so each `fitscript_<key>.jl` lists the full merged run spec, not only `key` plus a few overrides.

Values use `repr` (valid Julia literals) except:

- `method::AbstractString` is written verbatim (code snippet, e.g. `\"Tsit5()\"` or `\"(Tsit5(), true)\"`);
- `method` as a `Tsit5()` instance (or `(Tsit5(), bool)` for hierarchical) is shortened to `Tsit5()` / `(Tsit5(), true)`;
- `probfn` equal to `(prob_Gaussian, prob_Gaussian)` is written as module-qualified
  `(StochasticGene.prob_Gaussian, StochasticGene.prob_Gaussian)`.

`spec` should match the dict passed to [`write_run_spec_preset`](@ref). **`root`** is written exactly as
in `spec` (e.g. `"."`); no `abspath` / `expanduser` — path resolution is the user’s responsibility.
Keywords with value `nothing` are omitted.

This is used by [`_makeswarm_with_run_specs`](@ref) (coupled CSV and key-list [`makeswarmfiles`](@ref) paths).
"""
function write_fitscript_tracejoint_key(scriptpath::AbstractString, key::String, spec::Dict{Symbol,Any}; src::AbstractString="")
    sp = copy(spec)
    f = open(scriptpath, "w")
    write_prolog(f, src)
    write(f, "@time fit(; key=$(repr(key))")
    ks = sort!(collect(setdiff(keys(sp), (:key, :infolder, :traceinfo, :dttype))), by = x -> String(x))
    for sym in ks
        line = _fit_kw_value_for_fitscript_line(sym, sp[sym])
        line === nothing && continue
        write(f, ",\n    $(sym)=$(line)")
    end
    write(f, ")\n")
    close(f)
end

function _method_value_string_for_fitscript(v)::String
    v isa AbstractString && return String(v)
    if v isa Tuple && length(v) == 2 && v[2] isa Bool
        try
            if nameof(typeof(v[1])) === :Tsit5
                return v[2] ? "(Tsit5(), true)" : "(Tsit5(), false)"
            end
        catch
        end
    end
    try
        if nameof(typeof(v)) === :Tsit5
            return "Tsit5()"
        end
    catch
    end
    return repr(v)
end

function _probfn_value_string_for_fitscript(v)::String
    if v isa Tuple && length(v) == 2
        g1, g2 = v[1], v[2]
        if g1 isa Function && g2 isa Function && g1 === prob_Gaussian && g2 === prob_Gaussian
            return "(StochasticGene.prob_Gaussian, StochasticGene.prob_Gaussian)"
        end
    end
    return repr(v)
end

function _fit_kw_value_for_fitscript_line(k::Symbol, v)::Union{String,Nothing}
    v === nothing && return nothing
    if k === :method
        return _method_value_string_for_fitscript(v)
    end
    if k === :probfn
        return _probfn_value_string_for_fitscript(v)
    end
    return repr(v)
end

"""
`nchains` for the Julia `-p` flag: explicit `swarm_kw[:nchains]` wins; otherwise use the first
`nchains` found in `specs_by_key` (same order as `keys`) so the swarm matches
[`fit_coupled_default_spec`](@ref) / merged run specs (e.g. 16) without requiring a duplicate
`nchains=` in swarm-only kwargs.
"""
function _nchains_for_swarm_line(swarm_kw::Dict{Symbol,Any}, keys::Vector{String}, specs_by_key::Dict{String, Dict{Symbol,Any}})
    if haskey(swarm_kw, :nchains)
        return Int(swarm_kw[:nchains])
    end
    for k in keys
        sp = get(specs_by_key, k, nothing)
        sp === nothing && continue
        if haskey(sp, :nchains) && sp[:nchains] !== nothing
            return Int(sp[:nchains])
        end
    end
    return 2
end

function _makeswarm_with_run_specs(keys::Vector{String}, swarm_kw::Dict{Symbol,Any}, specs_by_key::Dict{String, Dict{Symbol,Any}})
    filedir = get(swarm_kw, :filedir, ".")
    !isempty(filedir) && !isdir(filedir) && mkpath(filedir)
    isempty(keys) && return
    nchains = _nchains_for_swarm_line(swarm_kw, keys, specs_by_key)
    nthreads = Int(get(swarm_kw, :nthreads, 1))
    swarmfile = get(swarm_kw, :swarmfile, "fit")
    juliafile = get(swarm_kw, :juliafile, "fitscript")
    project = get(swarm_kw, :project, "")
    sysimage = get(swarm_kw, :sysimage, "")
    src = get(swarm_kw, :src, "")
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    write_swarmfile_keys(joinpath(filedir, sfile), nchains, nthreads, juliafile, keys, project, sysimage)
    for k in keys
        haskey(specs_by_key, k) || throw(ArgumentError("missing run spec dict for key $(repr(k))"))
        scriptpath = joinpath(filedir, juliafile * "_" * sanitize_for_filename(k) * ".jl")
        write_fitscript_tracejoint_key(scriptpath, k, specs_by_key[k]; src=src)
    end
end

"""
    write_fitfile_key(fitfile, key::String; src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, ...)

Writes a script that runs `fit(; key=key, ...)` with optional overrides. Writable overrides: `String`, `Real`, `Bool`, and `Symbol`.
"""
function write_fitfile_key(fitfile, key::String; src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, warmupsteps=nothing, kwargs...)
    f = open(fitfile, "w")
    write_prolog(f, src)
    overrides = Dict{Symbol, String}()
    if resultfolder !== ""
        s = _format_fit_override(:resultfolder, resultfolder)
        s !== nothing && (overrides[:resultfolder] = s)
    end
    if root !== ""
        s = _format_fit_override(:root, root)
        s !== nothing && (overrides[:root] = s)
    end
    maxtime !== nothing && (s = _format_fit_override(:maxtime, maxtime); s !== nothing && (overrides[:maxtime] = s))
    samplesteps !== nothing && (s = _format_fit_override(:samplesteps, samplesteps); s !== nothing && (overrides[:samplesteps] = s))
    warmupsteps !== nothing && (s = _format_fit_override(:warmupsteps, warmupsteps); s !== nothing && (overrides[:warmupsteps] = s))
    for (k, v) in pairs(kwargs)
        s = _format_fit_override(Symbol(k), v)
        s !== nothing && (overrides[Symbol(k)] = s)
    end
    fit_args = ["key=$(repr(key))", values(overrides)...]
    write(f, "@time fit(; $(join(fit_args, ", ")))\n")
    close(f)
end

"""
    write_fitfile_genes(fitfile, nchains, datatype, datapath, cell, datacond, resultfolder, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)

Writes a fit script that takes the gene as ARGS[1] and calls `fit(nchains, datatype, ..., ARGS[1], cell, ...)`.
"""
function write_fitfile_genes(fitfile, nchains, datatype, datapath, cell, datacond, resultfolder, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor=1.0, trace_specs=[], dwell_specs=[]; kwargs...)
    s = '"'
    transitions = transitions isa AbstractVector && !(transitions isa Tuple) ? Tuple(transitions) : transitions
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $datapath, ARGS[1], $s$cell$s, $datacond, $s$resultfolder$s, $s$label$s, $fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $temp, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber, $yieldfactor, $trace_specs, $dwell_specs)" * _fit_positional_script_kw_suffix(; kwargs...) * "\n")
    close(f)
end

"""
    write_fitfile_coupled(fitfile, gene::String, nchains, ...)

Like `write_fitfile_genes` but with `gene` hardcoded in the script instead of `ARGS[1]`.
Used for coupled-model scripts where one gene is fixed per script.
`trace_specs` and `dwell_specs` are interpolated directly into the emitted `fit(...)` call.
"""
function write_fitfile_coupled(fitfile, gene::String, nchains, datatype, datapath, cell, datacond, resultfolder, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor=1.0, trace_specs=[], dwell_specs=[]; kwargs...)
    s = '"'
    transitions = transitions isa AbstractVector && !(transitions isa Tuple) ? Tuple(transitions) : transitions
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $datapath, $s$gene$s, $s$cell$s, $datacond, $s$resultfolder$s, $s$label$s, $fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $temp, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber, $yieldfactor, $trace_specs, $dwell_specs)" * _fit_positional_script_kw_suffix(; kwargs...) * "\n")
    close(f)
end

"""
    makeswarm_coupled(; gene, label, nchains, ..., trace_specs=[], dwell_specs=[], ...)

Write one fit script for a coupled model (gene hardcoded, not ARGS[1]) and append one
swarm line. Intended for coupled models where each label/coupling gets its own script.
`trace_specs` is written directly into the script for interval, per-trace window, and per-unit zero-centering.
"""
function makeswarm_coupled(; gene::String, label::String,
    nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript",
    filedir=".", project="", sysimage="", src="",
    datatype::String="tracejoint", datapath="", cell::String="HBEC",
    datacond=[], resultfolder::String="",
    fittedparam=Int[], fixedeffects=tuple(), transitions=tuple(), G=2, R=0, S=0, insertstep=1,
    coupling=tuple(), grid=nothing, root=".", maxtime=60.0, elongationtime=6.0,
    priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0,
    splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(),
    ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0,
    temp=1.0, temprna=1.0, burst=false,
    optimize=false, writesamples=false, method="Tsit5()", zeromedian=true,
    datacol=3, ejectnumber=1, yieldfactor::Float64=1.0,
    trace_specs=[], dwell_specs=[], kwargs...)
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    script_key = sanitize_for_filename(label)
    scriptfile = juliafile * "_" * script_key * ".jl"
    scriptpath = joinpath(filedir, scriptfile)
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    open(joinpath(filedir, sfile), "a") do f
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains scriptfile])
        elseif isempty(sysimage)
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains scriptfile])
        else
            writedlm(f, ["julia --project=$project --sysimage=$sysimage -t $nthreads -p" nchains scriptfile])
        end
    end
    write_fitfile_coupled(scriptpath, gene, nchains, datatype, datapath, cell,
        datacond, resultfolder, label, fittedparam, fixedeffects,
        transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime,
        priormean, nalleles, priorcv, onstates, decayrate, splicetype, probfn, noisepriors,
        hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp,
        temprna, burst, optimize, writesamples, method, src, zeromedian,
        datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs; kwargs...)
end

"""
    write_prolog(f, src)

Writes the prolog for the fitfile.

# Arguments
- `f`: The file handle to write to.
- `src`: The source directory.

# Description
This function writes the prolog for the fitfile. If `src` is empty, it assumes that `StochasticGene` is installed as a Julia package. Otherwise, it activates the specified project and uses `StochasticGene`.
"""
function write_prolog(f, src)
    s = '"'
    if isempty(src)
        write(f, "@everywhere using StochasticGene\n")
    else
        # write(f, "@everywhere using Pkg\n")
        # write(f, "@everywhere Pkg.activate($s$src$s)\n")
        write(f, "@everywhere using StochasticGene\n")
    end

end

"""
    create_label(label, datatype, datacond, cell)

Return the run label string used for outputs and for looking up prior runs under `resultfolder`.

# Arguments
- `label`: When non-empty, returned unchanged (after `datacond` is normalized for the vector method).
- `datatype`, `datacond`, `cell`: Used only when `label` is empty to build `datatype * "-" * cell * "_" * datacond`.

# Methods
- `create_label(label, datatype, datacond::String, cell)`
- `create_label(label, datatype, datacond::Vector, cell)`: joins `datacond` with `"-"`.

When `label` is empty, the built label is `datatype * "-" * cell * "_" * datacond`. Pass an explicit
`label` when filenames need extra disambiguators.
"""
function create_label(label, datatype, datacond::String, cell)
    if isempty(label)
        label = datatype * "-" * cell
        typeof(datacond) <: AbstractString && (label = label * "_" * datacond)
    end
    return label
end

function create_label(label, datatype, datacond::Vector, cell)
    create_label(label, datatype, join(datacond, "-"), cell)
end

"""
    create_modelstring(G, R, S, insertstep)

Creates a model string based on the given parameters.

# Arguments
- `G`: Parameter G.
- `R`: Parameter R.
- `S`: Parameter S.
- `insertstep`: The insert step.

# Returns
- A string representing the model.

# Description
This function generates a model string based on the provided parameters. If `R` is greater than 0, it adjusts `S` based on `R` and `insertstep` and returns a concatenated string of `G`, `R`, `S`, and `insertstep`. If `R` is not greater than 0, it returns `G` as the model string.
"""
function create_modelstring(G::Int, R, S, insertstep)
    if R > 0
        return "$G$R$S$insertstep"
    else
        return "$G"
    end
end

function create_modelstring(G::Tuple, R, S, insertstep)
    m = ""
    for i in eachindex(G)
        m *= "$(G[i])$(R[i])$(S[i])$(insertstep[i])"
    end
    return m
end

"""
    fix(folder)

Finds jobs that failed and writes a swarmfile for those genes.

# Arguments
- `folder`: The folder to search for failed jobs.

# Description
This function finds jobs that failed in the specified folder and writes a swarmfile for those genes.
"""
fix(folder) = writeruns(fixruns(findjobs(folder)))

function fix_filenames(folder, old="scRNA-ss-", new="scRNA-ss_")
    files = readdir(folder)
    for file in files
        if occursin(old, file)
            nfile = replace(file, old => new)
            mv(joinpath(folder, file), joinpath(folder, nfile), force=true)
        end
    end
end

"""
    rna_setup(root=".")

Sets up the folder system prior to use. Defaults to the current directory.

# Arguments
- `root`: The root directory for setting up the folder system. Defaults to the current directory.

# Description
This function sets up the necessary folder structure for RNA data processing, including creating directories for alleles, halflives, and test data if they do not already exist.
"""
function rna_setup(root=".")
    folder_setup(root)
    alleles = joinpath("data", "alleles")
    halflives = joinpath("data", "halflives")
    testdata = joinpath("data", "HCT116_testdata")

    if ~ispath(alleles)
        mkpath(alleles)
    end
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/alleles/CH12_alleles.csv", "$alleles/CH12_alleles.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/alleles/HCT116_alleles.csv", "$alleles/HCT116_alleles.csv")
    if ~ispath(halflives)
        mkpath(halflives)
    end
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/ESC_halflife.csv", "$halflives/ESC_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/CH12_halflife.csv", "$halflives/CH12_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/HCT116_halflife.csv", "$halflives/HCT116_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/OcaB_halflife.csv", "$halflives/OcaB_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/aB_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/CAST_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/FIBS_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/MAST_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/NK_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/TEC_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv", "$halflives/SKIN_halflife.csv")
    if ~ispath(testdata)
        mkpath(testdata)
    end
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/HCT116_testdata/CENPL_MOCK.txt", "$testdata/CENPL_MOCK.txt")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/HCT116_testdata/MYC_MOCK.txt", "$testdata/MYC_MOCK.txt")
end

"""
    folder_setup(root=".")

Creates data and results folders (if they do not already exist).

# Arguments
- `root`: The root directory for setting up the folder system. Defaults to the current directory.

# Description
This function creates the necessary data and results folders within the specified root directory. If the folders already exist, it does nothing.
"""
function folder_setup(root=".")
    data = joinpath(root, "data")
    results = joinpath(root, "results")
    if ~ispath(data)
        mkpath(data)
    end
    if ~ispath(results)
        mkpath(results)
    end
end


"""
    getbatches(genes, ngenes, batchsize)

Divides the list of genes into batches of a specified size.

# Arguments
- `genes`: A vector of gene names.
- `ngenes`: The total number of genes.
- `batchsize`: The size of each batch.

# Returns
- A vector of vectors, where each inner vector contains a batch of gene names.

# Description
This function divides the list of genes into batches of a specified size. The last batch may contain fewer genes if the total number of genes is not a multiple of the batch size.
"""
function getbatches(genes, ngenes, batchsize)
    nbatches = div(ngenes, batchsize)
    batches = Vector{Vector{String}}(undef, nbatches + 1)
    println(batchsize, " ", nbatches + 1)
    for i in 1:nbatches
        batches[i] = genes[batchsize*(i-1)+1:batchsize*(i)]
    end
    batches[end] = genes[batchsize*nbatches+1:end]
    return batches
end

"""
    checkgenes(root, conds, datapath, celltype::String, thresholdlow::Float64, thresholdhigh::Float64)

Checks and returns the genes that meet the specified thresholds across multiple conditions and data paths.

# Arguments
- `root`: The root directory.
- `conds`: A condition or a list of conditions to check.
- `datapath`: A data path or a list of data paths to check.
- `celltype`: The type of cell.
- `thresholdlow`: The lower threshold for filtering genes.
- `thresholdhigh`: The upper threshold for filtering genes.

# Returns
- `Vector{String}`: A vector of genes that meet the specified thresholds across the given conditions and data paths.

"""
function checkgenes(root, conds::Vector, datapath, celltype::String, thresholdlow::Float64, thresholdhigh::Float64)
    genes = Vector{Vector{String}}(undef, 0)
    typeof(conds) <: AbstractString && (conds = [conds])
    typeof(datapath) <: AbstractString && (datapath = [datapath])
    for d in datapath, c in conds
        push!(genes, checkgenes(root, c, d, celltype, thresholdlow, thresholdhigh))
    end
    geneset = genes[1]
    for g in genes
        genesest = convert(Vector{String}, intersect(geneset, g))
    end
    return genesest
end

function checkgenes(root, cond::String, datapath::String, cell::String, thresholdlow::Float64, thresholdhigh::Float64)
    if cell == "HBEC"
        return genes_hbec()
    elseif cell == ""
        return get_genes(cond, datapath)
    else
        datapath = folder_path(datapath, root, "data")
        genes = intersect(get_halflives(root, cell, thresholdlow, thresholdhigh), get_genes(cond, datapath))
        alleles = get_alleles(root, cell)
        if ~isnothing(alleles)
            return convert(Vector{String}, intersect(genes, alleles))
        else
            return convert(Vector{String}, genes)
        end
    end
end

"""
    folder_path(folder::String, root::String, folderatetype::String=""; make=false)

Returns the full path for a given folder, optionally creating the path if it does not exist.

# Arguments
- `folder`: A folder name.
- `root`: The root directory.
- `folderatetype`: Optional subfolder type (e.g. `"results"`); used to form `joinpath(root, folderatetype, folder)` when needed (default `""`).
- `make`: If `true`, create the path if it does not exist (default `false`).

# Returns
- `String`: The full path for the given folder.

If `folder` already begins with `folderatetype` as its first path component (e.g. `resultfolder="results/run-1"` with
`folderatetype="results"`), the resolved path is `joinpath(root, folder)` only — never `root/results/results/...`.
"""
function folder_path(folder::String, root::String, folderatetype::String=""; make=false)
    f = folder
    if ~ispath(folder) && ~isempty(folder)
        f = joinpath(root, folder)
        if ~ispath(f)
            # Only prepend folderatetype if folder doesn't already start with it (prevents results/results nesting)
            first_component = splitpath(folder)[1]
            if isempty(folderatetype) || first_component != folderatetype
                f = joinpath(root, folderatetype, folder)
            end
            if ~ispath(f) && ~make
                @warn "folder_path: no existing directory matched" folder=folder root=root folderatetype=folderatetype tried=f
            elseif ~ispath(f) && make
                mkpath(f)
            end
        end
    end
    f
end

"""
    folder_path(folder::Vector, root, foldertype)

Returns the full paths for a list of folders, optionally creating the paths if they do not exist.

# Arguments
- `folders`: A list of folder names.
- `root`: The root directory.
- `foldertype`: The type of folder (optional, default is an empty string).
- `make`: A boolean flag indicating whether to create the paths if they do not exist (optional, default is `false`).

# Returns
- `Vector{String}`: The full paths for the given folders.
"""
function folder_path(folder::Vector, root, foldertype)
    fv = folder
    for i in eachindex(fv)
        fv[i] = folder_path(fv[i], root, foldertype)
    end
    fv
end

"""
    get_halflives(root, cell, thresholdlow::Float64, thresholdhigh::Float64)

Retrieves the genes with halflives within the specified thresholds for a given cell type from the root directory.

# Arguments
- `root`: The root directory.
- `cell`: The type of cell.
- `thresholdlow`: The lower threshold for filtering halflives.
- `thresholdhigh`: The upper threshold for filtering halflives.

# Returns
- `Vector{String}`: A vector of genes that meet the specified halflife thresholds.

"""
function get_halflives(root, cell, thresholdlow::Float64, thresholdhigh::Float64)
    path = get_file(root, "data/halflives", cell, "csv")
    if isnothing(path)
        return nothing
    else
        get_halflives(path, thresholdlow, thresholdhigh)
    end
end

"""
    get_halflives(hlpath, thresholdlow::Float64, thresholdhigh::Float64)

Retrieves the genes with halflives within the specified thresholds from a given file path.

# Arguments
- `hlpath`: The file path to the halflives data.
- `thresholdlow`: The lower threshold for filtering halflives.
- `thresholdhigh`: The upper threshold for filtering halflives.

# Returns
- `Vector{String}`: A vector of genes that meet the specified halflife thresholds.
"""
function get_halflives(hlpath, thresholdlow::Float64, thresholdhigh::Float64)
    genes = Vector{String}(undef, 0)
    if ispath(hlpath)
        halflives = readdlm(hlpath, ',')
    else
        return nothing
    end
    for row in eachrow(halflives)
        if typeof(row[2]) <: Number
            if thresholdlow <= row[2] < thresholdhigh
                push!(genes, string(row[1]))
            end
        end
    end
    return genes
end


"""
    get_alleles(allelepath)
    get_alleles(root, cell)

Retrieves the alleles from a given file path.

# Arguments
- `root`: The root directory.
- `cell`: The type of cell.
- `allelepath`: The file path to the alleles data.

# Methods
- `get_alleles(root, cell)`: Retrieves the alleles for a given cell type from the root directory.
- `get_alleles(allelepath)`: Retrieves the alleles from a given file path.

# Returns
- `Vector{String}`: A vector of alleles if the file path is not `nothing`.
- `Nothing`: If the file path is `nothing`.


"""
function get_alleles(allelepath)
    if ~isnothing(allelepath)
        return readdlm(allelepath, ',')[2:end, 1]
    else
        return nothing
    end
end

get_alleles(root, cell) = get_alleles(get_file(root, "data/alleles", cell, "csv"))



"""
    get_file(root folder, filetype, suffix)
    get_file(folder, filetype, suffix)

Retrieves the file path for a file of a specified type and suffix within a given folder.

# Arguments
- `root`: The root directory.
- `folder`: The folder to search within.
- `filetype`: The type of file to search for (the prefix of the file name).
- `suffix`: The suffix of the file name to search for.

# Returns
- `String`: The full path to the file if found.
- `Nothing`: If the folder does not exist or the file is not found.
"""
get_file(root, folder, filetype, suffix) = get_file(joinpath(root, folder), filetype, suffix)


function get_file(folder, filetype, suffix)
    if ispath(folder)
        files = readdir(folder)
        for file in files
            name = split(file, "_")[1]
            if occursin(suffix, file) && name == filetype
                path = joinpath(folder, file)
                return path
            end
        end
    else
        println("$folder does not exist")
        return nothing
    end
end

"""
    findjobs(folder)

Finds and lists unique job identifiers from swarm files in the specified folder.

# Arguments
- `folder`: The folder to search for swarm files.

# Returns
- `Vector{String}`: A vector of unique job identifiers.

# Description
This function searches the specified folder for files with names starting with "swarm_". It extracts and returns a list of unique job identifiers from these file names.
"""
function findjobs(folder)
    files = readdir(folder)
    files = files[occursin.("swarm_", files)]
    for (i, file) in enumerate(files)
        files[i] = split(file, "_")[2]
    end
    unique(files)
end

"""
    fixruns(jobs, message="FAILED")

Identifies and processes jobs that failed based on a specified message.

# Arguments
- `jobs`: A vector of job identifiers.
- `message`: The message to search for in job histories to identify failures. Defaults to "FAILED".

# Returns
- `Vector{String}`: A vector of run identifiers that need to be fixed.

# Description
This function processes a list of job identifiers to find those that contain a specified failure message in their job history. It reads the job history, identifies the lines containing the failure message, and extracts the corresponding run identifiers. It returns a list of these run identifiers.

"""
function fixruns(jobs, message="FAILED")
    runlist = Vector{String}(undef, 0)
    for job in jobs
        if occursin(message, read(`jobhist $job`, String))
            swarmfile = findswarm(job)
            list = readdlm(swarmfile, ',')
            runs = chomp(read(pipeline(`jobhist $job`, `grep $message`), String))
            runs = split(runs, '\n')
            println(job)
            for run in runs
                linenumber = parse(Int, split(split(run, " ")[1], "_")[2]) + 1
                while linenumber < length(list)
                    a = String(list[linenumber])
                    linenumber += 1000
                    println(a)
                    push!(runlist, a)
                end
            end
        end
    end
    return runlist
end

"""
   writeruns(runs, outfile="fitfix.swarm")

Writes a list of runs to a specified output file.

# Arguments
- `runs`: A vector of run identifiers.
- `outfile`: The name of the output file. Defaults to "fitfix.swarm".

# Description
This function writes each run identifier from the provided list to the specified output file. Each run is written on a new line. If the output file already exists, it will be overwritten.

"""
function writeruns(runs, outfile="fitfix.swarm")

    f = open(outfile, "w")
    for run in runs
        writedlm(f, [run], quotes=false)
    end
    close(f)

end

"""
    findswarm(job)

Finds the swarm file associated with a given job.

# Arguments
- `job`: The job identifier.

# Returns
- `String`: The name of the swarm file associated with the job.

# Description
This function retrieves the swarm file name associated with a given job by reading the job history and searching for the swarm command. It returns the name of the swarm file.

"""
function findswarm(job)
    sc = "Swarm Command"
    line = read(pipeline(`jobhist $job`, `grep $sc`), String)
    list = split(line, " ")
    list[occursin.(".swarm", list)][1]
end

"""
    get_missing_genes(datapath::String, resultfolder::String, cell, filetype, label, cond, model, root=".")
    get_missing_genes(genes::Vector, resultfolder, filetype, label, cond, model)

Finds genes that are missing from the results folder based on the given parameters.

# Arguments
- `datapath`: The path to the data directory.
- `resultfolder`: The path to the results directory.
- `cell`: The cell type.
- `filetype`: The type of file to search for.
- `label`: The label to use for the search.
- `cond`: The condition to use for the search.
- `model`: The model to use for the search.
- `root`: The root directory. Defaults to the current directory.
- `genes`: A vector of gene identifiers.

# Returns
- `Vector{String}`: A vector of missing gene identifiers.

# Methods
- `get_missing_genes(datapath::String, resultfolder::String, cell, filetype, label, cond, model, root=".")`: Finds missing genes based on the specified parameters.
- `get_missing_genes(genes::Vector, resultfolder, filetype, label, cond, model)`: Finds missing genes based on a list of gene identifiers.

# Description
This function checks for genes that are missing from the results folder based on the specified parameters. It first retrieves the list of genes from the data directory and then compares it with the genes present in the results folder. It returns a list of missing gene identifiers.

"""
function get_missing_genes(datapath::String, resultfolder::String, cell, filetype, label, cond, model, root=".")
    genes = checkgenes(root, cond, datapath, cell, 0.0, 1e8)
    get_missing_genes(genes, folder_path(resultfolder, root, "results"), filetype, label, cond, model)
end

function get_missing_genes(genes::Vector, resultfolder, filetype, label, cond, model)
    genes1 = get_genes(resultfolder, filetype, label, cond, model)
    get_missing_genes(genes, genes1)
end

get_missing_genes(genes, genes1) = union(setdiff(genes1, genes), setdiff(genes, genes1))

"""
    scan_swarmfiles(jobid, folder=".")

Scans swarm files in the specified folder for a given job ID.

# Arguments
- `jobid`: The job identifier.
- `folder`: The folder to search for swarm files. Defaults to the current directory.

# Returns
- `Vector{String}`: A vector of gene identifiers found in the swarm files.

# Description
This function scans the swarm files in the specified folder for the given job ID. It reads the contents of the swarm files and extracts the gene identifiers associated with the job ID. It returns a list of these gene identifiers.
"""
function scan_swarmfiles(jobid, folder=".")
    if ~(typeof(jobid) <: String)
        jobid = string(jobid)
    end
    genes = Array{String,1}(undef, 0)
    files = readdir(folder)
    files = files[occursin.(jobid, files)]
    for file in files
        genes = vcat(genes, scan_swarmfile(file))
    end
    return genes
end

"""
    scan_swarmfile(file)

Scans a single swarm file for gene identifiers.

# Arguments
- `file`: The path to the swarm file.

# Returns
- `Vector{String}`: A vector of gene identifiers found in the swarm file.

# Description
This function reads the contents of a given swarm file and extracts the gene identifiers. It returns a list of these gene identifiers.
"""
function scan_swarmfile(file)
    genes = Array{String,1}(undef, 0)
    contents = readdlm(file, '\t')
    lines = contents[occursin.("[\"", string.(contents))]
    for line in lines
        push!(genes, split.(string(line), " ")[1])
    end
    return genes
end

"""
    scan_fitfile(file, folder=".")

Scans a fit file for gene identifiers.

# Arguments
- `file`: The name of the fit file.
- `folder`: The folder where the fit file is located. Defaults to the current directory.

# Returns
- `Vector{String}`: A vector of gene identifiers found in the fit file.

# Description
This function reads the contents of a given fit file and extracts the gene identifiers. It returns a list of these gene identifiers.
"""
function scan_fitfile(file, folder=".")
    genes = Array{String,1}(undef, 0)
    joinpath(folder, file)
    file = readdlm(file, '\t')
    for line in eachrow(file)
        push!(genes, line[4])
    end
    return genes
end

"""
    new_FISHfolder(newroot::String, oldroot::String, rep::String)

Creates a new folder structure for FISH data based on an existing structure.

# Arguments
- `newroot`: The root directory for the new folder structure.
- `oldroot`: The root directory of the existing folder structure.
- `rep`: The string to search for in directory names.

# Description
This function walks through the existing folder structure and creates a new folder structure in the specified new root directory. It searches for directories containing the specified string and replicates their structure in the new root directory.
"""
function new_FISHfolder(newroot::String, oldroot::String, rep::String)
    for (root, dirs, files) in walkdir(oldroot)
        for dir in dirs
            if occursin(rep, dir)
                oldtarget = joinpath(root, dir)
                newtarget = replace(oldtarget, oldroot => newroot)
                println(newtarget)
                if !ispath(newtarget)
                    mkpath(newtarget)
                end
                cp(oldtarget, newtarget, force=true)
            end
        end
    end
end

"""
    find_genes_to_rerun(datafolder::String, resultfolder::String, filetype::String, nlines_expected::Int; error_keywords=["nonconverged", "NaN", "median"], verbose=true)

Systematically finds genes that are missing or nonconverged in the results folder, based on the data folder gene list.

- `datafolder`: Path to the folder containing gene data files (e.g. data/U3AS4/WT-UTR/)
- `resultfolder`: Path to the folder containing measures files (e.g. results/WT-UTR/)
- `filetype`: e.g. "rnacount"
- `nlines_expected`: Number of lines expected in a complete measures file
- `error_keywords`: List of strings indicating nonconvergence or fallback (default: ["nonconverged", "NaN", "median"])
- `verbose`: If true, prints a summary

Returns a vector of gene names to rerun.
"""
function find_genes_to_rerun(datafolder::String, resultfolder::String, filetype::String, nlines_expected::Int; error_keywords=["nonconverged", "NaN", "median"], verbose=true)
    # Get all gene names from the data folder (assume <gene>.txt or <gene>_*.txt)
    datafiles = readdir(datafolder)
    genes = String[]
    for f in datafiles
        m = match(r"^([A-Za-z0-9_\-]+)", f)
        if m !== nothing
            gene = m.captures[1]
            push!(genes, gene)
        end
    end
    genes = unique(genes)

    # For each gene, check for measures file and its content
    to_rerun = String[]
    for gene in genes
        # Find measures file(s) for this gene
        pattern = r"measures_" * filetype * "-" * gene * ".*\\.txt"
        mfiles = filter(f -> occursin(pattern, f), readdir(resultfolder))
        if isempty(mfiles)
            push!(to_rerun, gene)
            if verbose
                println("Missing measures file for gene: $gene")
            end
            continue
        end
        # Check content of the first measures file (or all, if you want)
        file = joinpath(resultfolder, mfiles[1])
        lines = readlines(file)
        if length(lines) < nlines_expected || any(kw -> any(occursin(kw, l) for l in lines), error_keywords)
            push!(to_rerun, gene)
            if verbose
                println("Nonconverged or incomplete measures file for gene: $gene")
            end
        end
    end
    if verbose
        println("Total genes in data: ", length(genes))
        println("Genes to rerun: ", length(to_rerun))
    end
    return to_rerun
end

"""
    find_missing_genes(datafolder::String, resultfolder::String, filetype::String)

Returns a vector of gene names in the data folder that are missing a corresponding measures file in the results folder.
"""
function find_missing_genes(datafolder::String, resultfolder::String, filetype::String)
    datafiles = readdir(datafolder)
    genes = String[]
    for f in datafiles
        m = match(Regex("^([A-Za-z0-9_\\-]+)"), f)
        if m !== nothing
            gene = m.captures[1]
            push!(genes, gene)
        end
    end
    genes = unique(genes)
    missing = String[]
    for gene in genes
        pattern = "measures_" * filetype * "-" * gene * ".*\\.txt"
        mfiles = filter(f -> occursin(Regex(pattern), f), readdir(resultfolder))
        if isempty(mfiles)
            push!(missing, gene)
        end
    end
    return missing
end

"""
    find_nonconverged_genes(resultfolder::String, filetype::String, nlines_expected::Int; rhat_thresh=1.1, ess_min=100, geweke_thresh=2.0, verbose=true)

Checks all measures files in the results folder for nonconvergence based on R-hat, ESS, and Geweke diagnostics.
Prints out which diagnostic(s) failed for each nonconverged gene.
Returns a vector of nonconverged gene names.

NOTE: You must adapt the parsing code to your actual measures file format!
"""
function find_nonconverged_genes(resultfolder::String, filetype::String, nlines_expected::Int; rhat_thresh=1.1, ess_min=100, geweke_thresh=2.0, verbose=true)
    files = readdir(resultfolder)
    measure_files = filter(f -> startswith(f, "measures_" * filetype * "-"), files)
    nonconverged = String[]
    for file in measure_files
        pattern = "measures_" * filetype * "-([A-Za-z0-9_\\-]+)_"
        gene_match = match(Regex(pattern), file)
        gene = gene_match !== nothing ? gene_match.captures[1] : nothing
        lines = readlines(joinpath(resultfolder, file))
        # Placeholder: parse diagnostics from lines
        # You must adapt this to your file format!
        rhat = Float64[]  # fill with parsed values
        ess = Float64[]   # fill with parsed values
        geweke = Float64[] # fill with parsed values
        # Example: if lines 2,3,4 are rhat, ess, geweke (comma-separated)
        # rhat = parse.(Float64, split(lines[2], ","))
        # ess = parse.(Float64, split(lines[3], ","))
        # geweke = parse.(Float64, split(lines[4], ","))
        failed = false
        if length(lines) < nlines_expected
            failed = true
            if verbose && gene !== nothing
                println("$gene: incomplete file (only $(length(lines)) lines)")
            end
        end
        if any(rhat .> rhat_thresh)
            failed = true
            if verbose && gene !== nothing
                println("$gene: R-hat > $rhat_thresh")
            end
        end
        if any(ess .< ess_min)
            failed = true
            if verbose && gene !== nothing
                println("$gene: ESS < $ess_min")
            end
        end
        if any(abs.(geweke) .> geweke_thresh)
            failed = true
            if verbose && gene !== nothing
                println("$gene: |Geweke| > $geweke_thresh")
            end
        end
        if failed && gene !== nothing
            push!(nonconverged, gene)
        end
    end
    return nonconverged
end

"""
    find_highly_nonconverged_genes(measures_file::String; 
        rhat_thresh=1.2, ess_min=500, geweke_thresh=3.0, mcse_max=0.5, 
        accept_min=0.05, temp_val=1.0, verbose=true)

Scan an assembled measures file and return a vector of gene names that are highly nonconverged by relaxed criteria.
"""
function find_highly_nonconverged_genes(measures_file::String;
    rhat_thresh=1.2, ess_min=500, geweke_thresh=3.0, mcse_max=0.5,
    accept_min=0.05, temp_val=1.0, verbose=true)

    df = CSV.read(measures_file, DataFrame)
    bad_genes = String[]
    for row in eachrow(df)
        gene = row.Gene
        waic = row.WAIC
        rhat = row.Rhat
        ess = row.ESS
        geweke = row.Geweke
        mcse = row.MCSE
        accept = row.Acceptance
        temp = row.Temperature

        is_bad = false
        if ismissing(rhat) || isnan(rhat) || rhat > rhat_thresh
            is_bad = true
            verbose && println("$gene: Rhat = $rhat")
        end
        if ismissing(ess) || isnan(ess) || ess < ess_min
            is_bad = true
            verbose && println("$gene: ESS = $ess")
        end
        if ismissing(geweke) || isnan(geweke) || abs(geweke) > geweke_thresh
            is_bad = true
            verbose && println("$gene: Geweke = $geweke")
        end
        if ismissing(mcse) || isnan(mcse) || mcse > mcse_max
            is_bad = true
            verbose && println("$gene: MCSE = $mcse")
        end
        if ismissing(waic) || isnan(waic)
            is_bad = true
            verbose && println("$gene: WAIC is NaN")
        end
        if ismissing(accept) || isnan(accept) || accept < accept_min
            is_bad = true
            verbose && println("$gene: Acceptance = $accept")
        end
        if ismissing(temp) || isnan(temp) || temp != temp_val
            is_bad = true
            verbose && println("$gene: Temperature = $temp")
        end
        if is_bad
            push!(bad_genes, gene)
        end
    end
    return bad_genes
end
