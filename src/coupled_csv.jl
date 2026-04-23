# This file is part of StochasticGene.jl
#
# coupled_csv.jl — Coupled_models_to_test.csv or Coupled_models_free.csv → coupling connections and fit-spec pieces.
#
# Column naming (order-independent, searched by name pattern):
#   Model_name — run key (spaces → -)
#   enhancer_to_gene_1, enhancer_sign_1 — enhancer→gene connections for unit 1
#   enhancer_to_gene_2, enhancer_sign_2 — enhancer→gene connections for unit 2
#   gene_to_enhancer, gene_to_enhancer_sign — gene→enhancer connections
#   background_gene, background_gene_sign — genetic background connections
# Missing columns default to empty string (parse as :free sign mode).
#
# Legacy fixed-column format still supported via coupled_csv_cols parameter.
# Token codes: two-digit "st", "Rsumk", "Ranyk", legacy "Rk" (see csv_row_to_connections_simple).

"""
    Coupled_models CSV schema

This module implements parsing from coupled model spreadsheets (e.g., **Coupled_models_to_test.csv**
or **Coupled_models_free.csv**) into coupling connections and keyword dicts for 
[`build_coupled_fit_spec_from_csv_cells`](@ref), used by [`makeswarmfiles`](@ref) when `coupled_csv=true`.

# Column naming (flexible order and subset)

Columns are found by **name pattern** (case-insensitive, order-independent):

| Purpose | Column name patterns | Content |
|---------|----------------------|---------|
| Model key | `Model_name` | Run key; spaces → `-`. **Required.** |
| Enhancer→gene (unit 1) | `enhancer_to_gene_1`, `e1` | State token string. |
| Sign (unit 1) | `enhancer_sign_1`, `e1s` | Sign (`">0"` = activate, empty/`0`/`"free"` = free, else inhibit). |
| Enhancer→gene (unit 2) | `enhancer_to_gene_2`, `gene_to_enhancer`, `e2` | State token string. |
| Sign (unit 2) | `enhancer_sign_2`, `gene_to_enhancer_sign`, `e2s` | Sign. |
| Gene→enhancer (unit 1) | `background_gene`, `ge` | State token string. |
| Sign | `background_gene_sign`, `ges` | Sign. |

**Missing columns** default to empty string (which parse as `:free` mode).

# Token strings in cells

Comma-separated tokens are parsed per [`csv_row_to_connections_simple`](@ref):

- Two-digit **`st`** — G-state source `s` and target transition index `t`.
- **`Rsumk`**, **`Ranyk`** — sum or “any R step” coupling for elongation step `k`.
- Legacy **`Rk`** — may be rewritten to `Rsumk` / `Ranyk` via [`replace_csv_cell_legacy_r`](@ref) when emitting variants.

Sign cells use [`parse_coupling_sign_csv`](@ref):
- `">0"` → `:activate` (positive coupling)  
- empty string / `"0"` / `"free"` → `:free` (unconstrained coupling)
- anything else → `:inhibit` (negative coupling)

# Legacy fixed-column format

Older CSVs with 7 fixed columns use the `coupled_csv_cols` parameter (default `(2, 3, 4, 5, 6, 7)`).
The new flexible schema is recommended for clarity.

# See also

[`csv_row_to_connections_simple`](@ref), [`parse_coupling_sign_csv`](@ref), 
[`build_coupled_fit_spec_from_csv_cells`](@ref), [`makeswarmfiles`](@ref) docstring section *Coupled CSV*.
"""
const COUPLED_CSV_SCHEMA = true

const _COUPLED_CSV_DEFAULT_ECOLS = (2, 3, 4, 5, 6, 7)

"""Map CSV sign string to `:activate` or `:inhibit` (`>0` / `<0`)."""
function parse_coupling_sign_csv(s)::Symbol
    s_clean = strip(lowercase(string(s)))
    if s_clean == ">0"
        return :activate
    elseif s_clean == "<0"
        return :inhibit
    else
        return :free
    end
end

"""Default γ placeholder per mode (matches legacy `makescriptcoupled.jl`)."""
function default_coupling_gamma_csv(mode::Symbol)
    mode === :activate && return 0.1
    mode === :inhibit && return -0.1
    mode === :free && return 0.0
    return 0.0
end

function _legacy_r_token(tok::AbstractString)
    s = strip(string(tok))
    isempty(s) && return false
    startswith(s, "Rsum") && return false
    startswith(s, "Rany") && return false
    startswith(s, "R") || return false
    tail = s[2:end]
    !isempty(tail) && all(isdigit, tail)
end

"""
Replace legacy `Rk` tokens in a CSV cell with `Rsumk` or `Ranyk` (comma-separated tokens preserved).
"""
function replace_csv_cell_legacy_r(cell, kind::Symbol)
    s = strip(string(cell))
    (isempty(s) || s == "missing") && return cell
    parts = split(s, ',')
    out = map(parts) do p
        t = strip(p)
        if _legacy_r_token(t)
            k = t[2:end]
            (kind === :rsum) ? "Rsum$k" : "Rany$k"
        else
            t
        end
    end
    return join(out, ",")
end

function csv_row_has_legacy_r(e1, e2, ge)
    for cell in (e1, e2, ge)
        s = strip(string(cell))
        (isempty(s) || s == "missing") && continue
        for p in split(s, ',')
            _legacy_r_token(strip(p)) && return true
        end
    end
    return false
end

"""
    csv_row_to_connections_simple(e1_states, e1_sign, e2_states, e2_sign, ge_states, ge_sign, G, R; tie_rsum=true)

Parse six CSV cells into connection tuples `(β, s, α, t)`, default γs, tie groups (shared coupling indices),
and per-connection modes. See module docstring [Coupled model CSV format](@ref).
"""
function csv_row_to_connections_simple(e1_states, e1_sign, e2_states, e2_sign, ge_states, ge_sign, G::Tuple, R::Tuple; tie_rsum::Bool=true)
    connections = NTuple{4,Int}[]
    gammas = Float64[]
    tie_groups = Vector{Vector{Int}}()
    modes = Symbol[]

    function _parse_code(code, β, α)
        s = strip(string(code))
        isempty(s) && return NTuple{4,Int}[]
        Gβ = G[β]
        Rβ = R[β]
        if startswith(s, "Rany")
            tail = s[5:end]
            (isempty(tail) || !all(isdigit, tail)) && return NTuple{4,Int}[]
            t = parse(Int, tail)
            return [(β, Gβ + Rβ + 1, α, t)]
        elseif startswith(s, "Rsum")
            tail = s[5:end]
            (isempty(tail) || !all(isdigit, tail)) && return NTuple{4,Int}[]
            t = parse(Int, tail)
            return [(β, i, α, t) for i in Gβ+1:Gβ+Rβ]
        elseif startswith(s, "R")
            tail = s[2:end]
            (isempty(tail) || !all(isdigit, tail)) && return NTuple{4,Int}[]
            t = parse(Int, tail)
            return [(β, i, α, t) for i in Gβ+1:Gβ+Rβ]
        else
            length(s) == 2 || return NTuple{4,Int}[]
            isdigit(s[1]) && isdigit(s[2]) || return NTuple{4,Int}[]
            src = parse(Int, s[1])
            t = parse(Int, s[2])
            return [(β, src, α, t)]
        end
    end

    function _append_direction!(states, sign, β, α)
        s = strip(string(states))
        isempty(s) && return
        conns_dir = _parse_code(s, β, α)
        isempty(conns_dir) && return
        mode = parse_coupling_sign_csv(sign)
        γ = default_coupling_gamma_csv(mode)
        len_before = length(connections)
        append!(connections, conns_dir)
        append!(gammas, fill(γ, length(conns_dir)))
        append!(modes, fill(mode, length(conns_dir)))
        if startswith(s, "R") && !startswith(s, "Rany") && !startswith(s, "Rsum")
            push!(tie_groups, collect(len_before + 1:len_before + length(conns_dir)))
        elseif tie_rsum && startswith(s, "Rsum") && length(conns_dir) > 1
            push!(tie_groups, collect(len_before + 1:len_before + length(conns_dir)))
        end
    end

    _append_direction!(e1_states, e1_sign, 1, 2)
    _append_direction!(e2_states, e2_sign, 1, 2)
    _append_direction!(ge_states, ge_sign, 2, 1)

    return connections, gammas, tie_groups, modes
end

function _coupling_prior_mean_from_gammas(default_γ::Vector{Float64}, ncoupling::Int)
    isempty(default_γ) && return fill(0.0, ncoupling)
    length(default_γ) == ncoupling && return default_γ
    length(default_γ) > ncoupling && return default_γ[1:ncoupling]
    return vcat(default_γ, fill(default_γ[end], ncoupling - length(default_γ)))
end

"""
    build_coupled_fit_spec_from_csv_cells(e1, e1s, e2, e2s, ge, ges, key::AbstractString;
        G, R, S, insertstep, transitions, datapath, resultfolder, root,
        noisepriors, elongationtime, datacond, cell, hierarchical, interval, tracetime,
        zeromedian, probfn, method, trace_specs, initprior, tie_rsum, ...)

Build a `Dict{Symbol,Any}` suitable for [`write_run_spec_preset`](@ref) / `fit(; key=..., ...)`
from six CSV cells and shared tracejoint defaults.
"""
function build_coupled_fit_spec_from_csv_cells(
    e1, e1s, e2, e2s, ge, ges, key::AbstractString;
    G::Tuple,
    R::Tuple,
    S::Tuple,
    insertstep::Tuple,
    transitions::Tuple,
    datapath::AbstractString,
    resultfolder::AbstractString,
    root::AbstractString,
    gene::AbstractString="MYC",
    cell::AbstractString="HBEC",
    datacond::Union{AbstractVector{<:AbstractString},AbstractVector{String}}=["gene", "enhancer"],
    noisepriors=([0.0, 0.1, 0.5, 0.15], [0.0, 0.1, 0.9, 0.2]),
    elongationtime::Tuple=(20.0, 5.0),
    initprior::Float64=0.1,
    hierarchical::Bool=true,
    interval::Float64=5.0 / 3.0,
    tracetime::Float64=-1.0,
    zeromedian=Bool[true, true],
    probfn=(prob_Gaussian, prob_Gaussian),
    trace_specs=[],
    tie_rsum::Bool=true,
    maxtime=60.0,
    nchains::Int=16,
    samplesteps::Int=100_000,
    warmupsteps::Int=0,
    propcv::Float64=0.05,
    ratetype::AbstractString="median",
    datacol::Int=3,
    writesamples::Bool=false,
    prerun::Float64=0.0,
    splicetype::AbstractString="",
    decayrate::Float64=1.0,
)
    if noisepriors === nothing || noisepriors == [] || (noisepriors isa AbstractVector && isempty(noisepriors))
        noisepriors = ([0.0, 0.1, 0.5, 0.15], [0.0, 0.1, 0.9, 0.2])
    end
    if probfn === nothing || probfn isa Function
        probfn = (prob_Gaussian, prob_Gaussian)
    elseif !(probfn isa Tuple) || length(probfn) < 2
        probfn = (prob_Gaussian, prob_Gaussian)
    end
    conns, default_γ, tie_groups, modes = csv_row_to_connections_simple(e1, e1s, e2, e2s, ge, ges, G, R; tie_rsum=tie_rsum)
    isempty(conns) && return nothing
    ncoupling = length(conns)
    coupling = ((1, 2), conns, modes)
    coupling_prior_mean = _coupling_prior_mean_from_gammas(default_γ, ncoupling)

    # Base rate + noise priors (two units), then coupling block (matches makescriptcoupled template)
    priormean = Float64[fill(0.001, 2); fill(0.01, length(transitions[1]) - 2); initprior; fill(R[1] / elongationtime[1], R[1]); fill(0.05, max(0, S[1] - insertstep[1] + 1)); 0.03165055618995184; noisepriors[1]]
    priorcv = Float64[fill(1.0, length(transitions[1])); 1.0; fill(0.1, R[1]); fill(0.1, max(0, S[1] - insertstep[1] + 1)); 1.0; [0.1, 0.1, 0.1, 0.1]]
    for i in 2:length(R)
        priormean = vcat(priormean, Float64[fill(0.01, length(transitions[i])); initprior; fill(R[i] / elongationtime[i], R[i]); fill(0.05, max(0, S[i] - insertstep[i] + 1)); 1.0; noisepriors[i]])
        priorcv = vcat(priorcv, Float64[fill(1.0, length(transitions[i])); 1.0; fill(0.1, R[i]); fill(0.1, max(0, S[i] - insertstep[i] + 1)); 1.0; [0.1, 0.1, 0.01, 0.1]])
    end
    priormean = vcat(priormean, coupling_prior_mean)
    priorcv = vcat(priorcv, fill(10.0, ncoupling))

    hypercv_base = Float64[fill(1.0, length(transitions[1])); 1.0; fill(0.1, R[1] - 1); 1.0; fill(1.0, max(0, S[1] - insertstep[1] + 1)); 1.0; fill(0.25, 4)]
    for i in 2:length(R)
        hypercv_base = vcat(hypercv_base, Float64[fill(1.0, length(transitions[i])); 1.0; fill(0.1, R[i] - 1); 1.0; fill(1.0, max(0, S[i] - insertstep[i] + 1)); 1.0; fill(0.25, 4)])
    end
    h = if hierarchical
        priormean = vcat(priormean, vcat(hypercv_base, fill(1.0, ncoupling)))
        priorcv = vcat(priorcv, fill(2.0, length(priorcv)))
        (2, Int[], tuple())
    else
        tuple()
    end

    nr1 = num_rates(transitions[1], R[1], S[1], insertstep[1])
    n1 = nr1 + length(noisepriors[1])
    totalrates = 0
    for i in eachindex(R)
        totalrates += num_rates(transitions[i], R[i], S[i], insertstep[i]) + length(noisepriors[i])
    end

    target_inds = sort(unique([(α == 1 ? t : n1 + t) for (β, s, α, t) in conns]))
    tie_vectors = Vector{Vector{Int}}()
    tied_locals = Set{Int}()
    for g in tie_groups
        length(g) <= 1 && continue
        first_local = g[1]
        rest_locals = g[2:end]
        global_first = totalrates + first_local
        global_rest = [totalrates + i for i in rest_locals]
        push!(tie_vectors, vcat([global_first], global_rest))
        union!(tied_locals, rest_locals)
    end
    fixedeffects = Tuple(tie_vectors)
    all_locals = collect(1:ncoupling)
    fitted_coupling_locals = [i for i in all_locals if !(i in tied_locals)]
    coupling_fp = [totalrates + i for i in fitted_coupling_locals]
    fittedparam = Int[vcat(target_inds, coupling_fp)...]

    method = hierarchical ? (Tsit5(), true) : Tsit5()
    n_u = n_observed_trace_units(coupling)
    trace_specs_eff = if isempty(trace_specs)
        default_trace_specs_for_coupled((interval, 1.0, tracetime), zeromedian, n_u)
    else
        trace_specs
    end

    return Dict{Symbol,Any}(
        :key => string(key),
        :gene => gene,
        :cell => cell,
        :datacond => datacond,
        :datatype => "tracejoint",
        :datapath => datapath,
        :resultfolder => resultfolder,
        :root => root,
        :transitions => transitions,
        :G => G,
        :R => R,
        :S => S,
        :insertstep => insertstep,
        :coupling => coupling,
        :priormean => priormean,
        :priorcv => priorcv,
        :fittedparam => fittedparam,
        :fixedeffects => fixedeffects,
        :noisepriors => noisepriors,
        :hierarchical => h,
        :elongationtime => elongationtime,
        :zeromedian => zeromedian,
        :probfn => probfn,
        :method => method,
        :maxtime => maxtime,
        :nchains => nchains,
        :samplesteps => samplesteps,
        :warmupsteps => warmupsteps,
        :propcv => propcv,
        :ratetype => ratetype,
        :datacol => datacol,
        :writesamples => writesamples,
        :prerun => prerun,
        :splicetype => splicetype,
        :decayrate => decayrate,
        :trace_specs => trace_specs_eff,
        :label => "",
    )
end
