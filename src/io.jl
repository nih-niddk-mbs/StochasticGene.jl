# This file is part of StochasticGene.jl   

# io.jl
### Files for saving and reading mh results

# File handling types and structures

"""
Abstract type for file field definitions
"""
abstract type Fields end

"""
    Result_Fields

Structure for defining fields in result files

Fields:
- `name`: File name
- `label`: Run label
- `cond`: Condition identifier
- `gene`: Gene name
- `model`: Model type
- `nalleles`: Number of alleles
"""
struct Result_Fields <: Fields
    name::String
    label::String
    cond::String
    gene::String
    model::String
    nalleles::String
end

"""
    Summary_Fields

Structure for defining fields in summary files

Fields:
- `name`: File name
- `label`: Run label
- `cond`: Condition identifier
- `model`: Model type
"""
struct Summary_Fields <: Fields
    name::String
    label::String
    cond::String
    model::String
end

"""
    BurstMeasures

Structure for storing burst size statistics

Fields:
- `mean`: Mean burst size
- `std`: Standard deviation of burst size
- `median`: Median burst size
- `mad`: Median absolute deviation
- `quantiles`: Array of burst size quantiles
"""
struct BurstMeasures
    mean::Float64
    std::Float64
    median::Float64
    mad::Float64
    quantiles::Array
end

"""
    decompose_model(model::String)

Parse model string into component parameters

# Arguments
- `model::String`: Model identifier string

# Returns
- Tuple containing:
  - `G`: Number of gene states
  - `R`: Number of pre-RNA steps
  - `S`: Number of splice sites
  - `insertstep`: Reporter insertion step

# Notes
- Model string format: GRSI (e.g., "2121" for G=2, R=1, S=2, insertstep=1)
- For models with more than 4 digits, returns array of tuples
"""
function decompose_model(model::String)
    m = parse(Int, model)
    d = digits(m)
    l = length(d)
    if l > 4
        reshaped_arr = reverse(reshape(d, 4, div(l, 4)))
        return [Tuple(reshaped_arr[i, :]) for i in 1:4]
    else
        return reverse(d)
    end
end

"""
    decompose_cond(cond::String)

Decompose a condition string by checking for hyphens and returning appropriate format.

# Arguments
- `cond::String`: Condition string to decompose

# Returns
- If string contains hyphen: `Tuple{String, String}` - split on first hyphen
- If no hyphen: `String` - original string unchanged

# Examples
```julia
decompose_cond("enhancer-gene21")  # returns ("enhancer", "gene21")
decompose_cond("simple21")         # returns "simple21"
decompose_cond("no-digits")        # returns "no-digits"
decompose_cond("part1-part2-part3") # returns ("part1", "part2-part3")
```

# Notes
- Splits only on the first hyphen if multiple hyphens are present
- Preserves the original string if no hyphen is found
"""
function decompose_cond(cond::String)
    # Check if there's a hyphen in the base string
    if occursin("-", cond)
        # Split by hyphen and return as tuple
        return String.(split(cond, "-"))
    else
        # No hyphen, return the string as is
        return cond
    end
end

"""
    parse_filename(filename::String)

Parse a rate filename to extract all relevant model information and parameters.

# Arguments
- `filename::String`: Rate filename to parse

# Returns
- `datacond`: Data condition (string or tuple from decompose_cond)
- `transitions`: Model transitions structure
- `G`: Gene states (Int or Tuple)
- `R`: RNA states (Int or Tuple) 
- `S`: Splice states (Int or Tuple)
- `insertstep`: Insertion step (Int or Tuple)
- `hierarchical`: Bool indicating if hierarchical model
- `coupling_field`: String coupling field (empty for uncoupled models)

# Examples
```julia
# For uncoupled model
datacond, transitions, G, R, S, insertstep, hierarchical, coupling_field = 
    parse_filename("rates_gene_condition_3401_2.txt")

# For coupled/forced model with tracejoint
datacond, transitions, G, R, S, insertstep, hierarchical, coupling_field = 
    parse_filename("rates_gene_tracejoint_condition31_3401_2.txt")
```

# Notes
- Automatically detects hierarchical models from "-h" in label
- Extracts coupling field from condition for tracejoint models (2-char unidirectional or 4-char reciprocal)
- Uses decompose_cond to process the condition field
- Handles both uncoupled and coupled/forced model types
"""
function parse_filename(filename::String; hlabel="-h")
    parts = fields(filename)
    hierarchical = occursin(hlabel, parts.label)
    G, R, S, insertstep = decompose_model(parts.model)
    transitions = get_transitions(G, parts.label)
    if occursin("tracejoint", parts.label)
        # Coupling suffix: 2-char = unidirectional (e.g. "31"), 4-char = reciprocal (e.g. "3131")
        gene_pos = findfirst("gene", parts.cond)
        if gene_pos !== nothing
            coupling_field = parts.cond[gene_pos.stop+1:end]
            datacond = decompose_cond(parts.cond[1:gene_pos.stop])
        else
            coupling_field = length(parts.cond) >= 4 ? parts.cond[end-3:end] : (length(parts.cond) >= 2 ? parts.cond[end-1:end] : "")
            datacond = decompose_cond(length(parts.cond) >= 4 ? parts.cond[1:end-4] : (length(parts.cond) >= 2 ? parts.cond[1:end-2] : parts.cond))
        end
    else
        coupling_field = ""
        datacond = decompose_cond(parts.cond)
    end
    return datacond, transitions, G, R, S, insertstep, hierarchical, coupling_field
end

"""
    normalize_coupling_field(coupling_field::String)

Normalize a coupling-field string read from filenames or labels.

- **Legacy specs** (`"31"`, `"3131"`, `"R5"`) are returned unchanged.
- **Extended specs** that already contain `','` or `'|'` are returned unchanged.
- **Sanitized extended specs** (where `','` and `'|'` were replaced by `'-'` for filenames,
  e.g. `"24-33-33"` for `"24,33|33"`) are mapped back to a canonical extended form by
  treating all but the last `'-'` as `','` and the last `'-'` as `'|'`.

This keeps filenames shell-safe while letting the rest of the code work with the
intended coupling syntax.
"""
function normalize_coupling_field(coupling_field::String)
    # Empty / trivial
    isempty(coupling_field) && return coupling_field

    # Already extended: contains ',' or '|'
    (occursin(',', coupling_field) || occursin('|', coupling_field)) && return coupling_field

    # Legacy 2- or 4-character specs: digits and/or 'R' only
    if all(c -> (isdigit(c) || c == 'R'), coupling_field) &&
       (length(coupling_field) == 2 || length(coupling_field) == 4)
        return coupling_field
    end

    # Heuristic for sanitized extended specs: contains '-' but no ',' or '|'
    if occursin('-', coupling_field)
        chars = collect(coupling_field)
        dash_idxs = findall(==('-'), chars)
        nd = length(dash_idxs)
        if nd > 0
            for (i, idx) in enumerate(dash_idxs)
                chars[idx] = (i == nd) ? '|' : ','
            end
            return String(chars)
        end
    end

    # Fallback: return as-is
    return coupling_field
end

"""
    make_coupling(coupling_field::String, G, R)

Construct coupling structure from coupling field and model parameters.

# Arguments
- `coupling_field::String`: Coupling field (e.g., "31", "R5", "24,33|33")
- `G`: Number of gene states (Int or Tuple)
- `R`: Number of RNA states (Int or Tuple)

# Returns
- `Tuple`: Coupling structure in format ((1, 2), (sources...), (source_state...), (target_transition...), ncoupling[, coupling_ranges])

# Examples
```julia
# State 3 → State 1 coupling
make_coupling("31", 3, 4)  # returns coupling from state 3 to state 1

# All R states → State 5 coupling
make_coupling("R5", 3, 4)  # returns coupling from states 4,5,6,7 to state 5

# Extended multi-connection (models 9–12)
make_coupling("24,33|33", (3, 3), (0, 0))
```

# Notes
- If coupling_field starts with "R": all RNA states (G+1 to G+R) → target state
- Otherwise: single source state → target state
- Target state is always the last character of coupling_field
- Source state is first character for non-R couplings
- Extended format (models 9–12): contains ',' or '|', e.g. "24,33|33" for multiple (s,t) per direction.
  Sanitized filename forms like "24-33-33" are normalized back to this format internally.
"""
function make_coupling(coupling_field::String, G, R; coupling_ranges=nothing)
    coupling_field = normalize_coupling_field(coupling_field)
    if isempty(coupling_field)
        return tuple()
    end
    Gt = G isa Tuple ? G : (G, G)
    Rt = R isa Tuple ? R : (R, R)
    # Extended multi-connection format (models 9–12)
    if occursin(',', coupling_field) || occursin('|', coupling_field)
        c = make_coupling_extended(coupling_field, Gt, Rt; coupling_ranges=coupling_ranges)
        return c
    end
    if length(coupling_field) == 4
        c = make_coupling_reciprocal(coupling_field, Gt, Rt; coupling_ranges=coupling_ranges)
        return c
    end
    # Unidirectional (2-char or "R" + digit)
    if startswith(coupling_field, "R")
        source = collect(Gt[1]+1:Gt[1]+Rt[1])
    else
        source = parse(Int, coupling_field[1])
    end
    target = parse(Int, coupling_field[2:end])
    base = ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1)
    return coupling_ranges === nothing ? base : (base..., coupling_ranges)
end

"""
    make_coupling_extended(coupling_field::String, G, R)

Construct coupling structure from extended coupling field (models 9–12).
Format: "[dir1]|[dir2]" where each dir is comma-separated "st" pairs (s = digit or R, t = digit).
E.g. "24,33|33" = (1→2): (2,4),(3,3); (2→1): (3,3).
"""
function make_coupling_extended(coupling_field::String, G::Tuple, R::Tuple; coupling_ranges=nothing)
    parts = split(coupling_field, '|'; limit=2)
    # Parse direction 1→2 (target 2): first part
    dir_12 = strip(parts[1])
    sources_2 = Int[]
    source_state_2 = Union{Int,Vector{Int}}[]
    target_transition_2 = Int[]
    for st in split(dir_12, ',')
        st = strip(st)
        isempty(st) && continue
        s_char, t_char = st[1], st[2]
        s = s_char == 'R' ? collect(G[1]+1:G[1]+R[1]) : parse(Int, string(s_char))
        t = parse(Int, string(t_char))
        push!(sources_2, 1)
        push!(source_state_2, s)
        push!(target_transition_2, t)
    end
    # Parse direction 2→1 (target 1): second part if present
    sources_1 = Int[]
    source_state_1 = Union{Int,Vector{Int}}[]
    target_transition_1 = Int[]
    if length(parts) >= 2
        dir_21 = strip(parts[2])
        for st in split(dir_21, ',')
            st = strip(st)
            isempty(st) && continue
            s_char, t_char = st[1], st[2]
            s = s_char == 'R' ? collect(G[2]+1:G[2]+R[2]) : parse(Int, string(s_char))
            t = parse(Int, string(t_char))
            push!(sources_1, 2)
            push!(source_state_1, s)
            push!(target_transition_1, t)
        end
    end
    ncoupling = length(sources_1) + length(sources_2)
    sources = (tuple(sources_1...), tuple(sources_2...))
    source_state = (tuple(source_state_1...), tuple(source_state_2...))
    target_transition = (tuple(target_transition_1...), tuple(target_transition_2...))
    base = ((1, 2), sources, source_state, target_transition, ncoupling)
    return coupling_ranges === nothing ? base : (base..., coupling_ranges)
end

"""
    to_connections(coupling)

Normalize coupling 5-tuple to a flat list of connections in canonical order.

# Arguments
- `coupling`: 5-tuple (unit_model, sources, source_state, target_transition, ncoupling).
  Empty tuple returns empty list.

# Returns
- `Vector{Tuple}`: List of (β, α, s, t) with β = source unit, α = target unit,
  s = source state (Int or Vector{Int} for R states), t = target transition index.
  Order: for α in 1:n for k in 1:length(sources[α]) with β=sources[α][k], s=source_state[α][k], t=target_transition[α][k].

# Notes
- Use this order everywhere for γ (fit, simulator, make_mat_TC).
- Handles legacy (scalar source_state/target_transition per unit) and extended (tuples).
"""
function to_connections(coupling)
    isempty(coupling) && return Tuple{Int,Int,Union{Int,Vector{Int}},Int}[]
    unit_model, sources, source_state, target_transition, _ncoupling = coupling[1], coupling[2], coupling[3], coupling[4], coupling[5]
    n = length(unit_model)
    conns = Tuple{Int,Int,Union{Int,Vector{Int}},Int}[]
    for α in 1:n
        sα = sources[α]
        nconn = length(sα)
        nconn == 0 && continue
        ss = source_state[α]
        tt = target_transition[α]
        for k in 1:nconn
            β = sα[k]
            s = (ss isa Tuple || ss isa AbstractVector) ? ss[k] : ss
            t = (tt isa Tuple || tt isa AbstractVector) ? tt[k] : tt
            push!(conns, (β, α, s, t))
        end
    end
    conns
end

"""
    _connection_source_label(s)

Format source state `s` (Int or Vector{Int} for R states) for use in connection names.
"""
function _connection_source_label(s)
    if s isa Vector || s isa Tuple
        return "R"  # All R states (or multi-state block)
    end
    return "s" * string(s)
end

"""
    connection_name(β, α, s, t; unit_labels=nothing)

Return a short human-readable name for one coupling connection (β→α, state s, transition t).

# Arguments
- `β`: Source unit index (whose state is read).
- `α`: Target unit index (whose transition rate is modulated).
- `s`: Source state: `Int` (e.g. 3) or `Vector{Int}` for R states.
- `t`: Target transition index.
- `unit_labels::Union{Nothing,AbstractVector{String}}`: Optional labels for units (e.g. `["enhancer", "gene"]`).
  If provided, names use labels instead of indices for the "β→α" part.

# Returns
- `String`: e.g. `"2→1_s3t5"` (unit 2 state 3 → unit 1 transition 5) or `"gene→enhancer_s3t5"` if `unit_labels` given.

# Examples
```julia
connection_name(2, 1, 3, 5)                    # "2→1_s3t5"
connection_name(2, 1, 3, 5; unit_labels=["enhancer", "gene"])  # "gene→enhancer_s3t5"
connection_name(1, 2, [4,5,6], 5)              # "1→2_Rt5"
```
"""
function connection_name(β::Int, α::Int, s, t::Int; unit_labels=nothing)
    src = _connection_source_label(s)
    dir = if unit_labels !== nothing && length(unit_labels) >= max(β, α)
        string(unit_labels[β], "→", unit_labels[α])
    else
        string(β, "→", α)
    end
    return dir * "_" * src * "t" * string(t)
end

"""
    coupling_connection_names(coupling; unit_labels=nothing)

Return a vector of names for each coupling constant, in the same order as the rate vector
(i.e. `to_connections(coupling)` and the γ order in combined rate files).

Use these names when labeling columns, writing summaries, or plotting so that "first γ"
is unambiguously tied to the actual connection (e.g. 2→1_s3t5) instead of the 4-char
spec order (s1t1s2t2).

# Arguments
- `coupling`: 5-tuple (or 6-tuple with coupling_ranges) from `make_coupling` / `make_coupling_reciprocal`.
- `unit_labels::Union{Nothing,AbstractVector{String}}`: Optional (e.g. `["enhancer", "gene"]`) to use
  semantic names in the "β→α" part.

# Returns
- `Vector{String}`: One name per coupling parameter, e.g. `["2→1_s3t5", "1→2_s2t3"]` for a "2335" model
  (first γ = gene→enhancer, second γ = enhancer→gene).

# Examples
```julia
c = make_coupling_reciprocal("2335", (3, 3), (0, 0))
coupling_connection_names(c)   # ["2→1_s3t5", "1→2_s2t3"]
coupling_connection_names(c; unit_labels=["enhancer", "gene"])
# ["gene→enhancer_s3t5", "enhancer→gene_s2t3"]
```
"""
function coupling_connection_names(coupling; unit_labels=nothing)
    isempty(coupling) && return String[]
    conns = to_connections(coupling)
    return [connection_name(β, α, s, t; unit_labels=unit_labels) for (β, α, s, t) in conns]
end

"""
    coupling_parameter_labels(coupling; unit_labels=nothing)

Canonical names for coupling parameters, in the same order as the rate vector (and `to_connections`).

Use these for file headers and column labels (e.g. in combined rate files) instead of generic
"γ_1", "γ_2". Returns the same as `coupling_connection_names(coupling; unit_labels=unit_labels)`.

# Arguments
- `coupling`: 5-tuple (or 6-tuple with coupling_ranges) from `make_coupling` / `make_coupling_reciprocal`.
- `unit_labels`: Optional (e.g. `["enhancer", "gene"]`) for semantic unit names in labels.

# Returns
- `Vector{String}`: One name per coupling parameter, e.g. `["2→1_s3t5", "1→2_s2t3"]`.
"""
coupling_parameter_labels(coupling; unit_labels=nothing) = coupling_connection_names(coupling; unit_labels=unit_labels)

"""
    coupling_ranges(coupling) -> Vector{Symbol}

Return the range constraint for each coupling constant (:free, :activate, or :inhibit).
If coupling is a 5-tuple, returns fill(:free, ncoupling).
If coupling has a 6th element, it may be a single Symbol (apply to all) or a tuple/vector of length ncoupling.
"""
function coupling_ranges(coupling)
    isempty(coupling) && return Symbol[]
    ncoupling = coupling[5]
    if length(coupling) >= 6
        r = coupling[6]
        if r isa Symbol
            return fill(r, ncoupling)
        else
            return collect(r)
        end
    end
    return fill(:free, ncoupling)
end

"""
    make_coupling_reciprocal(coupling_field::String, G, R)

Construct reciprocal (bidirectional) coupling structure from 4-character coupling field.

# Arguments
- `coupling_field::String`: 4-character string "s1t1s2t2" where:
  - s1, t1: unit 1→2 direction (enhancer state s1 affects gene transition t1)
  - s2, t2: unit 2→1 direction (gene state s2 affects enhancer transition t2)
  - Source chars (s1, s2) can be "R" for all R states; target chars (t1, t2) are digits
- `G`: Tuple of gene states per unit
- `R`: Tuple of RNA states per unit

# Returns
- `Tuple`: Reciprocal coupling structure ((1, 2), (tuple(2), tuple(1)), (s2, s1), (t2, t1), 2)

# Examples
```julia
# Symmetric: both units in state 3 affect the other's transition 1
make_coupling_reciprocal("3131", (3, 3), (3, 3))

# Asymmetric: enhancer state 3→gene transition 1, gene R states→enhancer transition 5
make_coupling_reciprocal("31R5", (3, 3), (3, 3))
```
"""
function make_coupling_reciprocal(coupling_field::String, G, R; coupling_ranges=nothing)
    if length(coupling_field) != 4
        throw(ArgumentError("coupling_field for reciprocal must be 4 characters (s1t1s2t2), got \"$coupling_field\""))
    end
    # Parse unit 1→2: s1 (char1), t1 (char2)
    s1_char, t1_char = coupling_field[1], coupling_field[2]
    s1 = s1_char == 'R' ? collect(G[1]+1:G[1]+R[1]) : parse(Int, string(s1_char))
    t1 = parse(Int, string(t1_char))
    # Parse unit 2→1: s2 (char3), t2 (char4)
    s2_char, t2_char = coupling_field[3], coupling_field[4]
    s2 = s2_char == 'R' ? collect(G[2]+1:G[2]+R[2]) : parse(Int, string(s2_char))
    t2 = parse(Int, string(t2_char))
    # coupling: sources=(tuple(2), tuple(1)), source_state=(s2, s1), target_transition=(t2, t1)
    base = ((1, 2), (tuple(2), tuple(1)), (s2, s1), (t2, t1), 2)
    return coupling_ranges === nothing ? base : (base..., coupling_ranges)
end

# raterow_dict() = Dict([("ml", 1), ("mean", 2), ("median", 3), ("last", 4)])
# statrow_dict() = Dict([("mean", 1), ("SD", 2), ("median", 3), ("MAD", 4)])

"""
    write_dataframes(resultfolder::String, datapath::String; kwargs...)

Write and assemble model fitting results into CSV files

# Arguments
- `resultfolder::String`: Path to folder containing result files
- `datapath::String`: Path to folder containing input data
- `measure::Symbol=:AIC`: Model selection criterion
- `assemble::Bool=true`: Whether to assemble results into summary files
- `multicond::Bool=false`: Whether to handle multiple conditions
- `datatype::String="rna"`: Type of data ("rna", "trace", etc.)

# Returns
- Nothing, but writes:
  - Individual result files
  - Summary files (if assemble=true)
  - Winner files based on selected measure

# Notes
- Creates hierarchical directory structure for results
- Uses CSV format for better readability and compatibility
- Automatically handles multiple conditions if specified
"""
function write_dataframes(resultfolder::String, datapath::String; measure::Symbol=:AIC, assemble::Bool=true, multicond::Bool=false, datatype::String="rna")
    write_dataframes_only(resultfolder, datapath, assemble=assemble, multicond=multicond, datatype=datatype)
    write_winners(resultfolder, measure)
end

"""
    write_dataframes_only(resultfolder::String, datapath::String; kwargs...)

Write dataframes to CSV files without selecting winners.

# Arguments
- `resultfolder::String`: Path to folder containing result files
- `datapath::String`: Path to folder containing input data
- `assemble::Bool=true`: Whether to assemble results into summary files
- `multicond::Bool=false`: Whether to handle multiple conditions
- `datatype::String="rna"`: Type of data ("rna", "trace", etc.)

# Returns
- Nothing, but writes individual result files and summary files

# Notes
- This function only writes the dataframes, unlike write_dataframes which also selects winners
- Creates CSV files for each dataframe in the results
"""
function write_dataframes_only(resultfolder::String, datapath::String; assemble::Bool=true, multicond::Bool=false, datatype::String="rna")
    dfs = make_dataframes(resultfolder, datapath, assemble, multicond, datatype)
    for df in dfs
        for dff in dfs
            for dfff in dff
                csvfile = joinpath(resultfolder, dfff[1])
                CSV.write(csvfile, dfff[2])
            end
        end
    end
    nothing
end

"""
    write_winners(resultfolder, measure)

Write best performing model parameters based on selected measure

# Arguments
- `resultfolder`: Path to results folder
- `measure`: Model selection criterion (e.g., :AIC, :BIC, :WAIC)

# Returns
- Nothing, but writes:
  - CSV files containing best model parameters
  - Summary statistics for winning models

# Notes
- Automatically selects best model based on measure
- Handles multiple conditions if present
- Preserves original file structure
"""
function write_winners(resultfolder, measure)
    df = best_measure(resultfolder, measure)
    if ~isempty(df)
        for i in eachindex(df)
            csvfile = joinpath(resultfolder, df[i][1])
            CSV.write(csvfile, df[i][2])
        end
    end
    nothing
end

"""
    write_augmented(summaryfile::String, resultfolder::String)

Augment summary file with G=2 burst size, model predicted moments, and fit measures.

# Arguments
- `summaryfile::String`: Path to summary file to augment
- `resultfolder::String`: Path to folder containing result files

# Returns
- Nothing, but writes augmented summary file

# Notes
- Adds burst size statistics for G=2 models
- Includes model predicted moments
- Appends fit quality measures
- Automatically handles relative paths by prepending resultfolder if needed
"""
function write_augmented(summaryfile::String, resultfolder::String)
    if ~ispath(summaryfile)
        summaryfile = joinpath(resultfolder, summaryfile)
    end
    CSV.write(summaryfile, augment_dataframe(read_dataframe(summaryfile), resultfolder))
end

"""
    read_dataframe(csvfile::String)

Read a CSV file into a DataFrame.

# Arguments
- `csvfile::String`: Path to CSV file to read

# Returns
- `DataFrame`: DataFrame containing the CSV data

# Notes
- Uses CSV.File for robust CSV parsing
- Handles various CSV formats automatically
"""
read_dataframe(csvfile::String) = DataFrame(CSV.File(csvfile))

"""
    get_suffix(file::String)

Extract the base name and file extension from a filename.

# Arguments
- `file::String`: Filename to process

# Returns
- `Tuple{String, String}`: (base_name, extension)

# Examples
```julia
get_suffix("rates_gene_condition_3401_2.txt")  # returns ("rates_gene_condition_3401_2", "txt")
get_suffix("summary.csv")                      # returns ("summary", "csv")
```

# Notes
- Removes the last 4 characters (including the dot) to get the base name
- Returns the last 3 characters as the extension
"""
get_suffix(file::String) = chop(file, tail=4), last(file, 3)

"""
    fields(file::String)

Parse a filename to extract structured field information.

# Arguments
- `file::String`: Filename to parse

# Returns
- `Result_Fields` or `Summary_Fields`: Structured field information

# Notes
- For CSV files: expects 4 fields (name, label, cond, model)
- For TXT files: expects 6 fields (name, label, cond, gene, model, nalleles) or 5 fields (name, label, "", gene, model, nalleles)
- Throws ArgumentError for incorrect file name formats
- Does not account for CSV files with less than 4 fields
"""
function fields(file::String)
    file, suffix = get_suffix(file)
    v = split(file, "_")
    if suffix == "csv"
        if length(v) == 4
            s = Summary_Fields(v[1], v[2], v[3], v[4])
        else
            println(file)
            throw(ArgumentError("Incorrect file name format"))
        end
    else
        if length(v) == 6
            s = Result_Fields(v[1], v[2], v[3], v[4], v[5], v[6])
        elseif length(v) == 5
            s = Result_Fields(v[1], v[2], "", v[3], v[4], v[5])
        else
            println(file)
            throw(ArgumentError("Incorrect file name format"))
        end
    end
    return s
end

"""
    isratefile(folder::String)

Check if a folder contains rate files.

# Arguments
- `folder::String`: Path to folder to check

# Returns
- `Bool`: True if folder contains CSV files with "rates" in the name

# Notes
- Looks for files that contain both ".csv" and "rates" in their names
- Uses logical AND operation to ensure both conditions are met
"""
function isratefile(folder::String)
    files = readdir(folder)
    any(occursin.(".csv", files) .& occursin.("rates", files))
end

"""
    isfish(string::String)

Check if a string contains "FISH".

# Arguments
- `string::String`: String to check

# Returns
- `Bool`: True if string contains "FISH"

# Notes
- Case-sensitive search for "FISH" substring
"""
isfish(string::String) = occursin("FISH", string)

"""
    get_genes(file::String)

Extract gene names from a CSV file.

# Arguments
- `file::String`: Path to CSV file containing gene data

# Returns
- `Vector{String}`: Array of gene names from the first column

# Notes
- Assumes CSV format with header
- Returns the first column as gene names
"""
function get_genes(file::String)
    r, header = readdlm(file, ',', header=true)
    return r[:, 1]
end

"""
    get_genes(root, cond, datapath)

Get gene names from a data path with condition filtering.

# Arguments
- `root`: Root directory path
- `cond`: Condition string to filter files
- `datapath`: Path to data directory

# Returns
- `Vector{String}`: Array of gene names matching the condition

# Notes
- Combines root and datapath to form full path
- Delegates to get_genes(cond, datapath)
"""
get_genes(root, cond, datapath) = get_genes(cond, joinpath(root, datapath))

"""
    get_genes(cond, datapath)

Extract gene names from files in a directory that match a condition.

# Arguments
- `cond`: Condition string to match in filenames
- `datapath`: Path to directory containing gene files

# Returns
- `Vector{String}`: Array of gene names from matching files

# Notes
- Searches for files containing the condition string
- Extracts gene name from the first part of filename (before first underscore)
- Returns empty vector if no matching files found
"""
function get_genes(cond, datapath)
    genes = Vector{String}(undef, 0)
    files = readdir(datapath)
    for file in files
        if occursin(cond, file)
            push!(genes, split(file, "_")[1])
        end
    end
    return genes
end

"""
    get_genes(folder, type, label, cond, model)

Get gene names from files matching specific criteria.

# Arguments
- `folder`: Directory containing result files
- `type`: Type of result file (e.g., "rates", "measures")
- `label`: Label to match in filenames
- `cond`: Condition to match in filenames
- `model`: Model identifier to match in filenames

# Returns
- `Vector{String}`: Array of gene names from matching files

# Notes
- Uses get_files to find matching files
- Extracts gene name from each matching file
- Returns empty array if no matching files found
"""
function get_genes(folder, type, label, cond, model)
    genes = Array{String,1}(undef, 0)
    files = get_files(folder, type, label, cond, model)
    for file in files
        push!(genes, get_gene(file))
    end
    return genes
end

"""
    get_files(folder, resultname)

Get files matching a result name from a folder.

# Arguments
- `folder`: Directory to search
- `resultname`: Name of result type to match

# Returns
- `Vector{String}`: Array of matching filenames

# Notes
- Delegates to get_files with additional parameters
- Requires label, cond, and model parameters to be defined in scope
"""
get_files(folder, resultname) =
    get_files(folder::String, resultname, label, cond, model) = get_files(get_resultfiles(folder), resultname, label, cond, model)

"""
    file_indices(parts, resultname, label, cond, model)

Find indices of files matching specific criteria.

# Arguments
- `parts`: Vector of field structures from parsed filenames
- `resultname`: Name of result type to match
- `label`: Label to match
- `cond`: Condition to match
- `model`: Model identifier to match

# Returns
- `BitArray`: Boolean array indicating which files match all criteria

# Notes
- Uses logical AND to combine all matching conditions
- Checks for exact matches on name, label, and condition
- Uses occursin for model matching to allow partial matches
"""
file_indices(parts, resultname, label, cond, model) = (getfield.(parts, :name) .== resultname) .& (getfield.(parts, :label) .== label) .& (getfield.(parts, :cond) .== cond) .& occursin.(model, getfield.(parts, :model))

"""
    get_files(files::Vector, resultname, label, cond, model)

Filter files based on specific criteria.

# Arguments
- `files::Vector`: Array of filenames to filter
- `resultname`: Name of result type to match
- `label`: Label to match in filenames
- `cond`: Condition to match in filenames
- `model`: Model identifier to match in filenames

# Returns
- `Vector{String}`: Array of filenames matching all criteria

# Notes
- Parses each filename to extract fields
- Uses file_indices to find matching files
- Returns subset of input files that match all criteria
"""
function get_files(files::Vector, resultname, label, cond, model)
    parts = fields.(files)
    files[file_indices(parts, resultname, label, cond, model)]
    # files[(getfield.(parts, :name).==resultname).&(getfield.(parts, :label).==label).&(getfield.(parts, :cond).==cond).&(getfield.(parts, :model).==model)]
end

"""
    get_gene(file::String)

Extract gene name from a filename.

# Arguments
- `file::String`: Filename to parse

# Returns
- `String`: Gene name extracted from filename

# Notes
- Uses fields() to parse filename structure
- Returns the gene field from the parsed structure
"""
get_gene(file::String) = fields(file).gene

"""
    get_model(file::String)

Extract model identifier from a filename.

# Arguments
- `file::String`: Filename to parse

# Returns
- `String`: Model identifier extracted from filename

# Notes
- Uses fields() to parse filename structure
- Returns the model field from the parsed structure
"""
get_model(file::String) = fields(file).model

"""
    get_label(file::String)

Extract label from a filename.

# Arguments
- `file::String`: Filename to parse

# Returns
- `String`: Label extracted from filename

# Notes
- Uses fields() to parse filename structure
- Returns the label field from the parsed structure
"""
get_label(file::String) = fields(file).label

"""
    get_cond(file::String)

Extract condition from a filename.

# Arguments
- `file::String`: Filename to parse

# Returns
- `String`: Condition extracted from filename

# Notes
- Uses fields() to parse filename structure
- Returns the cond field from the parsed structure
"""
get_cond(file::String) = fields(file).cond

"""
    get_nalleles(file::String)

Extract number of alleles from a filename.

# Arguments
- `file::String`: Filename to parse

# Returns
- `String`: Number of alleles extracted from filename

# Notes
- Uses fields() to parse filename structure
- Returns the nalleles field from the parsed structure
"""
get_nalleles(file::String) = fields(file).nalleles

"""
    get_fields(parts::Vector{T}, field::Symbol) where {T<:Fields}

Extract unique values for a specific field from a vector of field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures
- `field::Symbol`: Field symbol to extract (e.g., :gene, :model, :label)

# Returns
- `Vector{String}`: Array of unique values for the specified field

# Notes
- Uses getfield to extract the specified field from each structure
- Returns unique values only
- Works with any field type that inherits from Fields
"""
get_fields(parts::Vector{T}, field::Symbol) where {T<:Fields} = unique(getfield.(parts, field))

"""
    get_models(parts::Vector{T}) where {T<:Fields}

Extract unique model identifiers from field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures

# Returns
- `Vector{String}`: Array of unique model identifiers

# Notes
- Delegates to get_fields with :model symbol
"""
get_models(parts::Vector{T}) where {T<:Fields} = get_fields(parts, :model)

"""
    get_genes(parts::Vector{T}) where {T<:Fields}

Extract unique gene names from field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures

# Returns
- `Vector{String}`: Array of unique gene names

# Notes
- Delegates to get_fields with :gene symbol
"""
get_genes(parts::Vector{T}) where {T<:Fields} = get_fields(parts, :gene)

"""
    get_conds(parts::Vector{T}) where {T<:Fields}

Extract unique conditions from field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures

# Returns
- `Vector{String}`: Array of unique conditions

# Notes
- Delegates to get_fields with :cond symbol
"""
get_conds(parts::Vector{T}) where {T<:Fields} = get_fields(parts, :cond)

"""
    get_labels(parts::Vector{T}) where {T<:Fields}

Extract unique labels from field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures

# Returns
- `Vector{String}`: Array of unique labels

# Notes
- Delegates to get_fields with :label symbol
"""
get_labels(parts::Vector{T}) where {T<:Fields} = get_fields(parts, :label)

"""
    get_names(parts::Vector{T}) where {T<:Fields}

Extract unique names from field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures

# Returns
- `Vector{String}`: Array of unique names

# Notes
- Delegates to get_fields with :name symbol
"""
get_names(parts::Vector{T}) where {T<:Fields} = get_fields(parts, :name)

"""
    get_nalleles(parts::Vector{T}) where {T<:Fields}

Extract unique allele counts from field structures.

# Arguments
- `parts::Vector{T}`: Vector of field structures

# Returns
- `Vector{String}`: Array of unique allele counts

# Notes
- Delegates to get_fields with :nalleles symbol
"""
get_nalleles(parts::Vector{T}) where {T<:Fields} = get_fields(parts, :nalleles)

"""
    get_resultfiles(folder::String)

Get all result files from a folder.

# Arguments
- `folder::String`: Directory to search

# Returns
- `Vector{String}`: Array of result filenames

# Notes
- Delegates to get_resultfiles with readdir(folder)
- Looks for files with both ".txt" extension and "_" in name
"""
get_resultfiles(folder::String) = get_resultfiles(readdir(folder))

"""
    get_resultfiles(files::Vector)

Filter files to find result files.

# Arguments
- `files::Vector`: Array of filenames to filter

# Returns
- `Vector{String}`: Array of result filenames

# Notes
- Looks for files with both ".txt" extension and "_" in name
- Uses logical AND to ensure both conditions are met
"""
get_resultfiles(files::Vector) = files[occursin.(".txt", files).&occursin.("_", files)]

"""
    get_resultfiles(folder::String, name)

Get result files matching a specific name pattern.

# Arguments
- `folder::String`: Directory to search
- `name`: Name pattern to match in filenames

# Returns
- `Vector{String}`: Array of matching result filenames

# Notes
- First gets all result files from the folder
- Then filters by the specified name pattern
"""
function get_resultfiles(folder::String, name)
    files = get_resultfiles(readdir(folder))
    files[occursin.(name, files)]
end

"""
    get_resultfile(type::String, infolder, label, gene, G, R, S, insertstep, nalleles)

Generate the full path for a result file.

# Arguments
- `type::String`: Type of result file (e.g., "rates", "measures")
- `infolder`: Directory containing the file
- `label`: Label for the file
- `gene`: Gene name
- `G`: Number of gene states
- `R`: Number of RNA steps (0 for simple models)
- `S`: Number of splice states
- `insertstep`: Insertion step
- `nalleles`: Number of alleles

# Returns
- `String`: Full path to the result file

# Notes
- Uses different filename generation for R=0 (simple models) vs R>0 (complex models)
- Joins the type prefix with the generated filename
"""
function get_resultfile(type::String, infolder, label, gene, G, R, S, insertstep, nalleles)
    if R == 0
        name = filename(label, gene, G, nalleles)
    else
        name = filename(label, gene, G, R, S, insertstep, nalleles)
    end
    joinpath(infolder, type * name)
end

"""
    get_measurefiles(folder::String)

Get all measure files from a folder.

# Arguments
- `folder::String`: Directory to search

# Returns
- `Vector{String}`: Array of measure filenames

# Notes
- Delegates to get_resultfiles with "measures" pattern
"""
get_measurefiles(folder::String) = get_resultfiles(folder, "measures")

"""
    get_summaryfiles(folder::String)

Get all summary files from a folder.

# Arguments
- `folder::String`: Directory to search

# Returns
- `Vector{String}`: Array of summary filenames

# Notes
- Delegates to get_summaryfiles with readdir(folder)
- Looks for files with both ".csv" extension and "_" in name
"""
get_summaryfiles(folder::String) = get_summaryfiles(readdir(folder))

"""
    get_summaryfiles(files::Vector)

Filter files to find summary files.

# Arguments
- `files::Vector`: Array of filenames to filter

# Returns
- `Vector{String}`: Array of summary filenames

# Notes
- Looks for files with both ".csv" extension and "_" in name
- Uses logical AND to ensure both conditions are met
"""
get_summaryfiles(files::Vector) = files[occursin.(".csv", files).&occursin.("_", files)]

"""
    get_summaryfiles(files::Vector, name)

Get summary files matching a specific name pattern.

# Arguments
- `files::Vector`: Array of filenames to filter
- `name`: Name pattern to match in filenames

# Returns
- `Vector{String}`: Array of matching summary filenames

# Notes
- First filters for summary files, then filters by name pattern
"""
get_summaryfiles(files::Vector, name) = files[occursin.(".csv", files).&occursin.(name, files)]

"""
    get_ratesummaryfiles(files::Vector)

Get rate summary files from a list of files.

# Arguments
- `files::Vector`: Array of filenames to filter

# Returns
- `Vector{String}`: Array of rate summary filenames

# Notes
- Delegates to get_summaryfiles with "rates" pattern
"""
get_ratesummaryfiles(files::Vector) = get_summaryfiles(files, "rates")

"""
    get_ratesummaryfiles(folder::String)

Get rate summary files from a folder.

# Arguments
- `folder::String`: Directory to search

# Returns
- `Vector{String}`: Array of rate summary filenames

# Notes
- Delegates to get_ratesummaryfiles with summary files from the folder
"""
get_ratesummaryfiles(folder::String) = get_ratesummaryfiles(get_summaryfiles(folder))

"""
    get_measuresummaryfiles(files::Vector)

Get measure summary files from a list of files.

# Arguments
- `files::Vector`: Array of filenames to filter

# Returns
- `Vector{String}`: Array of measure summary filenames

# Notes
- Delegates to get_summaryfiles with "measures" pattern
"""
get_measuresummaryfiles(files::Vector) = get_summaryfiles(files, "measures")

"""
    get_measuresummaryfiles(folder::String)

Get measure summary files from a folder.

# Arguments
- `folder::String`: Directory to search

# Returns
- `Vector{String}`: Array of measure summary filenames

# Notes
- Delegates to get_measuresummaryfiles with summary files from the folder
"""
get_measuresummaryfiles(folder::String) = get_measuresummaryfiles(get_summaryfiles(folder))

"""
    get_burstsummaryfiles(files::Vector)

Get burst summary files from a list of files.

# Arguments
- `files::Vector`: Array of filenames to filter

# Returns
- `Vector{String}`: Array of burst summary filenames

# Notes
- Delegates to get_summaryfiles with "burst" pattern
"""
get_burstsummaryfiles(files::Vector) = get_summaryfiles(files, "burst")

"""
    get_burstsummaryfiles(folder::String)

Get burst summary files from a folder.

# Arguments
- `folder::String`: Directory to search

# Returns
- `Vector{String}`: Array of burst summary filenames

# Notes
- Delegates to get_burstsummaryfiles with summary files from the folder
"""
get_burstsummaryfiles(folder::String) = get_burstsummaryfiles(get_summaryfiles(folder))

"""
    write_moments(outfile, genelist, cond, datapath, root)

Write gene expression moments (mean and variance) to a CSV file.

# Arguments
- `outfile`: Output file path for the moments data
- `genelist`: List of gene names to process
- `cond`: Condition identifier for the data
- `datapath`: Path to the data directory
- `root`: Root directory path

# Returns
- Nothing, but writes CSV file with columns: Gene, Expression Mean, Expression Variance

# Notes
- Creates a CSV file with header row
- For each gene, calculates mean and variance from RNA histogram data
- Uses get_histogram_rna, mean_histogram, and var_histogram functions
"""
function write_moments(outfile, genelist, cond, datapath, root)
    f = open(outfile, "w")
    writedlm(f, ["Gene" "Expression Mean" "Expression Variance"], ',')
    for gene in genelist
        h = get_histogram_rna(gene, cond, datapath, root)
        writedlm(f, [gene mean_histogram(h) var_histogram(h)], ',')
    end
    close(f)
end

"""
    write_histograms(resultfolder, ratefile, cell, datacond, G::Int, datapath::String, root, outfolder="histograms")

Write model-predicted histograms to files for each rate parameter.

# Arguments
- `resultfolder`: Directory containing result files
- `ratefile`: Name of rate file to process
- `cell`: Cell type identifier
- `datacond`: Data condition string (can contain hyphens for multiple conditions)
- `G::Int`: Number of gene states
- `datapath::String`: Path to input data directory
- `root`: Root directory path
- `outfolder`: Output folder name (default: "histograms")

# Returns
- Nothing, but writes histogram files for each rate parameter and condition

# Notes
- Reads rate parameters from CSV file with header
- Creates output folder if it doesn't exist
- Splits datacond on hyphens to handle multiple conditions
- For each rate parameter, generates histograms and writes to separate files
- File naming: "{gene}_{condition}.txt"
"""
function write_histograms(resultfolder, ratefile, cell, datacond, G::Int, datapath::String, root, outfolder="histograms")
    ratefile = joinpath(resultfolder, ratefile)
    rates, head = readdlm(ratefile, ',', header=true)
    outfolder = joinpath(resultfolder, outfolder)
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    cond = string.(split(datacond, "-"))
    for r in eachrow(rates)
        h = histograms(r, cell, cond, G, datapath, root)
        for i in eachindex(cond)
            f = open(joinpath(outfolder, string(r[1]) * cond[i] * ".txt"), "w")
            writedlm(f, h[i])
            close(f)
        end
    end
end

"""
    assemble_all(folder::String, multicond::Bool=false)

Assemble all result files in a folder into summary files.

# Arguments
- `folder::String`: Directory containing result files
- `multicond::Bool=false`: Whether to handle multiple conditions

# Returns
- Nothing, but creates summary files for rates, measures, stats, and burst sizes

# Notes
- Main entry point for assembling results
- Extracts all unique labels, conditions, and models from result files
- Calls assemble_all for each combination of label, condition, and model
- Only processes combinations that have rate files available
"""
function assemble_all(folder::String, multicond::Bool=false)
    files = get_resultfiles(folder)
    parts = fields.(files)
    labels = get_labels(parts)
    conds = get_conds(parts)
    models = get_models(parts)
    names = get_names(parts)
    # if isempty(fittedparams)
    #     fittedparams = collect(1:num_rates(models[1])-1)
    # end
    assemble_all(folder, files, labels, conds, models, names, multicond)
end

"""
    assemble_all(folder::String, files::Vector, labels::Vector, conds::Vector, models::Vector, names, multicond::Bool=false)

Assemble results for all combinations of labels, conditions, and models.

# Arguments
- `folder::String`: Directory containing result files
- `files::Vector`: Array of result filenames
- `labels::Vector`: Array of unique labels
- `conds::Vector`: Array of unique conditions
- `models::Vector`: Array of unique models
- `names`: Array of result file types
- `multicond::Bool=false`: Whether to handle multiple conditions

# Returns
- Nothing, but creates summary files for each valid combination

# Notes
- Iterates through all combinations of labels, conditions, and models
- Only processes combinations that have rate files available
- Calls assemble_all for each valid combination
"""
function assemble_all(folder::String, files::Vector, labels::Vector, conds::Vector, models::Vector, names, multicond::Bool=false)
    parts = fields.(files)
    for l in labels, c in conds, g in models
        any(file_indices(parts, "rates", l, c, g) .== 1) && assemble_all(folder, files, l, c, g, names, multicond)
    end
end

"""
    assemble_all(folder::String, files::Vector, label::String, cond::String, model::String, names, multicond::Bool=false)

Assemble all result types for a specific label, condition, and model combination.

# Arguments
- `folder::String`: Directory containing result files
- `files::Vector`: Array of result filenames
- `label::String`: Specific label to process
- `cond::String`: Specific condition to process
- `model::String`: Specific model to process
- `names`: Array of result file types
- `multicond::Bool=false`: Whether to handle multiple conditions

# Returns
- Nothing, but creates summary files for the specified combination

# Notes
- Assembles rates, measures, and stats for the combination
- Conditionally assembles burst sizes for non-simple models (model != "1")
- Conditionally assembles optimized results if available
- Uses the rate labels returned from assemble_rates for optimized results
"""
function assemble_all(folder::String, files::Vector, label::String, cond::String, model::String, names, multicond::Bool=false)
    labels = assemble_rates(folder, files, label, cond, model, multicond)
    assemble_measures(folder, files, label, cond, model)
    assemble_stats(folder, files, label, cond, model)
    if model != "1" && "burst" ∈ names
        assemble_burst_sizes(folder, files, label, cond, model)
    end
    if "optimized" ∈ names
        assemble_optimized(folder, files, label, cond, model, labels)
    end
end

"""
    assemble_files(folder::String, files::Vector, outfile::String, header, readfunction)

Assemble data from multiple files into a single output file.

# Arguments
- `folder::String`: Directory containing the input files
- `files::Vector`: Array of filenames to process
- `outfile::String`: Output file path
- `header`: Header row for the output file
- `readfunction`: Function to read data from each file

# Returns
- Nothing, but writes assembled data to outfile

# Notes
- Creates output file with specified header
- For each file, extracts gene name and reads data using readfunction
- Writes gene name and data as a row in the output file
- Skips processing if files array is empty
"""
function assemble_files(folder::String, files::Vector, outfile::String, header, readfunction)
    if ~isempty(files)
        f = open(outfile, "w")
        writedlm(f, header, ',')
        for file in files
            gene = get_gene(file)
            r = readfunction(joinpath(folder, file))
            writedlm(f, [gene r], ',')
        end
        close(f)
    end
end

"""
    assemble_rates(folder::String, files::Vector, label::String, cond::String, model::String, multicond::Bool=false, readfunction=readml)

Assemble rate parameters from multiple files into a single CSV file.

# Arguments
- `folder::String`: Directory containing result files
- `files::Vector`: Array of result filenames
- `label::String`: Label to match
- `cond::String`: Condition to match
- `model::String`: Model to match
- `multicond::Bool=false`: Whether to handle multiple conditions
- `readfunction=readml`: Function to read rate data (default: readml for maximum likelihood)

# Returns
- `Matrix`: Rate labels extracted from the first rate file

# Notes
- Creates output file named "rates_{label}_{cond}_{model}.csv"
- Filters files to get only rate files matching the criteria
- Extracts labels from the first rate file
- Uses split_conditions to handle multiple conditions if multicond=true
- Calls assemble_files to combine all rate data
"""
function assemble_rates(folder::String, files::Vector, label::String, cond::String, model::String, multicond::Bool=false, readfunction=readml)
    outfile = joinpath(folder, "rates_" * label * "_" * cond * "_" * model * ".csv")
    ratefiles = get_files(files, "rates", label, cond, model)
    labels = readdlm(joinpath(folder, ratefiles[1]), ',', header=true)[2]
    # header = ratelabels(model, split(cond, "-"))
    header = ratelabels(labels, split_conditions(cond, multicond))
    assemble_files(folder, ratefiles, outfile, header, readfunction)
    return labels
end

"""
    assemble_measures(folder::String, files, label::String, cond::String, model::String)

Assemble fit measures from multiple files into a single CSV file.

# Arguments
- `folder::String`: Directory containing result files
- `files`: Array of result filenames
- `label::String`: Label to match
- `cond::String`: Condition to match
- `model::String`: Model to match

# Returns
- Nothing, but writes measures file named "measures_{label}_{cond}_{model}.csv"

# Notes
- Creates output file with comprehensive header including all fit statistics
- Filters files to get only measure files matching the criteria
- For each file, extracts gene name, number of alleles, and fit measures
- Includes error handling to skip missing or empty files
- Writes warning messages for problematic files
- Header includes: Gene, Nalleles, LogMaxLikelihood, normalized_LL, n_obs, n_params, Deviance, WAIC, WAIC SE, AIC, Acceptance, Samples, Temperature, Rhat, ESS, Geweke, MCSE
"""
function assemble_measures(folder::String, files, label::String, cond::String, model::String)
    outfile = joinpath(folder, "measures_" * label * "_" * cond * "_" * model * ".csv")
    header = ["Gene" "Nalleles" "LogMaxLikelihood" "normalized_LL" "n_obs" "n_params" "Deviance" "WAIC" "WAIC SE" "AIC" "Acceptance" "Samples" "Temperature" "Rhat" "ESS" "Geweke" "MCSE"]
    files = get_files(files, "measures", label, cond, model)
    f = open(outfile, "w")
    writedlm(f, header, ',')
    for file in files
        try
            fullfile = joinpath(folder, file)
            if !isfile(fullfile) || isempty(read(fullfile))
                @warn "Skipping missing or empty file: $file"
                continue
            end
            gene = get_gene(file)
            nalleles = get_nalleles(file)
            r = readmeasures(fullfile)
            writedlm(f, [gene nalleles r], ',')
        catch e
            @warn "Skipping file due to error: $file" exception = (e, catch_backtrace())
            continue
        end
    end
    close(f)
end


"""
    assemble_measures_model(folder::String, files::Vector, outfile::String)

Assemble measure files into a single summary file with model comparison.

# Arguments
- `folder::String`: Directory containing measure files
- `files::Vector`: Array of measure filenames to process
- `outfile::String`: Output file path

# Returns
- Nothing, but writes assembled measures to outfile

# Notes
- Creates output file with header for model comparison
- For each file, extracts model name and fit measures
- Removes common prefixes and suffixes from filenames to get model names
- Header includes: Model, LogMaxLikelihood, normalized_LL, n_obs, n_params, Deviance, WAIC, WAIC_SE, AIC, Acceptance, Samples, Temperature, Rhat, ESS, Geweke, MCSE
"""
function assemble_measures_model(folder::String, files::Vector, outfile::String)
    header = ["Model" "LogMaxLikelihood" "normalized_LL" "n_obs" "n_params" "Deviance" "WAIC" "WAIC_SE" "AIC" "Acceptance" "Samples" "Temperature" "Rhat" "ESS" "Geweke" "MCSE"]
    f = open(outfile, "w")
    writedlm(f, header, ',')
    for file in files
        nalleles = get_nalleles(file)
        r = readmeasures(joinpath(folder, file))
        writedlm(f, [remove_string(file, "measures_trace-HBEC-", "_$nalleles.txt") r], ',')
    end
    close(f)
end
"""
    assemble_measures_model(folder::String, label::String, cond::String, gene::String)

Assemble measures for a specific gene across different models.

# Arguments
- `folder::String`: Directory containing result files
- `label::String`: Label to match
- `cond::String`: Condition to match
- `gene::String`: Gene name to match

# Returns
- Nothing, but writes measures file named "measures_{label}_{cond}_{gene}.csv"

# Notes
- Creates output file for comparing measures across models for a specific gene
- Filters files to get measure files matching label and condition
- Uses empty string for model to match all models
- Delegates to assemble_measures_model with files and output path
"""
function assemble_measures_model(folder::String, label::String, cond::String, gene::String)
    outfile = joinpath(folder, "measures_" * label * "_" * cond * "_" * gene * ".csv")
    files = get_files(get_resultfiles(folder), "measures", label, cond, "")
    assemble_measures_model(folder, files, outfile)
end
"""
    assemble_measures_model(folder::String)

Assemble all measure files in a folder into a single summary file.

# Arguments
- `folder::String`: Directory containing measure files

# Returns
- Nothing, but writes "measures.csv" file

# Notes
- Gets all measure files from the folder
- Creates output file named "measures.csv"
- Delegates to assemble_measures_model with files and output path
"""
function assemble_measures_model(folder::String)
    assemble_measures_model(folder, get_measurefiles(folder), joinpath(folder, "measures.csv"))
end
"""
    assemble_optimized(folder::String, files, label::String, cond::String, model::String, labels)

Assemble optimized parameter results into a single CSV file.

# Arguments
- `folder::String`: Directory containing result files
- `files`: Array of result filenames
- `label::String`: Label to match
- `cond::String`: Condition to match
- `model::String`: Model to match
- `labels`: Rate labels to use as header

# Returns
- Nothing, but writes optimized file named "optimized_{label}_{cond}_{model}.csv"

# Notes
- Creates output file with rate labels as header
- Filters files to get only optimized files matching the criteria
- Uses read_optimized function to read data from each file
- Calls assemble_files to combine all optimized data
"""
function assemble_optimized(folder::String, files, label::String, cond::String, model::String, labels)
    outfile = joinpath(folder, "optimized_" * label * "_" * cond * "_" * model * ".csv")
    assemble_files(folder, get_files(files, "optimized", label, cond, model), outfile, labels, read_optimized)
end

"""
    assemble_stats(folder::String, files, label::String, cond::String, model::String)

Assemble parameter statistics into a single CSV file.

# Arguments
- `folder::String`: Directory containing result files
- `files`: Array of result filenames
- `label::String`: Label to match
- `cond::String`: Condition to match
- `model::String`: Model to match

# Returns
- Nothing, but writes stats file named "stats_{label}_{cond}_{model}.csv"

# Notes
- Creates output file with statistical labels as header
- Filters files to get only param-stats files matching the criteria
- Extracts labels from the first stats file
- Uses statlabels to format the header with statistical measures
- Calls assemble_files to combine all statistical data
"""
function assemble_stats(folder::String, files, label::String, cond::String, model::String)
    outfile = joinpath(folder, "stats_" * label * "_" * cond * "_" * model * ".csv")
    statfiles = get_files(files, "param-stats", label, cond, model)
    labels = readdlm(joinpath(folder, statfiles[1]), ',', header=true)[2]
    assemble_files(folder, statfiles, outfile, statlabels(labels), readstats)
end

"""
    assemble_burst_sizes(folder, files, label, cond, model)

Assemble burst size statistics into a single CSV file.

# Arguments
- `folder`: Directory containing result files
- `files`: Array of result filenames
- `label`: Label to match
- `cond`: Condition to match
- `model`: Model to match

# Returns
- Nothing, but writes burst file named "burst_{label}_{cond}_{model}.csv"

# Notes
- Creates output file with burst statistics header
- Filters files to get only burst files matching the criteria
- Header includes: Gene, BurstMean, BurstSD, BurstMedian, BurstMAD
- Uses read_burst function to read burst size data
- Calls assemble_files to combine all burst data
"""
function assemble_burst_sizes(folder, files, label, cond, model)
    outfile = joinpath(folder, "burst_" * label * "_" * cond * "_" * model * ".csv")
    assemble_files(folder, get_files(files, "burst", label, cond, model), outfile, ["Gene" "BurstMean" "BurstSD" "BurstMedian" "BurstMAD"], read_burst)
end

"""
    assemble_info(folder, files, label, cond, model)

Assemble model information into a single CSV file.

# Arguments
- `folder`: Directory containing result files
- `files`: Array of result filenames
- `label`: Label to match
- `cond`: Condition to match
- `model`: Model to match

# Returns
- Nothing, but writes info file named "info_{label}_{cond}_{model}.csv"

# Notes
- Creates output file with information header
- Filters files to get only info files matching the criteria
- Header includes: Gene, Nalleles, Interval
- Uses read_info function to read model information
- Calls assemble_files to combine all info data
"""
function assemble_info(folder, files, label, cond, model)
    outfile = joinpath(folder, "info_" * label * "_" * cond * "_" * model * ".csv")
    assemble_files(folder, get_files(files, "info", label, cond, model), outfile, ["Gene" "Nalleles" "Interval"], read_info)
end

"""
    rlabels(model::String)

Generate rate labels for a simple gene model.

# Arguments
- `model::String`: Model identifier (number of gene states)

# Returns
- `Matrix{String}`: 1×n matrix of rate labels

# Notes
- Parses model string to get number of gene states (G)
- Creates labels for gene state transitions: Rate01, Rate10, Rate12, Rate21, etc.
- Adds "Eject" and "Decay" labels at the end
- For G states, creates 2*(G-1) transition labels plus 2 additional labels
- Example: model="3" creates ["Rate01" "Rate10" "Rate12" "Rate21" "Eject" "Decay"]
"""
function rlabels(model::String)
    G = parse(Int, model)
    n = G - 1
    Grates = Array{String,2}(undef, 1, 2 * n)
    for i = 0:n-1
        Grates[1, 2*i+1] = "Rate$i$(i+1)"
        Grates[1, 2*i+2] = "Rate$(i+1)$i"
    end
    return [Grates "Eject" "Decay"]
end

"""
    rlabels(model::String, conds::Vector)

Generate rate labels for a model with multiple conditions.

# Arguments
- `model::String`: Model identifier (number of gene states)
- `conds::Vector`: Vector of condition identifiers

# Returns
- `Matrix{String}`: Matrix of rate labels with condition suffixes

# Notes
- Uses rlabels(model) to get base rate labels
- For single condition, returns base labels unchanged
- For multiple conditions, appends condition suffix to each label
- Example: model="2", conds=["A", "B"] creates ["Rate01_A" "Rate10_A" "Eject_A" "Decay_A" "Rate01_B" "Rate10_B" "Eject_B" "Decay_B"]
"""
function rlabels(model::String, conds::Vector)
    nsets = length(conds)
    r = rlabels(model)
    if nsets == 1
        return r
    else
        rates = r .* conds[1]
        for i = 2:nsets
            rates = [rates r .* conds[i]]
        end
        return rates
    end
end

"""
    rlabels(labels::Matrix, conds::Vector)

Generate rate labels from existing labels with condition suffixes.

# Arguments
- `labels::Matrix`: Matrix of existing rate labels
- `conds::Vector`: Vector of condition identifiers

# Returns
- `Matrix{String}`: Matrix of rate labels with condition suffixes

# Notes
- Similar to rlabels(model::String, conds::Vector) but uses existing labels
- For single condition, returns labels unchanged
- For multiple conditions, appends condition suffix to each label
- Reshapes condition suffixes to match label dimensions
"""
function rlabels(labels::Matrix, conds::Vector)
    nsets = length(conds)
    r = labels
    if nsets == 1
        return r
    else
        rates = r .* conds[1]
        for i = 2:nsets
            rates = [rates r .* reshape(conds[i], 1, length(conds))]
        end
        return rates
    end
end

"""
    rlabels(model::String, conds, fittedparams)

Generate rate labels for fitted parameters only.

# Arguments
- `model::String`: Model identifier
- `conds`: Condition identifiers
- `fittedparams`: Indices of fitted parameters

# Returns
- `Matrix{String}`: Rate labels for fitted parameters only

# Notes
- Delegates to rlabels(model::String, conds::Vector)
- Returns only the columns specified by fittedparams
"""
rlabels(model::String, conds, fittedparams) = rlabels(model, conds)[1:1, fittedparams]

"""
    rlabels(labels::Matrix, conds, fittedparams)

Generate rate labels for fitted parameters from existing labels.

# Arguments
- `labels::Matrix`: Matrix of existing rate labels
- `conds`: Condition identifiers
- `fittedparams`: Indices of fitted parameters

# Returns
- `Matrix{String}`: Rate labels for fitted parameters only

# Notes
- Delegates to rlabels(labels::Matrix, conds::Vector)
- Returns only the columns specified by fittedparams
"""
rlabels(labels::Matrix, conds, fittedparams) = rlabels(labels, conds)[1:1, fittedparams]

"""
    ratelabels(labels::Matrix, conds)

Generate rate labels with "Gene" column header.

# Arguments
- `labels::Matrix`: Matrix of rate labels
- `conds`: Condition identifiers

# Returns
- `Matrix{String}`: Rate labels with "Gene" as first column

# Notes
- Delegates to rlabels(labels::Matrix, conds::Vector)
- Prepends "Gene" column to the result
"""
ratelabels(labels::Matrix, conds) = ["Gene" rlabels(labels, conds)]

"""
    ratelabels(labels::Matrix)

Generate rate labels with "Gene" column header.

# Arguments
- `labels::Matrix`: Matrix of rate labels

# Returns
- `Matrix{String}`: Rate labels with "Gene" as first column

# Notes
- Prepends "Gene" column to the input labels
"""
ratelabels(labels::Matrix) = ["Gene" labels]

"""
    rlabels(model::AbstractGMmodel)

Generate rate labels for an abstract gene model.

# Arguments
- `model::AbstractGMmodel`: Abstract gene model instance

# Returns
- `Matrix{String}`: Matrix of rate labels

# Notes
- Creates labels for each gene transition in model.Gtransitions
- Adds "Eject" and "Decay" labels
- For HMMReporter models, adds noise parameter labels
- Reshapes result to 1×n matrix
"""
function rlabels(model::AbstractGMmodel)
    labels = String[]
    for t in model.Gtransitions
        push!(labels, "Rate$(t[1])$(t[2])")
    end
    push!(labels, "Eject")
    push!(labels, "Decay")
    if typeof(model.reporter) == HMMReporter
        for i in 1:model.reporter.n
            push!(labels, "noiseparam$i")
        end
    end
    reshape(labels, 1, :)
end

"""
    rlabels_GRSM(transitions, R, S, reporter, unit=:"")

Generate rate labels for GRSM (Gene-RNA-Splice-Model) components.

# Arguments
- `transitions`: Gene transition structure
- `R`: Number of RNA steps
- `S`: Number of splice states
- `reporter`: Reporter type
- `unit`: Unit identifier for multi-unit models (default: "")

# Returns
- `Matrix{String}`: Matrix of GRSM rate labels

# Notes
- Creates labels for gene transitions with unit suffix
- Adds "Initiate" label for transcription initiation
- Adds "Rshift" labels for RNA processing steps (R-1 labels)
- Adds "Eject" label for RNA release
- Adds "Splice" labels for splicing steps (S labels)
- Adds "Decay" label for RNA degradation
- For HMMReporter, adds noise parameter labels
- Reshapes result to 1×n matrix
"""
function rlabels_GRSM(transitions, R, S, reporter, unit=:"")
    labels = String[]
    for t in transitions
        push!(labels, "Rate$(unit)_$(t[1])$(t[2])")
    end
    push!(labels, "Initiate$unit")
    for i in 1:R-1
        push!(labels, "Rshift$(unit)_$i")
    end
    push!(labels, "Eject$unit")
    for i in 1:S
        push!(labels, "Splice$(unit)_$i")
    end
    push!(labels, "Decay$unit")
    if typeof(reporter) == HMMReporter
        for i in 1:reporter.n
            push!(labels, "noiseparam$(unit)_$i")
        end
    end
    reshape(labels, 1, :)
end

"""
    rlabels_GRSM(model::AbstractGRSMmodel)

Generate rate labels for an abstract GRSM model.

# Arguments
- `model::AbstractGRSMmodel`: Abstract GRSM model instance

# Returns
- `Matrix{String}`: Matrix of GRSM rate labels

# Notes
- Handles coupled models with multiple units
- For multi-unit models, generates labels for each unit separately
- Adds coupling labels for coupled models
- For single-unit models, uses standard GRSM label generation
- Delegates to rlabels_GRSM for individual unit labels
"""
function rlabels_GRSM(model::AbstractGRSMmodel)
    if hastrait(model, :coupling)
        labels = Array{String}(undef, 1, 0)
        if length(model.G) > 1
            for i in eachindex(model.G)
                labels = hcat(labels, rlabels_GRSM(model.Gtransitions[i], model.R[i], model.S[i], model.reporter[i], i))
            end
        else
            labels = hcat(labels, rlabels_GRSM(model.Gtransitions, model.R, model.S, model.reporter))
        end
        cplabels = model.trait.coupling.labels
        if !isnothing(cplabels) && length(cplabels) == model.trait.coupling.ncoupling
            labels = hcat(labels, reshape(cplabels, 1, :))
        else
            for i in 1:model.trait.coupling.ncoupling
                labels = hcat(labels, ["Coupling_$i"])
            end
        end
    else
        labels = rlabels_GRSM(model.Gtransitions, model.R, model.S, model.reporter)
    end
    labels
end

"""
    rlabels_GRSM_grid(labels)

Add grid probability label to GRSM rate labels.

# Arguments
- `labels`: Matrix of existing GRSM rate labels

# Returns
- `Matrix{String}`: Rate labels with "GridProb" column added

# Notes
- Horizontally concatenates input labels with "GridProb" column
- Used for grid-based GRSM models
"""
function rlabels_GRSM_grid(labels)
    hcat(labels, "GridProb")
end

"""
    rlabels_GRSM_hierarchical(labelsin, model)

Generate hierarchical rate labels for GRSM models.

# Arguments
- `labelsin`: Matrix of base rate labels
- `model`: GRSM model with hierarchical trait

# Returns
- `Matrix{String}`: Hierarchical rate labels

# Notes
- Creates hyperparameter labels by prefixing "hyper_" to base labels
- Repeats hyperparameter labels for each hyperset (nhypersets times)
- Adds individual parameter labels for each individual (nindividuals times)
- Used for hierarchical GRSM models with shared and individual parameters
"""
function rlabels_GRSM_hierarchical(labelsin, model)
    labels = "hyper_" .* labelsin
    for i in 2:model.trait.hierarchical.nhypersets
        labels = hcat(labels, "hyper_" .* labelsin)
    end
    for i in 1:model.trait.hierarchical.nindividuals
        labels = hcat(labels, labelsin)
    end
    labels
end

"""
    rlabels(model::GRSMmodel)

Generate rate labels for a GRSM model with optional traits.

# Arguments
- `model::GRSMmodel`: GRSM model instance

# Returns
- `Matrix{String}`: Complete rate labels for the model

# Notes
- Starts with base GRSM labels from rlabels_GRSM
- Adds grid probability label if model has :grid trait
- Adds hierarchical labels if model has :hierarchical trait
- Handles multiple traits by chaining label modifications
"""
function rlabels(model::GRSMmodel)
    labels = rlabels_GRSM(model)
    if hastrait(model, :grid)
        labels = rlabels_GRSM_grid(labels)
    end
    if hastrait(model, :hierarchical)
        labels = rlabels_GRSM_hierarchical(labels, model)
    end
    labels
end


"""
    statlabels(model::String, conds, fittedparams)

Generate statistical labels for parameter statistics.

# Arguments
- `model::String`: Model identifier
- `conds`: Condition identifiers
- `fittedparams`: Indices of fitted parameters

# Returns
- `Matrix{String}`: Statistical labels with "Gene" column

# Notes
- Creates labels for statistical measures: _Mean, _SD, _Median, _MAD, _CI2.5, _CI97.5
- Uses rlabels to get base rate labels for fitted parameters
- Appends statistical suffixes to each rate label
- Prepends "Gene" column to the result
- Example: Rate01 becomes Rate01_Mean, Rate01_SD, Rate01_Median, etc.
"""
function statlabels(model::String, conds, fittedparams)
    label = ["_Mean", "_SD", "_Median", "_MAD", "_CI2.5", "_CI97.5"]
    Grates = rlabels(model, conds, fittedparams)
    rates = Matrix{String}(undef, 1, 0)
    for l in label
        rates = [rates Grates .* l]
    end
    return ["Gene" rates]
end

"""
    statlabels(labels::Matrix)

Generate statistical labels from existing rate labels.

# Arguments
- `labels::Matrix`: Matrix of existing rate labels

# Returns
- `Matrix{String}`: Statistical labels with "Gene" column

# Notes
- Creates labels for statistical measures: _Mean, _SD, _Median, _MAD, _CI2.5, _CI97.5
- Appends statistical suffixes to each input label
- Prepends "Gene" column to the result
- Similar to statlabels(model::String, conds, fittedparams) but uses existing labels
"""
function statlabels(labels::Matrix)
    label = ["_Mean", "_SD", "_Median", "_MAD", "_CI2.5", "_CI97.5"]
    rates = Matrix{String}(undef, 1, 0)
    for l in label
        rates = [rates labels .* l]
    end
    return ["Gene" rates]
end

"""
    optlabels(model::String, conds, fittedparams)

Generate labels for optimized parameter results.

# Arguments
- `model::String`: Model identifier
- `conds`: Condition identifiers
- `fittedparams`: Indices of fitted parameters

# Returns
- `Matrix{String}`: Optimized labels with "Gene", rate labels, "LL", and "Convergence"

# Notes
- Uses rlabels to get rate labels for fitted parameters
- Adds "LL" (log-likelihood) and "Convergence" columns
- Prepends "Gene" column to the result
- Used for optimized parameter output files
"""
optlabels(model::String, conds, fittedparams) = ["Gene" rlabels(model, conds, fittedparams) "LL" "Convergence"]

"""
    optlabels(labels::Matrix, conds, fittedparams)

Generate labels for optimized parameter results from existing labels.

# Arguments
- `labels::Matrix`: Matrix of existing rate labels
- `conds`: Condition identifiers
- `fittedparams`: Indices of fitted parameters

# Returns
- `Matrix{String}`: Optimized labels with "Gene", rate labels, "LL", and "Convergence"

# Notes
- Uses rlabels with existing labels to get rate labels for fitted parameters
- Adds "LL" (log-likelihood) and "Convergence" columns
- Prepends "Gene" column to the result
- Used for optimized parameter output files
"""
optlabels(labels::Matrix, conds, fittedparams) = ["Gene" rlabels(labels, conds) "LL" "Convergence"]

"""
    filename(label::String, gene::String, G::Int, R::Int, S::Int, insertstep::Int, nalleles::Int)

Generate filename for complex gene models.

# Arguments
- `label::String`: Label for the file
- `gene::String`: Gene name
- `G::Int`: Number of gene states
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep::Int`: Insertion step
- `nalleles::Int`: Number of alleles

# Returns
- `String`: Generated filename

# Notes
- Creates model string by concatenating G, R, S, and insertstep
- Delegates to filename(label, gene, model, nalleles)
- Used for complex models with RNA and splicing steps
"""
filename(label::String, gene::String, G::Int, R::Int, S::Int, insertstep::Int, nalleles::Int) = filename(label, gene, "$G" * "$R" * "$S" * "$insertstep", "$(nalleles)")

"""
    filename(label::String, gene, G::Int, nalleles::Int)

Generate filename for simple gene models.

# Arguments
- `label::String`: Label for the file
- `gene`: Gene name
- `G::Int`: Number of gene states
- `nalleles::Int`: Number of alleles

# Returns
- `String`: Generated filename

# Notes
- Creates model string from G only
- Delegates to filename(label, gene, model, nalleles)
- Used for simple models without RNA or splicing steps
"""
filename(label::String, gene, G::Int, nalleles::Int) = filename(label, gene, "$G", "$nalleles")

"""
    filename(label::String, gene::String, model::String, nalleles::String)

Generate filename with string allele count.

# Arguments
- `label::String`: Label for the file
- `gene::String`: Gene name
- `model::String`: Model identifier
- `nalleles::String`: Number of alleles as string

# Returns
- `String`: Generated filename

# Notes
- Creates filename in format: "_label_gene_model_nalleles.txt"
- Used when allele count is already a string
"""
filename(label::String, gene::String, model::String, nalleles::String) = "_" * label * "_" * gene * "_" * model * "_" * nalleles * ".txt"

"""
    filename(label::String, gene::String, model::String, nalleles::Int)

Generate filename with integer allele count.

# Arguments
- `label::String`: Label for the file
- `gene::String`: Gene name
- `model::String`: Model identifier
- `nalleles::Int`: Number of alleles as integer

# Returns
- `String`: Generated filename

# Notes
- Converts allele count to string
- Delegates to filename(label, gene, model, nalleles::String)
- Used when allele count is an integer
"""
filename(label::String, gene::String, model::String, nalleles::Int) = "_" * label * "_" * gene * "_" * model * "_" * "$nalleles" * ".txt"

"""
    filename(label, gene, G::Tuple, R, S, insertstep, nalleles)

Generate filename for multi-unit models.

# Arguments
- `label`: Label for the file
- `gene`: Gene name
- `G::Tuple`: Tuple of gene states for each unit
- `R`: RNA steps (can be tuple for multi-unit)
- `S`: Splice states (can be tuple for multi-unit)
- `insertstep`: Insertion steps (can be tuple for multi-unit)
- `nalleles`: Number of alleles

# Returns
- `String`: Generated filename

# Notes
- Uses create_modelstring to generate model string from tuples
- Delegates to filename(label, gene, model, nalleles)
- Used for multi-unit models with different parameters per unit
"""
function filename(label, gene, G::Tuple, R, S, insertstep, nalleles)
    # m = ""
    # for i in eachindex(G)
    #     m *= "$(G[i])$(R[i])$(S[i])$(insertstep[i])"
    # end
    m = create_modelstring(G::Tuple, R, S, insertstep)
    filename(label, gene, m, nalleles)
end


"""
    filename(data, model::AbstractGRSMmodel)
    filename(data, model::AbstractGMmodel)
    filename(data, model::GRSMmodel)

Generate output filenames for different model types.

# Arguments
- `data`: Data structure containing label and gene information
- `model`: Model instance (AbstractGRSMmodel, AbstractGMmodel, or GRSMmodel)

# Returns
- `String`: Generated filename

# Notes
- Extracts label and gene from data structure
- Extracts model parameters from model instance
- Delegates to appropriate filename function based on model type
- Used for generating consistent filenames across different model types
"""
filename(data, model::AbstractGRSMmodel) = filename(data.label, data.gene, model.G, model.R, model.S, model.insertstep, model.nalleles)
filename(data, model::AbstractGMmodel) = filename(data.label, data.gene, model.G, model.nalleles)
filename(data, model::GRSMmodel) = filename(data.label, data.gene, model.G, model.R, model.S, model.insertstep, model.nalleles)

"""
    writeall(path::String, fits::Fit, stats::Stats, measures::Measures, data, temp, model::AbstractGeneTransitionModel; optimized=0, burst=0, writesamples=false)

Write all model fitting results to files.

# Arguments
- `path::String`: Directory path for output files
- `fits::Fit`: Fit results structure
- `stats::Stats`: Statistical results structure
- `measures::Measures`: Fit measures structure
- `data`: Data structure used for fitting
- `temp`: Temperature parameter
- `model::AbstractGeneTransitionModel`: Model instance
- `optimized=0`: Optimized parameters (optional)
- `burst=0`: Burst size statistics (optional)
- `writesamples=false`: Whether to write sample data (optional)

# Returns
- Nothing, but writes multiple output files

# Notes
- Creates output directory if it doesn't exist
- Writes rates, measures, parameter stats, and info files
- For hierarchical GRSM models, writes shared parameters
- Conditionally writes optimized and burst files if provided
- Conditionally writes sample data if writesamples=true
- Uses consistent filename generation for all output files
"""
function writeall(path::String, fits::Fit, stats::Stats, measures::Measures, data, temp, model::AbstractGeneTransitionModel; optimized=0, burst=0, writesamples=false)
    if !isdir(path)
        mkpath(path)
    end
    name = filename(data, model)
    labels = rlabels(model)
    write_rates(joinpath(path, "rates" * name), fits, stats, model, labels)
    write_measures(joinpath(path, "measures" * name), fits, measures, deviance(fits, data, model), temp, data, model)
    write_param_stats(joinpath(path, "param-stats" * name), stats, model, labels)
    write_info(joinpath(path, "info" * name), fits, data, model, labels)
    if typeof(model) <: GRSMmodel && hastrait(model, :hierarchical)
        write_hierarchy(joinpath(path, "shared" * name), fits, stats, model, labels)
    end
    if optimized != 0
        write_optimized(joinpath(path, "optimized" * name), optimized)
    end
    if burst != 0
        write_burstsize(joinpath(path, "burst" * name), burst)
    end
    if writesamples
        write_array(joinpath(path, "ll_sampled_rates" * name), fits.ll)
        write_array(joinpath(path, "sampled_rates" * name), permutedims(inverse_transform_params(fits.param, model)))
    end
end

"""
    write_rates(file::String, fits::Fit, stats, model, labels)

Write rate parameters to a file.

# Arguments
- `file::String`: Output file path
- `fits::Fit`: Fit results structure
- `stats`: Statistical results structure
- `model`: Model instance
- `labels`: Rate labels for the file header

# Returns
- Nothing, but writes rate file with 4 rows:
  1. Maximum likelihood parameters
  2. Mean posterior parameters
  3. Median posterior parameters
  4. Last accepted sample parameters

# Notes
- Writes labels as header row
- Writes maximum likelihood parameters from fits.parml
- Writes mean posterior parameters from stats.meanparam
- Writes median posterior parameters from stats.medparam
- Writes last sample parameters from fits.param[:, end]
- Uses get_rates to extract rate parameters from parameter vectors
"""
function write_rates(file::String, fits::Fit, stats, model, labels)
    f = open(file, "w")

    writedlm(f, labels, ',')  # labels
    writedlm(f, [get_rates(fits.parml, model)], ',')  # max posterior
    writedlm(f, [get_rates(stats.meanparam, model, false)], ',')  # mean posterior
    writedlm(f, [get_rates(stats.medparam, model, false)], ',')  # median posterior
    writedlm(f, [get_rates(fits.param[:, end], model)], ',')  # last sample

    # if haskey(model.trait, :hierarchical)
    #     writedlm(f, rlabels(model), ',')  # labels
    #     writedlm(f, [get_rates_hierarchical(fits.parml, model)], ',')  # max posterior
    #     writedlm(f, [get_rates_hierarchical(stats.meanparam, model, false)], ',')  # mean posterior
    #     writedlm(f, [get_rates_hierarchical(stats.medparam, model, false)], ',')  # median posterior
    #     writedlm(f, [get_rates_hierarchical(fits.param[:, end], model)], ',')  # last sample
    # else
    #     writedlm(f, rlabels(model), ',')  # labels
    #     writedlm(f, [get_rates(fits.parml, model)], ',')  # max posterior
    #     writedlm(f, [get_rates(stats.meanparam, model, false)], ',')  # mean posterior
    #     writedlm(f, [get_rates(stats.medparam, model, false)], ',')  # median posterior
    #     writedlm(f, [get_rates(fits.param[:, end], model)], ',')  # last sample
    # end
    close(f)
end

"""
    write_rates(file::String, fits::Fit, stats, model)

Write rate parameters to a file using default labels.

# Arguments
- `file::String`: Output file path
- `fits::Fit`: Fit results structure
- `stats`: Statistical results structure
- `model`: Model instance

# Returns
- Nothing, but writes rate file

# Notes
- Delegates to write_rates with labels generated from rlabels(model)
- Uses model to generate appropriate rate labels
"""
function write_rates(file::String, fits::Fit, stats, model)
    write_rates(file, fits, stats, model, rlabels(model))
end

"""
    remove_rates(r, transitions, R, S, insertstep, nreporters, setnumber)

Remove rate parameters for a specific parameter set.

# Arguments
- `r`: Rate parameter matrix
- `transitions`: Gene transition structure
- `R`: Number of RNA steps
- `S`: Number of splice states
- `insertstep`: Insertion step
- `nreporters`: Number of reporter parameters
- `setnumber`: Set number to remove

# Returns
- `Vector{Int}`: Indices of parameters to keep

# Notes
- Calculates total number of parameters per set
- Identifies the range of parameters to remove for the specified set
- Returns indices of all parameters except those in the specified set
- Used for multi-set models where some parameter sets should be excluded
"""
function remove_rates(r, transitions, R, S, insertstep, nreporters, setnumber)
    n = num_rates(transitions, R, S, insertstep) + nreporters
    removeset = n*(setnumber-1)+1:n*setnumber
    setdiff(1:size(r, 2), removeset)
end


"""
    write_hierarchy(file::String, fits::Fit, stats, model::GRSMmodel, labels)

Write hierarchy parameters into a file for hierarchical models.

# Arguments
- `file::String`: Output file path
- `fits::Fit`: Fit results structure
- `stats`: Statistical results structure
- `model::GRSMmodel`: GRSM model with hierarchical trait
- `labels`: Rate labels for the file header

# Returns
- Nothing, but writes hierarchy file with shared parameters

# Notes
- Writes only the shared (hyper) parameters for hierarchical models
- Uses first nrates columns from labels where nrates is the number of shared parameters
- Writes 4 rows: max likelihood, mean, median, and last sample
- Used for hierarchical models to separate shared from individual parameters
"""
function write_hierarchy(file::String, fits::Fit, stats, model::GRSMmodel, labels)
    f = open(file, "w")
    writedlm(f, labels[1:1, 1:model.trait.hierarchical.nrates], ',')  # labels
    writedlm(f, [get_rates(fits.parml, model)[1:model.trait.hierarchical.nrates]], ',')  # max posterior
    writedlm(f, [get_rates(stats.meanparam, model, false)[1:model.trait.hierarchical.nrates]], ',')  # mean posterior
    writedlm(f, [get_rates(stats.medparam, model, false)[1:model.trait.hierarchical.nrates]], ',')  # median posterior
    writedlm(f, [get_rates(fits.param[:, end], model)[1:model.trait.hierarchical.nrates]], ',')  # last sample
    close(f)
end

"""
    write_hierarchy(file::String, fits::Fit, stats, model::GRSMmodel)

Write hierarchy parameters using default labels.

# Arguments
- `file::String`: Output file path
- `fits::Fit`: Fit results structure
- `stats`: Statistical results structure
- `model::GRSMmodel`: GRSM model with hierarchical trait

# Returns
- Nothing, but writes hierarchy file

# Notes
- Delegates to write_hierarchy with labels generated from rlabels(model)
- Uses model to generate appropriate rate labels
"""
function write_hierarchy(file::String, fits::Fit, stats, model::GRSMmodel)
    write_hierarchy(file, fits, stats, model, rlabels(model))
end

"""
    write_measures(file::String, fits::Fit, measures::Measures, dev, temp, data, model)

Write fit measures into a file.

# Arguments
- `file::String`: Output file path
- `fits::Fit`: Fit results structure
- `measures::Measures`: Fit measures structure
- `dev`: Deviance value
- `temp`: Temperature parameter
- `data`: Data structure used for fitting
- `model`: Model instance

# Returns
- Nothing, but writes measures file with comprehensive fit statistics

# Notes
- Calculates n_obs based on data type (histogram, RNA count, or trace data)
- Calculates n_params from model.fittedparam
- Calculates normalized log-likelihood
- Writes multiple rows with different statistics:
  - Row 1: LogMaxLikelihood, normalized_LL, n_obs, n_params, mean_ll, std_ll, quantiles, WAIC, AIC
  - Row 2: Deviance
  - Row 3: Acceptance rate and total samples
  - Row 4: Temperature
  - Row 5: Rhat values
  - Row 6: ESS values
  - Row 7: Geweke values
  - Row 8: MCSE values
  - Row 9: Maximum Rhat
  - Row 10: Minimum ESS
  - Row 11: Maximum Geweke
  - Row 12: Maximum MCSE
"""
function write_measures(file::String, fits::Fit, measures::Measures, dev, temp, data, model)
    f = open(file, "w")
    # Calculate n_obs based on data type
    if is_histogram_compatible(data)
        n_obs = sum(datahistogram(data))
    elseif typeof(data) <: RNACountData
        n_obs = length(data.countsRNA)
    elseif typeof(data) <: AbstractTraceData
        n_obs = length(data.trace[1])
    else
        n_obs = 1  # Default to 1 for other data types
    end

    # Calculate n_params from model
    n_params = length(model.fittedparam)

    # Calculate normalized LL
    normalized_ll = fits.llml / n_obs

    writedlm(f, [fits.llml normalized_ll Int(n_obs) Int(n_params) mean(fits.ll) std(fits.ll) quantile(fits.ll, [0.025; 0.5; 0.975])' measures.waic[1] measures.waic[2] aic(fits)], ',')
    writedlm(f, dev, ',')
    writedlm(f, [fits.accept fits.total], ',')
    writedlm(f, temp, ',')
    writedlm(f, measures.rhat', ',')
    writedlm(f, measures.ess', ',')
    writedlm(f, measures.geweke', ',')
    writedlm(f, measures.mcse', ',')
    writedlm(f, maximum(measures.rhat), ',')
    writedlm(f, minimum(measures.ess), ',')
    writedlm(f, maximum(measures.geweke), ',')
    writedlm(f, maximum(measures.mcse), ',')
    close(f)
end

"""
    write_param_stats(file, stats::Stats, model, labels)

Write parameter statistics into a file.

# Arguments
- `file`: Output file path
- `stats::Stats`: Statistical results structure
- `model`: Model instance
- `labels`: Parameter labels for the file header

# Returns
- Nothing, but writes parameter statistics file

# Notes
- Writes labels for fitted parameters only
- Writes multiple rows with different statistics:
  - Row 1: Mean parameters
  - Row 2: Standard deviation parameters
  - Row 3: Median parameters
  - Row 4: Median absolute deviation parameters
  - Row 5: Quantile parameters
  - Row 6: Correlation parameters
  - Row 7: Covariance parameters
  - Row 8: Log covariance parameters
"""
function write_param_stats(file, stats::Stats, model, labels)
    f = open(file, "w")
    writedlm(f, labels[1:1, model.fittedparam], ',')
    writedlm(f, stats.meanparam', ',')
    writedlm(f, stats.stdparam', ',')
    writedlm(f, stats.medparam', ',')
    writedlm(f, stats.madparam', ',')
    writedlm(f, stats.qparam, ',')
    writedlm(f, stats.corparam, ',')
    writedlm(f, stats.covparam, ',')
    writedlm(f, stats.covlogparam, ',')
    close(f)
end

"""
    write_param_stats(file, stats::Stats, model)

Write parameter statistics using default labels.

# Arguments
- `file`: Output file path
- `stats::Stats`: Statistical results structure
- `model`: Model instance

# Returns
- Nothing, but writes parameter statistics file

# Notes
- Delegates to write_param_stats with labels generated from rlabels(model)
- Uses model to generate appropriate parameter labels
"""
function write_param_stats(file, stats::Stats, model)
    write_param_stats(file, stats, model, rlabels(model))
end

"""
    write_optimized(file::String, optimized)

Write optimized parameter results to a file.

# Arguments
- `file::String`: Output file path
- `optimized`: Optimized results structure

# Returns
- Nothing, but writes optimized file with 3 rows:
  1. Optimized parameters (exponentiated)
  2. Minimum objective value
  3. Convergence status

# Notes
- Exponentiates the optimized parameters using exp()
- Writes the minimum objective value from optimization
- Writes the convergence status (true/false)
- Used for storing results from optimization algorithms
"""
function write_optimized(file::String, optimized)
    f = open(file, "w")
    writedlm(f, exp.(Optim.minimizer(optimized))', ',')
    writedlm(f, Optim.minimum(optimized), ',')
    writedlm(f, Optim.converged(optimized), ',')
    close(f)
end

"""
    write_burstsize(file::String, b::BurstMeasures)

Write burst size statistics to a file.

# Arguments
- `file::String`: Output file path
- `b::BurstMeasures`: Burst size statistics structure

# Returns
- Nothing, but writes burst file with 5 rows:
  1. Mean burst size
  2. Standard deviation of burst size
  3. Median burst size
  4. Median absolute deviation
  5. Burst size quantiles

# Notes
- Writes comprehensive burst size statistics
- Used for analyzing transcriptional bursting behavior
- Quantiles provide distribution information beyond mean/median
"""
function write_burstsize(file::String, b::BurstMeasures)
    f = open(file, "w")
    writedlm(f, b.mean, ',')
    writedlm(f, b.std, ',')
    writedlm(f, b.median, ',')
    writedlm(f, b.mad, ',')
    writedlm(f, b.quantiles, ',')
    close(f)
end

"""
    write_info(file::String, fits, data, model, labels)

Write model information and metadata to a file.

# Arguments
- `file::String`: Output file path
- `fits`: Fit results structure
- `data`: Data structure used for fitting
- `model`: Model instance
- `labels`: Parameter labels

# Returns
- Nothing, but writes info file with model metadata

# Notes
- Writes prior information for fitted parameters
- Handles hierarchical models differently (uses fittedpriors)
- Writes transformed prior means and standard deviations
- Writes model structure information (Gtransitions)
- Writes data metadata (label, gene, nalleles)
- For trace data, writes interval information
- Used for documenting model setup and priors
"""
function write_info(file::String, fits, data, model, labels)
    f = open(file, "w")
    if typeof(model) <: GRSMmodel && hasproperty(model.trait, :hierarchical)
        # writedlm(f, labels[1:1, 1:num_all_parameters(model)], ',')  # labels
        # println(model.trait) 
        writedlm(f, [" " labels[1:1, model.trait.hierarchical.fittedpriors]], ',')  # labels
        writedlm(f, ["prior mean" apply_transform(mean.(model.rateprior), model.transforms.f_inv[model.trait.hierarchical.fittedpriors])'], ',')
    else
        writedlm(f, [" " labels[1:1, model.fittedparam]], ',')  # labels
        # writedlm(f, ["prior mean" apply_transform(mean.(model.rateprior), model.transforms.f_inv[model.fittedparam])'], ',')
    end

    writedlm(f, ["transformed prior mean" reshape(mean.(model.rateprior), 1, :)], ',')
    writedlm(f, ["transformedprior std" reshape(std.(model.rateprior), 1, :)], ',')
    writedlm(f, ["Gtransitions" model.Gtransitions], ',')
    writedlm(f, ["label" data.label], ',')
    writedlm(f, ["gene" data.gene], ',')
    writedlm(f, ["nalleles" model.nalleles], ',')
    if typeof(data) <: AbstractTraceData
        writedlm(f, ["interval" data.interval], ',')
    end
    # writedlm(f, ["G state probability" residenceprob_G(get_rates(fits.parml, model), model.G)], ',')
    # writedlm(f, ["ON probability" onstate_prob(get_rates(fits.parml, model), model)], ',')
    close(f)
end

"""
    write_array(file::String, d::Array)

Write array data to a file.

# Arguments
- `file::String`: Output file path
- `d::Array`: Array data to write

# Returns
- Nothing, but writes array data to file

# Notes
- Uses writedlm to write array data
- Sets header=false to avoid automatic header generation
- Used for writing sample data, likelihood arrays, etc.
"""
write_array(file::String, d::Array) = writedlm(file, d, header=false)

"""
    get_row()

Get mapping from rate type names to row indices.

# Returns
- `Dict{String, Int}`: Mapping of rate types to row numbers

# Notes
- "ml" → row 1 (maximum likelihood)
- "mean" → row 2 (mean posterior)
- "median" → row 3 (median posterior)
- "last" → row 4 (last accepted sample)
- Used for selecting which row to read from rate files
"""
get_row() = Dict([("ml", 1); ("mean", 2); ("median", 3); ("last", 4)])

"""
    get_ratetype()

Get inverse mapping from row indices to rate type names.

# Returns
- `Dict{Int, String}`: Mapping of row numbers to rate types

# Notes
- Inverse of get_row() function
- Row 1 → "ml" (maximum likelihood)
- Row 2 → "mean" (mean posterior)
- Row 3 → "median" (median posterior)
- Row 4 → "last" (last accepted sample)
- Used for identifying rate types from row numbers
"""
get_ratetype() = invert_dict(get_row())

"""
    occursin_file(a, b, file::String)

Determine if string a or string b occurs in file (case insensitive).

# Arguments
- `a`: First string to search for
- `b`: Second string to search for
- `file::String`: Filename to search in

# Returns
- `Bool`: True if both strings are found in filename

# Notes
- Case-insensitive search using regex
- Excludes .DS_Store files automatically
- Uses word boundaries to avoid partial matches
- If a is empty, only searches for b
- If b is empty, only searches for a
- If both are empty, returns false
- Used for filtering files by gene and condition names
"""
function occursin_file(a, b, file::String)
    occursin(Regex("DS_Store", "i"), file) && return false

    function token_regex(token)
        isempty(token) && return ""
        # Match token not surrounded by alphanumerics (word boundary for filenames)
        return "(?<![A-Za-z0-9])" * escape_string(token) * "(?![A-Za-z0-9])"
    end

    ra = token_regex(a)
    rb = token_regex(b)

    if isempty(a)
        return occursin(Regex(rb, "i"), file)
    elseif isempty(b)
        return occursin(Regex(ra, "i"), file)
    else
        return occursin(Regex(ra, "i"), file) && occursin(Regex(rb, "i"), file)
    end
end

"""
    read_rna(filepath::String)

Read RNA count histogram from a file (first column = counts).

# Arguments
- `filepath`: Path to the histogram file

# Returns
- `Tuple{Int, Vector{Float64}}`: (histogram length, histogram data)

# Notes
- Truncates if longer than 300 elements (keeps 99th percentile, max 1000); pads with zeros if fewer than 4.
"""
function read_rna(filepath::String)
    h = readfile(filepath)[:, 1]
    if length(h) > 300
        h = truncate_histogram(h, 0.99, 1000)
    end
    if length(h) < 4
        h = vcat(h, zeros(4 - length(h)))
    end
    println("Histogram count: ", sum(h))
    return length(h), h
end

"""
    read_rna(gene, cond, datapath)

Read RNA histogram from a file `{gene}_{cond}.txt` in the given directory.

# Arguments
- `gene`: Gene name
- `cond`: Condition identifier (e.g. `""` for files like CANX_.txt)
- `datapath`: Path to directory containing the file

# Returns
- `Tuple{Int, Vector{Float64}}`: (histogram length, histogram data)

# Notes
- Constructs path as `joinpath(datapath, "{gene}_{cond}.txt")` and calls `read_rna(filepath)`.
"""
function read_rna(gene, cond, datapath)
    t = joinpath(datapath, "$gene" * "_" * "$cond.txt")
    read_rna(t)
end

"""
    read_rnacount(gene, cond, datapath)

Read RNA count data from a file.

# Arguments
- `gene`: Gene name
- `cond`: Condition identifier
- `datapath`: Path to data directory

# Returns
- `Tuple{Vector{Int}, Vector{Float64}, Int}`: (RNA counts, yield factors, histogram length)

# Notes
- Constructs filename as "{gene}_{cond}.txt"
- Reads first column as RNA counts (rounded to integers)
- Reads second column as yield factors
- Calculates histogram length as 99th percentile + 1 (minimum 1)
- Used for reading individual RNA count measurements
"""
function read_rnacount(gene, cond, datapath)
    t = joinpath(datapath, "$gene" * "_" * "$cond.txt")
    # println(t)
    # c = readfile(gene, cond, datapath)
    c = readfile(t)
    countsRNA = round.(Int, c[:, 1])
    yield = c[:, 2]  # Renamed from yieldfactor to yield
    nhist = round(Int, max(quantile(countsRNA, 0.99), 1) + 1)
    return countsRNA, yield, nhist
end

"""
    readfiles(gene::String, cond::String, datapath::Vector)

Read data from multiple files.

# Arguments
- `gene::String`: Gene name
- `cond::String`: Condition identifier
- `datapath::Vector`: Vector of data directory paths

# Returns
- `Vector{Vector}`: Array of data from each file

# Notes
- Reads data from each directory in datapath
- Uses readfile(gene, cond, path) for each path
- Returns vector of data arrays
- Used for reading data from multiple sources or replicates
"""
function readfiles(gene::String, cond::String, datapath::Vector)
    c = Vector{Vector}(undef, 0)
    for i in eachindex(datapath)
        push!(c, readfile(gene, cond, datapath[i]))
    end
    c
end

"""
    readfile(gene::AbstractString, cond::AbstractString, path::AbstractString)

Read file if name includes gene and condition.

# Arguments
- `gene::AbstractString`: Gene name to search for
- `cond::AbstractString`: Condition to search for
- `path::AbstractString`: Path to search in

# Returns
- `Matrix{Float64}`: Data from matching file

# Notes
- If path is a file, reads it directly
- Otherwise searches recursively through directory tree
- Uses occursin_file to find files containing both gene and condition
- Returns data from first matching file found
- Used for flexible file discovery in data directories
"""
function readfile(gene::AbstractString, cond::AbstractString, path::AbstractString)
    if isfile(path)
        return readfile(path)
    else
        for (root, dirs, files) in walkdir(path)
            for file in files
                target = joinpath(root, file)
                if occursin_file(gene, cond, target)
                    return readfile(target)
                end
            end
        end
    end
end

"""
    readfile(file::String)

Read file with floats accounting for delimiter and headers.

# Arguments
- `file::String`: File path to read

# Returns
- `Matrix{Float64}`: Numeric data from file

# Notes
- Automatically detects CSV files and uses appropriate delimiter
- Handles files with or without headers
- Converts string data to float64
- Skips header row if first row contains strings
- Falls back to comma delimiter if first row contains commas
- Throws ArgumentError if file doesn't exist
"""
function readfile(file::String)
    if isfile(file)
        if occursin("csv", file)
            c = readfile_csv(file)
        else
            c = readdlm(file)
            if typeof(c[1]) <: AbstractString && occursin(",", c[1])
                c = readdlm(file, ',')
            end
        end
        if typeof(c[1, 1]) <: AbstractString
            c = float.(c[2:end, :])
        end
        return c
    else
        throw(ArgumentError("File $file does not exist"))
    end
end

"""
    readfile_csv(file::String)

Read CSV file with proper handling of headers.

# Arguments
- `file::String`: CSV file path to read

# Returns
- `Matrix{Float64}`: Numeric data from CSV file

# Notes
- Uses comma delimiter for CSV files
- Converts string data to float64
- Skips header row if first row contains strings
- Used by readfile for CSV-specific handling
"""
function readfile_csv(file::String)
    c = readdlm(file, ',')
    if typeof(c[1, :]) <: AbstractString
        c = float.(c[2:end, :])
    end
    return c
end

"""
    readfile(file::String, col::Int)

read file and return given column
"""
readfile(file::String, col::Int) = readfile(file)[:, col]

"""
    read_cell_counts(file)

read cell counts from a file
"""
function read_cell_counts(file)
    a, h = readdlm(file, header=true)
    return a, h
end

"""
    validate_dwelltime_compat(dttype, onstates, dwellpath; datatype="rnadwelltime")

Check that dttype, onstates, and the dwell-time file path(s) are compatible.

- `length(onstates) == length(dttype)` (one onstate set per dwell type).
- Either one file per type (`length(dwellpath) == length(dttype)`), or the HJ-style layout:
  `length(dttype)==4`, `length(dwellpath)==2`, and dttype contains exactly two R-step types
  (no \"G\") and two G-step types (contain \"G\").

# Throws
- `ArgumentError` if any check fails.
"""
function validate_dwelltime_compat(dttype::Vector, onstates, dwellpath; datatype="rnadwelltime")
    ntypes = length(dttype)
    nfiles = length(dwellpath)
    if length(onstates) != ntypes
        throw(ArgumentError("dttype and onstates length mismatch: dttype has $ntypes entries, onstates has $(length(onstates)). Each dwell type (e.g. ON, OFF, ONG, OFFG) must have one onstate set."))
    end
    if nfiles == ntypes
        return nothing
    end
    if ntypes == 4 && nfiles == 2
        rstep = count(x -> !occursin("G", x), dttype)
        gstep = count(x -> occursin("G", x), dttype)
        if rstep != 2 || gstep != 2
            throw(ArgumentError("For 4 dttypes and 2 files, dttype must contain exactly two R-step types (e.g. ON, OFF) and two G-step types (e.g. ONG, OFFG). Got $rstep without 'G' and $gstep with 'G'."))
        end
        return nothing
    end
    throw(ArgumentError("dwellpath and dttype are incompatible: $nfiles file(s) for $ntypes dwell type(s). Use either one file per type, or 4 types with 2 files (first file = ON, OFF; second = ONG, OFFG)."))
end

"""
    read_dwelltimes(datapath::String)
    read_dwelltimes(datapath::Vector{String})

Read dwell-time file(s) and return (bins, DT) as vectors.

Column convention: first column is always bins; remaining columns are histograms in order.
- 2 columns: one (bins, DT) pair.
- 3 columns: two pairs (same bins for both; col 2 and col 3 are the two histograms).
- More than 3 columns: throws `ArgumentError`.

# Arguments
- `datapath`: Single file path (string) or vector of file paths.

# Returns
- `bins::Vector{Vector}`: One bin vector per histogram, in order.
- `DT::Vector{Vector}`: One dwell-time histogram per entry, in order.
"""
function read_dwelltimes(datapath::String)
    bins = Vector{Vector}(undef, 0)
    DT = Vector{Vector}(undef, 0)
    read_dwelltimes!(bins, DT, datapath)
    return bins, DT
end

function read_dwelltimes(datapath::Vector{String})
    bins = Vector{Vector}(undef, 0)
    DT = Vector{Vector}(undef, 0)
    for i in eachindex(datapath)
        read_dwelltimes!(bins, DT, datapath[i])
    end
    return bins, DT
end

"""
    read_dwelltimes!(bins, DT, datapath::String)

In-place: read a single dwell-time file and append (bins, DT) onto `bins` and `DT`.

First column is bins; remaining columns are histograms. 2 columns → append one (bins, DT)
pair; 3 columns → append same bins twice and the two histogram columns to DT; more than
3 columns → `ArgumentError`.

# Arguments
- `bins`: Vector to append bin vectors to.
- `DT`: Vector to append dwell-time histogram vectors to.
- `datapath`: Path to a single dwell-time file.

# Returns
- `(bins, DT)` (mutates the input vectors).
"""
function read_dwelltimes!(bins, DT, datapath::String)
    c = readfile(datapath)
    ncol = size(c, 2)
    if ncol == 2
        push!(bins, c[:, 1])
        push!(DT, c[:, 2])
    elseif ncol == 3
        bins_col = c[:, 1]
        push!(bins, bins_col)
        push!(bins, bins_col)
        push!(DT, c[:, 2])
        push!(DT, c[:, 3])
    else
        throw(ArgumentError("Dwell-time file must have 2 or 3 columns (bins + histogram(s)), got $ncol"))
    end
    return bins, DT
end

"""
    read_tracefile(path::String, start::Int, stop::Int, col=3)

Read intensity of a trace file into a vector

Arguments:
- `path`: folder holding trace files
- `start`: start index for trace
- `stop`: stop index for trace, -1 for end
- `col`: column to read


"""
function read_tracefile(path::String, start::Int, stop::Int, col=3)
    t = readfile(path)
    if size(t, 2) > 1 && col == 1
        println("Warning: trace file has multiple columns, but col=1.")
    end
    t = t[:, col]
    if stop < 0 || stop > length(t)
        (start <= length(t)) && return t[start:end]
    else
        (stop <= length(t)) && return t[start:stop]
    end
end

"""
    read_tracefiles(path::String, label::String, start::Int, stop::Int, col=3)

    Args: same as read_tracefile with
    label: string to match in file name


read in a set of trace files 
"""
function read_tracefiles(path::String, label::String, start::Int, stop::Int, col=3, uniquetrace::Bool=false)
    traces = Vector[]
    if isempty(path)
        return traces
    else
        for (root, dirs, files) in walkdir(path)
            for file in files
                if occursin(label, file) && ~occursin(".DS_Store", file)
                    t = read_tracefile(joinpath(root, file), start, stop, col)
                    !isempty(t) && push!(traces, t)
                end
            end
        end
        if uniquetrace
            set = sum.(traces)
            return traces[unique(i -> set[i], eachindex(set))]  # only return unique traces
        else
            return traces
        end
    end
end
"""
    read_tracefiles(path::String, label::Vector{String}, start::Int, stop::Int, col=3)

Args: same as read_tracefile with
label: vector of strings for traces to be fit simultaneously

read in joint trace files
"""
function read_tracefiles(path::String, label::Vector{String}, start::Int, stop::Int, col=3)
    l = length(label)
    traces = Matrix[]
    if isempty(path)
        return traces
    else
        # Get all files from the directory
        all_files = readdir(path)
        # Split files by label (e.g., all files with "enhancer", all files with "gene")
        files_by_label = split_files_by_label(all_files, label)
        # Check that we have files for all labels
        if any(isempty, files_by_label)
            return traces
        end
        # Sort each group to ensure consistent pairing by numeric prefix
        files_by_label = [sort(fs) for fs in files_by_label]
        # Pair files by index (assumes they're in the same order after sorting)
        n_pairs = minimum(length.(files_by_label))
        for i in 1:n_pairs
            tset = Vector{Vector}(undef, l)
            for j in 1:l
                tset[j] = read_tracefile(joinpath(path, files_by_label[j][i]), start, stop, col)
            end
            push!(traces, hcat(tset...))
        end
        return traces
    end
end

function read_tracefiles_unbalanced(path::String, label::Vector{String}, start::Int, stop::Int, col=3; backup_path::String="")
    l = length(label)
    traces = Matrix[]
    if isempty(path)
        return traces
    else
        for (root, dirs, files) in walkdir(path)
            files = split_files_by_label(readdir(path), label)
            nmin = minimum(length.(files))
            nmax = maximum(length.(files))
            imax = argmax(length.(files))
            imin = argmin(length.(files))
            # Pair up to the shortest list
            for i in 1:nmin
                tset = Vector{Vector}(undef, l)
                for j in 1:l
                    tset[j] = read_tracefile(joinpath(root, files[j][i]), start, stop, col)
                end
                push!(traces, hcat(tset...))
            end
            # Handle unpaired files in the longer list (only for l == 2)
            if l == 2 && nmax > nmin && !isempty(backup_path)
                for i in nmin+1:nmax
                    tset = Vector{Vector}(undef, l)
                    tset[imax] = read_tracefile(joinpath(root, files[imax][i]), start, stop, col)
                    # Find a file in backup directory with same label and sufficient length
                    backup_files = split_files_by_label(readdir(backup_path), [label[imin]])
                    if !isempty(backup_files[1])
                        # Find all backup files with length >= the unmatched file length
                        target_length = length(tset[imax])
                        suitable_backups = Vector{Vector{Float64}}()
                        for backup_file in backup_files[1]
                            backup_data = read_tracefile(joinpath(backup_path, backup_file), start, stop, col)
                            if length(backup_data) >= target_length
                                push!(suitable_backups, backup_data[1:target_length])  # Truncate if longer
                            end
                        end
                        if !isempty(suitable_backups)
                            # Randomly select one of the suitable backups
                            selected_backup = suitable_backups[rand(1:length(suitable_backups))]
                            # Add random noise with STD of 1% of the median
                            median_val = median(selected_backup)
                            noise_std = 0.01 * median_val
                            tset[imin] = selected_backup .+ noise_std .* randn(length(selected_backup))
                        else
                            # Fallback to random data if no suitable backup found
                            tset[imin] = randn(target_length)
                        end
                    else
                        # Fallback to random data if no backup files found
                        # Add random noise with STD of 1% of the median of the real trace
                        median_val = median(tset[imax])
                        noise_std = 0.01 * median_val
                        tset[imin] = randn(length(tset[imax])) .* noise_std
                    end
                    push!(traces, hcat(tset...))
                end
            end
            return traces
        end
    end
end

"""
    read_tracefiles(path, label, traceinfo::Tuple, col=3)

#Arguments: same as read_tracefiles with
- `traceinfo`: tuple of (dt, start, stop, ...)

read in a set of trace files from a folder
"""
function read_tracefiles(path, label, traceinfo::Tuple, col=3)
    start = max(round(Int, traceinfo[2] / traceinfo[1]), 1)
    stop = traceinfo[3] < 0 ? -1 : max(round(Int, traceinfo[3] / traceinfo[1]), 1)
    read_tracefiles(path, label, start, stop, col)
end
"""
    read_tracefiles_grid(path::String, label::String, start::Int, stop::Int, col=3)

Reads trace files from a specified directory that match a given label and extracts data from them.

# Arguments
- `path::String`: The directory path to search for trace files.
- `label::String`: The label to match in the filenames.
- `start::Int`: The starting index for reading the trace data.
- `stop::Int`: The stopping index for reading the trace data.
- `col::Int`: The column index to read from each trace file (default is 3).

# Returns
- `Vector`: A vector of unique traces read from the files.
"""

function read_tracefiles_grid(path, label, traceinfo)
    start = max(round(Int, traceinfo[2] / traceinfo[1]), 1)
    stop = traceinfo[3] < 0 ? -1 : max(round(Int, traceinfo[3] / traceinfo[1]), 1)
    read_tracefiles_grid(path, label, start, stop)
end

function read_tracefiles_grid(path::String, label::String, start::Int, stop::Int)
    traces = Matrix[]
    if isempty(path)
        return traces
    else
        for (root, dirs, files) in walkdir(path)
            for file in files
                if occursin(label, file) && ~occursin(".DS_Store", file)
                    t = read_tracefile_grid(joinpath(root, file), start, stop)
                    !isempty(t) && push!(traces, t)
                end
            end
        end
        set = sum.(traces)
        return traces[unique(i -> set[i], eachindex(set))]  # only return unique traces
    end
end

function read_tracefile_grid(path, start, stop, header=true)
    t = readdlm(path, ',', header=header)[1]
    if stop < 0
        (start <= length(t)) && return t[:, start:end]
    else
        (stop <= length(t)) && return t[:, start:stop]
    end
end




# function read_tracefiles_r(path::String, label::Vector{String}, start::Int, stop::Int, col=3)
#     l = length(label)
#     traces = [Vector[] for _ in 1:l]
#     if isempty(path)
#         return traces
#     else
#         for (root, dirs, files) in walkdir(path)
#             tset = Vector{Vector}(undef, l)
#             files = sort(readdir(path))
#             for file in files
#                 complete = true
#                 for i in eachindex(label)
#                     if occursin(label[i], file)
#                         tset[i] = read_tracefile(joinpath(root, file), start, stop, col)
#                     end
#                     complete &= isassigned(tset, i)
#                 end
#                 if complete
#                     for i in eachindex(label)
#                         push!(traces[i], tset[i])
#                     end
#                     tset = Vector{Vector}(undef, length(label))
#                 end
#             end
#         end
#         return traces
#     end
# end

# function read_tracefiles_v(path::String, label::Vector{String}, start::Int, stop::Int, col=3)
#     l = length(label)
#     traces = Vector[]
#     if isempty(path)
#         return traces
#     else
#         for (root, dirs, files) in walkdir(path)
#             tset = Vector{Vector}(undef, l)
#             files = sort(readdir(path))
#             for file in files
#                 complete = true
#                 for i in eachindex(label)
#                     if occursin(label[i], file)
#                         tset[i] = read_tracefile(joinpath(root, file), start, stop, col)
#                     end
#                     complete &= isassigned(tset, i)
#                 end
#                 if complete
#                     push!(traces, tset)
#                     tset = Vector{Vector}(undef, length(label))
#                 end
#             end
#         end
#         return traces
#     end
# end

"""
    fix_tracefiles(path::String)

reorder trace files into standard form
"""
function fix_tracefiles(path::String)
    for (root, dirs, files) in walkdir(path)
        for file in files
            target = joinpath(root, file)
            t = readdlm(target, header=true)
            writedlm(target, [t[1] t[1][:, 2]])
        end
    end
end

function read_allrates_csv(file::String, header::Bool)
    r = readdlm(file, ',', header=header)
    if header
        r = r[1]
    end
    return r
end

"""
    readrates(infolder::String, label::String, gene::String, G::Int, R::Int, S::Int, insertstep::Int, nalleles::Int, ratetype::String="median")

Read in rates from a file.

# Arguments
- `infolder::String`: Folder holding rate files.
- `label::String`: Label for rate files.
- `gene::String`: Gene name.
- `G::Int`: Number of G states.
- `R::Int`: Number of R steps.
- `S::Int`: Number of splice states.
- `insertstep::Int`: R step where splicing starts.
- `nalleles::Int`: Number of alleles.
- `ratetype::String`: Type of rate to read. Defaults to `"median"`.

# Returns
- `Array{Float64, 2}`: The read rates from the file.

# Examples
```julia
rates = readrates("data", "label", "gene", 3, 5, 2, 1, 2, "mean")

"""
function readrates(infolder, label, gene, G, R, S, insertstep, nalleles, ratetype="median")
    if R == 0
        name = filename(label, gene, G, nalleles)
    else
        name = filename(label, gene, G, R, S, insertstep, nalleles)
    end
    readrates(joinpath(infolder, "rates" * name), get_row(ratetype))
end


"""
readrates(file::String,row::Int)
readrates(file::String)

row
1       maximum likelihood
2       mean
3       median
4       last value of previous run
"""
readrates(file::String, row::Int) = readrow(file, row)

readrates(file::String) = readrates(file, 3)



"""
    get_row(ratetype)

"""
function get_row(ratetype)
    if ratetype == "ml"
        row = 1
    elseif ratetype == "mean"
        row = 2
    elseif ratetype == "median"
        row = 3
    elseif ratetype == "last"
        row = 4
    else
        row = 3
    end
    row
end

function readrow(file::String, row, delim=',')
    if isfile(file) && ~isempty(read(file))
        contents = readdlm(file, delim, header=false)
        if ~(typeof(contents[1]) <: Number)
            contents = readdlm(file, delim, header=true)[1]
        end
        if row <= size(contents, 1)
            m = contents[row, :]
            return m[.~isempty.(m)]
        else
            println("Row $row too large for file $file, returning median")
            return contents[3, :]
        end
    else
        println("File $file does not exist")
        return Float64[]
    end
end

function readrow_flip(file, row)
    m = readrow(file, row)
    reshape(m, 1, length(m))
end

function readmeasures(file::String)
    d = readdeviance(file)
    # w = readwaic(file)
    a = readaccept(file)
    t = readtemp(file)
    r = readrhat(file)
    e = readess(file)
    g = readgeweke(file)
    m = readmcse(file)
    ll = readrow(file, 1)
    # Order matches assemble_measures header:
    # LogMaxLikelihood, normalized_LL, n_obs, n_params, Deviance, WAIC, WAIC SE, AIC, Acceptance, Temperature, Rhat, ESS, Geweke, MCSE
    # writedlm(f, [fits.llml normalized_ll Int(n_obs) Int(n_params) mean(fits.ll) std(fits.ll) quantile(fits.ll, [0.025; 0.5; 0.975])' measures.waic[1] measures.waic[2] aic(fits)], ',')

    [ll[1] ll[2] ll[3] ll[4] d[1] ll[10] ll[11] ll[12] a[1] a[2] t[1] r[1] e[1] g[1] m[1]]
end

readdeviance(file::String) = readrow(file, 2)

# readwaic(file::String) = readrow(file, 1)

function readaccept(file::String)
    a = readrow(file, 3)
    a[1] / a[2], a[2]
end

readtemp(file::String) = readrow(file, 4)

readrhat(file::String) = readrow(file, 9)

readess(file::String) = readrow(file, 10)

readgeweke(file::String) = readrow(file, 11)

readmcse(file::String) = readrow(file, 12)

function readml(ratefile::String)
    m = readrow(ratefile, 1)
    reshape(m, 1, length(m))
end

function readmean(ratefile::String)
    m = readrow(ratefile, 2)
    reshape(m, 1, length(m))
end

function readsd(ratefile::String)
    m = readrow(ratefile, 2)
    reshape(m, 1, length(m))
end

function readstats(statfile::String)
    mean = readrow_flip(statfile, 1)
    sd = readrow_flip(statfile, 2)
    median = readrow_flip(statfile, 3)
    mad = readrow_flip(statfile, 4)
    credl = readrow_flip(statfile, 5)
    credh = readrow_flip(statfile, 7)
    [mean sd median mad credl credh]
end

function readmedian(statfile::String)
    m = readrow(statfile, 3)
    reshape(m, 1, length(m))
end

function readmad(statfile::String)
    m = readrow(file, 4)
    reshape(m, 1, length(m))
end

function read_corparam(file::String)
    c = readdlm(file, ',')
    n = length(c[1, :])
    # c[5+n:4+2*n,1:n]
    c[8:7+n, 1:n]
end

function read_covparam(file::String)
    c = readdlm(file, ',')
    read_covparam(c)
end

function read_covparam(c::Matrix)
    n = length(c[1, :])
    # c[5+2*n:4+3*n,1:n]
    c[8+n:7+2*n, 1:n]
end

function read_covlogparam2(file::String)
    c = readdlm(file, ',')
    # n = length(c[1, :])
    nrates = num_fitted_core_params(model)
    c[8+2*nrates+1:7+3*nrates+1, 1:nrates]
end

"""
    read_covlogparam(file::String)

Robustly reads the covariance matrix of the log parameters from a param-stats file.
Finds the last square block of size n x n, where n is the number of parameters.
"""
function read_covlogparam(file::String)
    if !isfile(file)
        println(file, " does not exist")
        return Float64[]
    end
    c = readdlm(file, ',')
    n = length(c[1, :])
    nrates = num_fitted_core_params(model)
    # Search for all n x n blocks in the file
    last_block_start = 0
    for i in 1:(size(c, 1)-nrates+1)
        block = c[i:i+nrates-1, 1:nrates]
        # Check if block is square and symmetric (optional, for extra safety)
        if size(block, 1) == nrates && size(block, 2) == nrates && isapprox(block, block', atol=1e-8)
            last_block_start = i
        end
    end
    if last_block_start == 0
        println("No n x n block found in ", file)
        return Float64[]
    end
    return c[last_block_start:last_block_start+nrates-1, 1:nrates]
end

read_crosscov(statfile::String) = read_crosscov(read_covparam(statfile))

function read_crosscov(C::Matrix)
    c = Float64[]
    N = size(C, 1)
    for i in 1:N
        for j in i+1:N
            push!(c, C[i, j])
        end
    end
    c
end

function read_burst(file::String)
    b = readdlm(file, ',')
    reshape(b[1:4], 1, 4)
end

function read_optimized(file::String)
    rates = readrow_flip(file, 1)
    ll = readrow(file, 2)
    conv = readrow(file, 3)
    [rates ll conv]
end

function read_info(file::String)
    info = readdlm(file, ',')
    info[1, :]
end

"""
    change_suffix(old, new, folder)

TBW
"""
function change_suffix(old, new, folder)
    for (root, dirs, files) in walkdir(folder)
        for file in files
            target = joinpath(root, file)
            if endswith(target, old)
                mv(target, replace(target, old => new))
            end
        end
    end
end

"""
    change_pattern(old, new, folder)

TBW
"""
function change_pattern(old, new, folder)
    for (root, dirs, files) in walkdir(folder)
        for file in files
            target = joinpath(root, file)
            if occursin(old, target)
                mv(target, replace(target, old => new))
                println(target)
            end
        end
    end
end
# Functions for saving and loading data and models

# """
# write_log(file,datafile,data,model)
# write all information necessary for rerunning
# """
# function save_data(file::String,data::TransientRNAData)
#     f = open(file,"w")
#     writedlm(f, [typeof(data)])
#     writedlm(f,[data.label])
#     writedlm(f,[data.gene])
#     writedlm(f,[data.nRNA])
#     writedlm(f,[data.time])
#     writedlm(f,[data.histRNA])
#     close(f)
# end

# function load_data(file::String, model::AbstractGMmodel)


# end

function save_model(file::String, model::AbstractGMmodel)
    f = open(file, "w")
    write(f, model.G)
    writedlm(f, model.nalleles)
    writedlm(f, model.ratepriors)
    writedlm(f, model.proposal)
    writedlm(f, model.fittedparam)
    writedlm(f, model.method)
    close(f)
end

# function load_model(file::String, model::AbstractGRSMmodel)

# end

function read_bottom_float_block(file::String)
    c = readdlm(file, ',')
    nrows, ncols = size(c)
    isnumlike(x) = (x isa Number) || (x isa AbstractString && occursin(r"^-?\d*\.?\d+([eE][\+\-]?\d+)?$", x)) || (x isa AbstractString && isempty(x))
    # Find the last row that contains only numbers or numeric strings
    last_row = nrows
    while last_row > 0
        row = c[last_row, :]
        if all(isnumlike, row)
            break
        end
        last_row -= 1
    end
    # Now, scan upwards to find the first row of the block
    first_row = last_row
    while first_row > 1
        row = c[first_row-1, :]
        if all(isnumlike, row)
            first_row -= 1
        else
            break
        end
    end
    block = c[first_row:last_row, :]
    # Convert all entries to Float64, set empty to NaN
    block_float = map(x -> (x isa Number) ? float(x) : (x isa AbstractString && !isempty(x) ? parse(Float64, x) : NaN), block)
    # Find the bottom row, count the number of valid floats (not NaN)
    last_valid_row = block_float[end, :]
    k = count(!isnan, last_valid_row)
    if k == 0 || k > size(block_float, 1) || k > size(block_float, 2)
        return Array{Float64}(undef, 0, 0)
    end
    # Return the bottom k rows and first k columns as a k x k square
    return block_float[end-k+1:end, 1:k]
end

"""
    write_track_files(df::DataFrame, col_pattern::String, base_name::String, condition::Function, output_dir::String=".")

Write data columns from a DataFrame into separate track files.

# Arguments
- `df`: DataFrame containing the data columns
- `col_pattern`: Pattern to match column names (e.g., "data" will match data1, data2, etc.)
- `base_name`: Base name for the output files (e.g., "track" will create "001_track.trk", "002_track.trk", etc.)
- `condition`: Function to apply to data before writing to file, e.g. x -> x > 0
- `output_dir`: Directory to write the files to (defaults to current directory)

# Returns
- Vector of created filenames

# Example
```julia
df = DataFrame(...)  # Your DataFrame with data1, data2, etc. columns
filenames = write_track_files(df, "data", "track", x -> x > 0, "output")
```
"""
function write_track_files(df::DataFrame, col_pattern::String, base_name::String, condition::Function, output_dir::String=".")
    # Create output directory if it doesn't exist
    mkpath(output_dir)

    # Find all columns matching the pattern
    data_cols = filter(name -> startswith(String(name), col_pattern), names(df))

    # Store created filenames
    created_files = String[]

    # Process each data column
    for (i, col) in enumerate(data_cols)
        # Get the data, filtering out missing values
        data = filter(!ismissing, df[!, col])

        # Skip if no valid data
        isempty(data) && continue

        # Create filename with padded number
        filename = joinpath(output_dir, string(lpad(i, 3, "0"), "_", base_name, ".trk"))
        push!(created_files, filename)

        # Write data to file using writedlm, applying the condition function
        writedlm(filename, condition.(data))
    end

    return created_files
end

# Helper function to split files into groups by label
function split_files_by_label(files::Vector{String}, labels::Vector{String})
    [filter(f -> occursin(label, f), files) for label in labels]
end

"""
    write_traces_to_files(traces::Vector{Matrix}, labels::Vector{String}, output_dir::String=".")

Write trace data from matrices to individual .trk files.

# Arguments
- `traces::Vector{Matrix}`: Vector of matrices from read_tracefiles_unbalanced
- `labels::Vector{String}`: Labels corresponding to each column
- `output_dir::String`: Directory to write files to (defaults to current directory)

# Returns
- Vector of created filenames

# Example
```julia
traces = read_tracefiles_unbalanced(path, labels, start, stop, col, backup_path=backup_path)
filenames = write_traces_to_files(traces, labels, "output_folder")
```
"""
function write_traces_to_files(traces::Vector{Matrix}, labels::Vector{String}, output_dir::String=".")
    # Create output directory if it doesn't exist
    mkpath(output_dir)

    # Store created filenames
    created_files = String[]

    # Process each matrix in traces
    for (trace_idx, trace_matrix) in enumerate(traces)
        # Process each column (label) in the matrix
        for (col_idx, label) in enumerate(labels)
            if col_idx <= size(trace_matrix, 2)  # Ensure column exists
                # Create filename with padded number and label
                filename = joinpath(output_dir, string(lpad(trace_idx, 3, "0"), "_", label, ".trk"))
                push!(created_files, filename)

                # Write the column data to file
                writedlm(filename, trace_matrix[:, col_idx])
            end
        end
    end

    return created_files
end

"""
    remove_string(str, st1)

Remove a substring from a string.

# Arguments
- `str`: Original string
- `st1`: Substring to remove

# Returns
- `String`: String with st1 removed

# Notes
- Uses replace function to remove all occurrences of st1
"""
remove_string(str, st1) = replace(str, st1 => "")

"""
    remove_string(str, str1, str2)

Remove two substrings from a string.

# Arguments
- `str`: Original string
- `str1`: First substring to remove
- `str2`: Second substring to remove

# Returns
- `String`: String with both str1 and str2 removed

# Notes
- Removes str1 first, then str2 from the result
- Uses chained replace operations
"""
remove_string(str, str1, str2) = replace(remove_string(str, str1), str2 => "")


