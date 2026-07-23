#!/usr/bin/env julia

# Static documentation audit for maintainer hygiene.
# This intentionally avoids `using StochasticGene` so it can run before package
# precompilation or dependency setup is healthy.

const ROOT = normpath(joinpath(@__DIR__, ".."))
const STRICT_DOCSTRINGS = "--strict-docstrings" in ARGS

relpath(parts...) = joinpath(ROOT, parts...)
readrepo(parts...) = read(relpath(parts...), String)

function collect_files(rootdir; exts=(".md", ".jl"))
    files = String[]
    isdir(rootdir) || return files
    for (dir, _, names) in walkdir(rootdir)
        occursin(joinpath("docs", "build"), dir) && continue
        for name in names
            any(endswith(name, ext) for ext in exts) || continue
            push!(files, joinpath(dir, name))
        end
    end
    sort!(files)
    return files
end

function strip_comment(line::AbstractString)
    i = findfirst('#', line)
    i === nothing && return String(line)
    return String(line[begin:prevind(line, i)])
end

function exported_symbols()
    text = readrepo("src", "StochasticGene.jl")
    symbols = Symbol[]
    in_export = false
    for raw in split(text, '\n')
        line = strip_comment(raw)
        stripped = strip(line)
        if startswith(stripped, "export")
            in_export = true
            stripped = strip(replace(stripped, r"^export\b" => ""))
        elseif in_export
            if isempty(stripped)
                in_export = false
                continue
            end
            if !startswith(raw, r"\s") && !startswith(stripped, ",")
                in_export = false
                continue
            end
        else
            continue
        end
        for m in eachmatch(r"[A-Za-z_][A-Za-z0-9_!]*", stripped)
            push!(symbols, Symbol(m.match))
        end
    end
    return unique(symbols)
end

function docstring_candidates(symbols)
    src = join(read.(collect_files(relpath("src"); exts=(".jl",)), String), "\n")
    documented = Set{Symbol}()
    for sym in symbols
        name = String(sym)
        pattern = Regex("\"\"\"[\\s\\S]{0,2500}\\b" * escape_string(name) * "\\b[\\s\\S]{0,2500}\"\"\"")
        occursin(pattern, src) && push!(documented, sym)
    end
    return documented
end

const DOCSTRING_ALLOWLIST = Set(Symbol[
    :CSV, :DataFrame, :Tsit5, :lsoda, :mean, :norm, :sparse,
    :INFERENCE_ADVI, :INFERENCE_CHOICES, :INFERENCE_MH, :INFERENCE_NUTS,
    :HMM_STACK_MH, :HMM_STACK_AD, :DEFAULT_CORRELATION_ALGORITHM,
    :COUPLING_MODE_RECIPROCAL_DEFAULT,
])

const STALE_DOC_PATTERNS = [
    r"\bsave_results\(",
    r"\bsetup_parallel\(",
    r"\bfit_hierarchical\(",
    r"\bcheck_workers\(",
    r"\bcleanup_parallel\(",
    r"\bplot_average_trace\(",
    r"\bplot_traces\(",
    r"\banalyze_bursts\(",
    r"\bplot_bursts\(",
    r"\banalyze_time_series\(",
    r"\bplot_time_series\(",
    r"\bcalculate_autocorrelation\(",
    r"\bplot_autocorrelation\(",
    r"\bcalculate_cross_correlation\(",
    r"\bplot_cross_correlation\(",
    r"\banalyze_transitions\(",
    r"\bplot_transition\(",
    r"\bcompare_models\(",
    r"\bplot_model_comparison\(",
    r"\bcalculate_correlation\(",
    r"\bcalculate_correlations\(",
    r"\bplot_correlation\(",
    r"\bplot_correlations\(",
    r"\bplot_reporter\(",
    r"fit_\*\.swarm",
    r"for f in fit_",
    r"batch command files",
    r"batchsize\s*=\s*(1000|4800)",
]

function stale_doc_hits()
    files = vcat(
        [relpath("README.md")],
        collect_files(relpath("docs", "src"); exts=(".md",)),
    )
    hits = Tuple{String,Int,String,String}[]
    for file in files
        basename(file) == "documentation_audit.md" && continue
        lines = split(read(file, String), '\n')
        for (i, line) in pairs(lines)
            for pat in STALE_DOC_PATTERNS
                if occursin(pat, line)
                    push!(hits, (relpath(ROOT, file), i, string(pat), String(strip(line))))
                end
            end
        end
    end
    return hits
end

function main()
    problems = String[]
    exports_path = relpath("src", "exports.jl")
    isfile(exports_path) && push!(problems, "retired src/exports.jl exists")

    symbols = exported_symbols()
    documented = docstring_candidates(symbols)
    missing_docs = sort!(setdiff(symbols, union(documented, DOCSTRING_ALLOWLIST)); by=String)

    hits = stale_doc_hits()
    if !isempty(hits)
        push!(problems, "stale documentation phrases found")
    end
    if STRICT_DOCSTRINGS && !isempty(missing_docs)
        push!(problems, "exported symbols without static docstring candidates")
    end

    println("Documentation audit")
    println("  exported symbols: ", length(symbols))
    println("  static docstring candidates: ", length(documented))
    println("  missing docstring candidates: ", length(missing_docs))
    if !isempty(missing_docs)
        shown = first(missing_docs, min(length(missing_docs), 30))
        println("  first missing candidates: ", join(String.(shown), ", "))
        length(missing_docs) > length(shown) && println("  ... ", length(missing_docs) - length(shown), " more")
    end

    if isempty(hits)
        println("  stale phrase check: ok")
    else
        println("  stale phrase check: ", length(hits), " hit(s)")
        for (file, line, pat, text) in hits
            println("    ", file, ":", line, " matches ", pat, " :: ", text)
        end
    end

    println("  retired exports file: ", isfile(exports_path) ? "present" : "absent")

    if isempty(problems)
        println("audit ok")
        return 0
    end
    println("audit failed: ", join(problems, "; "))
    return 1
end

exit(main())
