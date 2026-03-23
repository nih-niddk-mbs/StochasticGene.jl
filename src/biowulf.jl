# This file is part of StochasticGene.jl  

# biowulf.jl
# functions for use on the NIH Biowulf super computer


"""
struct ModelArgs

For passing model information to fit function from swarmfile (limited to numbers and strings (i.e. no vectors or tuples))

#Fields
`inlabel::String`
`label::String`
`G::Int`
`R::Int`
`S::Int`
`insertstep::Int`
`TransitionType::String`: type of model, e.g. "nstate", "KP", "cyclic"
`fixedeffects::String`: two numbers separated by a hyphen, e.g. "3-4", indicating parameters 3 and 4 are fixed to each other

"""

struct ModelArgs
    inlabel::String
    label::String
    G::Int
    R::Int
    S::Int
    insertstep::Int
    TransitionType::String
    fixedeffects::String
end

"""
    sanitize_for_filename(s::AbstractString)

Replace characters that are unsafe in shell filenames (e.g. `,` and `|`) with `-`.
Underscore is reserved for preset filename fields (see io.jl).
Use when building swarm/.jl filenames from label or coupling spec (e.g. "24,35|35").
"""
sanitize_for_filename(s::AbstractString) = replace(replace(string(s), "," => "-"), "|" => "-")

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
- Overrides (optional): `resultfolder`, `root`, `maxtime`, `samplesteps`, etc. are written into each
  script so `fit(; key=key, resultfolder=..., ...)` uses them. Only simple types (String, Number, Bool) are serialized.

# Examples
```julia
makeswarm(["33il", "44il"]; filedir="swarm", resultfolder="HCT116_test", root=".")
makeswarm(; key="33il", resultfolder="HCT116_test", maxtime=120.0)
```
"""
function makeswarm(keys::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript", filedir=".", project="", sysimage="", src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, kwargs...)
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    isempty(keys) && return
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    write_swarmfile_keys(joinpath(filedir, sfile), nchains, nthreads, juliafile, keys, project, sysimage)
    for k in keys
        scriptpath = joinpath(filedir, juliafile * "_" * sanitize_for_filename(k) * ".jl")
        write_fitfile_key(scriptpath, k; src=src, resultfolder=resultfolder, root=root, maxtime=maxtime, samplesteps=samplesteps, kwargs...)
    end
end

function makeswarm(; key::String, kwargs...)
    isempty(key) && error("makeswarm(; key=...) requires a non-empty key")
    makeswarm([key]; kwargs...)
end

"""
    makeswarm_genes(genes::Vector{String}; <keyword arguments>)

Write a swarm file and one shared fit script to run each gene in `genes`. Each swarm line runs the same
script with the gene as argument: `julia -t nthreads -p nchains fitscript_<label>_<model>.jl gene`.
The script calls `fit(...)` with `ARGS[1]` as the gene.

# Arguments
- `genes`: vector of gene names
- `batchsize=1000`: number of jobs per swarm file when `genes` is large
- `filedir="."`: directory to write swarm and script files
- Plus the same fit/model kwargs (datatype, dttype, datapath, cell, datacond,
  transitions, G, R, S, insertstep, etc.) and swarm options (nchains, nthreads, swarmfile, juliafile, project, src).

# Example
```julia
makeswarm_genes(["MYC", "SOX9"]; cell="HBEC", datatype="rna", datapath="data/", resultfolder="out")
```
"""
function makeswarm_genes(genes::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="rna", dttype=String[], datapath="", cell::String="HCT116", datacond="MOCK", traceinfo=(1.0, 1.0, -1, 1.0), infolder::String="HCT116_test", resultfolder::String="HCT116_test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, coupling=tuple(), TransitionType="", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian=false, datacol=3, ejectnumber=1, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[], filedir=".", project="", sysimage="")
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    modelstring = create_modelstring(G, R, S, insertstep)
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, TransitionType)
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
    write_fitfile_genes(joinpath(filedir, juliafile_full), nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)
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
    return nothing
end

"""
    write_fitfile_key(fitfile, key::String; src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, ...)

Writes a script that runs `fit(; key=key, ...)` with optional overrides. Only String, Real, and Bool overrides are written.
"""
function write_fitfile_key(fitfile, key::String; src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, kwargs...)
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
    for (k, v) in pairs(kwargs)
        s = _format_fit_override(Symbol(k), v)
        s !== nothing && (overrides[Symbol(k)] = s)
    end
    fit_args = ["key=$(repr(key))", values(overrides)...]
    write(f, "@time fit(; $(join(fit_args, ", ")))\n")
    close(f)
end

"""
    write_fitfile_genes(fitfile, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)

Writes a fit script that takes the gene as ARGS[1] and calls `fit(nchains, datatype, ..., ARGS[1], cell, ...)`.
"""
function write_fitfile_genes(fitfile, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor=1.0, trace_specs=[], dwell_specs=[])
    s = '"'
    transitions = transitions isa AbstractVector && !(transitions isa Tuple) ? Tuple(transitions) : transitions
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, ARGS[1], $s$cell$s, $datacond, $traceinfo, $s$infolder$s, $s$resultfolder$s, $s$inlabel$s, $s$label$s, $fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber, $yieldfactor, $trace_specs, $dwell_specs)")
    close(f)
end

"""
    write_fitfile_coupled(fitfile, gene::String, nchains, ...)

Like `write_fitfile_genes` but with `gene` hardcoded in the script instead of `ARGS[1]`.
Used for coupled-model scripts where one gene is fixed per script.
`trace_specs` and `dwell_specs` are interpolated directly — when non-empty they override
the legacy `traceinfo`/`zeromedian` inside `fit`.
"""
function write_fitfile_coupled(fitfile, gene::String, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber, yieldfactor=1.0, trace_specs=[], dwell_specs=[])
    s = '"'
    transitions = transitions isa AbstractVector && !(transitions isa Tuple) ? Tuple(transitions) : transitions
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, $s$gene$s, $s$cell$s, $datacond, $traceinfo, $s$infolder$s, $s$resultfolder$s, $s$inlabel$s, $s$label$s, $fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber, $yieldfactor, $trace_specs, $dwell_specs)\n")
    close(f)
end

"""
    makeswarm_coupled(; gene, inlabel, label, nchains, ..., trace_specs=[], dwell_specs=[], ...)

Write one fit script for a coupled model (gene hardcoded, not ARGS[1]) and append one
swarm line. Intended for coupled models where each label/coupling gets its own script.
`trace_specs` is written directly into the script so `fit` uses it instead of the legacy
`traceinfo`/`zeromedian` path for determining interval and per-unit zero-centering.
"""
function makeswarm_coupled(; gene::String, inlabel::String, label::String,
    nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript",
    filedir=".", project="", sysimage="", src="",
    datatype::String="tracejoint", dttype=String[], datapath="", cell::String="HBEC",
    datacond=[], traceinfo=(1.0, 1.0, -1.0), infolder::String="", resultfolder::String="",
    fittedparam=Int[], fixedeffects=tuple(), transitions=tuple(), G=2, R=0, S=0, insertstep=1,
    coupling=tuple(), grid=nothing, root=".", maxtime=60.0, elongationtime=6.0,
    priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0,
    splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(),
    ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0,
    annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false,
    optimize=false, writesamples=false, method="Tsit5()", zeromedian=false,
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
    write_fitfile_coupled(scriptpath, gene, nchains, datatype, dttype, datapath, cell,
        datacond, traceinfo, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects,
        transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime,
        priormean, nalleles, priorcv, onstates, decayrate, splicetype, probfn, noisepriors,
        hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp,
        tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian,
        datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)
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
    create_label(label, inlabel, datatype, datacond, cell, TransitionType)

Create a label string based on the provided parameters.

# Arguments
- `label`: Initial label string (can be empty).
- `inlabel`: Initial inlabel string (can be empty).
- `datatype`: Type of data.
- `datacond`: Condition of the data. Can be a `String` or a `Vector` of strings.
- `cell`: Cell type.
- `TransitionType`: Type of transition.

# Returns
- `label`: Generated label string.
- `inlabel`: Generated inlabel string.

# Methods
- `create_label(label, inlabel, datatype, datacond::String, cell, TransitionType)`: Handles `datacond` as a `String`.
- `create_label(label, inlabel, datatype, datacond::Vector, cell, TransitionType)`: Handles `datacond` as a `Vector` of strings.
"""
function create_label(label, inlabel, datatype, datacond::String, cell, TransitionType)
    if isempty(label)
        label = datatype * "-" * cell
        ~isempty(TransitionType) && (label = label * "-" * TransitionType)
        typeof(datacond) <: AbstractString && (label = label * "_" * datacond)
    end
    isempty(inlabel) && (inlabel = label)
    return label, inlabel
end

function create_label(label, inlabel, datatype, datacond::Vector, cell, TransitionType)
    create_label(label, inlabel, datatype, join(datacond, "-"), cell, TransitionType)
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
"""
function folder_path(folder::String, root::String, folderatetype::String=""; make=false)
    f = folder
    if ~ispath(folder) && ~isempty(folder)
        f = joinpath(root, folder)
        if ~ispath(f)
            f = joinpath(root, folderatetype, folder)
            if ~ispath(f) && ~make
                println("$folder not found")
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
