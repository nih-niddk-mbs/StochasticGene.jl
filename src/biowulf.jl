# This file is part of StochasticGene.jl  

# biowulf.jl
# functions for use on the NIH Biowulf super computer

"""
    sanitize_for_filename(s::AbstractString)

Replace characters that are unsafe in shell filenames (e.g. `,` and `|`) with `-`.
Underscore is reserved for preset filename fields (see io.jl).
Use when building swarm/.jl filenames from label or coupling spec (e.g. "24,35|35").
"""
sanitize_for_filename(s::AbstractString) = replace(replace(string(s), "," => "-"), "|" => "-")

function _scheduler_symbol(scheduler)
    s = Symbol(lowercase(string(scheduler)))
    s in (:swarm, :biowulf) && return :swarm
    s in (:slurm, :sbatch) && return :slurm
    s in (:parallel, :gnu_parallel, :gnuparallel) && return :parallel
    s in (:command, :commands, :bash, :serial, :none) && return :command
    throw(ArgumentError("unknown scheduler=$(repr(scheduler)); use :swarm, :slurm, :parallel, or :command"))
end

function _scheduler_script_base(commandfile::AbstractString, scheduler::Symbol)
    stem = splitext(basename(String(commandfile)))[1]
    stem = isempty(stem) ? "fit" : stem
    return "$(stem)_$(scheduler).sh"
end

"""
    write_scheduler_launcher(commandfile; scheduler=:swarm, jobs=16, cpus_per_task=1,
        mem="4G", time="02:00:00", jobname="sg-fit", launcherfile=nothing)

Optionally write a shell launcher next to a StochasticGene command file.

- `scheduler=:swarm` or `:command`: write no extra launcher; submit/run the command file yourself.
- `scheduler=:slurm`: write a Slurm array script, default `<commandfile stem>_slurm.sh`.
  The script contains `#SBATCH --array=1-N%jobs`, where `N` is the number of command lines.
- `scheduler=:parallel`: write a GNU Parallel script, default `<commandfile stem>_parallel.sh`.

The command file itself is unchanged: one `julia ... fitscript.jl` command per line.
"""
function write_scheduler_launcher(
    commandfile::AbstractString;
    scheduler=:swarm,
    jobs::Int=16,
    cpus_per_task::Int=1,
    mem::AbstractString="4G",
    time::AbstractString="02:00:00",
    jobname::AbstractString="sg-fit",
    launcherfile::Union{Nothing,AbstractString}=nothing,
)
    sched = _scheduler_symbol(scheduler)
    sched in (:swarm, :command) && return nothing
    cpath = String(commandfile)
    dir = dirname(cpath)
    cbase = basename(cpath)
    lfile = launcherfile === nothing ? _scheduler_script_base(cpath, sched) : String(launcherfile)
    lpath = joinpath(dir, lfile)
    ncommands = isfile(cpath) ? count(line -> !isempty(strip(line)), readlines(cpath)) : 0
    ncommands > 0 || throw(ArgumentError("cannot write scheduler launcher for empty command file: $(repr(cpath))"))
    open(lpath, "w") do io
        if sched == :slurm
            array_spec = jobs > 0 ? "1-$(ncommands)%$(jobs)" : "1-$(ncommands)"
            write(io, """
#!/bin/bash
#SBATCH --job-name=$(jobname)
#SBATCH --array=$(array_spec)
#SBATCH --cpus-per-task=$(cpus_per_task)
#SBATCH --mem=$(mem)
#SBATCH --time=$(time)
#SBATCH --output=logs/slurm_%A_%a.out
#SBATCH --error=logs/slurm_%A_%a.err

set -euo pipefail
mkdir -p logs

COMMAND_FILE="$(cbase)"
CMD=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" "\${COMMAND_FILE}")
echo "Running task \${SLURM_ARRAY_TASK_ID}: \${CMD}"
bash -lc "\${CMD}"
""")
        elseif sched == :parallel
            write(io, """
#!/usr/bin/env bash
set -euo pipefail
mkdir -p logs

COMMAND_FILE="$(cbase)"
JOBS="\${JOBS:-$(jobs)}"
parallel -j "\${JOBS}" --joblog "logs/$(splitext(cbase)[1]).joblog" --results "logs/$(splitext(cbase)[1])" < "\${COMMAND_FILE}"
""")
        end
    end
    try
        chmod(lpath, 0o755)
    catch
    end
    return lpath
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
- `scheduler=:swarm`: launcher style. `:swarm` preserves the traditional Biowulf command
  file only; `:slurm` also writes `<swarmfile>_slurm.sh`; `:parallel` also writes
  `<swarmfile>_parallel.sh`; `:command` writes only the command list.
- `scheduler_jobs=16`: max Slurm array jobs at once (`%scheduler_jobs`) or GNU Parallel jobs.
- `slurm_mem="4G"`, `slurm_time="02:00:00"`, `slurm_jobname="sg-fit"`:
  Slurm wrapper settings.
- Overrides (optional): `resultfolder`, `root`, `maxtime`, `samplesteps`, etc. are written into each
  script so `fit(; key=key, resultfolder=..., ...)` uses them. Only simple types (String, Number, Bool) are serialized.

# Examples
```julia
makeswarm(["33il", "44il"]; filedir="swarm", resultfolder="HCT116_test", root=".")
makeswarm(; key="33il", resultfolder="HCT116_test", maxtime=120.0)
makeswarm(["33il", "44il"]; filedir="slurm-jobs", scheduler=:slurm, scheduler_jobs=50)
```
"""
function makeswarm(keys::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript", filedir=".", project="", sysimage="", src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, scheduler=:swarm, scheduler_jobs::Int=16, slurm_mem::String="4G", slurm_time::String="02:00:00", slurm_jobname::String="sg-fit", kwargs...)
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    isempty(keys) && return
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    commandfile = joinpath(filedir, sfile)
    write_swarmfile_keys(commandfile, nchains, nthreads, juliafile, keys, project, sysimage)
    write_scheduler_launcher(commandfile; scheduler=scheduler, jobs=scheduler_jobs, cpus_per_task=max(1, nchains * Int(nthreads)), mem=slurm_mem, time=slurm_time, jobname=slurm_jobname)
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
    model_grid(; Gset=[2], Rset=[0], Sset=[0], insertset=[1])

Cartesian product of model dimensions for batch grids. Each element is a named tuple
`(G=..., R=..., S=..., insertstep=...)` suitable for `fit`/`makeswarm_genes` loops.
"""
function model_grid(; Gset=[2], Rset=[0], Sset=[0], insertset=[1])
    return [
        (; G=g, R=r, S=s, insertstep=ins)
        for g in Gset, r in Rset, s in Sset, ins in insertset
    ]
end

"""
    makeswarm_folder(resultfolder; root=".", filedir=".", nchains=2, nthreads=1, swarmfile="fit", juliafile="fitscript", ...)

Scan `joinpath(root, \"results\", resultfolder)` for `info_<key>.jld2`, collect each `<key>`, then write
[`write_fitfile_key`](@ref) scripts and a [`write_swarmfile_keys`](@ref) swarm under `filedir`.
Returns the discovered keys (sorted). Extra keywords are forwarded to `write_fitfile_key`.
Pass `scheduler=:slurm` or `scheduler=:parallel` to also write a portable launcher for the
generated command file.
"""
function makeswarm_folder(
    resultfolder::AbstractString;
    root::AbstractString=".",
    filedir::AbstractString=".",
    nchains::Int=2,
    nthreads::Int=1,
    swarmfile::String="fit",
    juliafile::String="fitscript",
    project::String="",
    sysimage::String="",
    src::String="",
    scheduler=:swarm,
    scheduler_jobs::Int=16,
    slurm_mem::String="4G",
    slurm_time::String="02:00:00",
    slurm_jobname::String="sg-fit",
    kwargs...,
)
    rdir = joinpath(root, "results", String(resultfolder))
    isdir(rdir) || throw(ArgumentError("results folder not found: $(repr(rdir))"))
    rx = r"^info_(.+)\.jld2$"
    keys = String[]
    for fn in readdir(rdir)
        m = match(rx, fn)
        m === nothing && continue
        push!(keys, String(m.captures[1]))
    end
    sort!(keys)
    isempty(keys) && return keys
    !isempty(filedir) && !isdir(filedir) && mkpath(filedir)
    sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
    commandfile = joinpath(filedir, sfile)
    write_swarmfile_keys(commandfile, nchains, nthreads, juliafile, keys, project, sysimage)
    write_scheduler_launcher(commandfile; scheduler=scheduler, jobs=scheduler_jobs, cpus_per_task=max(1, nchains * Int(nthreads)), mem=slurm_mem, time=slurm_time, jobname=slurm_jobname)
    for k in keys
        scriptpath = joinpath(filedir, juliafile * "_" * sanitize_for_filename(k) * ".jl")
        write_fitfile_key(scriptpath, k; src=src, kwargs...)
    end
    return keys
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
- Plus the same fit/model kwargs (datatype, `datapath`, cell, `datacond`, optional legacy `dttype` / `traceinfo`,
  `trace_specs`, `dwell_specs`, transitions, G, R, S, insertstep, etc.) and swarm options (nchains, nthreads, swarmfile, juliafile, project, src).
- `scheduler=:swarm`: output style. `:swarm` writes the command file only; `:slurm`
  also writes a submit-ready Slurm array script; `:parallel` also writes a GNU Parallel wrapper;
  `:command` writes only the command list.
- `scheduler_jobs=16`: max Slurm array jobs at once or GNU Parallel jobs.

Legacy separate input-folder routing (`infolder`) is retired; use `datapath`, `resultfolder`, `label`, and `root` like [`fit`](@ref).

# Example
```julia
makeswarm_genes(["MYC", "SOX9"]; cell="HBEC", datatype="rna", datapath="data/", resultfolder="out")
makeswarm_genes(; cell="HCT116", datatype="rna", datapath="data/HCT116_testdata", datacond="MOCK")
makeswarm_genes(; datapath="data/HCT116_testdata", datacond="MOCK", scheduler=:slurm, scheduler_jobs=100)
```
"""
function _genes_from_rna_datafolder(datapath::String, datacond; root::String=".")
    isempty(strip(datapath)) && throw(ArgumentError("datapath must be set when genes are inferred from a folder"))
    datadir = folder_path(datapath, root, "data")
    isdir(datadir) || throw(ArgumentError("data folder not found: $(repr(datadir))"))
    conds = datacond isa AbstractVector ? String.(datacond) : [String(datacond)]
    genesets = Vector{String}[]
    for cond in conds
        push!(genesets, unique(String.(get_genes(cond, datadir))))
    end
    isempty(genesets) && return String[]
    genes = genesets[1]
    for g in genesets[2:end]
        genes = intersect(genes, g)
    end
    sort!(unique!(genes))
    return genes
end

function makeswarm_genes(; datapath::String="", datacond="MOCK", root=".", cell::String="HCT116", kwargs...)
    genes = _genes_from_rna_datafolder(datapath, datacond; root=String(root))
    isempty(genes) && throw(ArgumentError("no genes found in datapath=$(repr(datapath)) matching datacond=$(repr(datacond))"))
    makeswarm_genes(genes; datapath=datapath, datacond=datacond, root=root, cell=cell, kwargs...)
end

function makeswarm_genes(genes::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="rna", dttype=String[], datapath="", cell::String="HCT116", datacond="MOCK", traceinfo=nothing, resultfolder::String="HCT116_test", label::String="",
    fittedparam::Vector=Int[], fixedeffects=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, coupling=tuple(), TransitionType="", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime=60.0, samplesteps::Int=1000000, warmupsteps=0, temp=1.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian=false, datacol=3, ejectnumber=1, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[], filedir=".", project="", sysimage="", scheduler=:swarm, scheduler_jobs::Int=16, slurm_mem::String="4G", slurm_time::String="02:00:00", slurm_jobname::String="sg-fit",
    inference_method=nothing, sample_stride=nothing, parallel=nothing, device=nothing, merge_max_memory=nothing, merge_max_gb=nothing, likelihood_executor=nothing, gradient_checkpoint_length=nothing, proposal_cv_rates=nothing, proposal_cv_noise=nothing, init_jitter=nothing, init_jitter_individual=nothing, init_jitter_noise=nothing)
    if !isempty(filedir) && !isdir(filedir)
        mkpath(filedir)
    end
    modelstring = create_modelstring(G, R, S, insertstep)
    label = create_label(label, datatype_label(datatype), datacond, cell; transition_type=TransitionType)
    label_safe = sanitize_for_filename(label)
    model_safe = sanitize_for_filename(modelstring)
    ngenes = length(genes)
    juliafile_full = juliafile * "_" * label_safe * "_" * model_safe * ".jl"
    commandfiles = String[]
    if ngenes > batchsize
        batches = getbatches(genes, ngenes, batchsize)
        for batch in eachindex(batches)
            sfile = (endswith(swarmfile, ".swarm") ? swarmfile[1:end-6] : swarmfile) * "_" * label_safe * "_" * model_safe * "_" * "$batch" * ".swarm"
            commandfile = joinpath(filedir, sfile)
            push!(commandfiles, commandfile)
            write_swarmfile(commandfile, nchains, nthreads, juliafile_full, batches[batch], project, sysimage)
            write_scheduler_launcher(commandfile; scheduler=scheduler, jobs=scheduler_jobs, cpus_per_task=max(1, nchains * Int(nthreads)), mem=slurm_mem, time=slurm_time, jobname=slurm_jobname)
        end
    else
        sfile = endswith(swarmfile, ".swarm") ? swarmfile : swarmfile * ".swarm"
        commandfile = joinpath(filedir, sfile)
        push!(commandfiles, commandfile)
        write_swarmfile(commandfile, nchains, nthreads, juliafile_full, genes, project, sysimage)
        write_scheduler_launcher(commandfile; scheduler=scheduler, jobs=scheduler_jobs, cpus_per_task=max(1, nchains * Int(nthreads)), mem=slurm_mem, time=slurm_time, jobname=slurm_jobname)
    end
    write_fitfile_genes(joinpath(filedir, juliafile_full), nchains, datatype, datapath, cell, datacond, resultfolder, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs; src=src, dttype=dttype, traceinfo=traceinfo,
        inference_method=inference_method, sample_stride=sample_stride, parallel=parallel, device=device, merge_max_memory=merge_max_memory, merge_max_gb=merge_max_gb, likelihood_executor=likelihood_executor, gradient_checkpoint_length=gradient_checkpoint_length, proposal_cv_rates=proposal_cv_rates, proposal_cv_noise=proposal_cv_noise, init_jitter=init_jitter, init_jitter_individual=init_jitter_individual, init_jitter_noise=init_jitter_noise)
    return (genes=genes, fitfile=joinpath(filedir, juliafile_full), commandfiles=commandfiles)
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
    v === prob_Gaussian && return "$k=prob_Gaussian"
    if occursin("Tsit5", string(typeof(v)))
        return "$k=StochasticGene.Tsit5()"
    end
    v isa AbstractString && return "$k=\"$(escape_string(String(v)))\""
    v isa Symbol && return "$k=$(repr(v))"
    v isa Bool && return "$k=$v"
    v isa Real && return "$k=$v"
    (v isa Tuple || v isa NamedTuple || v isa AbstractVector || v isa AbstractDict) && return "$k=$(repr(v))"
    return nothing
end

const _FITFILE_FULL_KEY_ORDER = Symbol[
    :key,
    :nchains, :datatype, :datapath, :gene, :cell, :datacond,
    :resultfolder, :label, :root,
    :G, :R, :S, :insertstep, :transitions, :coupling, :TransitionType,
    :fittedparam, :fixedeffects, :hierarchical, :shared_parameters, :priormean, :priorcv,
    :noisepriors, :nalleles, :onstates, :decayrate, :splicetype,
    :probfn, :method, :elongationtime, :grid,
    :trace_specs, :dwell_specs, :zeromedian, :datacol, :ejectnumber,
    :maxtime, :samplesteps, :warmupsteps, :sample_stride,
    :propcv, :ratetype, :writesamples, :temp, :temprna,
    :inference_method, :parallel, :device, :gradient,
    :n_samples, :n_adapts, :nuts_delta, :fd_ε, :max_depth,
    :nuts_verbose, :nuts_progress,
    :maxiter, :n_mc, :σ_floor, :init_s_raw,
    :advi_verbose, :advi_time_limit, :zygote_trace,
    :likelihood_executor, :gradient_checkpoint_length,
    :init_jitter, :init_jitter_individual, :init_jitter_noise,
    :propcv_rate, :propcv_noise, :propcv_levels,
    :proposal_cv_rates, :proposal_cv_noise, :proposal_cv_levels,
    :burst, :optimize, :prefer_legacy_ratefile, :init_ratefile,
    :rinit, :yieldfactor, :merge_max_memory, :merge_max_gb,
]

function _fitfile_full_key_order(spec::AbstractDict)
    seen = Set{Symbol}()
    ordered = Symbol[]
    for k in _FITFILE_FULL_KEY_ORDER
        haskey(spec, k) || continue
        push!(ordered, k)
        push!(seen, k)
    end
    rest = sort!(setdiff(Symbol.(collect(keys(spec))), collect(seen)); by=string)
    append!(ordered, rest)
    return ordered
end

function _format_fit_full_kw(k::Symbol, v)::String
    if k === :probfn
        return "$k=" * (v === prob_Gaussian ? "prob_Gaussian" : repr(v))
    elseif k === :method
        if v isa Tuple && length(v) == 2 && v[2] === true
            return "$k=(StochasticGene.Tsit5(), true)"
        elseif occursin("Tsit5", string(typeof(v)))
            return "$k=StochasticGene.Tsit5()"
        end
    end
    formatted = _format_fit_override(k, v)
    formatted !== nothing && return formatted
    return "$k=$(repr(v))"
end

"""
    write_fitfile_full(fitfile, spec; src="", key=nothing, skip_keys=(:infolder, :inlabel, :traceinfo, :dttype))

Write a fully formed multiline `@time fit(; ...)` script from a keyword spec.
Unlike [`write_fitfile_key`](@ref), this emits the full keyword surface into the
fitscript so the script can be inspected and edited without relying on a staged
`info_<key>.jld2` file.
"""
function write_fitfile_full(
    fitfile,
    spec::AbstractDict;
    src="",
    key=nothing,
    skip_keys=(:infolder, :inlabel, :traceinfo, :dttype),
)
    f = open(fitfile, "w")
    try
        write_prolog(f, src)
        write(f, "@time fit(;")
        order = _fitfile_full_key_order(spec)
        args = String[]
        skips = Set(Symbol.(skip_keys))
        for k in order
            k in skips && continue
            haskey(spec, k) || continue
            v = (k === :key && key !== nothing) ? key : spec[k]
            push!(args, _format_fit_full_kw(k, v))
        end
        for (i, arg) in enumerate(args)
            sep = i == length(args) ? "" : ","
            write(f, "\n    " * arg * sep)
        end
        write(f, "\n)\n")
    finally
        close(f)
    end
    return String(fitfile)
end

"""
    write_fitfile_key(fitfile, key::String; src="", resultfolder="", root=".", maxtime=nothing, samplesteps=nothing, ...)

Writes a script that runs `fit(; key=key, ...)` with optional overrides. Script-safe scalar,
symbol, tuple, named-tuple, and vector overrides are written.
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
    write_fitfile_genes(fitfile, nchains, datatype, datapath, cell, datacond, resultfolder, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs;
        src="", dttype=nothing, traceinfo=nothing)

Writes a fit script that takes the gene as `ARGS[1]` and calls the positional [`fit`](@ref) entry point:

`fit(nchains, datatype, datapath, gene, cell, datacond, resultfolder, label, fittedparam, …, trace_specs, dwell_specs)`,

optionally appending `; dttype=…` / `; traceinfo=…` for legacy loaders when those are set, and `; inference_method=…` when passed (forwarded as a trailing keyword on the `fit` call).
"""
function write_fitfile_genes(fitfile, nchains, datatype, datapath, cell, datacond, resultfolder, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, temp, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs;
    src="", dttype=nothing, traceinfo=nothing, inference_method=nothing, sample_stride=nothing, parallel=nothing, device=nothing, merge_max_memory=nothing, merge_max_gb=nothing, likelihood_executor=nothing, gradient_checkpoint_length=nothing, proposal_cv_rates=nothing, proposal_cv_noise=nothing, init_jitter=nothing, init_jitter_individual=nothing, init_jitter_noise=nothing)
    s = '"'
    transitions = transitions isa AbstractVector && !(transitions isa Tuple) ? Tuple(transitions) : transitions
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    datacond_q = datacond
    if datacond_q isa AbstractString
        datacond_q = "$s$datacond_q$s"
    end
    line = "@time fit($nchains, $s$datatype$s, $datapath, ARGS[1], $s$cell$s, $datacond_q, $s$resultfolder$s, $s$label$s, $fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $temp, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber, $yieldfactor, $trace_specs, $dwell_specs"
    extra = String[]
    if dttype !== nothing && !(dttype isa AbstractVector && isempty(dttype))
        push!(extra, "dttype=" * repr(dttype))
    end
    if traceinfo !== nothing
        push!(extra, "traceinfo=" * repr(traceinfo))
    end
    if inference_method !== nothing
        push!(extra, "inference_method=" * repr(inference_method))
    end
    for (name, value) in (
        (:sample_stride, sample_stride),
        (:parallel, parallel),
        (:device, device),
        (:merge_max_memory, merge_max_memory),
        (:merge_max_gb, merge_max_gb),
        (:likelihood_executor, likelihood_executor),
        (:gradient_checkpoint_length, gradient_checkpoint_length),
        (:proposal_cv_rates, proposal_cv_rates),
        (:proposal_cv_noise, proposal_cv_noise),
        (:init_jitter, init_jitter),
        (:init_jitter_individual, init_jitter_individual),
        (:init_jitter_noise, init_jitter_noise),
    )
        value === nothing && continue
        push!(extra, string(name) * "=" * repr(value))
    end
    if !isempty(extra)
        line *= "; " * join(extra, ", ")
    end
    line *= ")\n"
    write(f, line)
    close(f)
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
    batchsize > 0 || throw(ArgumentError("batchsize must be positive"))
    nbatches = cld(ngenes, batchsize)
    batches = Vector{Vector{String}}(undef, nbatches)
    for i in 1:nbatches
        first = batchsize * (i - 1) + 1
        last = min(batchsize * i, ngenes)
        batches[i] = genes[first:last]
    end
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
