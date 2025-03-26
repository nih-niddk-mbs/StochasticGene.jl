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
    makeswarm(;<keyword arguments>)

write swarm and fit files used on biowulf


#Arguments

- 'nthreads=1`: number of Julia threads per processesor, default = 1
- `swarmfile::String="fit"`: name of swarmfile to be executed by swarm
- `juliafile::String="fitscript`: name of file to be called by julia in swarmfile
- `src=""`: path to folder containing StochasticGene.jl/src (only necessary if StochasticGene not installed)

and all keyword arguments of function fit(; <keyword arguments> )

see fit

Note:: the keyword 'method' is handled slightly differently here than in the function fit.  In fit it is a function (i.e. numerical method function used in DifferentialEquations.jl) or a tuple
of the numerical method and a Bool for hierarchical models.  However, in biowulf.jl, the numerical method must be a String, i.e. use "lsoda()" for lsoda().  This is because when Julia writes
the function, it will parse the numerical method rather than just writing it.


"""
function makeswarm(; gene::String="", nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript", datatype::String="", dttype=String[], datapath="", cell::String="", datacond="", traceinfo=(1.0, 1.0, -1, 1.0), infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, coupling = tuple(), TransitionType="", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian::Bool=false, datacol=3, ejectnumber=1)
    modelstring = create_modelstring(G, R, S, insertstep)
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, TransitionType)
    juliafile = juliafile * "_" * label * "_" * "$modelstring" * ".jl"
    sfile = swarmfile * "_" * label * "_" * "$modelstring" * ".swarm"
    write_swarmfile(joinpath(root, sfile), nchains, nthreads, juliafile)
    write_fitfile(joinpath(root, juliafile), nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber)
end




"""
    makeswarm_genes(genes::Vector{String}; <keyword arguments> )

write a swarmfile and fit files to run all each gene in vector genes

# Arguments
- `genes`: vector of genes
- `batchsize=1000`: number of jobs per swarmfile, default = 1000

and all arguments in makeswarm(;<keyword arguments>)


    Examples

julia> genes = ["MYC","SOX9"]

julia> makeswarm(genes,cell="HBEC")
"""
function makeswarm(genes::Vector{String}; nchains::Int=2, nthreads=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="rna", dttype=String[], datapath="", cell::String="HCT116", datacond="MOCK", traceinfo=(1.0, 1.0, -1, 1.0), infolder::String="HCT116_test", resultfolder::String="HCT116_test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, coupling=tuple(), TransitionType="", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian::Bool=false, datacol=3, ejectnumber=1)
    modelstring = create_modelstring(G, R, S, insertstep)
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, TransitionType)
    ngenes = length(genes)
    println("number of genes: ", ngenes)
    juliafile = juliafile * "_" * label * "_" * "$modelstring" * ".jl"
    if ngenes > batchsize
        batches = getbatches(genes, ngenes, batchsize)
        for batch in eachindex(batches)
            sfile = swarmfile * "_" * label * "_" * "$modelstring" * "_" * "$batch" * ".swarm"
            write_swarmfile(sfile, nchains, nthreads, juliafile, batches[batch])
        end
    else
        sfile = swarmfile * "_" * label * "_" * "$modelstring" * ".swarm"
        write_swarmfile(joinpath(root, sfile), nchains, nthreads, juliafile, genes)
    end
    write_fitfile_genes(joinpath(root, juliafile), nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber)
end
"""
    makeswarm_genes(;<keyword arguments> )
 
#Arguments
    - `thresholdlow::Float=0`: lower threshold for halflife for genes to be fit
    - `threhsoldhigh::=Inf`: upper threshold

    and all keyword arguments in makeswarm(;<keyword arguments>)
"""
function makeswarm_genes(; nchains::Int=2, nthreads=1, swarmfile::String="fit", batchsize::Int=1000, juliafile::String="fitscript", thresholdlow::Float64=0.0, thresholdhigh::Float64=Inf, datatype::String="", dttype::Vector=String[], datapath="", cell::String="HBEC", datacond="", traceinfo=(1.0, 1.0, -1, 1.0), infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, coupling=tuple(), TransitionType="", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], priorcv::Float64=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian::Bool=false, datacol=3, ejectnumber=1)

    makeswarm(checkgenes(root, datacond, datapath, cell, thresholdlow, thresholdhigh), nchains=nchains, nthreads=nthreads, swarmfile=swarmfile, batchsize=batchsize, juliafile=juliafile, datatype=datatype, dttype=dttype, datapath=datapath, cell=cell, datacond=datacond, traceinfo=traceinfo, infolder=infolder, resultfolder=resultfolder, inlabel=inlabel, label=label,
        fittedparam=fittedparam, fixedeffects=fixedeffects, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, coupling=coupling, TransitionType=TransitionType, grid=grid, root=root, elongationtime=elongationtime, priormean=priormean, nalleles=nalleles, priorcv=priorcv, onstates=onstates, decayrate=decayrate, splicetype=splicetype, probfn=probfn, noisepriors=noisepriors, hierarchical=hierarchical, ratetype=ratetype,
        propcv=propcv, maxtime=maxtime, samplesteps=samplesteps, warmupsteps=warmupsteps, annealsteps=annealsteps, temp=temp, tempanneal=tempanneal, temprna=temprna, burst=burst, optimize=optimize, writesamples=writesamples, method=method, src=src, zeromedian=zeromedian, datacol=datacol, ejectnumber=ejectnumber)
end

"""
    makeswarm(models::Vector{ModelArgs}; <keyword arguments> )

creates a run for each model

#Arguments
- `models::Vector{ModelArgs}`: Vector of ModelArgs structures

and all keyword arguments in makeswarm(;<keyword arguments>)
"""
function makeswarm(models::Vector{ModelArgs}; gene="", nchains::Int=2, nthreads=1, swarmfile::String="fit", juliafile::String="fitscript", datatype::String="", dttype=String[], datapath="", cell::String="", datacond="", traceinfo=(1.0, 1.0, -1, 1.0), infolder::String="", resultfolder::String="test",
    fittedparam::Vector=Int[], grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], nalleles=1, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method="Tsit5()", src="", zeromedian::Bool=false, datacol=3, ejectnumber=1)
    juliafile = juliafile * "_" * gene * "_" * datacond * ".jl"
    sfile = swarmfile * "_" * gene * "_" * datacond * ".swarm"
    write_swarmfile(joinpath(root, sfile), nchains, nthreads, juliafile, datatype, datacond, cell, models)
    write_fitfile_models(joinpath(root, juliafile), nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, coupling, fittedparam,
        grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian, datacol, ejectnumber)
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
function write_swarmfile(sfile, nchains, nthreads, juliafile::String, project="")
    f = open(sfile, "w")
    if isempty(project)
        writedlm(f, ["julia -t $nthreads -p" nchains juliafile])
    else
        writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile])
    end
    close(f)
end

"""
    write_swarmfile(sfile, nchains, nthreads, juliafile::String, genes::Vector{String}, project="")

Writes a swarmfile for a vector of genes.

# Arguments
- `sfile`: The path to the swarmfile to be written.
- `nchains`: The number of chains to run.
- `nthreads`: The number of threads to use.
- `juliafile::String`: The path to the Julia file to be executed.
- `genes::Vector{String}`: A vector of gene names.
- `project`: The Julia project to use (optional, default is an empty string).

# Description
This function writes a swarmfile that specifies how to run a Julia file for each gene in the given vector of genes, with the specified number of chains and threads. If a project is specified, it includes the `--project` flag in the command.
"""
function write_swarmfile(sfile, nchains, nthreads, juliafile::String, genes::Vector{String}, project="")
    f = open(sfile, "w")
    for gene in genes
        gene = check_genename(gene, "(")
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains juliafile gene])
        else
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile gene])
        end
    end
    close(f)
end

"""
    write_swarmfile(sfile, nchains, nthreads, juliafile, datatype, datacond, cell, models::Vector{ModelArgs}, project="")

Writes a swarmfile for a vector of models.

# Arguments
- `sfile`: The path to the swarmfile to be written.
- `nchains`: The number of chains to run.
- `nthreads`: The number of threads to use.
- `juliafile`: The path to the Julia file to be executed.
- `datatype`: The type of data.
- `datacond`: The data condition.
- `cell`: The type of cell.
- `models::Vector{ModelArgs}`: A vector of models.
- `project`: The Julia project to use (optional, default is an empty string).

# Description
This function writes a swarmfile that specifies how to run a Julia file for each model in the given vector of models, with the specified number of chains and threads. If a project is specified, it includes the `--project` flag in the command.
"""
function write_swarmfile(sfile, nchains, nthreads, juliafile, datatype, datacond, cell, models::Vector{ModelArgs}, project="")
    f = open(sfile, "w")
    for model in models
        label, inlabel = create_label(model.label, model.inlabel, datatype, datacond, cell, model.TransitionType)
        if isempty(model.fixedeffects)
            fixedeffects = "1"
        else
            fixedeffects = model.fixedeffects
        end
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains juliafile inlabel label model.G model.R model.S model.insertstep model.TransitionType fixedeffects])
        else
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile inlabel label model.G model.R model.S model.insertstep model.TransitionType fixedeffects])
        end
    end
    close(f)
end


"""
    write_fitfile(fitfile, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src)

Writes a fitfile for the fit function parameters.

# Arguments
arguments are as in fit()

# Description
This function writes a fitfile that specifies how to run a fit for a given set of parameters. It includes the necessary prolog and handles the formatting of various parameters.
"""
function write_fitfile(fitfile, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian::Bool, datacol, ejectnumber)
    s = '"'
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, $s$gene$s, $s$cell$s, $datacond, $traceinfo, $s$infolder$s, $s$resultfolder$s, $s$inlabel$s, $s$label$s,$fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber)")
    close(f)
end

"""
    write_fitfile_genes(fitfile, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src)

Write fitfile for genes

# Arguments as in write_fitfile
"""
function write_fitfile_genes(fitfile, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian::Bool, datacol, ejectnumber)
    s = '"'
    # s3 = s * s * s
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, ARGS[1], $s$cell$s, $datacond, $traceinfo, $s$infolder$s, $s$resultfolder$s, $s$inlabel$s, $s$label$s, $fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber)")
    close(f)
end

"""
    write_fitfile(fitfile, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder,
    root, maxtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src)

write fitfile for models

# Arguments as in write_fitfile
"""
function write_fitfile_models(fitfile, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, coupling, fittedparam, grid,
    root, maxtime, elongationtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, src, zeromedian::Bool, datacol, ejectnumber)
    s = '"'
    f = open(fitfile, "w")
    write_prolog(f, src)
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")

    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, $s$gene$s, $s$cell$s, $datacond, $traceinfo, $s$infolder$s, $s$resultfolder$s, ARGS[1], ARGS[2], $fittedparam, ARGS[8], ARGS[3], ARGS[4], ARGS[5], ARGS[6], ARGS[7], $coupling, $grid, $s$root$s, $maxtime, $elongationtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noisepriors, $hierarchical, $s$ratetype$s, $propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples, $method, $zeromedian, $datacol, $ejectnumber)")
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
        write(f, "@everywhere using Pkg\n")
        write(f, "@everywhere Pkg.activate($s$src$s)\n")
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
    create_label(label, inlabel, datatype, join(datacond,"-"), cell, TransitionType)
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
function checkgenes(root, conds, datapath, celltype::String, thresholdlow::Float64, thresholdhigh::Float64)
    genes = Vector{Vector{String}}(undef, 0)
    typeof(conds) <: AbstractString && (conds = [conds])
    typeof(datapath) <: AbstractString && (datapath = [datapath])

    for d in datapath, c in conds
        push!(genes, checkgenes(root, c, d, celltype, thresholdlow, thresholdhigh))
    end
    geneset = genes[1]
    for g in genes
        genesest = intersect(geneset, g)
    end
    return geneset
end
"""
    checkgenes(root, conds, datapath, celltype::String, thresholdlow::Float64, thresholdhigh::Float64)

Checks and returns the genes that meet the specified thresholds across multiple conditions and data paths.

### Arguments
- `root`: The root directory.
- `cond`: The condition to check.
- `datapath`: The data path to check.
- `cell`: The type of cell.
- `thresholdlow`: The lower threshold for filtering genes.
- `thresholdhigh`: The upper threshold for filtering genes.

### Returns
- `Vector{String}`: A vector of genes that meet the specified thresholds for the given condition and data path.
"""
function checkgenes(root, cond::String, datapath::String, cell::String, thresholdlow::Float64, thresholdhigh::Float64)
    if cell == "HBEC"
        return genes_hbec()
    else
        datapath = folder_path(datapath, root, "data")
        genes = intersect(get_halflives(root, cell, thresholdlow, thresholdhigh), get_genes(cond, datapath))
        alleles = get_alleles(root, cell)
        if ~isnothing(alleles)
            return intersect(genes, alleles)
        else
            return genes
        end
    end
end

"""
    folder_path(folder::String, root::String, folderatetype::String=""; make=false)

Returns the full path for a given folder, optionally creating the path if it does not exist.

# Arguments
- `folder`: A folder name.
- `root`: The root directory.
- `foldertype`: The type of folder (optional, default is an empty string).
- `make`: A boolean flag indicating whether to create the path if it does not exist (optional, default is `false`).

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
            else
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
    get_halflives(path, thresholdlow, thresholdhigh)
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
    halflives = readdlm(hlpath, ',')
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
        println("folder $folder does not exist")
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
