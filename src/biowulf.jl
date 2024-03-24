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
`Gfamily::String`: type of model, e.g. "nstate", "KP", "Refract"
`fixedeffects::String`: two numbers separated by a hyphen, e.g. "3-4", indicating parameters 3 and 4 are fixed to each other

"""

struct ModelArgs
    inlabel::String
    label::String
    G::Int
    R::Int
    S::Int
    insertstep::Int
    Gfamily::String
    fixedeffects::String
end


"""
    makeswarm(genes::Vector{String}; <keyword arguments> )


write swarm and fit files used on biowulf
creates a run for each gene

#Arguments
- `genes`: vector of genes
- 'nthreads::Int=1`: number of Julia threads per processesor, default = 1
- `swarmfile::String="fit"`: name of swarmfile to be executed by swarm
- `batchsize=1000`: number of jobs per swarmfile, default = 1000
- `juliafile::String="fitscript`: name of file to be called by julia in swarmfile
- `thresholdlow::Float=0`: lower threshold for halflife for genes to be fit
- `threhsoldhigh::=Inf`: upper threshold
- `src=""`: path to folder containing StochasticGene.jl/src

and all keyword arguments of function fit(; <keyword arguments> )

see fit

Examples

julia> genes = ["MYC","SOX9"]

julia> makeswarm(genes,cell="HBEC")

"""
function makeswarm(genes::Vector{String}; nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="", dttype=String[], datapath="", cell::String="", datacond="", traceinfo=(1.0, 1.0, 0.65), nascent=(1, 2), infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, Gfamily="", root=".", priormean=Float64[], nalleles=2, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=0, hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")
    modelstring = create_modelstring(G, R, S, insertstep)
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, Gfamily)
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
    write_fitfile(joinpath(root, juliafile), nchains, datatype, dttype, datapath, cell, datacond, traceinfo, nascent, infolder, resultfolder, inlabel, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noiseparams, weightind, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)
end

"""
    makeswarm(models::Vector{ModelArgs}; <keyword arguments> )

creates a run for each model

#Arguments
- `models::Vector{ModelArgs}`: Vector of ModelArgs structures
- keyword arguments are same as above
"""
function makeswarm(models::Vector{ModelArgs}; gene="", nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", juliafile::String="fitscript", datatype::String="", dttype=String[], datapath="", cell::String="", datacond="", traceinfo=(1.0, 1.0, 0.65), nascent=(1, 2), infolder::String="", resultfolder::String="test",
    fittedparam::Vector=Int[], root=".", priormean=Float64[], nalleles=2, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=0, hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")
    juliafile = juliafile * "_" * gene * "_" * datacond * ".jl"
    sfile = swarmfile * "_" * gene * "_" * datacond * ".swarm"
    write_swarmfile(joinpath(root, sfile), nchains, nthreads, juliafile, datatype, datacond, cell, models)
    write_fitfile(joinpath(root, juliafile), nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, nascent, infolder, resultfolder,
        fittedparam, root, maxtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noiseparams, weightind, hierarchical, ratetype, propcv,
        samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)
end


function makeswarm(; nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", batchsize::Int=1000, juliafile::String="fitscript", thresholdlow::Float64=0.0, thresholdhigh::Float64=Inf, datatype::String="", dttype::Vector=String[], datapath="", cell::String="HBEC", datacond="", traceinfo=(1.0, 1.0, 0.65), nascent=(1, 2), infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, Gfamily="", root=".", priormean=Float64[], priorcv::Float64=10.0, nalleles=2, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=0, hierarchical=tuple(), ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")

    makeswarm(checkgenes(root, datacond, datapath, cell, thresholdlow, thresholdhigh), nchains=nchains, nthreads=nthreads, swarmfile=swarmfile, batchsize=batchsize, juliafile=juliafile, datatype=datatype, dttype=dttype, datapath=datapath, cell=cell, datacond=datacond, traceinfo=traceinfo, nascent=nascent, infolder=infolder, resultfolder=resultfolder, inlabel=inlabel, label=label,
        fittedparam=fittedparam, fixedeffects=fixedeffects, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, Gfamily=Gfamily, root=root, priormean=priormean, nalleles=nalleles, priorcv=priorcv, onstates=onstates, decayrate=decayrate, splicetype=splicetype, probfn=probfn, noiseparams=noiseparams, weightind=weightind, hierarchical=hierarchical, ratetype=ratetype,
        propcv=propcv, maxtime=maxtime, samplesteps=samplesteps, warmupsteps=warmupsteps, annealsteps=annealsteps, temp=temp, tempanneal=tempanneal, temprna=temprna, burst=burst, optimize=optimize, writesamples=writesamples, src=src)
end

"""
    write_swarmfile(sfile, nchains, nthreads, juliafile, genes::Vector)


"""
function write_swarmfile(sfile, nchains, nthreads, juliafile, genes::Vector{String},project="")
    f = open(sfile, "w")
    for gene in genes
        gene = check_genename(gene, "(")
        writedlm(f, ["julia -t $nthreads -p" nchains juliafile gene])
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains juliafile gene])
        else
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile gene])
        end
    end
    close(f)
end

"""
    write_swarmfile(sfile, nchains, nthreads, juliafile, datatype, datacond, cell, models::Vector{ModelArgs})


"""
function write_swarmfile(sfile, nchains, nthreads, juliafile, datatype, datacond, cell, models::Vector{ModelArgs},project="")
    f = open(sfile, "w")
    for model in models
        label, inlabel = create_label(model.label, model.inlabel, datatype, datacond, cell, model.Gfamily)
        if isempty(model.fixedeffects)
            fixedeffects = "1"
        else
            fixedeffects = model.fixedeffects
        end
        if isempty(project)
            writedlm(f, ["julia -t $nthreads -p" nchains juliafile inlabel label model.G model.R model.S model.insertstep model.Gfamily fixedeffects])
        else
            writedlm(f, ["julia --project=$project -t $nthreads -p" nchains juliafile inlabel label model.G model.R model.S model.insertstep model.Gfamily fixedeffects])
        end
    end
    close(f)
end

"""
    write_fitfile(fitfile, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, nascent, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noiseparams, weightind, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)

"""
function write_fitfile(fitfile, nchains, datatype, dttype, datapath, cell, datacond, traceinfo, nascent, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noiseparams, weightind, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)
    s = '"'
    # s3 = s * s * s
    f = open(fitfile, "w")
    if isempty(src)
        write(f, "@everywhere using StochasticGene\n")
    else
        write(f, "@everywhere include($s$src$s)\n")
        write(f, "@everywhere using .StochasticGene\n")
    end
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, ARGS[1], $s$cell$s, $datacond, $traceinfo, $nascent, $s$infolder$s, $s$resultfolder$s, $s$inlabel$s, $s$label$s,$fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $s$root$s, $maxtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noiseparams, $weightind, $hierarchical, $s$ratetype$s,$propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples)")
    close(f)
end

function write_fitfile(fitfile, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, nascent, infolder, resultfolder,
    fittedparam, root, maxtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noiseparams, weightind, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)
    s = '"'
    f = open(fitfile, "w")
    if isempty(src)
        write(f, "@everywhere using StochasticGene\n")
    else
        write(f, "@everywhere include($s$src$s)\n")
        write(f, "@everywhere using .StochasticGene\n")
    end
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")

    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, $s$gene$s, $s$cell$s, $datacond, $traceinfo, $nascent, $s$infolder$s, $s$resultfolder$s, ARGS[1], ARGS[2],$fittedparam, ARGS[8], ARGS[3], ARGS[4], ARGS[5], ARGS[6], ARGS[7], $s$root$s, $maxtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noiseparams, $weightind, $hierarchical, $s$ratetype$s,$propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples)")
    close(f)
end

"""
    create_label(label,inlabel,datacond,cell,Gfamily)

"""
function create_label(label, inlabel, datatype, datacond, cell, Gfamily)
    if isempty(label)
        label = datatype * "-" * cell
        ~isempty(Gfamily) && (label = label * "-" * Gfamily)
        typeof(datacond) <: AbstractString && (label = label * "_" * datacond)
    end
    isempty(inlabel) && (inlabel = label)
    return label, inlabel
end

"""
    create_modelstring(G,R,S,insertstep)

"""
function create_modelstring(G, R, S, insertstep)
    if R > 0
        if S > 0
            S = R - insertstep + 1
        end
        return "$G$R$S$insertstep"
    else
        return "$G"
    end
end

"""
    fix(folder)

    Finds jobs that failed and writes a swarmfile for those genes

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
    setup(rootfolder = "scRNA")

    Sets up the folder system prior to use
    Defaults to "scRNA"

"""
function rna_setup(root=".")
    folder_setup(root)
    alleles = joinpath(data, "alleles")
    halflives = joinpath(data, "halflives")
    testdata = joinpath(data, "HCT116_testdata")

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

creates data and results folders (if they do not already exit)
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
    checkgenes(root, cond::String, datapath::String, cell::String, thresholdlow::Float64, thresholdhigh::Float64)


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


"""
function folder_path(folder::String, root::String, folderatetype::String=""; make=false)
    f = folder
    if ~ispath(folder) && ~isempty(folder)
        f = joinpath(root, folder)
        if ~ispath(f)
            f = joinpath(root, folderatetype, folder)
            if ~ispath(f) && ~make
                throw("$folder not found")
            else
                mkpath(f)
            end
        end
    end
    f
end

"""
    folder_path(folder::Vector, root, foldertype)

TBW
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

TBW
"""
function get_halflives(root, cell, thresholdlow::Float64, thresholdhigh::Float64)
    path = get_file(root, "data/halflives", cell, "csv")
    get_halflives(path, thresholdlow, thresholdhigh)
end

"""
    get_halflives(hlpath, thresholdlow::Float64, thresholdhigh::Float64)

TBW
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
    get_alleles(root, cell) = get_alleles(get_file(root, "data/alleles", cell, "csv"))


TBW
"""
function get_alleles(allelepath)
    if ~isnothing(allelepath)
        return readdlm(allelepath, ',')[2:end, 1]
    else
        return nothing
    end
end

get_alleles(root, cell) = get_alleles(get_file(root, "data/alleles", cell, "csv"))

get_file(root, folder, filetype, suffix) = get_file(joinpath(root, folder), filetype, suffix)

"""
    get_file(folder, filetype, suffix)
    get_file(root, folder, filetype, suffix) = get_file(joinpath(root, folder), filetype, suffix)

TBW
"""
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

TBW
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

TBW
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

TBW
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

TBW
"""
function findswarm(job)
    sc = "Swarm Command"
    line = read(pipeline(`jobhist $job`, `grep $sc`), String)
    list = split(line, " ")
    list[occursin.(".swarm", list)][1]
end

"""
    get_missing_genes(datapath::String,folder::String,cell,type,label,cond,model)
"""

"""
    get_missing_genes(datapath::String, resultfolder::String, cell, filetype, label, cond, model, root=".")

TBW
"""
function get_missing_genes(datapath::String, resultfolder::String, cell, filetype, label, cond, model, root=".")
    genes = checkgenes(root, cond, datapath, cell, 0.0, 1e8)
    get_missing_genes(genes, folder_path(resultfolder, root, "results"), filetype, label, cond, model)
end

"""
    get_missing_genes(genes::Vector, resultfolder, filetype, label, cond, model)

TBW
"""
function get_missing_genes(genes::Vector, resultfolder, filetype, label, cond, model)
    genes1 = get_genes(resultfolder, filetype, label, cond, model)
    get_missing_genes(genes, genes1)
end

get_missing_genes(genes, genes1) = union(setdiff(genes1, genes), setdiff(genes, genes1))

"""
    scan_swarmfiles(jobid, folder=".")

TBW
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

TBW
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

TBW
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

TBW
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
