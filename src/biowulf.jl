# This file is part of StochasticGene.jl  

# biowulf.jl
# functions for use on the NIH Biowulf super computer


"""
    makeswarm(; nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", batchsize::Int=1000, juliafile::String="fitscript", thresholdlow::Float64=0.0, thresholdhigh::Float64=Inf, datatype::String="", dttype::Vector=String[], datapath="", cell::String="HBEC", datacond="", interval=1.0, nascent=0.5, infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, root=".", priormean=Float64[], priorcv::Float64=10.0, nalleles=2, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")
    fittedparam::Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1,2],[2,1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, root=".", nalleles=2, priormean=[],  priorcv::Float64=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false,src="")


returns swarmfile used on biowulf

Arguments
- `nchains`: number of MCMC chains = number of processors called by Julia, default = 2
- 'nthreads`: number of Julia threads per processesor, default = 1
- `swarmfile`: name of swarmfile to be executed by swarm
- `batchsize`: number of jobs per swarmfile, default = 1000
- `juliafile`: name of file to be called by julia in swarmfile
- `thresholdlow`: lower threshold for halflife for genes to be fit
- `threhsoldhigh`: upper threshold
- `datatype`: String that desecribes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracenascent", "tracerna"
- `dttype`
- `datapath`: path to data file or folder or array of files or folders
- `cell': cell type for halflives and allele numbers
- `datacond`: string or vector of strings describing data treatment condition, e.g. "WT", "DMSO" or ["DMSO","AUXIN"]
- `interval`: frame interval of intensity traces
- `nascent`: fraction of alleles exhibiting nascent rna
- `infolder`: result folder used for initial parameters
- `resultfolder`: folder for results of MCMC run
- `label`: label of output files produced
- `inlabel`: label of files used for initial conditions
- `fittedparam`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
- `fixedeffects`: tuple of vectors of rates that are fixed between control and treatment where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3
     (each vector in tuple is a fixed rate set)
- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G`: number of gene states
- `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S`: number of splice sites (set to 0 for classic telegraph models and R for GRS models)
- `insertstep`: R step where reporter is inserted
- `root`: name of root directory for project, e.g. "scRNA"
- `priormean`: mean of prior rate distribution
- 'priorcv`: coefficient of variation for the rate prior distributions, default is 10.
- `nalleles`: number of alleles, value in alleles folder will be used if it exists  
- `onstates`: vector of on or sojourn states
- `decayrate`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `splicetype`: RNA pathway for GRS models, (e.g. "offeject" =  spliced intron is not viable)
- `probfn`: probability function for hmm observation probability (e.g. prob_GaussianMixture)
- `noiseparams`: number of parameters of probfn
- `weightind`: parameter index of bias probability of mixtures, e.g. noiseparams=5, weightind=5 means last noise parameter is for mixture bias
- `hyperparams`: tuple of hyper parameters
- `ratetype`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
- `propcv`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `maxtime`: maximum wall time for run, default = 60 min
- `samplesteps`: number of MCMC sampling steps
- `warmupsteps`: number of MCMC warmup steps to find proposal distribution covariance
- `annealsteps`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `temp`: MCMC temperature
- `tempanneal`: annealing temperature
- `temprna`: reduce rna counts by temprna compared to dwell times
- `burst`: if true then compute burst frequency
- `optimize`: use optimizer to compute maximum likelihood value
- `writesamples`: write out MH samples if true, default is false
- `src`: path to folder containing StochasticGene.jl/src
   
"""
function makeswarm(; nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", batchsize::Int=1000, juliafile::String="fitscript", thresholdlow::Float64=0.0, thresholdhigh::Float64=Inf, datatype::String="", dttype::Vector=String[], datapath="", cell::String="HBEC", datacond="", interval=1.0, nascent=0.5, infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, root=".", priormean=Float64[], priorcv::Float64=10.0, nalleles=2, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, hyperparams=tuple(),ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")

    makeswarm(checkgenes(root, datacond, datapath, cell, thresholdlow, thresholdhigh), nchains=nchains, nthreads=nthreads, swarmfile=swarmfile, batchsize=batchsize, juliafile=juliafile, datatype=datatype, dttype=dttype, datapath=datapath, cell=cell, datacond=datacond, interval=interval, nascent=nascent, infolder=infolder, resultfolder=resultfolder, inlabel=inlabel, label=label,
        fittedparam=fittedparam, fixedeffects=fixedeffects, transitions=transitions, G=G, R=R, S=S, insertstep=insertstep, root=root, priormean=priormean, nalleles=nalleles, priorcv=priorcv, onstates=onstates, decayrate=decayrate, splicetype=splicetype, probfn=probfn, noiseparams=noiseparams, weightind=weightind, hyperparams=hyperparams,ratetype=ratetype,
        propcv=propcv, maxtime=maxtime, samplesteps=samplesteps, warmupsteps=warmupsteps, annealsteps=annealsteps, temp=temp, tempanneal=tempanneal, temprna=temprna, burst=burst, optimize=optimize, writesamples=writesamples, src=src)
end



"""
    makeswarm(genes::Vector; nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="", dttype=String[], datapath="", cell::String="", datacond="", interval=1.0, nascent=0.5, infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, root=".", priormean=Float64[], nalleles=2, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")
    fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, root=".", priormean=[], nalleles=2, priorcv::Float64=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false,src="")


"""
function makeswarm(genes::Vector; nchains::Int=2, nthreads::Int=1, swarmfile::String="fit", batchsize=1000, juliafile::String="fitscript", datatype::String="", dttype=String[], datapath="", cell::String="", datacond="", interval=1.0, nascent=0.5, infolder::String="", resultfolder::String="test", inlabel::String="", label::String="",
    fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, root=".", priormean=Float64[], nalleles=2, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5,   hyperparams=tuple(),ratetype="median",
    propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")
    if R > 0
        if S > 0
            S = R
        end
        modelstring = "$G$R$S$insertstep"
    else
        modelstring = "$G"
    end

    if isempty(label)
        label = datatype * "-" * cell
        typeof(datacond) <: AbstractString && (label = label * "_" * datacond)
    end
    isempty(inlabel) && (inlabel = label)

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
    write_fitfile(joinpath(root, juliafile), nchains, datatype, dttype, datapath, cell, datacond, interval, nascent, infolder, resultfolder, inlabel, label,
        fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, nalleles, priorcv, onstates,
        decayrate, splicetype, probfn, noiseparams, weightind, hyperparams,ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)
end
"""
    write_fitfile(fitfile, nchains, datatype, dttype, datapath, cell, datacond, interval, nascent, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noiseparams, weightind, ratetype,
    propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples)

"""
function write_fitfile(fitfile, nchains, datatype, dttype, datapath, cell, datacond, interval, nascent, infolder, resultfolder, inlabel, label,
    fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, nalleles, priorcv, onstates,
    decayrate, splicetype, probfn, noiseparams, weightind, hyperparams,ratetype,
    propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, src)
    # function write_fitfile(fitfile, nchains, cell, datacond, G, R, S, transitions, insertstep, onstates, maxtime, fittedparam, fixedeffects, infolder, resultfolder, datapath, datatype, inlabel, label, nsets, transient, samplesteps, warmupsteps, annealsteps, temp, tempanneal, root, cv, priorcv, decayrate, burst, nalleles, optimize, splicetype, ratetype, writesamples)
    s = '"'
    s3 = s * s * s
    f = open(fitfile, "w")
    if isempty(src)
        write(f, "@everywhere using StochasticGene\n")
    else
        write(f, "@everywhere include($s$src$s)\n")
        write(f, "@everywhere using .StochasticGene\n")
    end
    typeof(datapath) <: AbstractString && (datapath = "$s$datapath$s")
    typeof(datacond) <: AbstractString && (datacond = "$s$datacond$s")
    # write(f, "@time fit($nchains,ARGS[1],$s$cell$s,$fittedparam,$fixedeffects,$transitions,$s$datacond$s,$G,$R,$S,$insertstep,$maxtime,$s$infolder$s,$s$resultfolder$s,$s$datapath$s,$s$datatype$s,$s$inlabel$s,$s$label$s,$nsets,$cv,$transient,$samplesteps,$warmupsteps,$annealsteps,$temp,$tempanneal,$s$root$s,$priorcv,$decayrate,$burst,$nalleles,$optimize,$s$splicetype$s,$s$ratetype$s,$writesamples)\n")
    write(f, "@time fit($nchains, $s$datatype$s, $dttype, $datapath, ARGS[1], $s$cell$s, $datacond, $interval, $nascent, $s$infolder$s, $s$resultfolder$s, $s$inlabel$s, $s$label$s,$fittedparam, $fixedeffects, $transitions, $G, $R, $S, $insertstep, $s$root$s, $maxtime, $priormean, $priorcv, $nalleles, $onstates, $decayrate, $s$splicetype$s, $probfn, $noiseparams, $weightind, $hyperparams, $s$ratetype$s,$propcv, $samplesteps, $warmupsteps, $annealsteps, $temp, $tempanneal, $temprna, $burst, $optimize, $writesamples)")
    close(f)
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
    write_swarmfile(sfile, nchains, nthreads, juliafile, genes::Vector)


"""
function write_swarmfile(sfile, nchains, nthreads, juliafile, genes::Vector)
    f = open(sfile, "w")
    for gene in genes
        gene = check_genename(gene, "(")
        writedlm(f, ["julia -t $nthreads -p" nchains juliafile gene])
        # writedlm(f,["julia -p" nchains juliafile nchains gene cell cond G maxtime infolder resultfolder datapath fish inlabel label nsets runcycle transient fittedparam fixedeffects])
    end
    close(f)
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
    if ~ispath(folder)
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
