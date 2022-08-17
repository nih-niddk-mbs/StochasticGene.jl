# biowulf.jl
# functions for use on the NIH Biowulf super computer

"""
    function makeswarm(;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel="label",timestamp="",nsets=1,datafolder::String="HCT116_testdata",fish= false,thresholdlow::Float64=0.,thresholdhigh::Float64=1e8,conds::String="MOCK",resultfolder::String= "fit_result",infolder=resultfolder,batchsize=1000,maxtime = 60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=40000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv = 0.02,yieldprior=0.05,priorcv= 10.)

    function makeswarm(genes::Vector;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel="label",timestamp="",nsets=1,datafolder::String="HCT116_testdata",fish=false,conds::String="MOCK",resultfolder::String="fit_result",infolder=resultfolder,batchsize=1000,maxtime=60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=40000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv=0.02,yieldprior=0.05,priorcv=10.)
        if label == "label"

    Arguments
    - `G`: number of gene states
    - `cell': cell type for halflives and allele numbers
    - `infolder`: name of folder for initial parameters
    - `swarmfile`: name of swarmfile to be executed by swarm
    - `label`: label of output files produced
    - `inlabel`: label of files used for initial conditions
    - `timestamp`: label for time of sample (e.g. T120)
    - `nsets`: number of histograms to be fit (e.g. one for wild type and one for perturbation)
    - `datafolder`: folder holding histograms, if two folders use `-` (hyphen) to separate, e.g.  "data\folder1-data\folder2"
    - `datatype`: String that desecribes data file type, e.g. "scRNA", "fish", "genetrap", "RNA-counts", (presumes folder structure)
    - `thresholdlow`: lower threshold for halflife for genes to be fit
    - `threhsoldhigh`: upper threshold
    - `conds`: string describing conditions to be fit with `-` to separate if two conditions, e.g. "WT-AUXIN"
    - `result`: folder for results
    - `batchsize`: number of jobs per swarmfile, default = 1000
    - `maxtime`: maximum wall time for run, default = 2 hrs
    - `nchains`: number of MCMC chains = number of processors called by Julia, default = 2
    - 'nthreads`: number of Julia threads per processesor, default = 1
    - `transient::Bool`: true means fit transient model (T0, T30, T120)
    - `fittedparam`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
    - `fixedeffects`: tuple of vectors of rates that are fixed between control and treatment where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3
         (each vector in tuple is a fixed rate set)
    - `juliafile`: name of file to be called by julia in swarmfile
    - `root`: name of root directory for project, e.g. "scRNA\"
    - `samplesteps`: number of MCMC sampling steps
    - `warmupsteps`: number of MCMC warmup steps to find proposal distribution covariance
    - `annealsteps`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
    - `temp`: MCMC temperature
    - `tempanneal`: annealing temperature
    - `cv`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
    - `yieldprior`: prior for yield, default = .05, (set to 1 for FISH if using scRNA data format for FISH data)
    - 'priorcv`: coefficient of variation for the rate prior distributions, default is 10.


    returns swarmfile that calls a julia file that is executed on biowulf

"""

function makeswarm(;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel="label",timestamp="",nsets=1,datafolder::String="HCT116_testdata",datatype="scRNA",thresholdlow::Float64=0.,thresholdhigh::Float64=1e8,conds::String="MOCK",resultfolder::String= "fit_result",infolder=resultfolder,batchsize=1000,maxtime = 60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=40000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv = 0.02,yieldprior=0.05,priorcv= 10.,decayrate=-1.)
    if occursin.("-",conds)
        cond = string.(split(conds,"-"))
    else
        cond = conds
    end
    genes = checkgenes(root,cond,datafolder,cell,thresholdlow,thresholdhigh)
    makeswarm(genes,G=G,cell=cell,infolder=infolder,swarmfile=swarmfile,label=label,inlabel=inlabel,timestamp=timestamp,nsets=nsets,datafolder=datafolder,datatype=datatype,conds=conds,resultfolder=resultfolder,batchsize=batchsize,maxtime=maxtime,nchains=nchains,nthreads=nthreads,transient=transient,fittedparam=fittedparam,fixedeffects=fixedeffects,juliafile=juliafile,root=root,samplesteps=samplesteps,warmupsteps=warmupsteps,annealsteps=annealsteps,temp=temp,tempanneal=tempanneal,cv=cv,yieldprior=yieldprior,priorcv=priorcv,decayrate=decayrate)
end

function makeswarm(genes::Vector;G::Int=2,cell="HCT116",swarmfile::String="fit",label="label",inlabel="label",timestamp="",nsets=1,datafolder::String="HCT116_testdata",datatype="scRNA",conds::String="MOCK",resultfolder::String="fit_result",infolder=resultfolder,batchsize=1000,maxtime=60.,nchains::Int=2,nthreads::Int=1,transient::Bool=false,fittedparam=collect(1:2*G-1),fixedeffects=(),juliafile::String="fitscript",root=".",samplesteps::Int=40000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,cv=0.02,yieldprior=0.05,priorcv=10.,decayrate=-1.)
    if label == "label"
        if fish
            label = "FISH-ss"
        else
            label = "scRNA-ss"
        end
        label = label * "-" * cell
        if ~isempty(timestamp)
            label = label * "-" * timestamp
        end
        label = label * "_" * conds
        if inlabel == "label"
            inlabel = label
        end
    end
    ngenes = length(genes)
    println("number of genes: ",ngenes)
    juliafile = juliafile * "_" * label * "_" * "$G" * ".jl"
    if ngenes > batchsize
        batches = getbatches(genes,ngenes,batchsize)
        for batch in eachindex(batches)
            sfile = swarmfile * "_" * label * "_" * "$G" * "_" * "$batch" * ".swarm"
            write_swarmfile(sfile,nchains,nthreads,juliafile,batches[batch])
        end
    else
        sfile = swarmfile * "_" * label * "_" * "$G" * ".swarm"
        write_swarmfile(sfile,nchains,nthreads,juliafile,genes)
    end
    write_fitfile(juliafile,nchains,cell,conds,G,float(maxtime),fittedparam,fixedeffects,infolder,resultfolder,datafolder,datatype,inlabel,label,nsets,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root,cv,yieldprior,priorcv,decayrate)
end

"""
    fix(folder)

    Finds jobs that failed and writes a swarmfile for those genes

"""
fix(folder) = writeruns(fixruns(findjobs(folder)))

function fix_filenames(folder,old="scRNA-ss-",new="scRNA-ss_")
    files = readdir(folder)
    for file in files
        if occursin(old,file)
            nfile = replace(file,old => new)
            mv(joinpath(folder,file),joinpath(folder,nfile),force=true)
        end
    end
end

"""
    setup(rootfolder = "scRNA")

    Sets up the folder system prior to use
    Defaults to "scRNA"

"""
function rna_setup(root = ".")

    data = joinpath(root,"data")
    results = joinpath(root,"results")
    alleles = joinpath(data,"alleles")
    halflives = joinpath(data,"halflives")
    testdata = joinpath(data,"HCT116_testdata")

    if ~ispath(data)
        mkpath(data)
    end
    if ~ispath(results)
        mkpath(results)
    end
    if ~ispath(alleles)
        mkpath(alleles)
    end
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/alleles/CH12_alleles.csv","$alleles/CH12_alleles.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/alleles/HCT116_alleles.csv","$alleles/HCT116_alleles.csv")
    if ~ispath(halflives)
        mkpath(halflives)
    end
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/ESC_halflife.csv","$halflives/ESC_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/CH12_halflife.csv","$halflives/CH12_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/HCT116_halflife.csv","$halflives/HCT116_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/OcaB_halflife.csv","$halflives/OcaB_halflife.csv")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/halflives/aB_halflife.csv","$halflives/aB_halflife.csv")
    if ~ispath(testdata)
        mkpath(testdata)
    end
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/HCT116_testdata/CENPL_MOCK.txt","$testdata/CENPL_MOCK.txt")
    Downloads.download("https://raw.githubusercontent.com/nih-niddk-mbs/StochasticGene.jl/master/data/HCT116_testdata/MYC_MOCK.txt","$testdata/MYC_MOCK.txt")

end

"""
make_fitfile(fitfile,fittedparam,fixedeffects)

make the file the swarm file calls to execute julia code

"""
function write_fitfile(fitfile,nchains,cell,datacond,G,maxtime,fittedparam,fixedeffects,infolder,resultfolder,datafolder,datatype,inlabel,label,nsets,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root,cv,yieldprior,priorcv,decayrate)
        f = open(fitfile,"w")
        s =   '"'
        write(f,"@everywhere using StochasticGene\n")
        write(f,"@time fit_rna($nchains,ARGS[1],$s$cell$s,$fittedparam,$fixedeffects,$s$datacond$s,$G,$maxtime,$s$infolder$s,$s$resultfolder$s,$s$datafolder$s,$datatype,$s$inlabel$s,$s$label$s,$nsets,$cv,$transient,$samplesteps,$warmupsteps,$annealsteps,$temp,$tempanneal,$s$root$s,$yieldprior,$priorcv,$decayrate)\n")
        close(f)
end

function getbatches(genes,ngenes,batchsize)
    nbatches = div(ngenes,batchsize)
    batches = Vector{Vector{String}}(undef,nbatches+1)
    println(batchsize," ",nbatches+1)
    for i in 1:nbatches
        batches[i] = genes[batchsize*(i-1)+1:batchsize*(i)]
    end
    batches[end] = genes[batchsize*nbatches+1:end]
    return batches
end

function write_swarmfile(sfile,nchains,nthreads,juliafile,genes::Vector)
    f = open(sfile,"w")
    for gene in genes
        gene = check_genename(gene,"(")
        writedlm(f,["julia -t $nthreads -p" nchains juliafile gene])
        # writedlm(f,["julia -p" nchains juliafile nchains gene cell cond G maxtime infolder resultfolder datafolder fish inlabel label nsets runcycle transient fittedparam fixedeffects])
    end
    close(f)
end

function checkgenes(root,conds::Vector,datafolder,celltype::String,thresholdlow::Float64,thresholdhigh::Float64)
    genes = Vector{Vector{String}}(undef,2)
    if occursin.("-",datafolder)
        datafolder = string.(split(datafolder,"-"))
        for i in 1:2
            genes[i] = checkgenes(root,conds[i],datafolder[i],celltype,thresholdlow,thresholdhigh)
        end
    else
        for i in 1:2
            genes[i] = checkgenes(root,conds[i],datafolder,celltype,thresholdlow,thresholdhigh)
        end
    end
    intersect(genes[1],genes[2])
end

function checkgenes(root,cond::String,datafolder,cell::String,thresholdlow::Float64,thresholdhigh::Float64)
    datafolder = folder_path(datafolder,root,"data")
    genes = intersect(get_halflives(root,cell,thresholdlow,thresholdhigh), get_genes(cond,datafolder))
    alleles = get_alleles(root,cell)
    if ~isnothing(alleles)
        return intersect(genes,alleles)
    else
        return genes
    end
end

function folder_path(folder::String,root::String,foldertype::String)
    f = folder
    if ~ispath(folder)
        f = joinpath(root,folder)
        if ~ispath(f)
            f = joinpath(root,foldertype,folder)
            if ~ispath(f)
                throw("$folder not found")
            end
        end
    end
    f
end

function get_halflives(root,cell,thresholdlow::Float64,thresholdhigh::Float64)
    path = get_file(root,"data/halflives",cell,"csv")
    get_halflives(path,thresholdlow,thresholdhigh)
end

function get_halflives(hlpath,thresholdlow::Float64,thresholdhigh::Float64)
    genes = Vector{String}(undef,0)
    halflives = readdlm(hlpath,',')
    for row in eachrow(halflives)
        if typeof(row[2]) <: Number
            if thresholdlow <= row[2] < thresholdhigh
                push!(genes,string(row[1]))
            end
        end
    end
    return genes
end

get_alleles(root,cell) = get_alleles(get_file(root,"data/alleles",cell,"csv"))

function get_alleles(allelepath)
    if ~isnothing(allelepath)
        return readdlm(allelepath,',')[2:end,1]
    else
        return nothing
    end
end

get_file(root,folder,filetype,suffix) = get_file(joinpath(root,folder),filetype,suffix)

function get_file(folder,filetype,suffix)
    if ispath(folder)
        files = readdir(folder)
        for file in files
            name = split(file,"_")[1]
            if occursin(suffix,file) && name == filetype
                path = joinpath(folder,file)
                return path
            end
        end
    else
        throw("folder $folder does not exist")
    end
    nothing
end

function findjobs(folder)
    files = readdir(folder)
    files = files[occursin.("swarm_",files)]
    for (i,file) in enumerate(files)
        files[i] = split(file,"_")[2]
    end
    unique(files)
end

function fixruns(jobs,message="FAILED")
    runlist = Vector{String}(undef,0)
    for job in jobs
        if occursin(message,read(`jobhist $job`,String))
            swarmfile = findswarm(job)
            list = readdlm(swarmfile,',')
            runs =  chomp(read(pipeline(`jobhist $job`, `grep $message`),String))
            runs = split(runs,'\n')
            println(job)
            for run in runs
                linenumber = parse(Int,split(split(run," ")[1],"_")[2]) + 1
                while linenumber < length(list)
                    a = String(list[linenumber])
                    linenumber += 1000
                    println(a)
                    push!(runlist,a)
                end
            end
        end
    end
    return runlist
end

function writeruns(runs,outfile="fitfix.swarm")

    f = open(outfile,"w")
    for run in runs
        writedlm(f,[run],quotes=false)
    end
    close(f)

end

function findswarm(job)
    sc = "Swarm Command"
    line = read(pipeline(`jobhist $job`, `grep $sc`),String)
    list = split(line," ")
    list[occursin.(".swarm",list)][1]
end

"""
    get_missing_genes(datafolder::String,folder::String,cell,type,label,cond,model)
"""

function get_missing_genes(datafolder::String,resultfolder::String,cell,filetype,label,cond,model,root=".")
    genes = checkgenes(root,cond,datafolder,cell,0.,1e8)
    get_missing_genes(genes,folder_path(resultfolder,root,"results"),filetype,label,cond,model)
end

function get_missing_genes(genes::Vector,resultfolder,filetype,label,cond,model)
    genes1=get_genes(resultfolder,filetype,label,cond,model)
    get_missing_genes(genes,genes1)
end

get_missing_genes(genes,genes1) = union(setdiff(genes1,genes),setdiff(genes,genes1))

function scan_swarmfiles(jobid,folder=".")
    if ~(typeof(jobid) <: String)
        jobid = string(jobid)
    end
    genes = Array{String,1}(undef,0)
    files = readdir(folder)
    files = files[occursin.(jobid,files)]
    for file in files
        genes = vcat(genes,scan_swarmfile(file))
    end
    return genes
end

function scan_swarmfile(file)
    genes = Array{String,1}(undef,0)
    contents = readdlm(file,'\t')
    lines = contents[occursin.("[\"",string.(contents))]
    for line in lines
        push!(genes,split.(string(line)," ")[1])
    end
    return genes
end

function scan_fitfile(file,folder=".")
    genes = Array{String,1}(undef,0)
    joinpath(folder,file)
    file = readdlm(file,'\t')
    for line in eachrow(file)
        push!(genes,line[4])
    end
    return genes
end
