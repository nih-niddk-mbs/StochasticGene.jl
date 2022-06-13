# rna.jl
#
# Fit G state models (generalized telegraph models) to RNA abundance data
# For single cell RNA (scRNA) technical noise is included as a yieldfactor
# single molecule FISH (smFISH) is treated as loss less
#

"""
    fit_rna(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,fish::Bool,runcycle::Bool,inlabel,label,nsets::Int,cv=0.,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/",yieldprior = 0.05,ejectprior = 1.0)
    fit_rna(nchains::Int,data::AbstractRNAData,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,fish::Bool,runcycle,inlabel,label,nsets,cv=0.,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/",yieldprior = 0.05,ejectprior = 1.0)


Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

# Arguments
- `nchains`: number of MCMC chains
- `gene`: gene name
- `cell`: cell type
- `datacond`: condition, if more than one condition use vector of strings e.g. ["DMSO","AUXIN"]
- `maxtime`: float maximum time for entire run
- `infolder`: folder pointing to results used as initial conditions
- `resultfolder`: folder where results go
- `datafolder`: folder for data, string or array of strings
- `inlabel`: name of input files (not including gene name but including condition)
- `label`: = name of output files
- `nsets`: int number of rate sets
- `runcycle`: if true, cycle through all parameters sequentially in MCMC
- `cyclesteps`: int number of steps per cycle
- `samplesteps`: int number of samples
- `warmupsteps`: int number of warmup steps
- `annealsteps`: in number of annealing steps
- `temp`: MCMC temperature
- `tempanneal`: starting temperature for annealing
- `root`: root folder of data and Results folders
- `fittedparam`: vector of rate indices,  indices of parameters to be fit (input as string of ints separated by "-")
- `fixedeffects`: (tuple of vectors of rate indices) string indicating which rate is fixed, e.g. "eject"
- `data`: data structure

"""

function fit_rna(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,fish::Bool,runcycle::Bool,inlabel,label,nsets::Int,cv=0.,transient::Bool=false,samplesteps::Int=100000,warmupsteps=0,annealsteps=0,temp=1.,tempanneal=100.,root = ".",yieldprior::Float64=0.05,priorcv::Float64=10.,decayrate=-1.)
    gene = check_genename(gene,"[")
    datafolder = folder_path(datafolder,root,"data")
    if occursin("-",datafolder)
        datafolder = string.(split(datafolder,"-"))
    end
    if occursin("-",datacond)
        datacond = string.(split(datacond,"-"))
    end
    if transient
        data = data_rna(gene,datacond,datafolder,fish,label,["T0","T30","T120"],[0.,30.,120.])
    else
        data = data_rna(gene,datacond,datafolder,fish,label)
    end
    fit_rna(nchains,data,gene,cell,fittedparam,fixedeffects,datacond,G,maxtime,infolder,resultfolder,datafolder,fish,runcycle,inlabel,label,nsets,cv,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root,yieldprior,priorcv,decayrate)
end
function fit_rna(nchains::Int,data::AbstractRNAData,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,fish::Bool,runcycle,inlabel,label,nsets,cv=0.,transient::Bool=false,samplesteps::Int=100000,warmupsteps=0,annealsteps=0,temp=1.,tempanneal=100.,root = ".",yieldprior::Float64=0.05,priorcv::Float64=10.,decayrate=-1.)
    println(now())
    printinfo(gene,G,datacond,datafolder,infolder,resultfolder,maxtime)
    println("size of histogram: ",data.nRNA)

    resultfolder = joinpath("results",resultfolder)
    infolder = joinpath("results",infolder)

    if fish
        yieldprior = 1.
    end

    model = model_rna(gene,cell,G,cv,fittedparam,fixedeffects,inlabel,infolder,nsets,root,data,yieldprior,decayrate,LogNormal,priorcv,true)
    options = MHOptions(samplesteps,0,warmupsteps,annealsteps,maxtime*.9,temp,tempanneal)
    if runcycle
        model = cycle(nchains,fish,fixedeffects,model,data,options)
    end
    print_ll(data,model)
    fit,stats,measures = run_mh(data,model,options,nchains);
    optimized = 0
    try
        optimized = Optim.optimize(x -> loglikelihood(x,data,model)[1],fit.parml,LBFGS())
    catch
        @warn "Optimizer failed"
    end
    finalize(data,model,fit,stats,measures,temp,resultfolder,optimized,root)
    println(now())
    nothing
end

"""
cycle(nchains,fish,fixedeffects,model,data,options)

run_mh by cycling through individual parameters sequentially
"""
function cycle(nchains,fish,fixedeffects,model,data,options)
    maxtime = options.maxtime/10
    options = MHOptions(100,0,0,0,maxtime,options.temp,options.tempanneal)
    print_ll(data,model,"pre-cycle ll: ")
    t0 = time()
    nsets = length(data.nRNA)
    r = model.rates
    nalleles = model.nalleles
    G = model.G
    fittedparam = model.fittedparam
    cv = 0.02
    rateprior = model.rateprior
    while (time() - t0 < maxtime)
        for i in eachindex(fittedparam)
            model = model_rna(r,[rateprior[i]],G,nalleles,cv,[fittedparam[i]],fixedeffects,0)
            fit,_,_ = run_mh(data,model,options,nchains);
            r = get_rates(fit.parml,model)
        end
    end
    return model_rna(r,rateprior,G,nalleles,cv,fittedparam,fixedeffects,0)
end
"""
check_genename(gene,p1)

Check genename for p1
if p1 = "[" change to "("
(since swarm cannot parse "(")

"""
function check_genename(gene,p1)
    if occursin(p1,gene)
        if p1 == "["
            gene = replace(gene,"[" => "(")
            gene = replace(gene,"]" => ")")
        elseif p1 == "("
            gene = replace(gene,"(" => "]")
            gene = replace(gene,")" => "]")
        end
    end
    return gene
end


"""
print_ll(param,data,model,message="initial ll:")

compute and print initial loglikelihood
"""
function print_ll(param,data,model,message)
    ll,_ = loglikelihood(param,data,model)
    println(message,ll)
end
function print_ll(data,model,message="initial ll: ")
    ll,_ = loglikelihood(get_param(model),data,model)
    println(message,ll)
end

"""
printinfo(gene,G,cond,datafolder,infolder,resultfolder,maxtime)

print out run information
"""
function printinfo(gene,G,cond,datafolder,infolder,resultfolder,maxtime)
    println("Gene: ",gene," G: ",G," Treatment:  ",cond)
    println("data: ",datafolder)
    println("in: ", infolder," out: ",resultfolder)
    println("maxtime: ",maxtime)
end

"""
finalize(data,model,fit,stats,waic,temp,resultfolder,root)

write out run results and print out final loglikelihood and deviance
"""
function finalize(data,model,fit,stats,measures,temp,resultfolder,optimized,root)
    writefolder = joinpath(root,resultfolder)
    writeall(writefolder,fit,stats,measures,data,temp,model,optimized)
    println("final max ll: ",fit.llml)
    print_ll(vec(stats.meanparam),data,model,"mean ll: ")
    println("Mean fitted rates: ",stats.meanparam[:,1])
    println("ML rates: ",fit.parml)
    println("Acceptance: ",fit.accept,"/",fit.total)
    println("Deviance: ",deviance(fit,data,model))
    println("rhat: ",maximum(measures.rhat))
    if optimized != 0
        println("Optimized ML: ",Optim.minimum(optimized))
        println("Optimized rates: ",Optim.minimizer(optimized))
    end
end

#Prepare data structures
"""
    data_rna(gene::String,cond::String,datafolder,fish,label)
    data_rna(gene::String,cond::String,datafolder::String,fish::Bool,label,sets::Vector,time::Vector)
    data_rna(path,time,gene,time)
    data_rna(path,time,gene)

Load data structure
"""

function data_rna(gene::String,cond::String,datafolder::String,fish::Bool,label)
    if cond == "null"
        cond = ""
    end
    datafile = fish ? FISHpath(gene,cond,datafolder) : scRNApath(gene,cond,datafolder)
    data_rna(datafile,label,gene,fish)
end
function data_rna(gene::String,cond::Array,datafolder::String,fish::Bool,label)
    datafile = Array{String,1}(undef,length(cond))
    for i in eachindex(cond)
        datafile[i] = fish ? FISHpath(gene,cond[i],datafolder) : scRNApath(gene,cond[i],datafolder)
    end
    data_rna(datafile,label,gene,fish)
end
function data_rna(gene::String,cond::Array,datafolder::Array,fish::Bool,label)
    datafile = Array{String,1}(undef,length(cond))
    for i in eachindex(cond)
        datafile[i] = fish ? FISHpath(gene,cond[i],datafolder) : scRNApath(gene,cond[i],datafolder)
    end
    println(datafile)
    data_rna(datafile,label,gene,fish)
end
function data_rna(gene::String,cond::String,datafolder::String,fish::Bool,label,sets::Vector,time::Vector)
    if cond == "null"
        cond = ""
    end
    datafile =[]
    for set in sets
        folder = joinpath(datafolder,set)
        path = fish ? FISHpath(gene,cond,datafolder) : scRNApath(gene,cond,datafolder)
        datafile = vcat(datafile,path)
    end
    data_rna(datafile,label,times,gene,false)
end
function data_rna(path,name,time,gene::String,fish::Bool)
    len,h = histograms_rna(path,gene,fish)
    TransientRNAData(name,gene,len,time,h)
end
function data_rna(path,name,gene::String,fish::Bool)
    len,h = histograms_rna(path,gene,fish)
    RNAData(name,gene,len,h)
end
function data_rna(path,name,gene::String,fish::Array{Bool,1})
    len,h = histograms_rna(path,gene,fish)
    RNAMixedData(name,gene,len,fish,h)
end
"""
histograms_rna(path,gene)

prepare mRNA histograms for gene given data folder path
"""
function histograms_rna(path::Array,gene::String,fish::Bool)
    n = length(path)
    h = Array{Array,1}(undef,n)
    lengths = Array{Int,1}(undef,n)
    for i in eachindex(path)
        lengths[i], h[i] = histograms_rna(path[i],gene,fish)
    end
    return lengths,h
end
function histograms_rna(path::String,gene::String,fish::Bool)
    if fish
        h = read_fish(path,gene,.98)
    else
        h = read_scrna(path,.99)
        if length(h) == 0
            throw("data not found")
        end
    end
    return length(h),h
end
function histograms_rna(path::Array,gene::String,fish::Array{Bool,1})
    n = length(path)
    h = Array{Array,1}(undef,n)
    lengths = Array{Int,1}(undef,n)
    for i in eachindex(path)
        lengths[i], h[i] = histograms_rna(path[i],gene,fish[i])
    end
    return lengths,h
end

# Prepare model structures
"""
model_rna(gene,cell,G,fish,fittedparam,fixedeffects,inlabel,infolder,nsets,root,data,verbose=true)
model_rna(r::Vector,d::Vector,G::Int,nalleles::Int,propcv,fittedparam,fixedeffects,fish::Bool,method=0)
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,yieldprior,method)
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,method)
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,noisepriors,method)

make model structure

"""

function model_rna(gene::String,cell::String,G::Int,cv,fittedparam,fixedeffects,inlabel,infolder,nsets,root,data,yieldprior,decayrate::Float64,fprior=LogNormal,priorcv=10.,verbose=true)
    if decayrate < 0
        decayrate = get_decay(gene,cell,root)
    end
    nalleles = alleles(gene,cell,root)
    if verbose
        println("alleles: ",nalleles)
        if decayrate < 0
            throw("decayrate < 0")
        else
            println("decay rate: ",decayrate)
        end
    end
    if cv <= 0
        cv = getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root,verbose)
    end
    if G == 1
        ejectrate = mean_histogram(data.histRNA) * decayrate/nalleles * yieldprior
    else
        ejectrate = yieldprior
    end
    r = getr(gene,G,nalleles,decayrate,ejectrate,inlabel,infolder,nsets,root,verbose)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,ejectrate,fprior,priorcv)
    return model
end


function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam,fixedeffects,decayprior,ejectprior,f=LogNormal,cv=10.)
    d = prior_rna(r,G,nsets,fittedparam,decayprior,ejectprior,f,cv)
    model_rna(r::Vector,d,G::Int,nalleles,propcv,fittedparam,fixedeffects)
end
function model_rna(r::Vector,d,G::Int,nalleles,propcv,fittedparam,fixedeffects,method=0)
    if length(fixedeffects) > 0
        model = GMfixedeffectsmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,fixedeffects,method)
    else
        model = GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    end
    return model
end

# prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,ejectprior::Float64,f=LogNormal,cv::Float64 = 10.)
# function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior::Float64,method::Int)
#     d = prior_rna(r,G,nsets,fittedparam,decayprior,LogNormal)
#     GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
# end
# function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior::Float64,yieldprior::Float64,method::Int)
#     d = prior_rna(r,G,nsets,fittedparam,decayprior,yieldprior)
#     GMlossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
# end
# function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,fixedeffects::Tuple,decayprior::Float64,method::Int)
#     d = prior_rna(r,G,nsets,fittedparam,decayprior)
#     GMfixedeffectsodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,fixedeffects,method)
# end
# function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,fixedeffects::Tuple,decayprior::Float64,yieldprior::Float64,method::Int)
#     # propcv = proposal_cv_rna(propcv,fittedparam)
#     d = prior_rna(r,G,nsets,fittedparam,decayprior,yieldprior)
#     GMfixedeffectslossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,fixedeffects,method)
# end
function model_rna(r::Vector,G::Int,nalleles::Int,propcv,fittedparam::Array,decayprior::Float64,ejectprior,noisepriors::Array,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G::Int,1,fittedparam,decayprior,ejectprior,noisepriors)
    if method == 1
        GMrescaledmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    else
        GMmultimodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    end
end
function model_delay_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior,ejectprior,delayprior)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G,nsets,fittedparam,decayprior,ejectprior,delayprior)
    GMdelaymodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),Int64}(G,nalleles,r,d,propcv,fittedparam,1)
end

"""
transient_rna(nchains,gene::String,fittedparam,cond::Vector,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,samplesteps::Int=40000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

make structures for transient model fit
"""
function transient_rna(nchains,gene::String,cell,fittedparam,cond::Vector,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,samplesteps::Int=40000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    data = make_data(gene,cond,G,datafolder,label,nsets,root)
    model = make_model(gene,cell,G,fittedparam,inlabel,infolder,nsets,root,data)
    param,_ = initial_proposal(model)
    return param, data, model
end

"""
getr(gene,G,nalleles,decayrate,ejectrate,inlabel,infolder,nsets::Int,root,verbose)

"""
function getr(gene,G,nalleles,decayrate,ejectrate,inlabel,infolder,nsets::Int,root,verbose)
    r = getr(gene,G,nalleles,inlabel,infolder,root,verbose)
    if ~isnothing(r)
        if length(r) == 2*G*nsets + 1
            for n in nsets
                r[2*G*n-1] *= clamp(r[2*G*nsets + 1],eps(Float64),1-eps(Float64))
            end
            r = r[1:2*G*nsets]
        end
        if length(r) == 2*G*nsets
            if verbose
                println("init rates: ",r)
            end
            return r
        end
    end
    println("No r")
    setr(G,nsets,decayrate,ejectrate)
end

function getr(gene,G,nalleles,inlabel,infolder,root,verbose)
    ratefile = path_Gmodel("rates",gene,G,nalleles,inlabel,infolder,root)
    if verbose
        println("rate file: ",ratefile)
    end
    if isfile(ratefile)
        r = readrates(ratefile,2)
    else
        return nothing
    end
end
function getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root,verbose = true)
    paramfile = path_Gmodel("param-stats",gene,G,nalleles,inlabel,infolder,root)
    if isfile(paramfile)
        cv = read_covlogparam(paramfile)
        cv = float.(cv)
        if ~ isposdef(cv) || size(cv)[1] != length(fittedparam)
            cv = .02
        end
    else
        cv = .02
    end
    if verbose
        println("cv: ",cv)
    end
    return cv
end
"""
    get_decay(gene::String,cell::String,root::String,col::Int=2)
    get_decay(gene::String,path::String,col::Int)

    Get decay rate for gene and cell

"""
function get_decay(gene::String,cell::String,root::String,col::Int=2)
    path = get_file(root,"data/halflives",cell,"csv")
    if isnothing(path)
        println(gene," has no decay time")
        return -1.
    else
        get_decay(gene,path,col)
    end
end
function get_decay(gene::String,path::String,col::Int)
    a = nothing
    in = readdlm(path,',')
    ind = findfirst(in[:,1] .== gene)
    if ~isnothing(ind)
        a = in[ind,col]
    end
    get_decay(a,gene)
end
function get_decay(a,gene::String)
    if typeof(a) <: Number
        return get_decay(float(a))
    else
        println(gene," has no decay time")
        return -1.
    end
end
get_decay(a::Float64) = log(2)/a/60.

"""
    alleles(gene::String,cell::String,root::String,col::Int=3)
    alleles(gene::String,path::String,col::Int=3)

    Get allele number for gene and cell
"""
function alleles(gene::String,cell::String,root::String;nalleles::Int=2,col::Int=3)
    path = get_file(root,"data/alleles",cell,"csv")
    if isnothing(path)
        return 2
    else
        alleles(gene,path,nalleles=nalleles,col=col)
    end
end

function alleles(gene::String,path::String;nalleles::Int=2,col::Int=3)
    a = nothing
    in,h = readdlm(path,',',header=true)
    ind = findfirst(in[:,1] .== gene)
    if isnothing(ind)
        return nalleles
    else
        a = in[ind,col]
        return isnothing(a) ? nalleles : Int(a)
    end
end

"""
prior_rna(r::Vector,G::Int,nsets::Int,propcv,fittedparam::Array,decayprior,yieldprior)
prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,noisepriors::Array)

prepare prior distribution
r[mod(1:2*G,nsets)] = rates for each model set  (i.e. rates for each set are stacked)
r[2*G*nsets + 1] == yield factor (i.e remaining after loss due to technical noise)

prior for multithresholded smFISH data
r[1:2G] = model rates
r[2G+1] = additive noise mean
r[2G + 1 + 1:length(noisepriors)] = remaining fraction after thresholding (i.e. yield)
"""

function prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,ejectprior::Float64,f=LogNormal,cv::Float64 = 10.)
        if length(r) == 2*G*nsets
            rm,rcv = setrate(G,nsets,decayprior,ejectprior,cv)
            return distribution_array(rm[fittedparam],rcv[fittedparam],f)
        else
            throw("rates have wrong length")
        end
end
function prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,ejectprior::Float64,noisepriors::Array,f=LogNormal,cv::Float64 = 10.)
        if length(r) == 2*G * nsets + length(noisepriors)
            rm,rcv = setrate(G,nsets,decayprior,ejectprior,noisepriors,cv)
            return distribution_array(rm[fittedparam],rcv[fittedparam],f)
        else
            throw("rates have wrong length")
        end
end

"""
    setrate(G::Int,nsets::Int,decayrate::Float64,ejectrate::Float64,cv::Float64 = 10.)
    setrate(G::Int,nsets::Int,decayrate::Float64,ejectrate::Float64,noisepriors::Array,cv::Float64 = 10.)

    Set mean and cv of rates

"""

function setr(G,nsets,decayrate,ejectrate)
    r = setrate(G,nsets,decayrate,ejectrate)[1]
    println("init rates: ",r)
    return r
end

function setrate(G::Int,nsets::Int,decayrate::Float64,ejectrate::Float64,cv::Float64 = 10.)
    rm = Vector{Float64}(undef,0)
    rc = similar(rm)
    for i in 1:nsets
        rm = vcat(rm,[repeat([0.01;0.1],(G-1));ejectrate;decayrate])
        rc = vcat(rc,[cv*ones(2*(G-1));cv;0.01])
    end
    return rm,rc
end

function setrate(G::Int,nsets::Int,decayrate::Float64,ejectrate::Float64,noisepriors::Array,cv::Float64 = 10.)
    rm,rc = setrate(G,nsets,decayrate,ejectrate,cv)
    for nm in noisepriors
        rm = vcat(rm,nm)
        rc = vcat(rc,.25)
    end
    return rm,rc
end

"""
proposal_cv_rna(cv)
proposal_cv_rna(propcv,fittedparam)

set proposal cv to a vector or matrix
"""
proposal_cv_rna(cv) = typeof(cv) == Float64 ? propcv*ones(length(fittedparam)) : propcv

proposal_cv_rna(propcv,fittedparam) = typeof(propcv) == Float64 ? propcv*ones(length(fittedparam)) : propcv

"""
fittedparam_rna(G,nsets,loss)
fittedparam_rna(G,nsets)

select all parameters to be fitted except decay rates
yieldfactor is last parameter in loss models
"""
function fittedparam_rna(G,nsets,loss)
    fp = fittedparam_rna(G,nsets)
    if loss == 1
        return vcat(fp,nsets*2*G+1)
    else
        return fp
    end
end
function fittedparam_rna(G,nsets)
    nrates = 2*G  # total number of rates in GM model
    k = nrates - 1  # adjust all rate parameters except decay time
    fp = Array{Int,1}(undef,nsets*k)
    for j in 0:nsets-1
        for i in 1:k
            fp[k*j + i] = nrates*j + i
        end
    end
    return fp
end


"""
rescale_rate_rna(r,G,decayrate::Float64)

set new decay rate and rescale
transition rates such that steady state distribution is the same
"""
function rescale_rate_rna(r,G,decayrate::Float64)
    rnew = copy(r)
    if mod(length(r),2*G) == 0
        rnew *= decayrate/r[2*G]
    else
        stride = fld(length(r),fld(length(r),2*G))
        for i in 0:stride:length(r)-1
            rnew[i+1:i+2*G] *= decayrate/r[2*G]
        end
    end
    return rnew
end


# Read in data and construct histograms
"""
read_scrna(filename::String,yield::Float64=.99,nhistmax::Int=1000)
Construct mRNA count per cell histogram array of a gene
"""
function read_scrna(filename::String,threshold::Float64=.99,nhistmax::Int=500)
    if isfile(filename) && filesize(filename) > 0
        x = readdlm(filename)[:,1]
        x = truncate_histogram(x,threshold,nhistmax)
        if x == 0
            dataFISH = Array{Int,1}(undef,0)
        else
            dataFISH = x
        end
        return dataFISH
    else
        # println("data file not found")
        return Array{Int,1}(undef,0)
    end
end

"""
read_fish(path,gene,threshold)
Read in FISH data from 7timepoint type folders

"""
function read_fish(path::String,cond::String,threshold::Float64=.98,maxlength = 1800)
    xr = zeros(maxlength)
    lx = 0
    for (root,dirs,files) in walkdir(path)
        for file in files
            target = joinpath(root, file)
            if occursin(cond,target) && occursin("cellular",target)
                # println(target)
                x1 = readdlm(target)[:,1]
                x1 = truncate_histogram(x1,threshold,maxlength)
                lx = length(x1)
                # println(lx)
                xr[1:min(lx,maxlength)] += x1[1:min(lx,maxlength)]
            end
        end
    end
    return truncate_histogram(xr,1.0,maxlength)
end

function read_fish(path1::String,cond1::String,path2::String,cond2::String,threshold::Float64=.98)
    x1 = read_fish(path1,cond1,threshold)
    x2 = read_fish(path2,cond2,threshold)
    combine_histogram(x1,x2)
end

# functions to build paths to RNA/FISH data and results
"""
scRNApath(gene::String,cond::String,datapath::String,root::String)

generate path to scRNA data for given gene and condition cond
    data files must have format genename_cond.txt
"""
function scRNApath(gene::String,cond::String,datapath::String,root::String)
    datapath = joinpath(root,datapath)
    if cond == ""
        joinpath(datapath,gene * ".txt")
    else
        joinpath(datapath,gene * "_" * cond * ".txt")
    end
end
scRNApath(gene,cond,datapath) = joinpath(datapath,gene * "_" * cond * ".txt")

"""
FISHpath(gene,cond,datapath,root)

generate path to FISH data
"""
# FISHpath(gene,cond,datapath,root) = joinpath(joinpath(joinpath(root,datapath),gene),cond)
FISHpath(gene,cond,datapath,root) = joinpath(root,datapath,gene,cond)
FISHpath(gene,cond,datapath) = joinpath(datapath,gene,cond)

"""
ratepath_Gmodel(gene::String,cond::String,G::Int,nalleles::Int,label,folder,root)

"""
function ratepath_Gmodel(gene::String,cond::String,G::Int,nalleles::Int,label,folder,root)
    path_Gmodel("rates",gene,G,nalleles,label * "-" * cond,folder,root)
end

"""
path_Gmodel(type,gene::String,G::Int,nalleles::Int,label::String,folder,root)
"""
function path_Gmodel(type,gene::String,G::Int,nalleles::Int,label::String,folder,root)
    filelabel = label  * "_" * gene *  "_" * "$G" * "_" * "$nalleles" * ".txt"
    ratefile = type * "_" * filelabel
    joinpath(root,folder,ratefile)
end

"""
stats_rna(genes::Vector,conds,datapath,threshold=.99)
"""
function stats_rna(genes::Vector,conds,datapath,threshold=.99)
    g = Vector{String}(undef,0)
    z = Vector{Float64}(undef,0)
    h1 = Vector{Float64}(undef,0)
    h2 = Vector{Float64}(undef,0)
    for gene in genes
        t,m1,m2 = expression_rna(gene,conds,datapath,threshold)
        push!(g,gene)
        push!(z,t)
        push!(h1,m2)
        push!(h2,m1)
    end
    return g,z,h1,h2
end

"""
expression_rna(gene,conds::Vector,folder::String,threshold=.99)
"""
function expression_rna(gene,conds::Vector,folder::String,threshold=.99)
    h = Vector{Vector}(undef,2)
    for i in eachindex(conds)
        datapath = scRNApath(gene,conds[i],folder)
        h[i] = read_scrna(datapath,threshold)
    end
    if length(h[1]) > 0 && length(h[2]) > 0
        return log_2sample(h[1],h[2]), mean_histogram(h[1]), mean_histogram(h[2])
    else
        return 0,0,0
    end
end

"""
expression_rna(gene,conds::Vector,folder::Vector,threshold=.99)
"""
function expression_rna(gene,conds::Vector,folder::Vector,threshold=.99)
    h = Vector{Vector}(undef,2)
    for i in eachindex(conds)
        datapath = scRNApath(gene,conds[i],folder[i])
        h[i] = read_scrna(datapath,threshold)
    end
    if length(h[1]) > 0 && length(h[2]) > 0
        return log_2sample(h[1],h[2]), mean_histogram(h[1]), mean_histogram(h[2])
    else
        return 0,0,0
    end
end

"""
expression_rna(genes::Vector,cond::String,folder,threshold=.99)

"""
function expression_rna(genes::Vector,cond::String,folder,threshold=.99)
    h1 = Array{Any,2}(undef,0,2)
    for gene in genes
        datapath = scRNApath(gene,cond,folder)
        h = read_scrna(datapath,threshold)
        if length(h) > 0
            h1 = vcat(h1,[gene mean_histogram(h)])
        end
    end
    return h1
end

"""
rna_fish(gene,cond,fishfolder,rnafolder,yield,root)

output RNA histogram and downsampled FISH histogram with loss
"""
function rna_fish(gene,cond,fishfolder,rnafolder,yield)
    datarna = histograms_rna(scRNApath(gene,cond,rnafolder),gene,false)
    f=reduce_fish(gene,cond,datarna[1],fishfolder,yield)
    return datarna[2],f
end

"""
reduce_fish(gene,cond,nhist,fishfolder,yield,root)

sample fish histogram with loss (1-yield)
"""
function reduce_fish(gene,cond,nhist,fishfolder,yield)
    fish = histograms_rna(FISHpath(gene,cond,fishfolder),gene,true)
    technical_loss(fish[2],yield,nhist)
    fish[2]
end

function make_histograms(folder,file,label)
    if ~ispath(folder)
        mkpath(folder)
    end
    a,h =readdlm(file,',',header=true)
    for r in eachrow(a)
        f= open("$folder/$(r[1])_$label.txt","w")
        a = r[2:end]
        a = a[(a .!= "") .& (a.!= "NA")]
        h = make_histogram(Int.(a) .+ 1)
        writedlm(f,h)
        close(f)
    end
end

function make_histogram(r)
    nhist = maximum(r)
    h = zeros(Int,nhist+1)
    for c in r
        h[c] += 1
    end
    h
end
