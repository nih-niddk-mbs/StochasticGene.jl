"""
rna.jl

Fit G state models (generalized telegraph models) to RNA abundance data
For single cell RNA (scRNA) technical noise is included as a yieldfactor
single molecule FISH (smFISH) is treated as loss less
"""

"""
    fit_rna(nchains::Int,gene::String,cell::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=40000,cyclesteps=0,warmupsteps=0,annealsteps=0,temp=100.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fit_rna(nchains::Int,gene::String,cell::String,fittedparam::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fit_rna(nchains::Int,gene::String,cell::String,fittedparam::String,fixedeffects::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fit_rna(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fit_rna(nchains::Int,data,gene::String,decayrate::Float64,nalleles::Int,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

# Arguments
- `nchains`: number of MCMC chains
- `gene`: gene name
- `cell`: cell type
- `datacond`: condition, if more than one condition is used enter as a single string separated by hyphen, e.g. "DMSO-AUXIN"
- `maxtime`: float maximum time for entire run
- `infolder`: folder pointing to results used as initial conditions
- `resultfolder`: folder where results go
- `datafolder`: folder for data
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

"""
# function fit_rna(nchains::Int,gene::String,cell::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,datafish,inlabel,label,nsets,runcycle::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
#     fittedparam = make_fittedparam(G,nsets)
#     fit_rna(nchains,gene,cell,fittedparam,(),datacond,G,maxtime,infolder,resultfolder,datafolder,datafish,inlabel,label,nsets,runcycle,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root)
# end
# function fit_rna(nchains::Int,gene::String,cell::String,fittedparam::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,datafish::Bool,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
#     fittedparam = parse.(Int,split(fittedparam,"-"))
#     fit_rna(nchains,gene,cell,fittedparam,(),datacond,G,maxtime,infolder,resultfolder,datafolder,datafish,inlabel,label,nsets,transient,samplesteps,cyclesteps,warmupsteps,annealsteps,temp,tempanneal,root)
# end
# function fit_rna(nchains::Int,gene::String,cell::String,fittedparam::String,fixedeffects::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,datafish::Bool,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
#     fittedparam = parse.(Int,split(fittedparam,"-"))
#     fixedeffects = make_fixedeffects(fixedeffects,G,nsets)
#     fit_rna(nchains,gene,cell,fittedparam,fixedeffects,datacond,G,maxtime,infolder,resultfolder,datafolder,datafish,inlabel,label,nsets,transient,samplesteps,cyclesteps,warmupsteps,annealsteps,temp,tempanneal,root)
# end
function fit_rna(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,datafish::Bool,inlabel,label,nsets::Int,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

    gene = check_genename(gene,"[")
    # println(gene)
    # decayrate = get_decay(gene,cell,root)
    # nalleles = alleles(gene,cell,root)
    # println("alleles: ",nalleles)
    # if decayrate < 0
    #     throw("decayrate < 0")
    # else
    #     println("decay rate: ",decayrate)
    # end
    if occursin("-",datafolder)
        datafolder = string.(split(datafolder,"-"))
    end
    if occursin("-",datacond)
        datacond = string.(split(datacond,"-"))
    end
    if transient
        data = make_data(gene,datacond,datafolder,datafish,label,root,["T0","T30","T120"],[0.,30.,120.])
    else
        data = make_data(gene,datacond,datafolder,datafish,label,root)
    end
    fit_rna(nchains,data,gene,cell,fittedparam,fixedeffects,datacond,G,maxtime,infolder,resultfolder,datafolder,datafish,inlabel,label,nsets,transient,samplesteps,cyclesteps,warmupsteps,annealsteps,temp,tempanneal,root)
end
function fit_rna(nchains::Int,data::AbstractRNAData,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,datafish::Bool,inlabel,label,nsets,transient::Bool=false,samplesteps::Int=100000,cyclesteps=0,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    println(now())
    printinfo(gene,G,datacond,datafolder,infolder,resultfolder,maxtime)
    println(data.nRNA)
    yieldprior = .05
    # r = getr(gene,G,nalleles,decayrate,inlabel,infolder,nsets,root,data,datafish)
    # cv = getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root)
    # if cv != 0.02
    #     warmupsteps = 0
    # end
    if cyclesteps > 0
        cv = .02
        annealsteps = 0
        warmupsteps = div(samplesteps,5)
    end
    # options,model,param = make_structures(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
    # model = make_model(r,G,nalleles,datafish,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior)
    model = make_model(gene,cell,G,datafish,fittedparam,fixedeffects,inlabel,infolder,nsets,root,data)
    options = StochasticGene.MHOptions(samplesteps,cyclesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
    param = StochasticGene.get_param(model)

    print_ll(param,data,model)
    fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);
    finalize(data,model,fit,stats,waic,temp,resultfolder,root)
    println(now())
    nothing
end

# """
# make_structures(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
#
# return structures options, model, and param
# """
# function make_structures(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
#     if length(fixedeffects) > 0
#         model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,0)
#     else
#         model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,decayrate,yieldprior,0)
#     end
#     options = StochasticGene.MHOptions(samplesteps,cyclesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
#     param = StochasticGene.get_param(model)
#     return options,model,param
# end

"""
make_model(gene,cell,G,datafish,fittedparam,fixedeffects,inlabel,infolder,nsets,root,data,verbose=true)

make model structure

"""
function make_model(gene,cell,G,datafish,fittedparam,fixedeffects,inlabel,infolder,nsets,root,data,verbose=true)
    decayrate = get_decay(gene,cell,root)
    nalleles = alleles(gene,cell,root)
    if verbose
        println("alleles: ",nalleles)
        if decayrate < 0
            throw("decayrate < 0")
        else
            println("decay rate: ",decayrate)
        end
    end
    r = getr(gene,G,nalleles,decayrate,inlabel,infolder,nsets,root,data,verbose)
    cv = getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root,verbose)
    if datafish
        model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,decayrate,0)
    else
        yieldprior = 0.05
        if length(fixedeffects) > 0
            model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,0)
        else
            model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,decayrate,yieldprior,0)
        end
    end
    return model
end

make_model(gene,r::Vector,G,fittedparam,nsets,root) = StochasticGene.model_rna(r,G,alleles(gene,cell,root),nsets,.02,fittedparam,r[2*G],.1,0)


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
initial_ll(param,data,model,message="initial ll:")

compute and print initial loglikelihood
"""
function print_ll(param,data,model,message="initial ll: ")
    ll,_ = StochasticGene.loglikelihood(param,data,model)
    println(message,ll)
end

"""
printinfo(gene,G,cond,datafolder,infolder,resultfolder,maxtime)

print out run information
"""
function printinfo(gene,G,cond,datafolder,infolder,resultfolder,maxtime)
    println(gene," ",G," ",cond)
    println(datafolder)
    println("in: ", infolder," out: ",resultfolder)
    println(maxtime)
end

"""
finalize(data,model,fit,stats,waic,temp,resultfolder,root)

write out run results and print out final loglikelihood and deviance
"""
function finalize(data,model,fit,stats,waic,temp,resultfolder,root)
    writefile = joinpath(root,resultfolder)
    StochasticGene.writeall(writefile,fit,stats,waic,data,temp,model)
    println("final max ll: ",fit.llml)
    # println(stats.meanparam)
    print_ll(stats.meanparam,data,model,"mean ll: ")
    println(fit.accept," ",fit.total)
    println("Deviance: ",StochasticGene.deviance(fit,data,model))
    println(stats.meanparam)
end

"""
transient_rna(nchains,gene::String,fittedparam,cond::Vector,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,samplesteps::Int=40000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

make structures for transient model fit
"""
function transient_rna(nchains,gene::String,cell,fittedparam,cond::Vector,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,samplesteps::Int=40000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    data = make_data(gene,cond,G,datafolder,label,nsets,root)
    model = make_model(gene,cell,G,fittedparam,inlabel,infolder,nsets,root,data)
    param,_ = StochasticGene.initial_proposal(model)
    return param, data, model
end

"""
steadystate_rna(gene::String,fittedparam,cond,G::Int,folder::String,datafolder,label,nsets,root)
steadystate_rna(r::Vector,gene::String,fittedparam,cond,G::Int,datafolder::String,label,nsets,root)

make structures for steady state rna fit

"""
function steadystate_rna(gene::String,fittedparam,cond,G::Int,folder::String,datafolder,label,nsets,root)
    datacond = string.(split(cond,"-"))
    data = make_data(gene,datacond,datafolder,label,root)
    model = make_model(gene,G,fittedparam,label * "_" * cond,folder,nsets,root,data,false)
    param,_ = StochasticGene.initial_proposal(model)
    return param, data, model
end

function steadystate_rna(r::Vector,gene::String,fittedparam,cond,G::Int,datafolder::String,label,nsets,root)
    datacond = string.(split(cond,"-"))
    data = make_data(gene,datacond,datafolder,label,root)
    model = make_model(gene,r,G,fittedparam,nsets,root)
    param,_ = StochasticGene.initial_proposal(model)
    return param, data, model
end


# function make_fixedeffects(fixedeffects,G,nsets)
#     if fixedeffects == "eject"
#         fixedeffects = ([2*G*(nsets-1) + 2*G-1],[2*G-1])
#     elseif fixedeffects == "off"
#         fixedeffects = ([2*G*(nsets-1) + 2*G-2],[2*G-2])
#     elseif fixedeffects == "on"
#         fixedeffects = ([2*G*(nsets-1) + 2*G-3],[2*G-3])
#     elseif fixedeffects == "offeject" || fixedeffects == "ejectoff"
#         fixedeffects = ([2*G*(nsets-1) + 2*G-2,2*G*(nsets-1) + 2*G-1],[2*G-2,2*G-1])
#     elseif fixedeffects == "oneject" || fixedeffects == "ejecton"
#         fixedeffects = ([2*G*(nsets-1) + 2*G-3,2*G*(nsets-1) + 2*G-1],[2*G-3,2*G-1])
#     elseif fixedeffects == "onoff" || fixedeffects == "offon"
#         fixedeffects = ([2*G*(nsets-1) + 2*G-3,2*G*(nsets-1) + 2*G-2],[2*G-3,2*G-2])
#     end
#     return fixedeffects
# end
# function make_fittedparam(G::Int,nsets)
#
#     if nsets == 1
#         if G == 3
#             fittedparam = [1,2,3,4,5,7]
#         elseif G == 2
#             fittedparam = [1,2,3,5]
#         elseif G == 1
#             fittedparam = [1,3]
#         end
#     else
#         if G == 3
#             # fittedparam = [1,2,3,4,5,7,8,9,10,11,13]
#             # fittedparam = [7,8,9,10,11,13]
#             fittedparam = [7,8,9,10,11]
#         elseif G == 2
#             # fittedparam = [1,2,3,5,6,7,9]
#             # fittedparam = [5,6,7,9]
#             fittedparam = [5,6,7]
#         elseif G == 1
#             fittedparam = [3,5]
#         end
#     end
#     return fittedparam
# end

function make_data(gene::String,cond::String,datafolder,fish,label,root)
    if cond == "null"
        cond = ""
    end
    datafile = fish ? StochasticGene.FISHpath(gene,cond,datafolder,root) : StochasticGene.scRNApath(gene,cond,datafolder,root)
    StochasticGene.data_rna(datafile,label,gene,fish)
end

function make_data(gene::String,cond::Array,datafolder::String,fish,label,root)
    datafile = Array{String,1}(undef,length(cond))
    for i in eachindex(cond)
        datafile[i] = fish ? StochasticGene.FISHpath(gene,cond[i],datafolder,root) : StochasticGene.scRNApath(gene,cond[i],datafolder,root)
    end
    StochasticGene.data_rna(datafile,label,gene,fish)
end

function make_data(gene::String,cond::Array,datafolder::Array,fish,label,root)
    datafile = Array{String,1}(undef,length(cond))
    for i in eachindex(cond)
        datafile[i] = fish ? StochasticGene.FISHpath(gene,cond[i],datafolder,root) : StochasticGene.scRNApath(gene,cond[i],datafolder,root)
    end
    println(datafile)
    StochasticGene.data_rna(datafile,label,gene,fish)
end

function make_data(gene::String,cond::String,datafolder,fish,label,root,sets::Vector,time::Vector)
    if cond == "null"
        cond = ""
    end
    datafile =[]
    for set in sets
        folder = joinpath(datafolder,set)
        path = fish ? StochasticGene.FISHpath(gene,cond,datafolder,root) : StochasticGene.scRNApath(gene,cond,datafolder,root)
        datafile = vcat(datafile,path)
    end
    StochasticGene.data_rna(datafile,label,times,gene,false)
end

function getr(gene,G,nalleles,decayrate,inlabel,infolder,nsets::Int,root,data,fish::Bool,verbose=true)
    ratefile = StochasticGene.path_Gmodel("rates",gene,G,nalleles,inlabel,infolder,root)
    if verbose
        println(ratefile)
    end
    if isfile(ratefile)
        r = StochasticGene.readrates(ratefile,2)
        if ~fish
            r[end] = clamp(r[end],eps(Float64),1-eps(Float64))
            if length(r) == 2*G*nsets + 1
                if verbose
                    println(r)
                end
                return r
            end
        else
            if length(r) == 2*G*nsets
                if verbose
                    println(r)
                end
                return r
            end
        end
    end
    println("No r")
    setr(G,decayrate,nsets,data,fish)
end

function setr(G,decayrate,nsets,data,fish)
    if G == 2
        r = [0.015,0.015,0.5,.01,1.]*decayrate/.01
    elseif G == 3
        r = [0.015,.2,.2,0.015,1.5,.01,1.]*decayrate/.01
    elseif G == 1
        if typeof(data.nRNA) <: Vector
            r = Array{Float64,1}(undef,0)
            for hist in data.histRNA
                mu = StochasticGene.mean_histogram(hist)
                r = vcat([10*mu,1.],r)
            end
            r *= decayrate
            r = [r;1.]
            nsets = 1
        else
            mu=StochasticGene.mean_histogram(data.histRNA)
            r = [10*mu,1.,1.]*decayrate
        end
    end
    if nsets > 1
        r = [r[1:end-1];r]
    end
    if fish
        r = r[1:end-1]
    else
        r[end] = .05
    end
    println(r)
    return r
end

function getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root,verbose = true)
    paramfile = StochasticGene.path_Gmodel("param_stats",gene,G,nalleles,inlabel,infolder,root)
    if isfile(paramfile)
        cv = StochasticGene.read_covlogparam(paramfile)
        cv = float.(cv)
        if ~StochasticGene.isposdef(cv) || size(cv)[1] != length(fittedparam)
            cv = .02
        end
    else
        cv = .02
    end
    if verbose
        println(cv)
    end
    return cv
end

function get_decay(gene,cell,root,col=2)
    folder = joinpath(root,"data/halflives")
    files = readdir(folder)
    a = nothing
    for file in files
        if occursin(cell,file)
            path = joinpath(folder,file)
            in = readdlm(path,',')
            ind = findfirst(in[:,1] .== gene)
            if ~isnothing(ind)
                a = in[ind,col]
            end
        end
    end
    get_decay(a,gene)
end
get_decay(a::Float64) = log(2)/a/60.
function get_decay(a,gene)
    if typeof(a) <: Number
        return get_decay(a)
    else
        println(gene," has no decay time")
        return -1.
    end
end

function alleles(gene,cell,root,col=3)
    folder = joinpath(root,"data/alleles")
    files = readdir(folder)
    a = nothing
    for file in files
        if occursin(cell,file)
            path = joinpath(folder,file)
            in = readdlm(path)
            ind = findfirst(in[:,1] .== gene)
            if ~isnothing(ind)
                a = in[ind,col]
            end
        end
    end
    if isnothing(a)
        return 2
    else
        return Int(a)
    end
end

"""
rna_fish(gene,cond,fishfolder,rnafolder,yield,root)

output RNA histogram and downsampled FISH histogram with loss
"""
function rna_fish(gene,cond,fishfolder,rnafolder,yield,root)
    datarna = StochasticGene.histograms_rna(StochasticGene.scRNApath(gene,cond,rnafolder,root),gene,false)
    f=reduce_fish(gene,cond,datarna[1],fishfolder,yield,root)
    return datarna[2],f
end

"""
reduce_fish(gene,cond,nhist,fishfolder,yield,root)

sample fish histogram with loss (1-yield)
"""
function reduce_fish(gene,cond,nhist,fishfolder,yield,root)
    datafish = StochasticGene.histograms_rna(StochasticGene.FISHpath(gene,cond,fishfolder,root),gene,true)
    StochasticGene.technical_loss(datafish[2],yield,nhist)
end

# Load data, model, and option structures
# """
# transient_rna(nsets::Int,control::String,treatment::String,time::Float64,gene::String,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
# Fit transient G model to time dependent mRNA data
# """
#
# function transient_rna(path,name::String,time,gene::String,nsets::Int,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,time,gene,false)
#     model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,yieldprior,method)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return data,model,options
# end
# function transient_fish(path,name::String,time,gene::String,r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,time,gene,true)
#     model,options = transient_fish(r,decayprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,method,warmupsteps,annealsteps)
#     return data,model,options
# end
# function transient_fish(path,name::String,time,gene::String,r::Vector,decayprior::Float64,delayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,time,gene,true)
#     model,options = transient_fish(r,decayprior,delayprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,method,warmupsteps,annealsteps)
#     return data,model,options
# end
# function transient_fish(r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     model = model_rna(r,G,nalleles,2,cv,fittedparam,decayprior,method)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return model,options
# end
# function transient_fish(r::Vector,decayprior::Float64,delayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     model = model_delay_rna(r,G,nalleles,2,cv,fittedparam,decayprior,delayprior)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return model,options
# end
# function transient_rnafish(path,name::String,time,gene::String,nsets::Int,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,time,gene)
#     model,options = transient_rnafish(r,decayprior,yieldprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,method,warmupsteps,annealsteps)
#     return data,model,options
# end
# function transient_rnafish(r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
#     model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,yieldprior,method)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return model,options
# end
#
# """
# steadystate_rna(nsets::Int,file::String,gene::String,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
# Fit G model to steady state data
# """
# function steadystate_rna(path,name::String,gene::String,nsets,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,gene,false)
#     model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,yieldprior,0)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return data,model,options
# end
# function steadystate_rna(path,name::String,gene::String,nsets,r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,gene,false)
#     model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,0)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return data,model,options
# end
# function steadystate_rnafish(path,name::String,gene::String,fish::Array,r::Vector,decayprior::Float64,noisepriors::Vector,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,method=1,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,gene,fish)
#     model = model_rna(r,G,nalleles,cv,fittedparam,decayprior,noisepriors,method)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return data,model,options
# end
# function thresholds_fish(path,name::String,gene::String,r::Vector,decayprior::Float64,noisepriors::Vector,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
#     data = data_rna(path,name,gene,true)
#     model,options = thresholds_fish(r,decayprior,noisepriors,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,warmupsteps,annealsteps)
#     return data,model,options
# end
# function thresholds_fish(r::Vector,decayprior::Float64,noisepriors::Vector,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
#     model = model_rna(r,G,nalleles,cv,fittedparam,decayprior,noisepriors,0)
#     options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,temp)
#     return model,options
# end

#Prepare data structures
"""
data_rna(path,time,gene,time)
data_rna(path,time,gene)

Load data structure
"""
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
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,yieldprior,method)
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,method)
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,noisepriors,method)

Load model structure
"""
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior::Float64,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G::Int,nsets,fittedparam,decayprior)
    GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
end
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior::Float64,yieldprior::Float64,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G,nsets,fittedparam,decayprior,yieldprior)
    GMlossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
end
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,randomeffects::Tuple,decayprior::Float64,yieldprior::Float64,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G,nsets,fittedparam,decayprior,yieldprior)
    GMfixedeffectslossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,randomeffects,method)
end
function model_rna(r::Vector,G::Int,nalleles::Int,propcv,fittedparam::Array,decayprior::Float64,noisepriors::Array,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G::Int,1,fittedparam,decayprior,noisepriors)
    if method == 1
        GMrescaledmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    else
        GMmultimodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    end
end
function model_delay_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior,delayprior)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G,nsets,fittedparam,decayprior,delayprior)
    GMdelaymodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),Int64}(G,nalleles,r,d,propcv,fittedparam,1)
end

"""
prior_rna(r::Vector,G::Int,nsets::Int,propcv,fittedparam::Array,decayprior,yieldprior)

prepare prior distribution
r[mod(1:2*G,nsets)] = rates for each model set  (i.e. rates for each set are stacked)
r[2*G*nsets + 1] == yield factor (i.e remaining after loss due to technical noise)
"""
function prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,yieldprior::Float64,f=Gamma)
        ind = 2*G * nsets + 1
        if length(r) == ind
            rm,rcv = setpriorrate(G,nsets,decayprior,yieldprior)
            if in(ind,fittedparam)
                return distributionBeta_array(rm[fittedparam],rcv[fittedparam],findfirst(ind.==fittedparam),f)
            else
                return distribution_array(rm[fittedparam],rcv[fittedparam],f)
            end
        else
            throw("rates have wrong length")
        end
end
function prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,f=Gamma)
        if length(r) == 2*G*nsets
            rm,rcv = setpriorrate(G,nsets,decayprior)
            return distribution_array(rm[fittedparam],rcv[fittedparam],f)
        else
            throw("rates have wrong length")
        end
end
"""
prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,noisepriors::Array)

prior for multithresholded smFISH data
r[1:2G] = model rates
r[2G+1] = additive noise mean
r[2G + 1 + 1:length(noisepriors)] = remaining fraction after thresholding (i.e. yield)
"""
function prior_rna(r::Vector,G::Int,nsets::Int,fittedparam::Array,decayprior::Float64,noisepriors::Array)
        if length(r) == 2*G * nsets + length(noisepriors)
            rm,rcv = setpriorrate(G,nsets,decayprior,noisepriors)
            return distribution_array(rm[fittedparam],rcv[fittedparam])
        else
            throw("rates have wrong length")
        end
end
"""
setpriorrate(G::Int,nsets::Int,decayrate::Float64,yieldfactor::Float64)
setpriorrate(G::Int,nsets::Int,decayrate::Float64)
setpriorrate(G::Int,nsets::Int,decayrate::Float64,noisepriors::Array)

Set prior distribution for mean and cv of rates
"""
function setpriorrate(G::Int,nsets::Int,decayrate::Float64,yieldfactor::Float64)
    rm,rcv = setpriorrate(G,nsets,decayrate)
    return [rm;yieldfactor],[rcv;.25]
end
function setpriorrate(G::Int,nsets::Int,decayrate::Float64)
    r0 = [.01*ones(2*(G-1));.1;decayrate]
    rc = [.1*ones(2*(G-1));.01;0.01]
    rm = r0
    rcv = rc
    for i in 2:nsets
        rm = vcat(rm,r0)
        rcv = vcat(rcv,rc)
    end
    return rm,rcv
end
function setpriorrate(G::Int,nsets::Int,decayrate::Float64,noisepriors::Array)
    rm,rcv = setpriorrate(G,nsets,decayrate)
    for nm in noisepriors
        rm = vcat(rm,nm)
        rcv = vcat(rcv,.25)
    end
    return rm,rcv
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
function read_fish(path::String,cond::String,threshold::Float64=.98)
    xr = zeros(1000)
    lx = 0
    for (root,dirs,files) in walkdir(path)
        for file in files
            target = joinpath(root, file)
            if occursin(cond,target) && occursin("cellular",target)
                # println(target)
                x1 = readdlm(target)[:,1]
                x1 = truncate_histogram(x1,threshold,1000)
                lx = length(x1)
                # println(lx)
                xr[1:min(lx,1000)] += x1[1:min(lx,1000)]
            end
        end
    end
    return truncate_histogram(xr,1.0,1000)
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
FISHpath(gene,cond,datapath,root) = joinpath(joinpath(joinpath(root,datapath),gene),cond)

"""
ratepath_Gmodel(gene::String,cond::String,G::Int,nalleles::Int,label,folder,root)

"""
function ratepath_Gmodel(gene::String,cond::String,G::Int,nalleles::Int,label,folder,root)
    path_Gmodel("rates",gene,G,nalleles,label * "_" * cond,folder,root)
end

"""
path_Gmodel(type,gene::String,G::Int,nalleles::Int,label::String,folder,root)
"""
function path_Gmodel(type,gene::String,G::Int,nalleles::Int,label::String,folder,root)
    filelabel = label  * "_" * gene *  "_" * "$G" * "_" * "$nalleles" * ".txt"
    ratefile = type * "_" * filelabel
    joinpath(root, joinpath(folder,ratefile))
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
