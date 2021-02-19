"""
rna.jl

Fit G state models (generalized telegraph models) to RNA abundance data
For single cell RNA (scRNA) technical noise is included as a yieldfactor
single molecule FISH (smFISH) is treated as loss less
"""


# Load data, model, and option structures
"""
transient_rna(nsets::Int,control::String,treatment::String,time::Float64,gene::String,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
Fit transient G model to time dependent mRNA data
"""
function transient_rna(control::String,treatment::String,name::String,time::Float64,gene::String,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna([control;treatment],name,time,gene)
    model = model_rna(r,G,nalleles,2,cv,fittedparam,decayprior,yieldprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function transient_rna(path,name::String,time,gene::String,nsets::Int,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,time,gene)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,yieldprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function transient_fish(path,name::String,time,gene::String,r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,time,gene,true)
    model,options = transient_fish(r,decayprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,method,warmupsteps,annealsteps)
    return data,model,options
end
function transient_fish(path,name::String,time,gene::String,r::Vector,decayprior::Float64,delayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,time,gene,true)
    model,options = transient_fish(r,decayprior,delayprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,method,warmupsteps,annealsteps)
    return data,model,options
end
function transient_fish(r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    model = model_rna(r,G,nalleles,2,cv,fittedparam,decayprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return model,options
end
function transient_fish(r::Vector,decayprior::Float64,delayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    model = model_delay_rna(r,G,nalleles,2,cv,fittedparam,decayprior,delayprior)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return model,options
end
function transient_rnafish(path,name::String,time,gene::String,nsets::Int,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,time,gene)
    model,options = transient_rnafish(r,decayprior,yieldprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,method,warmupsteps,annealsteps)
    return data,model,options
end
function transient_rnafish(r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,yieldprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return model,options
end

"""
steadystate_rna(nsets::Int,file::String,gene::String,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
Fit G model to steady state data
"""
function steadystate_rna(path,name::String,gene::String,nsets,r::Vector,decayprior::Float64,yieldprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,gene,false)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,yieldprior,0)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function steadystate_rna(path,name::String,gene::String,nsets,r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,gene,false)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,0)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function steadystate_rnafish(path,name::String,gene::String,fish::Array,r::Vector,decayprior::Float64,noisepriors::Vector,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,method=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,gene,fish)
    model = model_rna(r,G,nalleles,cv,fittedparam,decayprior,noisepriors,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function thresholds_fish(path,name::String,gene::String,r::Vector,decayprior::Float64,noisepriors::Vector,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,gene,true)
    model,options = thresholds_fish(r,decayprior,noisepriors,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,warmupsteps,annealsteps)
    return data,model,options
end
function thresholds_fish(r::Vector,decayprior::Float64,noisepriors::Vector,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
    model = model_rna(r,G,nalleles,cv,fittedparam,decayprior,noisepriors,0)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return model,options
end

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
prepare mRNA histograms
"""
function histograms_rna(path::Array,gene::String,fish::Bool)
    n = length(path)
    h = Array{Array,1}(undef,n)
    lengths = Array{Int,1}(undef,n)
    for i in eachindex(path)
        lengths[i], h[i] = histograms_rna(path[i],gene,fish[i])
    end
    return lengths,h
end
function histograms_rna(path::String,gene::String,fish::Bool)
    if fish
        h = read_fish(path,gene,.98)
    else
        h = read_scrna(path,.99)
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
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior::Float64,yieldprior::Float64,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G,nsets,fittedparam,decayprior,yieldprior)
    GMlossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
end
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior::Float64,method::Int)
    # propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r,G::Int,nsets,fittedparam,decayprior)
    GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
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
function setpriorrate(G)
Set prior distribution for mean and cv of rates
"""
function setpriorrate(G::Int,nsets::Int,decayrate::Float64,yieldfactor::Float64)
    rm,rcv = setpriorrate(G,nsets,decayrate)
    return [rm;yieldfactor],[rcv;.25]
end
function setpriorrate(G::Int,nsets::Int,decayrate::Float64)
    r0 = [.01*ones(2*(G-1));.1;decayrate]
    rc = [.25*ones(2*(G-1));.1;0.01]
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
proposal_cv_rna(propcv)
set propcv to a vector or matrix
"""
proposal_cv_rna(cv) = typeof(cv) == Float64 ? propcv*ones(length(fittedparam)) : propcv

proposal_cv_rna(propcv,fittedparam) = typeof(propcv) == Float64 ? propcv*ones(length(fittedparam)) : propcv

"""
fittedparam_rna(G,nsets,loss)
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
function read_scrna(filename::String,threshold::Float64=.98,nhistmax::Int=300)
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
        throw("data file not found")
        # return Array{Int,1}(undef,0)
    end
end

"""
read_fish(path,gene,threshold)
Read in FISH data from 7timepoint type folders

"""
function read_fish(path::String,cond::String,threshold::Float64=.98)
    xr = zeros(1000)
    for (root,dirs,files) in walkdir(path)
        for file in files
            target = joinpath(root, file)
            if occursin(cond,target) && occursin("cellular",target)
                println(target)
                x1 = readdlm(target)[:,1]
                x1 = truncate_histogram(x1,threshold,1000)
                lx = length(x1)
                xr[1:min(lx,1000)] += x1[1:min(lx,1000)]
            end
        end
    end
    return xr
end

function read_fish(path1::String,cond1::String,path2::String,cond2::String,threshold::Float64=.98)
    x1 = read_fish(path1,cond1,threshold)
    x2 = read_fish(path2,cond2,threshold)
    combine_histogram(x1,x2)
end


"""
plot_histogram_rna()

functions to plot data and model predicted histograms

"""
function plot_histogram_rna(gene::String,datapaths::Array,modelfile::String,time=[0.;30.;120.],fittedparam=[7;8;9;10;11])
    r = readrates(modelfile,1)
    data,model,_ = transient_fish(datapaths,"",time,gene,r,1.,3,2,fittedparam,1.,1.,10)
    h=likelihoodarray(r[fittedparam],data,model)
    figure(gene)
    for i in eachindex(h)
        plot(h[i])
        plot(normalize_histogram(data.histRNA[i]))
    end
    return h
end

function plot_histogram_rna(gene,cond,G,nalleles,label,folder,root)
    datapath = scRNApath(gene,cond)
    h = read_scrna(datapath,.99)
    ratepath = ratepath_Gmodel(gene,cond,G,nalleles,label,folder,root)
    println(ratepath)
    r = readrates(ratepath)
    println(r)
    plot(steady_state(r[1:2*G],r[end],G-1,length(h),nalleles))
    plot(normalize_histogram(h))
end

# functions to build paths to data and results

function scRNApath(gene,cond,folder = "data/HCT116_T0_T30_T120_Histograms/T120",root= "/Users/carsonc/Box/scrna/")
    datapath = joinpath(root,folder)
    joinpath(datapath,gene * "_" * cond * ".txt")
end

function ratepath_Gmodel(gene::String,cond::String,G::Int,nalleles::Int,label="scRNA_T120_ss_",folder="Results/2021-01-19" ,root="/Users/carsonc/Box/scrna/")
    path_Gmodel("rates",gene,G,nalleles,label * cond,folder,root)
end

function path_Gmodel(type,gene::String,G::Int,nalleles::Int,label::String,folder,root)
    filelabel = label  * "_" * gene *  "_" * "$G" * "_" * "$nalleles" * ".txt"
    ratefile = type * "_" * filelabel
    joinpath(root, joinpath(folder,ratefile))
end

function rna_stats(genes,conds)
    h = Vector{Vector}(undef,2)
    z = Vector{Float64}(undef,0)
    for gene in genes
        for i in eachindex(conds)
            datapath = scRNApath(gene,conds[i])
            h[i] = read_scrna(datapath,.99)
        end
        # push!(z,tstat_2sample(h[1],h[2]))
        # push!(z,2*(mean_histogram(h[1])-mean_histogram(h[2]))/(mean_histogram(h[1])+mean_histogram(h[2])))
        push!(z,(mean_histogram(h[1])-mean_histogram(h[2]))/(mean_histogram(h[1])))
    end
    z
end
