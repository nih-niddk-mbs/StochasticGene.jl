"""
rna.jl

Fit G state models (generalized telegraph models) to RNA abundance data
For single cell RNA (scRNA) technical noise is included as a lossfactor
single molecule FISH (smFISH) is treated as loss less
"""

# Functions necessary for metropolis_hastings.jl
"""
datahistogram(data)
Return the RNA histogram data as one vector
"""
function datahistogram(data::TransientRNAData)
    v = data.histRNA[1]
    for i in 2:length(data.histRNA)
        v = vcat(v,data.histRNA[i])
    end
    return v
end
datahistogram(data::RNAData) = data.histRNA
"""
logprior(param,model::AbstractGMmodel)
compute log of the prior
"""
function logprior(param,model::AbstractGMmodel)
    d = model.rateprior
    p=0
    for i in eachindex(d)
        p -= logpdf(d[i],param[i])
    end
    return p
end
"""
likelihoodfn(param,data,model)
model likelihoodfn
"""
function likelihoodfn(param,data::RNAData,model::GMmodel)
    r = get_rates(param,model)
    n = model.G-1
    steady_state(r[1:2*n+2],n,data.nRNA,model.nalleles)
end
function likelihoodfn(param,data::RNAData,model::GMLossmodel)
    r = get_rates(param,model)
    lossfactor = r[end]
    n = model.G-1
    steady_state(r[1:2*n+2],lossfactor,n,data.nRNA,model.nalleles)
end
function likelihoodfn(param,data::TransientRNAData,model::AbstractGMmodel)
    h1,h2 = likelihoodtuple(param,data,model)
    return [h1;h2]
end
function likelihoodfn(param,data::TransientRNAData,model::AbstractGMmodel)
    h = likelihoodarray(param,data,model)
    hconcat = Array{Float64,1}(undef,0)
    for h in h
        hconcat = vcat(hconcat,h)
    end
    return hconcat
end

# Load data, model, and option structures
"""
transient_rna(nsets::Int,control::String,treatment::String,time::Float64,gene::String,r::Vector,decayprior::Float64,lossprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
Fit transient G model to time dependent mRNA data
"""
function transient_rna(control::String,treatment::String,name::String,time::Float64,gene::String,r::Vector,decayprior::Float64,lossprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna([control;treatment],name,time,gene)
    model = model_rna(r,G,nalleles,2,cv,fittedparam,decayprior,lossprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function transient_rna(path,name::String,time,gene::String,nsets::Int,r::Vector,decayprior::Float64,lossprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,time,gene)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,lossprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function transient_fish(path,name::String,time,gene::String,r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp::Float64=10.,method::Int=1,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,time,gene,true)
    model = model_rna(r,G,nalleles,2,cv,fittedparam,decayprior,method)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
"""
steadystate_rna(nsets::Int,file::String,gene::String,r::Vector,decayprior::Float64,lossprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
Fit G model to steady state data
"""
function steadystate_rna(path,name::String,gene::String,nsets,r::Vector,decayprior::Float64,lossprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,gene)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,lossprior,0)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end
function steadystate_rna(path,name::String,gene::String,nsets,r::Vector,decayprior::Float64,G::Int,nalleles::Int,fittedparam::Vector,cv,maxtime::Float64,samplesteps::Int,temp=10.,warmupsteps=0,annealsteps=0)
    data = data_rna(path,name,gene,true)
    model = model_rna(r,G,nalleles,nsets,cv,fittedparam,decayprior,0)
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end

#Prepare data structures
"""
data_rna(path,time,gene,time)
data_rna(path,time,gene)
Load data structure
"""
function data_rna(path,name,time,gene::String,fish::Bool=false)
    len,h = histograms_rna(path,gene,fish)
    TransientRNAData(name,gene,len,time,h)
end
function data_rna(path,name,gene::String,fish::Bool=false)
    len,h = histograms_rna(path,gene,fish)
    RNAData(name,gene,len,h)
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
        lengths[i], h[i] = histograms_rna(path[i],gene,fish)
    end
    return lengths,h
end
function histograms_rna(path::String,gene::String,fish::Bool)
    if fish
        h = read_fish(path,gene)
    else
        h = read_scrna(joinpath(path , gene * txtstr))
    end
    return length(h),h
end
# Prepare model structures
"""
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,lossprior,method)
model_rna(r,G,nalleles,nsets,propcv,fittedparam,decayprior,method)
Load model structure
"""
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior,lossprior,method)
    propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r::Vector,G::Int,nsets::Int,decayprior::Float64,fittedparam::Array,lossprior::Float64)
    GMLossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
end
function model_rna(r::Vector,G::Int,nalleles::Int,nsets::Int,propcv,fittedparam::Array,decayprior,method)
    propcv = proposal_cv_rna(propcv,fittedparam)
    d = prior_rna(r::Vector,G::Int,nsets::Int,decayprior::Float64,fittedparam::Array)
    GMmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
end

"""
proposal_cv_rna(propcv)
set propcv to a vector or matrix
"""
proposal_cv_rna(cv) = typeof(cv) == Float64 ? propcv*ones(length(fittedparam)) : propcv
"""
prior_rna(r::Vector,G::Int,nsets::Int,propcv,fittedparam::Array,decayprior,lossprior)
compute prior distribution
"""
function prior_rna(r::Vector,G::Int,nsets::Int,decayprior::Float64,fittedparam::Array,lossprior::Float64)
        if length(r) == 2*G * nsets + 1
            rm,rcv = setpriorrate(G,nsets,decayprior,lossprior)
            return priorLogNormal(rm[fittedparam],rcv[fittedparam])
        else
            throw("rates have wrong length")
        end
end
function prior_rna(r::Vector,G::Int,nsets::Int,decayprior::Float64,fittedparam::Array)
        if length(r) == 2*G*nsets
            rm,rcv = setpriorrate(G,nsets,decayprior)
            return priorLogNormal(rm[fittedparam],rcv[fittedparam])
        else
            throw("rates have wrong length")
        end
end
"""
proposal_cv_rna(propcv)
set propcv to a vector or matrix
"""
proposal_cv_rna(propcv,fittedparam) = typeof(propcv) == Float64 ? propcv*ones(length(fittedparam)) : propcv
"""
prior_rna(r::Vector,G::Int,nsets::Int,propcv,fittedparam::Array,decayprior,lossprior)
compute prior distribution
"""
function prior_rna(r::Vector,G::Int,nsets::Int,decayprior::Float64,lossprior::Float64,fittedparam::Array)
        if length(r) == 2*G * nsets + 1
            rm,rcv = setpriorrate(G,nsets,decayprior,lossprior)
            return priorLogNormal(rm[fittedparam],rcv[fittedparam])
        else
            throw("rates have wrong length")
        end
end
function prior_rna(r::Vector,G::Int,nsets::Int,decayprior::Float64,fittedparam::Array)
        if length(r) == 2*G*nsets
            rm,rcv = setpriorrate(G,nsets,decayprior)
            return priorLogNormal(rm[fittedparam],rcv[fittedparam])
        else
            throw("rates have wrong length")
        end
end

"""
fittedparam_rna(G,nsets,loss)
select all parameters to be fitted except decay rates
lossfactor is last parameter in loss models
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
get_rates(param,model::AbstractGMmodel)
replace fitted rates with new values and return
"""
function get_rates(param,model::AbstractGMmodel)
    r = copy(model.rates)
    r[model.fittedparam] = param
    return r
end
"""
priorLogNormal(r,cv,G,R)
LogNormal Prior distribution
"""
function priorLogNormal(param,cv)
    sigma = sigmalognormal(cv)
    d = []
    for i in eachindex(param)
        push!(d,Distributions.LogNormal(log(param[i]),sigma[i]))
    end
    return d
end

"""
function setpriorrate(G)
Set prior distribution for mean and cv of rates
"""
function setpriorrate(G,nsets,decayrate,lossrate)
    rm,rcv = setpriorrate(G,nsets,decayrate)
    return [rm;.1],[rcv;.025]
end
function setpriorrate(G,nsets,decayrate)
    r0 = [.05*ones(2*(G-1));.4;decayrate]
    rc = [ones(2*(G-1));.1;0.2]
    rm = r0
    rcv = rc
    for i in 2:nsets
        rm = vcat(rm,r0)
        rcv = vcat(rcv,rc)
    end
    return rm,rcv
end


"""
likelihoodarray(param,data,model::AbstractGmodel)

Compute time dependent GM model likelihoods
first set of parameters gives the initial histogram
2nd set gives the new parameters at time 0
data.histRNA holds array of histograms for time points given by data.time
transient computes the time evolution of the histogram
model.method=1 specifies finite difference solution otherwise use eigendecomposition solution,
"""
function likelihoodarray(param,data::TransientRNAData,model::GMLossmodel)
    lossfactor = get_rates(param,model)[end]
    h = likelihoodarray(param,data::TransientRNAData,model,maximum(data.nRNA))
    noise_convolve!(h,lossfactor)
    trim(h,data.nRNA)
end
function likelihoodarray(param,data::TransientRNAData,model::GMmodel)
    h=likelihoodarray(param,data::TransientRNAData,model::GMmodel,maximum(data.nRNA))
    trim(h,data.nRNA)
end
function likelihoodarray(param,data::TransientRNAData,model::AbstractGMmodel,nRNAmax)
    nRNA = data.nRNA
    r = get_rates(param,model)
    G = model.G
    h0 = steady_state_full(r[1:2*G],G-1,nRNAmax)
    transient(data.time,r[2*G+1:4*G],G-1,model.nalleles,h0,model.method)
end

function trim(h::Array,nh::Array)
    for i in eachindex(h)
        h[i] = h[i][1:nh[i]]
    end
    return h
end

# Read in data and construct histograms
"""
readRNA(filename::String,yield::Float64=.99,nhistmax::Int=1000)
Construct mRNA count per cell histogram array of a gene
"""
function read_scrna(filename::String,yield::Float64=.99,nhistmax::Int=1000)
    if isfile(filename) && filesize(filename) > 0
        x = readdlm(filename)[:,1]
        x = round.(Int,x)
        # x = x[x[:,1] .!= "",:] .* counts[i]
        x = truncate_histogram(x,yield,nhistmax)
        if x == 0
            dataFISH = Array{Int,1}(undef,0)
        else
            dataFISH = x
        end
        return dataFISH
    else
        return Array{Int,1}(undef,0)
    end
end

"""
readRNA(control::String,treatment::String,yield::Float64=.99,nhistmax::Int=1000)
Construct 2D array of mRNA count per cell histograms for control and treatment of a gene
"""
function read_scrna(control::String,treatment::String,yield::Float64=.95,nhistmax::Int=1000)
    x = readRNA_scrna(control)[:,1]
    y = readRNA_scrna(treatment)[:,1]
    n = min(length(x),length(y))
    z = Array{Array{Int,1},1}(undef,2)
    z[1] = x[1:n]
    z[2] = y[1:n]
    return z
end

"""
readRNAFISH(scRNAfolder::String,FISHfolder::String,genein::String,cond::String)
Create a 2D data array of mRNA count/cell histograms for FISH and scRNA for same gene
"""
function read_fish_scrna(scRNAfolder::String,FISHfolder::String,genein::String,cond::String)
    histfile = "/cellular RNA histogram.csv"
    data = Array{Array{Int,1},1}(undef,2)
    scRNAfile = scRNAfolder * genein * ".txt"
    data[1] = readRNA_scrna(scRNAfile)
    rep = ["rep1/","rep2/"]
    x = Array{Array{Float64,1},1}(undef,2)
    for r in eachindex(rep)
        repfolder = FISHfolder * genein * "/" * cond * "/" * rep[r]
        FISHfiles = readdir(repfolder)
        FISH = FISHfiles[.~occursin.(".",FISHfiles)]
        xr = zeros(1000)
        for folder in FISH
            x1 = readdlm(repfolder * folder * histfile)[:,1]
            lx = length(x1)
            xr[1:min(lx,1000)] += x1[1:min(lx,1000)]
        end
        xr /= length(FISH)
        xr = round.(Int,xr)
        x[r] = truncate_histogram(xr,.99,1000)
    end
    l = min(length(x[1]),length(x[2]))
    data[2] = x[1][1:l] + x[2][1:l]
    return data
end
"""
read_fish(path,gene,threshold)
Read in FISH data from 7timepoint type folders

"""
function read_fish(path::String,gene::String,threshold=.99)
    folder = joinpath(path,"rep1" * "/" * gene)
    xr = zeros(1000)
    for (root,dirs,files) in walkdir(folder)
        for file in files
            x1 = readdlm(joinpath(root,file))[:,1]
            lx = length(x1)
            xr[1:min(lx,1000)] += x1[1:min(lx,1000)]
        end
    end
    truncate_histogram(xr,threshold,1000)
end
"""
read_rates_scrna(infile,rinchar,inpath)
read in saved rates
"""
function read_rates_scrna(infile::String,rinchar::String,inpath="/Users/carsonc/Box/scrna/Results/")
    infile = inpath * infile
    if isfile(infile) && ~isempty(read(infile))
        readdlm(infile)[:,1]
    else
        return 0
    end
end

function read_rates_scrna(infile::String,rinchar::String,gene::String,inpath="/Users/carsonc/Box/scrna/Results/")
    if rinchar == "rml" || rinchar == "rlast"
        rskip = 1
        rstart = 1
    elseif rinchar == "rmean"
        rskip = 2
        rstart = 1
    elseif rinchar == "rquant"
        rskip = 3
        rstart = 2
    end
    if isfile(infile) && ~isempty(read(infile))
        rall = readdlm(infile)
        rind = findfirst(rall[:,1] .== gene)
        if ~isnothing(rind)
            return convert(Array{Float64,1},rall[rind,rstart:rskip:end])
        else
            return 0
        end
    else
        return 0
    end
end
