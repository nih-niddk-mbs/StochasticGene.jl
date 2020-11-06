"""
rna.jl

Fit G state models (generalized telegraph models) to RNA abundance data
For single cell RNA (scRNA) technical noise is included as a lossfactor
single molecule FISH (smFISH) is treated as loss less
"""

function rna_transient(control,treatment,time,gene::String,r::Vector,G::Int,nalleles::Int,cv::Float64,maxtime::Float64,samplesteps::Int,temp=1000,annealsteps=0,warmupsteps=0)
    data = data_rna(control,treatment,time,gene)
    nsets = 2
    model = model_rna(r,G,nalleles,nsets,cv*ones(nsets*(2*G-1)+1))
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end

function rna_steadystate(control::String,gene::String,r::Vector,G::Int,nalleles::Int,cv::Float64,maxtime::Float64,samplesteps::Int,temp=1,annealsteps=0,warmupsteps=0)
    data = data_rna(control,gene)
    nsets = 1
    model = model_rna(r,G,nalleles,nsets,cv*ones(nsets*(2*G-1)+1))
    options = MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp)
    return data,model,options
end

function data_rna(control::String,treatment::String,time,gene::String)
    h = Array{Array,1}(undef,2)
    h[1] = readRNA_rna(control * gene * txtstr)
    h[2] = readRNA_rna(treatment * gene * txtstr)
    TransientRNAData(gene,[length(h[1]);length(h[2])],time,h)
end

function data_rna(control::Array,treatment::Array,time,gene::String)
    h = Array{Array,1}(undef,2)
    h[1] = readRNA_rna(control * gene * txtstr)
    h[2] = readRNA_rna(treatment * gene * txtstr)
    TransientRNAData(gene,[length(h[1]);length(h[2])],time,h)
end

function data_rna(control::String,gene::String)
    h = readRNA_rna(control * gene * txtstr)
    RNAData(gene,length(h),h)
end

function model_rna(r::Vector,G::Int,nalleles,nsets,propcv,fittedparam=collect(1:2*G-1),method=0,Gprior=(.01,.5),Sprior=(.1,10.),Rcv=1e-2)
    if
    rm,rcv = setpriorrate(G,nsets)
    fittedparam = fittedparam_rna(G,1,nsets)
    d = priordistributionLogNormal(rm[fittedparam],rcv[fittedparam],G)
    if r == 2*G * nsets + 1
        return GMLossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    elseif r == 2*G * nsets
        return GMLossmodel{typeof(r),typeof(d),typeof(propcv),typeof(fittedparam),typeof(method)}(G,nalleles,r,d,propcv,fittedparam,method)
    end


end

function fittedparam_rna(G,lossfactor,nsets)
    nrates = 2*G  # total number of rates in GM model
    k = nrates - 1  # adjust all rate parameters except decay time
    fp = Array{Int,1}(undef,nsets*k+1)
    for j in 0:nsets-1
        for i in 1:k
            fp[k*j + i] = nrates*j + i
        end
    end
    fp[nsets*k+1] = nsets*nrates+1
    return fp
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

function get_rates(param,model::AbstractGMmodel)
    r = copy(model.rates)
    r[model.fittedparam] = param
    return r
end

function logprior(param,model::AbstractGMmodel)
    d = model.rateprior
    p=0
    for i in eachindex(d)
        p -= logpdf(d[i],param[i])
    end
    return p
end

"""
priordistributionLogNormal(r,cv,G,R)
LogNormal Prior distribution
"""
function priorLogNormal(r,cv,G)
    sigma = sigmalognormal(cv)
    d = []
    for i in 1:2*G-1
        push!(d,Distributions.LogNormal(log(r[i]),sigma[j]))
    end
    return d
end

"""
function setpriorrate(G)
Set prior distribution for mean and cv of rates
"""
function setpriorrate(G,nsets,decayrate=log(2)/(60 * 16))
    r0 = [.01019*ones(2*(G-1));.0511;decayrate]
    rc = [ones(2*(G-1));1.;0.]
    rm = r0
    rcv = rc
    for i in 2:nsets
        rm = vcat(rm,r0)
        rcv = vcat(rcv,rc)
    end
    return [rm;.2],[rcv;.25]
end

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

function likelihoodtuple(param,data::TransientRNAData,model::GMLossmodel)
    r = get_rates(param,model)
    lossfactor = r[end]
    n = model.G-1
    # nRNA = data.nRNA[1] > data.nRNA[2] ? data.nRNA[1] : data.nRNA[2]
    nRNA = data.nRNA
    h0 = steady_state_full(r[1:2*n+2],n,nRNA[1])
    h1 = transient(data.time,r[1:2*n+2],lossfactor,n,nRNA[1],model.nalleles,h0)
    h0 = steady_state_full(r[1:2*n+2],n,nRNA[2])
    h2 = transient(data.time,r[2*n+3:4*n+4],lossfactor,n,nRNA[2],model.nalleles,h0)
    return h1, h2
end

function get_rates(param,model::AbstractGMmodel)
    r = copy(model.rates)
    r[model.fittedparam] = param
    return r
end

"""
readRNA(filename::String,yield::Float64=.99,nhistmax::Int=1000)
Construct mRNA count per cell histogram array of a gene
"""
# Construct mRNA count/cell histogram for a gene in filename
function readRNA_scrna(filename::String,yield::Float64=.9999,nhistmax::Int=1000)
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
function readRNA_scrna(control::String,treatment::String,yield::Float64=.9999,nhistmax::Int=1000)
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
function readRNAFISH_scrna(scRNAfolder::String,FISHfolder::String,genein::String,cond::String)
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
        x[r] = truncate_histogram(xr,.9999,1000)
    end
    l = min(length(x[1]),length(x[2]))
    data[2] = x[1][1:l] + x[2][1:l]
    return data
end

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
        # if isempty(x)
        #     println(genes[gene]," ",reffect," no prior")
        #     r = [rmm;rmm]
        # else
        #     if length(x) == 2*nparam
        #         println(gene, " ",reffect)
        #     elseif length(x) == nparam
        #         r = [x;x]
        #         println(gene, " ",reffect)
        #     else
        #         r[:,1] = [rmm;rmm]
        #         println(gene," ",reffect," rate mismatch")
        #     end
        # end
        # return r
    else
        return 0
    end
end

"""
datahistogram(data)

Return the RNA histogram data as one vector
"""
function datahistogram(data::TransientRNAData)
    v = [data.histRNA[1]]
    for i in 2:length(data.histRNA)
        v = vcat(v,data.histRNA[i])
    end
    return v
end
datahistogram(data::RNAData) = data.histRNA
