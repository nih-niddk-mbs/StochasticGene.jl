# functions to run StochasticGene for Bayesian parameter estimation of stochastic Markov models of gene transcription

using Dates
using DelimitedFiles

"""
fit_rna(nchains::Int,gene::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=40000,warmupsteps=0,annealsteps=0,temp=100.,tempanneal=100.,root = "/home/carsonc/scrna/")

Fit steady state or transient GM model to RNA data for a single gene

Arguments:
nchains = number of MCMC chains
gene = gene name
datacond = condition, if more than one condition is used enter as a single string separated by underscore, e.g. "DMSO_AUXIN"
maxtime = float maximum time for entire run
infolder = folder pointing to results used as initial conditions
resultfolder = folder where results go
datafolder = folder for data
inlabel = name of input files (not including gene name but including condition)
label = name of output files
nsets = int number of rate sets
runcycle = if true, cycle through all parameters sequentially in MCMC
samplesteps = int number of samples
warmupsteps = int number of warmup steps
annealsteps = in number of annealing steps
temp = MCMC temperature
tempanneal = starting temperature for annealing
root = root folder of data and Results folders
fittedparam = indices of parameters to be fit (input as string of ints separated by "-")

"""
# function fit_rna(nchains::Int,gene::String,G::Int,data::StochasticGene.HistogramData,maxtime::Float64,nsets,fittedparam,infolder,resultfolder,datafolder,inlabel,runcycle::Bool,params::Tuple,root)
function fit_rna(nchains::Int,gene::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fittedparam = make_fittedparam(G,nsets)
    fit_rna(nchains,gene,fittedparam,[],datacond,G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root)
    nothing
end

function fit_rna(nchains::Int,gene::String,fittedparam::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fittedparam = parse.(Int,split(fittedparam,"-"))
    fit_rna(nchains,gene,fittedparam,[],datacond,G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root)
    nothing
end

function fit_rna(nchains::Int,gene::String,fittedparam::String,fixedeffects::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    fittedparam = parse.(Int,split(fittedparam,"-"))
    fixedeffects = make_fixedeffects(fixedeffects)
    fit_rna(nchains,gene,fittedparam,fixedeffects,datacond,G,maxtime,infolder,resultfolder,datafolder,inlabel,label,nsets,runcycle,transient,samplesteps,warmupsteps,annealsteps,temp,tempanneal,root)
    nothing
end

function fit_rna(nchains::Int,gene::String,fittedparam::Vector,fixedeffects::Tuple,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    println(now())
    datacond = string.(split(datacond,"-"))
    if transient
        data = make_data(gene,datacond,datafolder,label,root,["T0","T30","T120"],[0.,30.,120.])
    else
        data = make_data(gene,datacond,datafolder,label,root)
    end
    initialize(gene,G,datacond,datafolder,infolder,resultfolder,maxtime)
    decayrate = decay(root,gene)
    if decayrate < 0
        throw("error")
    else
        println(decayrate)
    end
    nalleles = alleles(root,gene)
    yieldprior = .1
    r = getr(gene,G,nalleles,decayrate,inlabel,infolder,nsets,root,data)
    cv = getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root)
    if cv != 0.02
        warmupsteps = 0
    end
    if runcycle
        maxtime /= 2
        cv = .02
        r = cycle(nchains,data,r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,maxtime,temp,tempanneal)
        annealsteps = 0
        warmupsteps = div(samplesteps,5)
    end

    # options = StochasticGene.MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
    # model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,0)
    # param,_ = StochasticGene.initial_proposal(model)
    options,model,param = get_structures(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
    initial_ll(param,data,model)
    fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);
    finalize(data,model,fit,stats,waic,temp,resultfolder,root)
    println(now())
    nothing
end

function fit_fish(nchains::Int,gene::String,datacond,G::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder,inlabel,label,nsets,runcycle::Bool=false,transient::Bool=false,samplesteps::Int=100000,warmupsteps=20000,annealsteps=0,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")

end

function rna_fish(gene,cond,fishfolder,rnafolder,yield,root)
    datarna = StochasticGene.histograms_rna(StochasticGene.scRNApath(gene,cond,rnafolder,root),gene,false)
    f=reduce_fish(gene,cond,datarna[1],fishfolder,yield,root)
    return datarna[2],f
end

function reduce_fish(gene,cond,nhist,fishfolder,yield,root)
    datafish = StochasticGene.histograms_rna(StochasticGene.FISHpath(gene,cond,fishfolder,root),gene,true)
    StochasticGene.technical_loss(datafish[2],yield,nhist)
end

function get_structures(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
    if length(fixedeffects) > 0
        model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,0)
    else
        model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,0)
    end
    options = StochasticGene.MHOptions(samplesteps,annealsteps,warmupsteps,maxtime,temp,tempanneal)
    param,_ = StochasticGene.initial_proposal(model)
    return options,model,param
end

function cycle(nchains,data,r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,maxtime,temp,tempanneal)
    # options = StochasticGene.MHOptions(100,0,0,maxtime/10,temp,tempanneal)
    # model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,decayrate,yieldprior,0)
    # param,_ = StochasticGene.initial_proposal(model)
    options,model,param = get_structures(r,G,nalleles,nsets,cv,fittedparam,fixedeffects,decayrate,yieldprior,100,0,0,maxtime/10,temp,tempanneal)
    initial_ll(param,data,model,"pre-cycle ll:")
    t0 = time()
    while (time() - t0 < maxtime)
        for fp in fittedparam
            model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,[fp],decayrate,yieldprior,0)
            fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);
            r = StochasticGene.get_rates(fit.parml,model)
        end
    end
    return r
end

function initial_ll(param,data,model,message="initial ll:")
    ll,_ = StochasticGene.loglikelihood(param,data,model)
    println(message,ll)
end

function initialize(gene,G,cond,datafolder,infolder,resultfolder,maxtime)
    println(gene," ",G," ",cond)
    println(datafolder)
    println("in: ", infolder," out: ",resultfolder)
    println(maxtime)
end

function finalize(data,model,fit,stats,waic,temp,resultfolder,root)
    writefile = joinpath(root,resultfolder)
    StochasticGene.writeall(writefile,fit,stats,waic,data,temp,model)
    println("final ll: ",fit.llml)
    println(fit.accept," ",fit.total)
    println("Deviance: ",StochasticGene.deviance(fit,data,model))
    println(stats.meanparam)
end

function transient_rna(nchains,gene::String,fittedparam,cond::Vector,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,samplesteps::Int=40000,warmupsteps=20000,annealsteps=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
    data = make_data(gene,cond,G,datafolder,label,nsets,root)
    model = make_model(gene,G,fittedparam,inlabel,infolder,nsets,root,data)
    param,_ = StochasticGene.initial_proposal(model)
    return param, data, model
end

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

function make_model(gene,G,fittedparam,inlabel,infolder,nsets,root,data,verbose=true)
    decayrate = decay(root,gene)
    if verbose
        println(decayrate)
        if decayrate < 0
            throw("error")
        end
    end
    nalleles = alleles(root,gene)
    yieldprior = .1
    r = getr(gene,G,nalleles,decayrate,inlabel,infolder,nsets,root,data,verbose)
    cv = getcv(gene,G,nalleles,fittedparam,inlabel,infolder,root,verbose)
    model = StochasticGene.model_rna(r,G,nalleles,nsets,cv,fittedparam,decayrate,yieldprior,0)
end

make_model(gene,r::Vector,G,fittedparam,nsets,root) = StochasticGene.model_rna(r,G,alleles(root,gene),nsets,.02,fittedparam,r[2*G],.1,0)

function make_fixedeffects(fixedeffects)
    if fixedeffects == "eject"
        fixedeffects = ([2*G*(nsets-1) + 2*G-1],[2*G-1])
    elseif fixedeffects == "off"
        fixedeffects = ([2*G*(nsets-1) + 2*G-2],[2*G-2])
    elseif fixedeffects == "on"
        fixedeffects = ([2*G*(nsets-1) + 2*G-3],[2*G-3])
    elseif fixedeffects == "offeject" || fixedeffects == "offeject"
        fixedeffects = ([2*G*(nsets-1) + 2*G-2,2*G*(nsets-1) + 2*G-1],[2*G-2,2*G-1])
    elseif fixedeffects == "oneject" || fixedeffects == "ejecton"
        fixedeffects = ([2*G*(nsets-1) + 2*G-3,2*G*(nsets-1) + 2*G-1],[2*G-3,2*G-1])
    elseif fixedeffects == "onoff" || fixedeffects == "offon"
        fixedeffects = ([2*G*(nsets-1) + 2*G-3,2*G*(nsets-1) + 2*G-2],[2*G-3,2*G-2])
    end
    return fixedeffects
end
function make_fittedparam(G::Int,nsets)

    if nsets == 1
        if G == 3
            fittedparam = [1,2,3,4,5,7]
        elseif G == 2
            fittedparam = [1,2,3,5]
        elseif G == 1
            fittedparam = [1,3]
        end
    else
        if G == 3
            # fittedparam = [1,2,3,4,5,7,8,9,10,11,13]
            # fittedparam = [7,8,9,10,11,13]
            fittedparam = [7,8,9,10,11]
        elseif G == 2
            # fittedparam = [1,2,3,5,6,7,9]
            # fittedparam = [5,6,7,9]
            fittedparam = [5,6,7]
        elseif G == 1
            fittedparam = [3,5]
        end
    end
    return fittedparam
end

function make_data(gene::String,cond::String,datafolder,label,root)
    if cond == "null"
        cond = ""
    end
    datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
    StochasticGene.data_rna(datafile,label,gene,false)
end

function make_data(gene::String,cond::Array,datafolder,label,root)
    datafile = Array{String,1}(undef,length(cond))
    for i in eachindex(cond)
        datafile[i] = StochasticGene.scRNApath(gene,cond[i],datafolder,root)
    end
    StochasticGene.data_rna(datafile,label,gene,false)
end

function make_data(gene::String,cond::String,datafolder,label,root,sets::Vector,time::Vector)
    if cond == "null"
        cond = ""
    end
    datafile =[]
    for set in sets
        folder = joinpath(datafolder,set)
        datafile = vcat(datafile,StochasticGene.scRNApath(gene,cond,folder,root))
    end
    StochasticGene.data_rna(datafile,label,times,gene,false)
end


function getr(gene,G,nalleles,decayrate,inlabel,infolder,nsets::Int,root,data,verbose=true)
    ratefile = StochasticGene.path_Gmodel("rates",gene,G,nalleles,inlabel,infolder,root)
    if verbose
        println(ratefile)
    end
    if isfile(ratefile)
        r = StochasticGene.readrates(ratefile,2)
        r[end] = clamp(r[end],eps(Float64),1-eps(Float64))
        if length(r) == 2*G*nsets + 1
            if verbose
                println(r)
            end
            return r
        end
    end
    println("No r")
    setr(G,decayrate,nsets,data)
end


function setr(G,decayrate,nsets,data)
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
    r[end] = .1
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

function decay(root::String,gene,file="data/HCT116_all_cells_histograms_and_half_lives_March_2021/Starved_Rep7_half_lives.csv",col=2)
    path = joinpath(root,file)
    if isfile(path)
        in = readdlm(path,',')
        a = in[findfirst(in[:,1] .== gene),col]
        return decay(a)
    else
        println(gene," has no decay time")
        return -1.
    end
end

decay(a::Float64) = log(2)/a/60.

function decay(a,gene)
    if typeof(a) <: Number
        return decay(a)
    else
        println(gene," has no decay time")
        return -1.
    end
end

function alleles(root,gene,file="data/HCT116_alleles_number.txt")
    in = readdlm(joinpath(root,"data/HCT116_alleles_number.txt"))
    in[findfirst(in[:,1] .== gene),3]
end

function repair_measures(resultfolder,measurefile,ratefile::String,cond,n::Int,datafolder,root)
    front,back = split(measurefile,".")
    fix = front * "fix" * "." * back
    measurepath = joinpath(resultfolder,measurefile)
    ratepath = joinpath(resultfolder,ratefile)
    rates,_ = readdlm(ratepath,',',header=true)
    println(length(rates[:,1]))
    measures,head = readdlm(measurepath,',',header=true)
    f = open(joinpath(resultfolder,fix),"w")
    writedlm(f,[head "Deviance fixed" "LogML fixed" "AIC fixed"],',')
    for i in eachindex(rates[:,1])
        h,hd = histograms(rates[i,:],cond,n,datafolder,root)
        d = StochasticGene.deviance(h,hd)
        ll = StochasticGene.crossentropy(h,hd)
        a = 4*(n+1) + 2*ll
        writedlm(f,[measures[i:i,:] d ll a],',')
    end
    close(f)
end

function repair_measures(outfile,file,fittedparam,cond,G::Int,folder::String,datafolder,label,nsets,root)
    contents,header = readdlm(file,',',header=true)
    f = open(outfile,"w")
    writedlm(f,[header],',')
    for row in eachrow(contents)
        row[3],row[6],row[2] = measures(string(row[1]),fittedparam,cond,G,folder,datafolder,label,nsets,root)
        writedlm(f, [row],',')
    end
    close(f)
end

function measures(gene::String,fittedparam,cond,G::Int,folder::String,datafolder,label,nsets,root)
    param,data,model = steadystate_rna(gene,fittedparam,cond,G,folder,datafolder,label,nsets,root)
    ll = StochasticGene.loglikelihood(param,data,model)[1]
    return ll, StochasticGene.aic(length(fittedparam),ll), StochasticGene.deviance(data,model)
end

function deviance(r,cond,n,datafolder,root)
    h,hd = histograms(r,cond,n,datafolder,root)
    StochasticGene.deviance(h,hd)
end


function compute_deviance(outfile,ratefile::String,cond,n,datafolder,root)
    f = open(outfile,"w")
    rates,head = readdlm(ratefile,',',header=true)
    for r in eachrow(rates)
        d=deviance(r,cond,n,datafolder,root)
        writedlm(f,[gene d],',')
    end
    close(f)
end

function write_histograms(outfolder,ratefile::String,cond,n,datafolder,root)
    rates,head = readdlm(ratefile,',',header=true)
    for r in eachrow(rates)
        h,hd = histograms(r,cond,n,datafolder,root)
        f = open(joinpath(outfolder,r[1] * ".txt"),"w")
        writedlm(f,h)
        close(f)
    end
end

function histograms(r,cond,n,datafolder,root)
    gene = String(r[1])
    datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
    hd = StochasticGene.read_scrna(datafile)
    h = StochasticGene.steady_state(r[2:2*n+3],r[end],n,length(hd),alleles(root,gene))
    return h,hd
end

function write_histograms(outfolder,ratefile,fittedparam,datacond,G::Int,datafolder::String,label,nsets,root)
    rates,head = readdlm(ratefile,',',header=true)
    if ~isdir(outfolder)
        mkpath(outfolder)
    end
    cond = string.(split(datacond,"-"))
    for r in eachrow(rates)
        h = histograms(r,fittedparam,datacond,G,datafolder,label,nsets,root)
        for i in eachindex(cond)
            f = open(joinpath(outfolder,string(r[1]) * cond[i] * ".txt"),"w")
            writedlm(f,h[i])
            close(f)
        end
    end
end

function histograms(rin,fittedparam,cond,G::Int,datafolder,label,nsets,root)
    gene = string(rin[1])
    r = float.(rin[2:end])
    param,data,model = steadystate_rna(r,gene,fittedparam,cond,G,datafolder,label,nsets,root)
    StochasticGene.likelihoodarray(r,data,model)
end

function write_burst_stats(outfile,infile::String,G::String,folder,cond,root)
    condarray = split(cond,"-")
    g = parse(Int,G)
    lr = 2*g
    lc = 2*g-1
    freq = Array{Float64,1}(undef,2*length(condarray))
    burst = similar(freq)
    f = open(outfile,"w")
    contents,head = readdlm(infile,',',header=true)
    label = Array{String,1}(undef,0)
    for c in condarray
        label = vcat(label, "Freq " * c, "sd","Burst Size " * c, "sd")
    end
    writedlm(f,["gene" reshape(label,1,length(label))],',')
    for r in eachrow(contents)
        gene = String(r[1])
        rates = r[2:end]
        rdecay = decay(root,gene)
        cov = StochasticGene.read_covparam(joinpath(folder,getfile("param_stats",gene,G,folder,cond)[1]))
        # mu = StochasticGene.readmean(joinpath(folder,getfile("param_stats",gene,G,folder,cond)[1]))
        for i in eachindex(condarray)
            j = i-1
            freq[2*i-1], freq[2*i] = frequency(rates[1+lr*(i-1)],sqrt(cov[1+lc*j,1+lc*j]),rdecay)
            burst[2*i-1], burst[2*i] = burstsize(rates[3+lr*j],rates[2+lr*j],cov[3+lc*j,3+lc*j],cov[2+lc*j,2+lc*j],cov[2+lc*j,3+lc*j])
        end
        writedlm(f,[gene freq[1] freq[2] burst[1] burst[2] freq[3] freq[4] burst[3] burst[4]],',')
    end
    close(f)
end

frequency(ron,sd,rdecay) = (ron/rdecay, sd/rdecay)

function burstsize(reject::Float64,roff,covee,covoo,coveo::Float64)
        v = StochasticGene.var_ratio(reject,roff,covee,covoo,coveo)
        return reject/roff, sqrt(v)
end

function ratestats(gene,G,folder,cond)
    filestats=joinpath(folder,getfile("param_stats",gene,G,folder,cond)[1])
    filerates = joinpath(folder,getratefile(gene,G,folder,cond)[1])
    rates = StochasticGene.readrates(filerates)
    cov = StochasticGene.read_covparam(filestats)
    # mu = StochasticGene.readmean(filestats)
    # sd = StochasticGene.readsd(filestats)
    return rates, cov
end

meanofftime(gene::String,infile,n,method,root) = sum(1 .- offtime(gene,infile,n,method,root))

function meanofftime(r::Vector,n::Int,method::Int)
    if n == 1
        return 1/r[1]
    else
        return sum(1 .- offtime(r,n,method))
    end
end

function offtime(r::Vector,n::Int,method::Int)
    _,_,TI = StochasticGene.mat_G_DT(r,n)
    vals,_ = StochasticGene.eig_decompose(TI)
    minval = min(minimum(abs.(vals[vals.!=0])),.2)
    StochasticGene.offtimeCDF(collect(1.:5/minval),r,n,TI,method)
end

function offtime(gene::String,infile,n,method,root)
    contents,head = readdlm(infile,',',header=true)
    r = float.(contents[gene .== contents[:,1],2:end-1])[1,:]
    offtime(r,n,method)

end

function write_moments(outfile,genelist,cond,datafolder,root)
    f = open(outfile,"w")
    writedlm(f,["Gene" "Expression Mean" "Expression Variance"],',')
    for gene in genelist
        datafile = StochasticGene.scRNApath(gene,cond,datafolder,root)
        # data = StochasticGene.data_rna(datafile,label,gene,false)
        len,h = StochasticGene.histograms_rna(datafile,gene,false)
        # h,hd = histograms(r,cond,n,datafolder,root)
        writedlm(f,[gene StochasticGene.mean_histogram(h) StochasticGene.var_histogram(h)],',')
    end
    close(f)
end

function join_files(file1::String,file2::String,outfile::String,addlabel::Bool=true)
    contents1,head1 = readdlm(file1,',',header=true)   # model G=2
    contents2,head2 = readdlm(file2,',',header=true)   # model G=3
    f = open(outfile,"w")
    if addlabel
        header = vcat(String.(head1[2:end]) .* "_G2",String.(head2[2:end]) .* "_G3")
    else
        header = vcat(String.(head1[2:end]),String.(head2[2:end]))
    end
    header = reshape(permutedims(header),(1,length(head1)+length(head2)-2))
    header = hcat(head1[1],header)
    println(header)
    writedlm(f,header,',')
    for row in 1:size(contents1,1)
        if contents1[row,1] == contents2[row,1]
            contents = hcat(contents1[row:row,2:end],contents2[row:row,2:end])
            contents = reshape(permutedims(contents),(1,size(contents1,2)+size(contents2,2)-2))
            contents = hcat(contents1[row,1],contents)
            writedlm(f,contents,',')
        end
    end
    close(f)
end

function join_files(models::Array,files::Array,outfile::String,addlabel::Bool=true)
    m = length(files)
    contents = Array{Array,1}(undef,m)
    headers = Array{Array,1}(undef,m)
    len = 0
    for i in 1:m
        contents[i],headers[i] = readdlm(files[i],',',header=true)
        len += length(headers[i][2:end])
    end
    f = open(outfile,"w")
    header = Array{String,1}(undef,0)
    for i in 1:m
        if addlabel
            header = vcat(header,String.(headers[i][2:end]) .* "_G$(models[i])")
        else
            header = vcat(header,String.(headers[i][2:end]))
        end
    end
    header = reshape(permutedims(header),(1,len))
    header = hcat(headers[1][1],header)
    println(header)
    writedlm(f,header,',')
    for row in 1:size(contents[1],1)
        content = contents[1][row:row,2:end]
        for i in 1:m-1
            if contents[i][row,1] == contents[i+1][row,1]
                content = hcat(content,contents[i+1][row:row,2:end])
                # content = reshape(permutedims(content),(1,len))
            end
        end
        content = hcat(contents[1][row:row,1],content)
        writedlm(f,[content],',')
    end
    close(f)
end

function best_AIC(outfile,infile)
    contents,head = readdlm(infile,',',header=true)
    head = vec(head)
    ind = occursin.("AIC",string.(head)) .& .~ occursin.("WAIC",string.(head))
    f = open(outfile,"w")
    writedlm(f,["Gene" "Model 1" "Model 2" "Model 3" "Winning Model"],',')
    for row in eachrow(contents)
        writedlm(f,[row[1] reshape(row[ind],1,3) argmin(float.(row[ind]))],',')
    end
    close(f)
end

function sample_non1_genes(infile,n)
    contents,head = readdlm(infile,',',header=true)
    list = Array{String,1}(undef,0)
    for c in eachrow(contents)
        if c[5] != 1
            push!(list,c[1])
        end
    end
    a = sample(list,n)
end





function best_model(file::String)
    contents,head = readdlm(file,',',header=true)
    f = open(file,"w")
    head = hcat(head,"Winning Model")
    writedlm(f,head,',')
    for row in eachrow(contents)
        if abs(row[11] - row[4]) > row[5] + row[12]
            if row[11] < row[4]
                writedlm(f, hcat(permutedims(row),3),',')
            else
                writedlm(f, hcat(permutedims(row),2),',')
            end
        else
            if row[13] < row[6]
                writedlm(f, hcat(permutedims(row),3),',')
            else
                writedlm(f, hcat(permutedims(row),2),',')
            end
        end
    end
    close(f)
end

function bestmodel(measures2,measures3)
    m2,head = readdlm(infile,',',header=true)
    m3,head = readdlm(infile,',',header=true)
end


function add_best_burst(filein,fileout,filemodel2,filemodel3)
    contents,head = readdlm(filein,',',header=true)
    burst2,head2 = readdlm(filemodel2,',',header=true)
    burst3,head3 = readdlm(filemodel3,',',header=true)
    f = open(fileout,"w")
    head = hcat(head,["mean off period" "bust size"])
    writedlm(f,head,',')
    for row in eachrow(contents)
        if Int(row[end]) == 2
            writedlm(f, hcat(permutedims(row),permutedims(burst2[findfirst(burst2[:,1] .== row[1]),2:3])),',')
        else
            writedlm(f, hcat(permutedims(row),permutedims(burst3[findfirst(burst3[:,1] .== row[1]),2:3])),',')
        end
    end
    close(f)
end

function add_best_occupancy(filein,fileout,filemodel2,filemodel3)
    contents,head = readdlm(filein,',',header=true)
    occupancy2,head2 = readdlm(filemodel2,',',header=true)
    occupancy3,head3 = readdlm(filemodel3,',',header=true)
    f = open(fileout,"w")
    head = hcat(head,["Off -2" "Off -1" "On State" ])
    writedlm(f,head,',')
    for row in eachrow(contents)
        if Int(row[end-2]) == 2
            writedlm(f, hcat(permutedims(row),hcat("NA",permutedims(occupancy2[findfirst(occupancy2[:,1] .== row[1]),2:end]))),',')
        else
            writedlm(f, hcat(permutedims(row),permutedims(occupancy3[findfirst(occupancy3[:,1] .== row[1]),2:end])),',')
        end
    end
    close(f)
end


# precompile(fit_rna,(String,String,Int,Float64,String,String,String,String,String,Int,String,Bool))
# fit_rna(nchains,gene::String,cond::String,G::Int,maxtime::Float64,infolder::String,resultfolder,datafolder,inlabel,label,nsets,runcycle::Bool=false,sample::Int=40000,warmup=20000,anneal=100000,temp=1.,tempanneal=100.,root = "/home/carsonc/scrna/")
