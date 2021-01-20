@everywhere include("/home/carsonc/StochasticGene/src/StochasticGene.jl")

using Dates
using DelimitedFiles


function fit_rna(gene)
    nchains = div(nworkers(),2)
    #Model parameters
    fit_rna(gene,"AUXIN",nchains)
    fit_rna(gene,"DMSO",nchains)
end

function fit_rna(gene,cond,nchains)
    println(now())

    G = 2

    label = "scRNA_T120_ss_" * cond

    # Paths
    root = "/home/carsonc/scrna/"
    # root = "/Users/carsonc/Box/scrna"
    resultfolder = "Results/2021-01-19"

    datafolder = "data/datanew/T120"
    datapath = joinpath(root,datafolder)
    datafile = joinpath(datapath,gene * "_" * cond * ".txt")

    decayrate = decay(root,gene)
    nalleles = alleles(root,gene)
    yieldprior = .1

    fittedparam = [1,2,3,5]

    # Metropolis-Hastings parameters
    maxtime = 3600*2.
    # maxtime = 10.
    samplesteps = 100000
    warmupsteps = 10000
    temp = 10.

    infolder = "Results/2021-01-19"
    filelabel = label * cond * "_" * gene *  "_" * "$G" * "_" * "$nalleles" * ".txt"
    ratefile = "rates" * filelabel
    paramfile = "param_stats" * filelabel
    ratefile = joinpath(root, joinpath(infolder,ratefile))
    paramfile = joinpath(root, joinpath(infolder,paramfile))

    if isfile(ratefile)
        r = StochasticGene.read_rates(ratefile,1)
    else
        r = [0.015,0.015,1.5,decayrate,0.1]
    end
    if isfile(paramfile)
        cv = StochasticGene.read_covlogparam(paramfile)
        if ~StochasticGene.isposdef(cv)
            cv = .02
        end
    else
        cv = .02
    end

    println(cv)
    println(r)

    println(gene," ",G)
    data,model,options = StochasticGene.steadystate_rna(datafile,label,gene,1,r,decayrate,yieldprior,G,nalleles,fittedparam,cv,maxtime,samplesteps,temp,warmupsteps)

    @time fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);

    writefile = joinpath(root,resultfolder)
    StochasticGene.write_all(writefile,fit,stats,waic,data,model)

    println(fit.accept," ",fit.total)
    println(now())
end

function decay(root,gene)
    in = readdlm(joinpath(root,"data/HCT_HEK_Half_life_top_below_1h.csv"),',')
    in[findfirst(in[:,end] .== gene),12]
end

function alleles(root,gene)
    in = readdlm(joinpath(root,"data/HCT116_alleles_number.txt"))
    in[findfirst(in[:,1] .== gene),3]
end
