@everywhere include("StochasticGene.jl")

using Dates


function fit_rna()
    println(now())

    runname = "IAA_FISH_Transient_Delay"

    # Paths
    root = "/home/carsonc/scrna/"
    root = "/Users/carsonc/Box/scrna"
    resultfolder = "2020-12-17"
    datafolder = "smRNA_FISH_HCT_7timepoints"
    ratefolder = "2020-12-16/rates_IAA_FISH_Transient_Delay_Myc_3_2.txt"
    ratefolder = joinpath(root,"Results/" * ratefolder)
    path = joinpath(root,"data/" * datafolder)
    datafiles = [joinpath(path,"T0_IAA");joinpath(path,"T30_IAA");joinpath(path,"T120_IAA")]

    # Data parameters
    gene = "Myc"
    halflife = .52
    times = [0.;30.;120.]

    #Model parameters
    rMYC = StochasticGene.read_rates(ratefolder,1)
    decayprior = log(2)/(halflife*60)
    delayprior = .05
    G = 3
    nalleles = 2
    cv = .01
    fitted = [7;8;9;10;11;13]

    # Metropolis-Hastings parameters
    maxtime = 360.
    samplesteps = 10000
    warmupsteps = 1000
    temp = 2.
    method = 1

    println(gene," ",G)
    data,model,options = StochasticGene.transient_fish(datafiles,name,times,gene,rMYC,decayprior,delayprior,G,nalleles,fitted,cv,maxtime,samplesteps,temp,method,warmupsteps);

    nchains = nworkers()
    @time fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);

    writefile = joinpath(root,"Results/" * resultfolder)
    StochasticGene.write_all(writefile,fit,stats,waic,data,model)

    println(fit.accept," ",fit.total)
    println(now())
end
