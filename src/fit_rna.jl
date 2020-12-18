@everywhere include("StochasticGene.jl")

using Dates


function fit_rna()
    println(now())

    runname = "IAA_FISH_Transient_Delay"

    # Paths
    root = "/home/carsonc/scrna/"
    # root = "/Users/carsonc/Box/scrna"
    resultfolder = "Results/2020-12-17"
    datafolder = "smRNA_FISH_HCT_7timepoints"
    datapath = joinpath(root,"data/" * datafolder)
    datafiles = [joinpath(datapath,"T0_IAA");joinpath(datapath,"T30_IAA");joinpath(datapath,"T120_IAA")]

    infolder = "Results/2020-12-17"
    ratefile = "rates_IAA_FISH_Transient_Delay_Myc_3_2.txt"
    paramfile = "param_stats_IAA_FISH_Transient_Delay_Myc_3_2.txt"

    ratefile = joinpath(root, joinpath(infolder,ratefile))
    paramfile = joinpath(root, joinpath(infolder,paramfile))

    # Data parameters
    gene = "Myc"
    halflife = .52
    times = [0.;30.;120.]

    #Model parameters
    rMYC = StochasticGene.read_rates(ratefile,1)
    decayprior = log(2)/(halflife*60)
    delayprior = .05
    G = 3
    nalleles = 2
    cv = StochasticGene.read_covlogparam(paramfile)
    if ~isposdef(cv)
        cv = .01
    end
    fitted = [7;8;9;10;11;13]

    # Metropolis-Hastings parameters
    maxtime = 3600*12.
    samplesteps = 10000
    warmupsteps = 0
    temp = 10.
    method = 1

    println(gene," ",G)
    data,model,options = StochasticGene.transient_fish(datafiles,runname,times,gene,rMYC,decayprior,delayprior,G,nalleles,fitted,cv,maxtime,samplesteps,temp,method,warmupsteps);

    nchains = nworkers()
    @time fit,stats,waic = StochasticGene.run_mh(data,model,options,nchains);

    writefile = joinpath(root,resultfolder)
    StochasticGene.write_all(writefile,fit,stats,waic,data,model)

    println(fit.accept," ",fit.total)
    println(now())
end
