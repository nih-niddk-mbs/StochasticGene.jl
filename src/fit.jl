# fit.jl
#
# Fit GRS models (generalized telegraph models) to RNA abundance and live cell imaging data
#

"""
fit(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,transitions::Tuple,datacond,G::Int,R::Int,S::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder::String,datatype::String,inlabel::String,label::String,nsets::Int,cv=0.,transient::Bool=false,samplesteps::Int=1000000,warmupsteps=0,annealsteps=0,temp=1.,tempanneal=100.,root = ".",priorcv::Float64=10.,decayrate=-1.,burst=true,nalleles=2,optimize=true,rnatype="",rtype="median",writesamples=false)

Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

# Arguments
- `nchains`: number of MCMC chains
- `gene`: gene name
- `cell`: cell type
- `fittedparam`: vector of rate indices,  indices of parameters to be fit (input as string of ints separated by "-")
- `fixedeffects`: (tuple of vectors of rate indices) string indicating which rate is fixed, e.g. "eject"
- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `datacond`: condition, if more than one condition use vector of strings e.g. ["DMSO","AUXIN"]
- `G`: number of gene states
- `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S`: number of splice sites (set to 0 for classic telegraph models and R for GRS models)
- `insertstep`: R step where reporter is first observed
- `datatype`: data type, e.g. genetrap, scRNA, smFISH
- `maxtime`: float maximum time for entire run
- `infolder`: folder pointing to results used as initial conditions
- `resultfolder`: folder where results go
- `datafolder`: folder for data, string or array of strings
- `inlabel`: name of input files (not including gene name but including condition)
- `label`: = name of output files
- `nsets`: int number of rate sets
- `cv`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `transient::Bool`: true means fit a time dependent transient model (T0, T30, T120)
- `samplesteps`: int number of samples
- `warmupsteps`: int number of warmup steps
- `annealsteps`: in number of annealing steps
- `temp`: MCMC temperature
- `tempanneal`: starting temperature for annealing
- `root`: root folder of data and Results folders
- 'priorcv`: coefficient of variation for the rate prior distributions, default is 10.
- `decayrate`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `burst`: if true then compute burst frequency
- `nalleles`: number of alleles, value in alleles folder will be used if it exists
- `optimize`: use optimizer to compute maximum likelihood value
- `rnatype`: switch used for GRS models, choices include "", "offeject", "offdecay"
- `rtype`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
- `writesamples`: write out MH samples if true, default is false
- `data`: data structure

"""

function fit(nchains::Int, gene::String, cell::String, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, datacond, G::Int, R::Int, S::Int, insertstep::Int, maxtime::Float64, infolder::String, resultfolder::String, datafolder::String, datatype::String, inlabel::String, label::String, nsets::Int, propcv=0.0, transient::Bool=false, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, root=".", priorcv::Float64=10.0, decayrate=-1.0, burst=false, nalleles=2, optimize=false, rnatype="", rtype="ml", writesamples=false, onstates=Int[], tempfish=1.0, tracedata=true)
    println(now())
    gene = check_genename(gene, "[")
    printinfo(gene, G, R, S, insertstep, datacond, datafolder, infolder, resultfolder, maxtime)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    infolder = folder_path(infolder, root, "results")
    if datatype == "genetrap"
        data, model = genetrap(root, gene, transitions, G, R, S, insertstep, nalleles, rnatype, fittedparam, infolder, resultfolder, label, "ml", tempfish, priorcv, propcv, onstates, tracedata)
    elseif datatype == "scRNA" || datatype == "fish"
        datafolder = folder_path(datafolder, root, "data")
        if occursin("-", datafolder)
            datafolder = string.(split(datafolder, "-"))
        end
        if datatype == "fish"
            fish = true
            yieldprior = 1.0
        else
            fish = false
            yieldprior = 0.05
        end
        if occursin("-", datacond)
            datacond = string.(split(datacond, "-"))
        end
        if transient
            data = data_rna(gene, datacond, datafolder, fish, label, ["T0", "T30", "T120"], [0.0, 30.0, 120.0])
        else
            data = data_rna(gene, datacond, datafolder, fish, label)
        end
        model = model_rna(data, gene, cell, G, propcv, fittedparam, fixedeffects, transitions, inlabel, infolder, nsets, root, yieldprior, decayrate, Normal, priorcv, true)
    end
    println("size of histogram: ", data.nRNA)
    options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
    fit(nchains, data, model, options, temp, resultfolder, burst, optimize, writesamples, root)
end

"""
fit(nchains,data,model,options,temp,resultfolder,burst,optimize,writesamples,root)

"""

function fit(nchains, data, model, options, temp, resultfolder, burst, optimize, writesamples, root)
    print_ll(data, model)
    fit, stats, measures = run_mh(data, model, options, nchains)
    optimized = 0
    if optimize
        try
            optimized = Optim.optimize(x -> lossfn(x, data, model), fit.parml, LBFGS())
        catch
            @warn "Optimizer failed"
        end
    end
    if burst
        bs = burstsize(fit, model)
    else
        bs = 0
    end
    finalize(data, model, fit, stats, measures, temp, resultfolder, optimized, bs, writesamples, root)
    println(now())
    get_rates(transform_rates(stats.medparam, model), model)
end

"""
lossfn(x,data,model)

Compute loss function

"""
lossfn(x, data, model) = loglikelihood(x, data, model)[1]

"""
burstsize(fit,model::AbstractGMmodel)

Compute burstsize and stats using MCMC chain

"""
function burstsize(fit, model::AbstractGMmodel)
    if model.G > 1
        b = Float64[]
        for p in eachcol(fit.param)
            r = get_rates(p, model)
            push!(b, r[2*model.G-1] / r[2*model.G-2])
        end
        return BurstMeasures(mean(b), std(b), median(b), mad(b), quantile(b, [0.025; 0.5; 0.975]))
    else
        return 0
    end
end
function burstsize(fit::Fit, model::GRSMmodel)
    if model.G > 1
        b = Float64[]
        L = size(fit.param, 2)
        rho = 100 / L
        println(rho)
        for p in eachcol(fit.param)
            r = get_rates(p, model)
            if rand() < rho
                push!(b, burstsize(r, model))
            end
        end
        return BurstMeasures(mean(b), std(b), median(b), mad(b), quantile(b, [0.025; 0.5; 0.975]))
    else
        return 0
    end
end

burstsize(r, model::AbstractGRSMmodel) = burstsize(r, model.R, length(model.Gtransitions))

function burstsize(r, R, ntransitions)
    total = min(Int(div(r[ntransitions+1], r[ntransitions])) * 2, 400)
    M = make_mat_M(make_components_M([(2, 1)], 2, R, total, r[end], ""), r)
    # M = make_components_M(transitions, G, R, nhist, decay, rnatype)
    nT = 2 * 2^R
    L = nT * total
    S0 = zeros(L)
    S0[2] = 1.0
    s = time_evolve_diff([1, 10 / minimum(r[ntransitions:ntransitions+R+1])], M, S0)
    mean_histogram(s[2, collect(1:nT:L)])
end

"""
check_genename(gene,p1)

Check genename for p1
if p1 = "[" change to "("
(since swarm cannot parse "(")

"""
function check_genename(gene, p1)
    if occursin(p1, gene)
        if p1 == "["
            gene = replace(gene, "[" => "(")
            gene = replace(gene, "]" => ")")
        elseif p1 == "("
            gene = replace(gene, "(" => "]")
            gene = replace(gene, ")" => "]")
        end
    end
    return gene
end

"""
print_ll(param,data,model,message="initial ll:")

compute and print initial loglikelihood
"""
function print_ll(param, data, model, message)
    ll, _ = loglikelihood(param, data, model)
    println(message, ll)
end
function print_ll(data, model, message="initial ll: ")
    ll, _ = loglikelihood(get_param(model), data, model)
    println(message, ll)
end

"""
printinfo(gene,G,cond,datafolder,infolder,resultfolder,maxtime)

print out run information
"""
function printinfo(gene, G, cond, datafolder, infolder, resultfolder, maxtime)
    println("Gene: ", gene, " G: ", G, " Treatment:  ", cond)
    println("data: ", datafolder)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime)
end

function printinfo(gene, G, R, S, insertstep, cond, datafolder, infolder, resultfolder, maxtime)
    if R == 0
        printinfo(gene, G, cond, datafolder, infolder, resultfolder, maxtime)
    else
        println("Gene: ", gene, " G R S insertstep: ", G, R, S, insertstep)
        println("in: ", infolder, " out: ", resultfolder)
        println("maxtime: ", maxtime)
    end
end

"""
    finalize(data,model,fit,stats,measures,temp,resultfolder,optimized,burst,writesamples,root)

write out run results and print out final loglikelihood and deviance
"""
function finalize(data, model, fit, stats, measures, temp, resultfolder, optimized, burst, writesamples, root)
    writefolder = joinpath(root, resultfolder)
    writeall(writefolder, fit, stats, measures, data, temp, model, optimized=optimized, burst=burst, writesamples=writesamples)
    println("final max ll: ", fit.llml)
    print_ll(transform_rates(vec(stats.medparam), model), data, model, "median ll: ")
    println("Median fitted rates: ", stats.medparam[:, 1])
    println("ML rates: ", inverse_transform_rates(fit.parml, model))
    println("Acceptance: ", fit.accept, "/", fit.total)
    if typeof(data) <: AbstractHistogramData
        println("Deviance: ", deviance(fit, data, model))
    end
    println("rhat: ", maximum(measures.rhat))
    if optimized != 0
        println("Optimized ML: ", Optim.minimum(optimized))
        println("Optimized rates: ", exp.(Optim.minimizer(optimized)))
    end
end

function getratefile_genetrap(infolder::String, rtype::String, gene::String, label, G, R, S, insertstep, nalleles, rnatype::String)
    model = R == 0 ? "$G" : "$G$R$S$insertstep"
    file = "rates" * "_" * label * "_" * gene * "_" * model * "_" * "$(nalleles)" * "$rnatype" * ".txt"
    joinpath(infolder, file)
end

"""
getr(gene,G,nalleles,decayrate,ejectrate,inlabel,infolder,nsets::Int,root,verbose)

"""
function getr(gene, G, nalleles, decayrate, ejectrate, inlabel, infolder, nsets::Int, root, verbose)
    r = getr(gene, G, nalleles, inlabel, infolder, root, verbose)
    if ~isnothing(r)
        if length(r) == 2 * G * nsets + 1
            for n in nsets
                r[2*G*n-1] *= clamp(r[2*G*nsets+1], eps(Float64), 1 - eps(Float64))
            end
            return r[1:2*G*nsets]
        end
        if length(r) == 2 * G * nsets
            if verbose
                println("init rates: ", r)
            end
            return r
        end
    end
    println("No r")
    setr(G, nsets, decayrate, ejectrate)
end

function getr(gene, G, nalleles, inlabel, infolder, root, verbose)
    ratefile = path_Gmodel("rates", gene, G, nalleles, inlabel, infolder, root)
    if verbose
        println("rate file: ", ratefile)
    end
    if isfile(ratefile)
        return readrates(ratefile, 3)
    else
        return nothing
    end
end
function getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root, verbose=true)
    paramfile = path_Gmodel("param-stats", gene, G, nalleles, inlabel, infolder, root)
    if isfile(paramfile)
        cv = read_covlogparam(paramfile)
        cv = float.(cv)
        if ~isposdef(cv) || size(cv)[1] != length(fittedparam)
            cv = 0.02
        end
    else
        cv = 0.02
    end
    if verbose
        println("cv: ", cv)
    end
    return cv
end
"""
    get_decay(gene::String,cell::String,root::String,col::Int=2)
    get_decay(gene::String,path::String,col::Int)

    Get decay rate for gene and cell

"""
function get_decay(gene::String, cell::String, root::String, col::Int=2)
    path = get_file(root, "data/halflives", cell, "csv")
    if isnothing(path)
        println(gene, " has no decay time")
        return -1.0
    else
        get_decay(gene, path, col)
    end
end
function get_decay(gene::String, path::String, col::Int)
    a = nothing
    in = readdlm(path, ',')
    ind = findfirst(in[:, 1] .== gene)
    if ~isnothing(ind)
        a = in[ind, col]
    end
    get_decay(a, gene)
end
function get_decay(a, gene::String)
    if typeof(a) <: Number
        return get_decay(float(a))
    else
        println(gene, " has no decay time")
        return -1.0
    end
end
get_decay(a::Float64) = log(2) / a / 60.0

"""
    alleles(gene::String,cell::String,root::String,col::Int=3)
    alleles(gene::String,path::String,col::Int=3)

    Get allele number for gene and cell
"""
function alleles(gene::String, cell::String, root::String; nalleles::Int=2, col::Int=3)
    path = get_file(root, "data/alleles", cell, "csv")
    if isnothing(path)
        return 2
    else
        alleles(gene, path, nalleles=nalleles, col=col)
    end
end

function alleles(gene::String, path::String; nalleles::Int=2, col::Int=3)
    a = nothing
    in, h = readdlm(path, ',', header=true)
    ind = findfirst(in[:, 1] .== gene)
    if isnothing(ind)
        return nalleles
    else
        a = in[ind, col]
        return isnothing(a) ? nalleles : Int(a)
    end
end
