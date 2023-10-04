# fit.jl
#
# Fit GRS models (generalized telegraph models) to RNA abundance and live cell imaging data
#

"""
fit(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,transitions::Tuple,datacond,G::Int,R::Int,S::Int,maxtime::Float64,infolder::String,resultfolder::String,datafolder::String,datatype::String,inlabel::String,label::String,nsets::Int,cv=0.,transient::Bool=false,samplesteps::Int=1000000,warmupsteps=0,annealsteps=0,temp=1.,tempanneal=100.,root = ".",priorcv::Float64=10.,decayrate=-1.,burst=true,nalleles=2,optimize=true,splicetype="",ratetype="median",writesamples=false)

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
- `splicetype`: switch used for GRS models, choices include "", "offeject", "offdecay"
- `ratetype`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
- `writesamples`: write out MH samples if true, default is false
- `data`: data structure

"""

# function fit(nchains::Int, gene::String, cell::String, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, datacond, G::Int, R::Int, S::Int, insertstep::Int, maxtime::Float64, infolder::String, resultfolder::String, datafolder::String, datatype::String, inlabel::String, label::String, nsets::Int, propcv=0.0, transient::Bool=false, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, root=".", priorcv::Float64=10.0, decayrate=-1.0, burst=false, nalleles=2, optimize=false, splicetype="", ratetype="ml", writesamples=false, onstates=Int[], tempfish=1.0, tracedata=true)
#     println(now())
#     gene = check_genename(gene, "[")
#     printinfo(gene, G, R, S, insertstep, datacond, datafolder, infolder, resultfolder, maxtime)
#     resultfolder = folder_path(resultfolder, root, "results", make=true)
#     infolder = folder_path(infolder, root, "results")
#     if datatype == "genetrap"
#         data, model = genetrap(root, gene, transitions, G, R, S, insertstep, nalleles, splicetype, fittedparam, infolder, label, "ml", tempfish, priorcv, propcv, onstates, tracedata)
#     elseif datatype == "scRNA" || datatype == "fish"
#         datafolder = folder_path(datafolder, root, "data")
#         if occursin("-", datafolder)
#             datafolder = string.(split(datafolder, "-"))
#         end
#         if datatype == "fish"
#             fish = true
#             yieldprior = 1.0
#         else
#             fish = false
#             yieldprior = 0.05
#         end
#         if occursin("-", datacond)
#             datacond = string.(split(datacond, "-"))
#         end
#         if transient
#             data = data_rna(gene, datacond, datafolder, fish, label, ["T0", "T30", "T120"], [0.0, 30.0, 120.0])
#         else
#             data = data_rna(gene, datacond, datafolder, fish, label)
#         end
#         model = model_rna(data, gene, cell, G, propcv, fittedparam, fixedeffects, transitions, inlabel, infolder, nsets, root, yieldprior, decayrate, Normal, priorcv, true)
#     end
#     println("size of histogram: ", data.nRNA)
#     options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
#     fit(nchains, data, model, options, joinpath(root, resultfolder), burst, optimize, writesamples)
# end

function fit(nchains::Int, datatype::String, dttype, datafolder, gene::String, cell::String, datacond::String, interval, nascent,infolder::String, resultfolder::String, inlabel::String, label::String,
    fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, root=".", rmean=[],nalleles=2, priorcv::Float64=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
     propcv=0.01, maxtime::Float64=10.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, tempfish=1.0, burst=false, optimize=false, writesamples=false)
    if decayrate < 0
        decayrate = get_decay(gene, cell, root)
    end
    println(now())
    gene = check_genename(gene, "[")
    printinfo(gene, G, R, S, insertstep, datacond, datafolder, infolder, resultfolder, maxtime)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    infolder = folder_path(infolder, root, "results")
    datafolder = folder_path(datafolder, root, "data")
    data = load_data(datatype, dttype, datafolder, label, gene, datacond, interval, tempfish, nascent)
    r = readrates(infolder, label, gene, G, R, S, insertstep, nalleles, ratetype)
    r = prior_ratemean(transitions, R, S, insertstep, decayrate, noiseparams, weightind)
    model = load_model(data, r, rmean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind)
    options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)
end

function readrates(infolder, label, gene, G, R, S, insertstep, nalleles, ratetype="median")
    if R == 0
        name = filename(label, gene, G, nalleles)
    else
        name = filename(label, gene, G, R, S, insertstep, nalleles)
    end
    readrates(joinpath(infolder, "rates" * name), get_row(ratetype))
end

datatype_dict() = Dict("rna" => 1, "rnaoffon" => 2, "rnadwelltimes" => 3, "rnatrace" => 4, "dwelltimes" => 5, "trace" => 5, "tracenascent" => 6)


"""
    load_data(datatype, dttype, datafolder, label, gene, datacond, interval, tempfish, nascent)

return data structure
"""
function load_data(datatype, dttype, datafolder, label, gene, datacond, interval, tempfish, nascent)
    if datatype == "rna"
        len, h = read_rna(gene, datacond, datafolder)
        return RNAData(label, gene, len, h)
    elseif datatype == "rnaoffon"
        len, h = read_rna(gene, datacond, datafolder[1])
        h = div.(h, tempfish)
        LC = readfile(gene, datacond, datafolder[2])
        return RNALiveCellData(label, gene, len, h, LC[:, 1], LC[:, 2], LC[:, 3])
    elseif datatype == "rnadwelltimes"
        len, h = read_rna(gene, datacond, datafolder[1])
        h = div.(h, tempfish)
        LC = readfiles(gene, datacond, datafolders[2:end])
        bins, DT = read_dwelltimes(datafolders)
        return RNADwellTimeData(label, gene, len, h, bins, DT, dttype)
    elseif datatype == "trace"
        trace = read_tracefiles(datafolder, datacond)
        return TraceData("trace", gene, interval, trace)
    elseif datatype == "tracenascent"
        trace = read_tracefiles(datafolder, datacond)
        return TraceNascentData(label, gene, interval, trace, nascent)
    elseif datatype == "rnatrace"
        len, h = histograms_rna(datafolder[1], gene, fish)
        traces = read_tracefiles(datafolder[2], datacond)
        return TraceRNAData(label, gene, interval, traces, len, h)
    end
end


"""
    load_model(data, r, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind)

    trace_model(r::Vector, transitions::Tuple, G, R, S, insertstep, fittedparam; noiseparams=5, probfn=prob_GaussianMixture, weightind=5, propcv=0.05, priorprob=Normal, 
    priormean=[fill(.1, num_rates(transitions, R, S, insertstep)); fill(100,noiseparams-1);.9], priorcv=[fill(10, num_rates(transitions, R, S, insertstep)); fill(.5, noiseparams)], fixedeffects=tuple(), onstates::Vector=[G],genetrap=false,nhist=20,nalleles=2)

"""
function load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind)
    if typeof(data) <: AbstractRNAData
        components = make_components_M(transitions, G, 0, data.nRNA, decayrate, splicetype)
    elseif typeof(data) <: AbstractTraceData
        reporter = ReporterComponents(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions, R, S, insertstep) + weightind)
        components = make_components_T(transitions, G, R, S, insertstep, splicetype)
    elseif typeof(data) <: AbstractTraceHistogramData
        reporter = ReporterComponents(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions, R, S, insertstep) + weightind)
        components = make_components_MT(transitions, G, R, S, insertstep, data.nRNA, decayrate, splicetype)
    elseif typeof(data) <: AbstractHistogramData
        for i in eachindex(onstates)
            if isempty(onstates[i])
                onstates[i] = on_states(G, R, S, insertstep)
            end
        end
        reporter = onstates
        components = make_components_MTAI(transitions, G, R, S, insertstep, onstates, data.nRNA, decayrate)
    end
    priord = prior_distribution(rm, transitions, R, S, insertstep, fittedparam, decayrate, priorcv, noiseparams, weightind)
    if propcv < 0
        propcv = getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root)
    end
    load_model(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

"""
    load_model(r,fittedparam,fixedeffects,transitions,G,R,S,insertstep,priord,components,reporter)


"""
function load_model(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    if R == 0
        if isempty(fixedeffects)
            return GMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(r, transitions, G, nalleles, priord, propcv, fittedparam, method, components, reporter)
        else
            return GMfixedeffectsmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components)}(r, transitions, G, nalleles, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        end
    else
        if isempty(fixedeffects)
            return GRSMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, method, components, reporter)
        else
            return GRSMfixedeffectsmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        end
    end
end



"""
    model_prior(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, decayrate, priorcv=10.0, noiseparams, weightind)

TBW
"""
function prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, decayrate, priorcv, noiseparams, weightind)
    if isempty(rm)
        rm = prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind)
    end
    rcv = fill(priorcv, length(rm))
    rcv[num_rates(transitions, R, S, insertstep)] = .1
    distribution_array(log.(rm[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
end

function prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind)
    if S > 0
        S = R
    end
    ntransitions = length(transitions)
    nrates = num_rates(transitions, R, S, insertstep)
    rm = fill(0.1, nrates)
    rm[collect(1:ntransitions)] .= 0.01
    rm[nrates] = decayrate
    if noiseparams > 0
        rm = [rm; fill(100., noiseparams)]
        rm[nrates+weightind] = 0.9
    end
    rm
end

"""
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)

TBW
"""
function fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)
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
    finalize(data, model, fit, stats, measures, options.temp, resultfolder, optimized, bs, writesamples)
    println(now())
    get_rates(stats.medparam, model, false)
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
    # M = make_components_M(transitions, G, R, nhist, decay, splicetype)
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
printinfo(gene,G,datacond,datafolder,infolder,resultfolder,maxtime)

print out run information
"""
function printinfo(gene, G, datacond, datafolder, infolder, resultfolder, maxtime)
    println("Gene: ", gene, " G: ", G, " Treatment:  ", datacond)
    println("data: ", datafolder)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime)
end

function printinfo(gene, G, R, S, insertstep, datacond, datafolder, infolder, resultfolder, maxtime)
    if R == 0
        printinfo(gene, G, datacond, datafolder, infolder, resultfolder, maxtime)
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
function finalize(data, model, fit, stats, measures, temp, writefolder, optimized, burst, writesamples)
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
    if uppercase(cell) == "HBEC"
        if uppercase(gene) ∈ uppercase.(genes_gt())
            # return log(2.0) / (60 .* halflife_gt()[gene])
            return get_decay(halflife_gt()[gene])
        else
            println(gene, " has no decay time")
            return -1.0
        end
    else
        path = get_file(root, "data/halflives", cell, "csv")
        if isnothing(path)
            println(gene, " has no decay time")
            return -1.0
        else
            return get_decay(gene, path, col)
        end
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
