# This file is part of StochasticGene.jl   

# fit.jl
#
# Fit GRS models (generalized telegraph models) to RNA abundance and live cell imaging data
#
"""
HBEC gene information
"""
genes_hbec() = ["CANX"; "DNAJC5"; "ERRFI1"; "KPNB1"; "MYH9"; "Rab7a"; "RHOA"; "RPAP3"; "Sec16A"; "SLC2A1"]
genelength_hbec() = Dict([("Sec16A", 42960); ("SLC2A1", 33802); ("ERRFI1", 14615); ("RHOA", 52948); ("KPNB1", 33730); ("MYH9", 106741); ("DNAJC5", 40930); ("CANX", 32710); ("Rab7a", 88663); ("RPAP3", 44130); ("RAB7A", 88663); ("SEC16A", 42960)])
MS2end_hbec() = Dict([("Sec16A", 5220); ("SLC2A1", 26001); ("ERRFI1", 5324); ("RHOA", 51109); ("KPNB1", 24000); ("MYH9", 71998); ("DNAJC5", 14857); ("CANX", 4861); ("Rab7a", 83257); ("RPAP3", 38610);("SEC16A", 5220); ("RAB7A", 83257)])
halflife_hbec() = Dict([("CANX", 50.0), ("DNAJC5", 5.0), ("ERRFI1", 1.35), ("KPNB1", 9.0), ("MYH9", 10.0), ("Rab7a", 50.0), ("RHOA", 50.0), ("RPAP3", 7.5), ("Sec16A", 8.0), ("SLC2A1", 5.0), ("RAB7A", 50.0), ("SEC16A", 8.0)])


"""
    fit(nchains::Int, datatype::String, dttype, datapath, gene::String, cell::String, datacond::String, interval, nascent, infolder::String, resultfolder::String, inlabel::String, label::String,
    fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, root=".", maxtime::Float64=60.0, priormean::Vector=Float64[], nalleles=2, priorcv::Float64=10.0, onstates=Int[], 
    decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
    propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false)

fit(nchains::Int,gene::String,cell::String,fittedparam::Vector,fixedeffects::Tuple,transitions::Tuple,datacond,G::Int,R::Int,S::Int,maxtime::Float64,infolder::String,resultfolder::String,datapath::String,datatype::String,inlabel::String,label::String,nsets::Int,cv=0.,transient::Bool=false,samplesteps::Int=1000000,warmupsteps=0,annealsteps=0,temp=1.,tempanneal=100.,root = ".",priorcv::Float64=10.,decayrate=-1.,burst=true,nalleles=2,optimize=true,splicetype="",ratetype="median",writesamples=false)

Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

# Arguments
- `nchains`: number of MCMC chains
- `datatype`: choices "rna", "rnaonoff", "rnadwelltime", "trace", "tracenascent", "tracerna"
- `ddtype`: Vector of dwell time types, e.g. "ON", "OFF"
- `datapath`: path to file or folder for data, string or array of strings
- `gene`: gene name
- `cell`: cell type
- `datacond`: condition, if more than one condition use vector of strings e.g. ["DMSO","AUXIN"]
- `interval`: frame interval for traces
- `nascent`: fraction of alleles exhibiting nascent rna
- `infolder`: folder pointing to results used as initial conditions
- `resultfolder`: folder for results
- `inlabel`: name of input files (not including gene name but including condition)
- `label`: = name of output files
- `fittedparam`: vector of rate indices,  indices of parameters to be fit (input as string of ints separated by "-")
- `fixedeffects`: (tuple of vectors of rate indices) string indicating which rate is fixed, e.g. "eject"
- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G`: number of gene states
- `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S`: number of splice sites (set to 0 for classic telegraph models and R for GRS models)
- `insertstep`: R step where reporter is first observed
- `root`: root folder of data and Results folders
- `maxtime`: float maximum time for entire run
- `priormean`: Vector of prior rate means
- `nalleles`: number of alleles, value in alleles folder will be used if it exists
- 'priorcv`: coefficient of variation for the rate prior distributions, default is 10.
- `onstates`: vector of sojourn or on states
- `decayrate`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `splicetype`: switch used for GRS models, choices include "", "offeject", "offdecay"
- `probfn`: observation (noise) probability distribution for trace data
- `noiseparams`: number of noise distribution parameters
- `weightind`: noise parameter index of the first bias weight parameter for probfn mixture distributions (e.g. Gaussian Mixture)
- `ratetype`: which rate for initial condition, choices are "ml", "mean", "median", or "last"
- `propcv`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `samplesteps`: int number of samples
- `warmupsteps`: int number of warmup steps
- `annealsteps`: in number of annealing steps
- `temp`: MCMC temperature
- `tempanneal`: starting temperature for annealing
- `temprna`: temperature for scRNA distribution compared to dwell time distributions (reduces mRNA cell count by 1/temprna)
- `burst`: if true then compute burst frequency
- `optimize`: use optimizer to compute maximum likelihood value
- `writesamples`: write out MH samples if true, default is false


"""
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene::String, cell::String, datacond::String, interval, nascent, infolder::String, resultfolder::String, inlabel::String, label::String,
    fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, root=".", maxtime::Float64=60.0,priormean=Float64[], nalleles=2, priorcv::Float64=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, ratetype="median",
    propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false)
    println(now())
    gene = check_genename(gene, "[")
    printinfo(gene, G, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    infolder = folder_path(infolder, root, "results")
    datapath = folder_path(datapath, root, "data")
    data = load_data(datatype, dttype, datapath, label, gene, datacond, interval, temprna, nascent)
    ~occursin("trace", lowercase(datatype)) && (noiseparams = 0)
    decayrate < 0 && (decayrate = get_decay(gene, cell, root))
    isempty(priormean) && (priormean = prior_ratemean(transitions, R, S, insertstep, decayrate, noiseparams, weightind))
    isempty(fittedparam) && (fittedparam = collect(1:num_rates(transitions,R,S,insertstep)-1))
    r = readrates(infolder, inlabel, gene, G, R, S, insertstep, nalleles, ratetype)
    isempty(r) && (r = priormean)
    println(r)
    model = load_model(data, r, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind)
    options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
    # return data, model, options
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)
end

"""
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)


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
    load_data(datatype, dttype, datapath, label, gene, datacond, interval, temprna, nascent)

return data structure
"""
function load_data(datatype, dttype, datapath, label, gene, datacond, interval, temprna, nascent)
    if datatype == "rna"
        len, h = read_rna(gene, datacond, datapath)
        return RNAData(label, gene, len, h)
    elseif datatype == "rnaonoff"
        len, h = read_rna(gene, datacond, datapath[1])
        h = div.(h, temprna)
        LC = readfile(gene, datacond, datapath[2])
        return RNAOnOffData(label, gene, len, h, LC[:, 1], LC[:, 2], LC[:, 3])
    elseif datatype == "rnadwelltime"
        len, h = read_rna(gene, datacond, datapath[1])
        h = div.(h, temprna)
        bins, DT = read_dwelltimes(datapaths[2:end])
        return RNADwellTimeData(label, gene, len, h, bins, DT, dttype)
    elseif datatype == "trace"
        trace = read_tracefiles(datapath, datacond)
        return TraceData("trace", gene, interval, trace)
    elseif datatype == "tracenascent"
        trace = read_tracefiles(datapath, datacond)
        return TraceNascentData(label, gene, interval, trace, nascent)
    elseif datatype == "tracerna"
        len, h = read_rna(gene, datacond, datapath[1])
        traces = read_tracefiles(datapath[2], datacond)
        return TraceRNAData(label, gene, interval, traces, len, h)
    else
        throw("$datatype not included")
    end
end


"""
    load_model(data, r, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind)

return model structure
"""
function load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind)
    if typeof(data) <: AbstractRNAData
        reporter = onstates
        components = make_components_M(transitions, G, 0, data.nRNA, decayrate, splicetype)
    elseif typeof(data) <: AbstractTraceData
        reporter = ReporterComponents(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions, R, S, insertstep) + weightind)
        components = make_components_T(transitions, G, R, S, insertstep, splicetype)
    elseif typeof(data) <: AbstractTraceHistogramData
        reporter = ReporterComponents(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions, R, S, insertstep) + weightind)
        components = make_components_MT(transitions, G, R, S, insertstep, data.nRNA, decayrate, splicetype)
    elseif typeof(data) <: AbstractHistogramData
        if isempty(onstates)
            onstates = on_states(G, R, S, insertstep)
        else
            for i in eachindex(onstates)
                if isempty(onstates[i])
                    onstates[i] = on_states(G, R, S, insertstep)
                end
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
    load_model(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)

"""
function load_model(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    if R == 0
        if isempty(fixedeffects)
            return GMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, nalleles, priord, propcv, fittedparam, method, components, reporter)
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
    prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, decayrate, priorcv, noiseparams, weightind)


"""
function prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, decayrate, priorcv, noiseparams, weightind)
    if isempty(rm)
        rm = prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind)
    end
    rcv = fill(priorcv, length(rm))
    rcv[num_rates(transitions, R, S, insertstep)] = 0.1
    distribution_array(log.(rm[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
end

"""
    prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind)


"""
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
        rm = [rm; fill(100.0, noiseparams)]
        rm[nrates+weightind] = 0.9
    end
    rm
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
printinfo(gene,G,datacond,datapath,infolder,resultfolder,maxtime)

print out run information
"""
function printinfo(gene, G, datacond, datapath, infolder, resultfolder, maxtime)
    println("Gene: ", gene, " G: ", G, " Treatment:  ", datacond)
    println("data: ", datapath)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime)
end

function printinfo(gene, G, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime)
    if R == 0
        printinfo(gene, G, datacond, datapath, infolder, resultfolder, maxtime)
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
    writeall(writefolder, fit, stats, measures, data, temp, model, optimized=optimized, burst=burst, writesamples=writesamples)
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
        if uppercase(gene) ∈ uppercase.(genes_hbec())
            # return log(2.0) / (60 .* halflife_hbec()[gene])
            return get_decay(halflife_hbec()[gene])
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
