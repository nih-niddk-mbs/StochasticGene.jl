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
MS2end_hbec() = Dict([("Sec16A", 5220); ("SLC2A1", 26001); ("ERRFI1", 5324); ("RHOA", 51109); ("KPNB1", 24000); ("MYH9", 71998); ("DNAJC5", 14857); ("CANX", 4861); ("Rab7a", 83257); ("RPAP3", 38610); ("SEC16A", 5220); ("RAB7A", 83257)])
halflife_hbec() = Dict([("CANX", 50.0), ("DNAJC5", 5.0), ("ERRFI1", 1.35), ("KPNB1", 9.0), ("MYH9", 10.0), ("Rab7a", 50.0), ("RHOA", 50.0), ("RPAP3", 7.5), ("Sec16A", 8.0), ("SLC2A1", 5.0), ("RAB7A", 50.0), ("SEC16A", 8.0)])


"""
    fit(; <keyword arguments> )

Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

#Arguments
- `nchains::Int=2`: number of MCMC chains = number of processors called by Julia, default = 2
- `datatype::String=""`: String that desecribes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracenascent", "tracerna"
- `dttype=String[]`
- `datapath=""`: path to data file or folder or array of files or folders
- `cell::String=""': cell type for halflives and allele numbers
- `datacond=""`: string or vector of strings describing data treatment condition, e.g. "WT", "DMSO" or ["DMSO","AUXIN"]
- `interval=[1.0, 1.]`: vector of frame interval of intensity traces and transient time
- `nascent=[1,2,.7]`: vector of number of spots, total number of locations (e.g. number of cells times number of alleles/cell), and fraction of busting traces
- `traceinfo=tuple()`: tuple of trace information (transient::Float64,onfraction::Float64,background::tracetype)
- `infolder::String=""`: result folder used for initial parameters
- `resultfolder::String=test`: folder for results of MCMC run
- `label::String=""`: label of output files produced
- `inlabel::String=""`: label of files used for initial conditions
- `fittedparam::Vector=Int[]`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]
- `fixedeffects::Tuple=tuple()`: tuple of vectors of rates that are fixed where first index is fit and others are fixed to first, e.g. ([3,8],) means  index 8 is fixed to index 3
     (only first parameter should be included in fixedeffects)
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G::Int=2`: number of gene states
- `R::Int=0`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S::Int=0`: number of splice sites (set to 0 for classic telegraph models and R for GRS models)
- `insertstep::Int=1`: R step where reporter is inserted
- `Gfamily=""`: String describing type of G transition model, e.g. "3state", "KP" (kinetic proofreading), "Refractory"
- `root="."`: name of root directory for project, e.g. "scRNA"
- `priormean=Float64[]`: mean of prior rate distribution
- 'priorcv=10.`: coefficient of variation for the rate prior distributions, default is 10.
- `nalleles=2`: number of alleles, value in alleles folder will be used if it exists  
- `onstates=Int[]`: vector of on or sojourn states
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `splicetype=""`: RNA pathway for GRS models, (e.g. "offeject" =  spliced intron is not viable)
- `probfn=prob_GaussianMixture`: probability function for hmm observation probability (e.g. prob_GaussianMixture)
- `noiseparams=5`: number of parameters of probfn
- `weightind=5`: parameter index of bias probability of mixtures, e.g. noiseparams=5, weightind=5 means last noise parameter is for mixture bias
- `hierarchical=tuple()`: tuple of hierchical model parameters
- `ratetype="median"`: which rate to use for initial condition, choices are "ml", "mean", "median", or "last"
- `propcv=0.01`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `maxtime=Float64=60.`: maximum wall time for run, default = 60 min
- `samplesteps::Int=1000000`: number of MCMC sampling steps
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `annealsteps=0`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: annealing temperature
- `temprna=1.`: reduce rna counts by temprna compared to dwell times
- `burst=false`: if true then compute burst frequency
- `optimize=false`: use optimizer to compute maximum likelihood value
- `writesamples=false`: write out MH samples if true, default is false

"""

function fit(; nchains::Int=2, datatype::String="rna", dttype=String[], datapath="HCT116_testdata/", gene="MYC", cell::String="HCT116", datacond="MOCK", interval=1.0, nascent=[1, 2], infolder::String="HCT116_test", resultfolder::String="HCT116_test", inlabel::String="", label::String="", fittedparam::Vector=Int[], fixedeffects::Tuple=tuple(), transitions::Tuple=([1, 2], [2, 1]), G::Int=2, R::Int=0, S::Int=0, insertstep::Int=1, Gfamily="", root=".", priormean=Float64[], nalleles=2, priorcv=10.0, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, hierarchical=tuple(), ratetype="median", propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, src="")
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, Gfamily)
    fit(nchains, datatype, dttype, datapath, gene, cell, datacond, interval, nascent, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, root, maxtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noiseparams, weightind, hierarchical, ratetype,
        propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples)

end
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene::String, cell::String, datacond::String, interval, nascent, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, root=".", maxtime::Float64=60.0, priormean=Float64[], priorcv=10.0, nalleles=2, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_GaussianMixture, noiseparams=5, weightind=5, hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false)
    println(now())
    gene = check_genename(gene, "[")
    if S > 0 && S ≠ R
        println("Setting S = R")
        S = R
    end
    printinfo(gene, G, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    infolder = folder_path(infolder, root, "results")
    datapath = folder_path(datapath, root, "data")
    data = load_data(datatype, dttype, datapath, label, gene, datacond, interval, temprna, nascent)
    ~occursin("trace", lowercase(datatype)) && (noiseparams = 0)
    decayrate < 0 && (decayrate = get_decay(gene, cell, root))
    if isempty(priormean)
        if isempty(hierarchical)
            priormean = prior_ratemean(transitions, R, S, insertstep, decayrate, noiseparams, weightind)
        else
            priormean = prior_ratemean(transitions, R, S, insertstep, decayrate, noiseparams, weightind, length(data.trace), hierarchical[1])
        end
    end
    isempty(fittedparam) && (fittedparam = default_fittedparams(datatype, transitions, R, S, insertstep, noiseparams))
    r = readrates(infolder, inlabel, gene, G, R, S, insertstep, nalleles, ratetype)
    if isempty(r)
        r = priormean
        println("No prior")
    end
    if isempty(hierarchical)
        println(r)
    else
        println(length(r))
    end
    model = load_model(data, r, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind, hierarchical)
    options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)
end

"""
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)


"""
function fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)
    print_ll(data, model)
    fits, stats, measures = run_mh(data, model, options, nchains)
    optimized = 0
    if optimize
        try
            optimized = Optim.optimize(x -> lossfn(x, data, model), fits.parml, LBFGS())
        catch
            @warn "Optimizer failed"
        end
    end
    if burst
        bs = burstsize(fits, model)
    else
        bs = 0
    end
    finalize(data, model, fits, stats, measures, options.temp, resultfolder, optimized, bs, writesamples)
    println(now())
    # get_rates(stats.medparam, model, false)
    return fits, stats, measures, data, model, options
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
        bins, DT = read_dwelltimes(datapath[2:end])
        return RNADwellTimeData(label, gene, len, h, bins, DT, dttype)
    elseif datatype == "trace"
        trace = read_tracefiles(datapath[1], datacond, round(Int, interval[2] / interval[1]))
        background = read_tracefiles(datapath[2], datacond, round(Int, interval[2] / interval[1]))
        weight = (1 - nascent[3]) / nascent[3] * length(trace)
        return TraceData(label, gene, interval, (trace, background, weight))
    elseif datatype == "tracenascent"
        trace = read_tracefiles(datapath, datacond, round(Int, interval[2] / interval[1]))
        return TraceNascentData(label, gene, interval, trace, nascent)
    elseif datatype == "tracerna"
        len, h = read_rna(gene, datacond, datapath[1])
        traces = read_tracefiles(datapath[2], datacond, round(Int, interval[2] / interval[1]))
        return TraceRNAData(label, gene, interval, traces, len, h)
    else
        throw("$datatype not included")
    end
end

"""
    load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind, hierarchical)

return model structure
"""
function load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noiseparams, weightind, hierarchical)
    if typeof(data) <: AbstractRNAData
        reporter = onstates
        components = make_components_M(transitions, G, 0, data.nRNA, decayrate, splicetype)
    elseif typeof(data) <: AbstractTraceData
        reporter = HMMReporter(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions, R, S, insertstep) + weightind)
        components = make_components_T(transitions, G, R, S, insertstep, splicetype)
    elseif typeof(data) <: AbstractTraceHistogramData
        reporter = HMMReporter(noiseparams, num_reporters_per_state(G, R, S, insertstep), probfn, num_rates(transitions, R, S, insertstep) + weightind)
        components = make_components_MT(transitions, G, R, S, insertstep, data.nRNA, decayrate, splicetype)
    elseif typeof(data) <: AbstractHistogramData
        if isempty(onstates)
            onstates = on_states(G, R, S, insertstep)
        else
            for i in eachindex(onstates)
                if isempty(onstates[i])
                    onstates[i] = on_states(G, R, S, insertstep)
                end
                onstates[i] = Int64.(onstates[i])
            end
        end
        reporter = onstates
        if typeof(data) == RNADwellTimeData
            if length(onstates) == length(data.DTtypes)
                components = make_components_MTD(transitions, G, R, S, insertstep, onstates, data.DTtypes, data.nRNA, decayrate, splicetype)
            else
                throw("length of onstates and data.DTtypes not the same")
            end
        else
            components = make_components_MTAI(transitions, G, R, S, insertstep, onstates, data.nRNA, decayrate, splicetype)
        end
    end
    if isempty(hierarchical)
        priord = prior_distribution(rm, transitions, R, S, insertstep, fittedparam, decayrate, priorcv, noiseparams, weightind)
    else
        nsets = hierarchical[1]
        nrates = num_rates(transitions, R, S, insertstep) + noiseparams
        priord = prior_distribution(rm[1:nsets*nrates], transitions, R, S, insertstep, make_fitted(fittedparam, nsets, [], nrates, 0), decayrate, priorcv, noiseparams, weightind)
    end
    if propcv < 0
        propcv = getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root)
    end
    load_model(data, r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, 1, components, reporter, hierarchical)
end

"""
    load_model(data, r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter, hierarchical)

"""
function load_model(data, r, transitions::Tuple, G::Int, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter, hierarchical)
    if isempty(hierarchical)
        if R == 0
            return GMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, nalleles, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        else
            return GRSMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        end
    else
        nrates = num_rates(transitions, R, S, insertstep) + reporter.n
        nindividuals = length(data.trace)
        pool = Pool(hierarchical[1], length(fittedparam), nrates, length(hierarchical[2]), nindividuals)
        fittedparam = make_fitted(fittedparam, hierarchical[1], hierarchical[2], nrates, nindividuals)
        fixedeffects = make_fixed(fixedeffects, hierarchical[3], nrates, nindividuals)
        return GRSMhierarchicalmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, pool, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    end
end

"""
    make_fitted(fittedparams,N)


"""
function make_fitted(fittedpool, npool, fittedindividual, nrates, nindividuals)
    f = copy(fittedpool)
    for i in 1:npool-1
        append!(f, fittedpool .+ i * nrates)
    end
    for i in 1:nindividuals
        append!(f, fittedindividual .+ (i + npool - 1) * nrates)
    end
    f
end

function make_fixed(fixedpool, fixedindividual, nrates, nindividuals)
    fixed = Vector{Int}[]
    for f in fixedpool
        push!(fixed, f)
    end
    for h in fixedindividual
        push!(fixed, [h + i * nrates for i in 0:nindividuals-1])
    end
    tuple(fixed...)
end

"""
    prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, decayrate, priorcv, noiseparams, weightind)


"""
function prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, decayrate, priorcv, noiseparams, weightind)
    if isempty(rm)
        rm = prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind)
    end
    if priorcv isa Number
        rcv = fill(priorcv, length(rm))
        rcv[num_rates(transitions, R, S, insertstep)] = 0.1
    end
    if length(rcv) == length(rm)
        return distribution_array(log.(rm[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
    else
        throw("priorcv not the same length as prior mean")
    end
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

function prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind, nindividuals, npools, cv=1.0)
    rm = prior_ratemean(transitions::Tuple, R::Int, S::Int, insertstep, decayrate, noiseparams, weightind)
    r = copy(rm)
    if npools == 2
        append!(r, cv .* rm)
    end
    for i in 1:nindividuals
        append!(r, rm)
    end
    r
end
"""
    default_fittedparams(datatype, transitions, R, S, insertstep, noiseparams)

create vector of fittedparams that includes all rates except the decay time
"""
function default_fittedparams(datatype, transitions, R, S, insertstep, noiseparams)
    n = num_rates(transitions, R, S, insertstep)
    fittedparam = collect(1:n-1)
    if occursin("trace", datatype)
        fittedparam = vcat(fittedparam, collect(n+1:n+noiseparams))
    end
    fittedparam
end

"""
lossfn(x,data,model)

Compute loss function

"""
lossfn(x, data, model) = loglikelihood(x, data, model)[1]

"""
burstsize(fits,model::AbstractGMmodel)

Compute burstsize and stats using MCMC chain

"""
function burstsize(fits, model::AbstractGMmodel)
    if model.G > 1
        b = Float64[]
        for p in eachcol(fits.param)
            r = get_rates(p, model)
            push!(b, r[2*model.G-1] / r[2*model.G-2])
        end
        return BurstMeasures(mean(b), std(b), median(b), mad(b), quantile(b, [0.025; 0.5; 0.975]))
    else
        return 0
    end
end
function burstsize(fits::Fit, model::GRSMmodel)
    if model.G > 1
        b = Float64[]
        L = size(fits.param, 2)
        rho = 100 / L
        println(rho)
        for p in eachcol(fits.param)
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
    finalize(data,model,fits,stats,measures,temp,resultfolder,optimized,burst,writesamples,root)

write out run results and print out final loglikelihood and deviance
"""
function finalize(data, model, fits, stats, measures, temp, writefolder, optimized, burst, writesamples)
    println("final max ll: ", fits.llml)
    print_ll(transform_rates(vec(stats.medparam), model), data, model, "median ll: ")
    println("Median fitted rates: ", stats.medparam[:, 1])
    println("ML rates: ", inverse_transform_rates(fits.parml, model))
    println("Acceptance: ", fits.accept, "/", fits.total)
    if typeof(data) <: AbstractHistogramData
        println("Deviance: ", deviance(fits, data, model))
    end
    println("rhat: ", maximum(measures.rhat))
    if optimized != 0
        println("Optimized ML: ", Optim.minimum(optimized))
        println("Optimized rates: ", exp.(Optim.minimizer(optimized)))
    end
    writeall(writefolder, fits, stats, measures, data, temp, model, optimized=optimized, burst=burst, writesamples=writesamples)
end



# """
# getr(gene,G,nalleles,decayrate,ejectrate,inlabel,infolder,nsets::Int,root,verbose)

# """
# function getr(gene, G, nalleles, decayrate, ejectrate, inlabel, infolder, nsets::Int, root, verbose)
#     r = getr(gene, G, nalleles, inlabel, infolder, root, verbose)
#     if ~isnothing(r)
#         if length(r) == 2 * G * nsets + 1
#             for n in nsets
#                 r[2*G*n-1] *= clamp(r[2*G*nsets+1], eps(Float64), 1 - eps(Float64))
#             end
#             return r[1:2*G*nsets]
#         end
#         if length(r) == 2 * G * nsets
#             if verbose
#                 println("init rates: ", r)
#             end
#             return r
#         end
#     end
#     println("No r")
#     setr(G, nsets, decayrate, ejectrate)
# end

# function getr(gene, G, nalleles, inlabel, infolder, root, verbose)
#     ratefile = path_Gmodel("rates", gene, G, nalleles, inlabel, infolder, root)
#     if verbose
#         println("rate file: ", ratefile)
#     end
#     if isfile(ratefile)
#         return readrates(ratefile, 3)
#     else
#         return nothing
#     end
# end

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
            return 1.0
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
        return 1.0
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
