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
    get_transitions(G, TransitionType)

Create transitions tuple
"""
function get_transitions(G::Int, TransitionType)
    typeof(G) <: AbstractString && (G = parse(Int, G))
    if G == 2
        return ([1, 2], [2, 1])
    elseif G == 3
        if occursin("KP", TransitionType)
            return ([1, 2], [2, 1], [2, 3], [3, 1])
        elseif occursin("cyclic", TransitionType)
            return ([1, 2], [2, 3], [3, 1])
        else
            return ([1, 2], [2, 1], [2, 3], [3, 2])
        end
    elseif G == 4
        if occursin("KP", TransitionType)
            return ([1, 2], [2, 1], [2, 3], [3, 2], [3, 4], [4, 2])
        elseif occursin("cyclic", TransitionType)
            return ([1, 2], [2, 3], [3, 4], [4, 1])
        else
            return ([1, 2], [2, 1], [2, 3], [3, 2], [3, 4], [4, 3])
        end
    else
        throw("transition type unknown")
    end
end

"""
    get_transitions(G::Tuple, TransitionType)

TBW
"""
function get_transitions(G::Tuple, TransitionType)
    Tuple([get_transitions(G, TransitionType) for G in G])
end

function make_coupling(source::UnitRange{Int64}=1:3, target::UnitRange{Int64}=1:3)
    coupling = []
    for s in source, t in target
        push!(coupling, ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1))
    end
    return coupling
end

"""
    set_decayrate(decayrate::Float64, gene, cell, root)

TBW
"""
function set_decayrate(decayrate::Float64, gene, cell, root)
    if decayrate < 0
        return get_decay(gene, cell, root)
    else
        return decayrate
    end
end

"""
    set_decayrate(decayrate::Tuple, gene, cell, root)

TBW
"""
function set_decayrate(decayrate::Tuple, gene, cell, root)
    for i in eachindex(decayrate)
        if typeof(gene) < Tuple
            decayrate[i] = set_decayrate(decayrate[i], gene[i], cell, root)
        else
            decayrate[i] = set_decayrate(decayrate[i], gene, cell, root)
        end
    end
    return decayrate
end

"""
    reset_S(S::Int, R::Int, insertstep::Int)

Set S to R - insertstep + 1 if greater than zero
"""
function reset_S(S::Int, R::Int, insertstep::Int)
    if S > R - insertstep + 1
        S = R - insertstep + 1
        println("Setting S to ", S)
    end
    return S
end

"""
    reset_S(S::Tuple, R::Tuple, insertstep::Tuple)

"""
function reset_S(S::Tuple, R::Tuple, insertstep::Tuple)
    S = collect(S)
    for i in eachindex(S)
        if S[i] > R[i] - insertstep[i] + 1
            S[i] = R[i] - insertstep[i] + 1
            println("Setting S[$i] to ", S[i])
        end
    end
    return Tuple(S)
end

"""
    reset_nalleles(nalleles, coupling)

TBW
"""
function reset_nalleles(nalleles, coupling)
    if !isempty(coupling)
        nalleles = 1
        println("Setting nalleles to 1")
    end
    return nalleles
end

"""
    fit(; <keyword arguments> )

Fit steady state or transient GM model to RNA data for a single gene, write the result (through function finalize), and return nothing.

For coupled transcribing units, arguments transitions, G, R, S, insertstep, and trace become tuples of the single unit type, e.g. If two types of transcription models are desired with G= 2 and G=3 then
then G = (2,3). 

#Arguments
- `annealsteps=0`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `burst=false`: if true then compute burst frequency
- `cell::String=""`: cell type for halflives and allele numbers
- `coupling=tuple()`: if nonempty, a 4-tuple where elements are 
    1. tuple of model indices corresponding to each unit, e.g. (1, 1, 2) means that unit 1 and 2 use model 1 and unit 3 uses model 2
    2. tuple of vectors indicating source units for each unit, e.g. ([2,3], [1], Int[]) means unit 1 is influenced by source units 2 and 3, unit 2 is influenced by unit 1 and unit 3 is uninfluenced.
    3. source states: tuple of vectors of strings, e.g. (["G3","R1"], []) means that model 1 influences other units whenever it is in G state 3 or R step 1 (if a number is not included (e.g. (R,0)) then all the Gstates or R steps are included), 
        while model 2 does not influence any other unit
    4. target transitions:tuple, e.g. ([], 4) means that model 1 is not influenced by any source while model 2 is influenced by sources at transition 4. Transitions
        are number consecutively by order of the transition rates. So for a G=2 model, transition 1 is the G1 to G2 transition and transition 3 is the initiation transition
    5. Int indicating number of coupling parameters
- `datatype::String=""`: String that describes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracerna", "tracejoint", "tracegrid"
- `datacond=""`: string or vector of strings describing data, e.g. "WT", "DMSO" or ["DMSO","AUXIN"], ["gene","enhancer"]
- `datapath=""`: path to data file or folder or array of files or folders
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `dttype=String[]`: dwelltime types, choices are "OFF", "ON", for R states and "OFFG", "ONG" for G states
- `elongationtime=6.0`: average time for elongation, vector of times for coupled model
- `fittedparam::Vector=Int[]`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7]  (applies to shared rates for hierarchical models, fitted hyper parameters are specified by individual fittedparams)
- `fixedeffects::String`: if "fixed" is included after a hyphen, then fixedeffects Tuple will be created such that R transitions are fixed to be identical
- `fixedeffects::Tuple=tuple()`: tuple of vectors of rates that are fixed where first index is fit and others are fixed to first, e.g. ([3,8],) means index 8 is fixed to index 3 (only first parameter should be included in fittedparam) (applies to shared rates for hierarchical models)
- `gene::String="MYC"`: gene name
- `grid=nothing`: Int number of grid points for grid model
- `G=2`: number of gene states, for coupled models G, R, S, and insertstep are vectors (vector for coupled models)
- `hierarchical=tuple()`: empty tuple for nonhierarchical model; 3-tuple for hierarchical: hierarchical=(number of hyper parameter sets::Int, individual fittedparams::Vector, individual fixedeffects::Tuple),
    for hierarchical models the keywords `fittedparam` and `fixedeffects` pertain to shared rates.  rates are given by a single vector that can be reshaped into a matrix where the columns correspond to the model rates and noise params, the first nhyper rows pertain to the shared and hyper parameter rates (whether fit or not), 
    usually the first row is the shared and mean hyper parameters and the 2nd are the standard deviations, the rest of the rows are the individual rates and noise params
- `infolder::String=""`: result folder used for initial parameters
- `inlabel::String=""`: label of files used for initial conditions
- `insertstep=1`: R step where reporter is inserted
- `label::String=""`: label of output files produced
- `maxtime=Float64=60.`: maximum wall time for run, default = 60 min
- `method=Tsit5()`: DifferentialEquations.jl numerical method (e.g. Tsit5(), lsoda(),...); use a tuple for hierarchical models: method = tuple(method, Bool) = (numerical method (currently not used), true if transition rates are shared)
- `nalleles=1`: number of alleles, value in alleles folder will be used if it exists, always set to 1 if coupling nonempty (multiple alleles handled as extra identical units)
- `nchains::Int=2`: number of MCMC chains = number of processors called by Julia, default = 2
- `noisepriors=[]`: priors of observation noise (use empty set if not fitting traces), superseded if priormean is set
- `onstates::Vector{Int}=Int[]`: vector of on or sojourn states, e.g. [[2,3],Int[]], use empty vector for R states, do not use Int[] for R=0 models
- `optimize=false`: use optimizer to compute maximum likelihood value
- `priormean=Float64[]`: mean rates of prior distribution (must set priors for all rates including those that are not fitted)
- `priorcv=10.`: (vector or number) coefficient of variation(s) for the rate prior distributions, default is 10.
- `probfn=prob_Gaussian`: probability function for HMM observation probability (i.e., noise distribution)
- `propcv=0.01`: coefficient of variation (mean/std) of proposal distribution, if cv <= 0. then cv from previous run will be used
- `resultfolder::String=test`: folder for results of MCMC run
- `R=0`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `root="."`: name of root directory for project, e.g. "scRNA"
- `samplesteps::Int=1000000`: number of MCMC sampling steps
- `S=0`: number of splice sites (set to 0 for classic telegraph models and R - insertstep + 1 for GRS models)
- `splicetype=""`: RNA pathway for GRS models, (e.g., "offeject" = spliced intron is not viable)
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: annealing temperature
- `temprna=1.`: reduce RNA counts by temprna compared to dwell times
- `traceinfo=(1.0, 1., -1, 1.)`: 4-tuple = (frame interval of intensity traces, starting frame time in minutes, ending frame time (use -1 for last index), fraction of observed active traces); for simultaneous joint traces, the fraction of active traces is a vector of the active fractions for each trace, e.g. (1.0, 1., -1, [.5, .7]) 
- `TransitionType=""`: String describing G transition type, e.g. "3state", "KP" (kinetic proofreading), "cyclic", or if hierarchical, coupled
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2-state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3-state kinetic proofreading model
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `writesamples=false`: write out MH samples if true, default is false

Example:

If you are in the folder where data/HCT116_testdata is installed, then you can fit the mock RNA histogram running 4 mcmc chains with

bash> julia -p 4

julia> fits, stats, measures, data, model, options = fit(nchains = 4)

"""
function fit(; nchains::Int=2, datatype::String="rna", dttype=String[], datapath="HCT116_testdata/", gene="MYC", cell="HCT116", datacond="MOCK", traceinfo=(1.0, 1, -1, 1.0), infolder::String="HCT116_test", resultfolder::String="HCT116_test", inlabel::String="", label::String="", fittedparam=Int[], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, coupling=tuple(), TransitionType="nstate", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, maxtime::Float64=60.0, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, TransitionType)
    if typeof(fixedeffects) <: AbstractString
        fixedeffects, fittedparam = make_fixedfitted(datatype, fixedeffects, transitions, R, S, insertstep, length(noisepriors), coupling, grid)
        println(transitions)
        println(fixedeffects)
        println(fittedparam)
    end
    fit(nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype,
        propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method)
end

"""
    fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene::String, cell::String, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fixedeffects::String, G::String, R::String, S::String, insertstep::String, TransitionType="", root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())

"""
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene::String, cell::String, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fixedeffects::String, G::String, R::String, S::String, insertstep::String, TransitionType="", grid=nothing, root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())
    transitions = get_transitions(G, TransitionType)
    fixedeffects, fittedparam = make_fixedfitted(datatype, fixedeffects, transitions, parse(Int, R), parse(Int, S), parse(Int, insertstep), length(noisepriors), coupling, grid)
    println(transitions)
    println(fixedeffects)
    println(fittedparam)
    fit(nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects, transitions, parse(Int, G), parse(Int, R), parse(Int, S), parse(Int, insertstep), tuple(), root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype,
        propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method)
end

"""
    fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling=tuple(), root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())


"""
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling=tuple(), grid=nothing, root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())
    S = reset_S(S, R, insertstep)
    nalleles = reset_nalleles(nalleles, coupling)
    fit(readrates(folder_path(infolder, root, "results"), inlabel, gene, G, R, S, insertstep, nalleles, ratetype), nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method)
end

"""
    fit(rinit, nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())

"""
function fit(rinit, nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5())
    println(now())
    printinfo(gene, G, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    data, model, options = make_structures(rinit, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, method)
    fit(nchains, data, model, options, resultfolder, burst, optimize, writesamples)
end

function make_structures(rinit, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime::Float64=60.0, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, method=Tsit5())
    gene = check_genename(gene, "[")
    S = reset_S(S, R, insertstep)
    nalleles = reset_nalleles(nalleles, coupling)
    infolder = folder_path(infolder, root, "results")
    datapath = folder_path(datapath, root, "data")
    data = load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna)
    decayrate = set_decayrate(decayrate, gene, cell, root)
    priormean = set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]))
    fittedparam = set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noisepriors, coupling, grid)
    model = load_model(data, rinit, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid)
    options = MHOptions(samplesteps, warmupsteps, annealsteps, maxtime, temp, tempanneal)
    return data, model, options
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
    load_data_trace(datapath, label, gene, datacond, traceinfo, datatype)

TBW
"""
function load_data_trace(datapath, label, gene, datacond, traceinfo, datatype)
    if typeof(datapath) <: String
        trace = read_tracefiles(datapath, datacond, traceinfo)
    else
        trace = read_tracefiles(datapath[1], datacond, traceinfo)
    end
    (length(trace) == 0) && throw("No traces")
    println(length(trace))
    println(datapath)
    println(datacond)
    println(traceinfo)
    weight = (1 - traceinfo[4]) / traceinfo[4]
    nframes = round(Int, mean(length.(trace)))  #mean number of frames of all traces
    if length(traceinfo) > 4
        background = traceinfo[5][1] .+ randn(nframes) .* traceinfo[5][2]
    else
        background = Vector[]
    end
    if datatype == "trace"
        return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo[1], (trace, background, weight, nframes))
    elseif datatype == "tracerna"
        len, h = read_rna(gene, datacond, datapath[2])
        return TraceRNAData(label, gene, traceinfo[1], (trace, background, weight, nframes), len, h)
    end
end

"""
    load_data_tracejoint(datapath, label, gene, datacond, traceinfo)

data structure for joint traces
"""
function load_data_tracejoint(datapath, label, gene, datacond, traceinfo)
    trace = read_tracefiles(datapath, datacond, traceinfo)
    weight = Float64[]
    for f in traceinfo[4]
        # push!(weight, (1 - f) / f * length(trace))
        push!(weight, 1 - f)
    end
    nframes = traceinfo[3] < 0 ? floor(Int, (720 - traceinfo[2] + traceinfo[1]) / traceinfo[1]) : floor(Int, (traceinfo[3] - traceinfo[2] + traceinfo[1]) / traceinfo[1])
    return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo[1], (trace, Vector[], weight, nframes))
end

"""
    load_data_grid(datapath, label, gene, datacond, traceinfo)

TBW
"""
function load_data_tracegrid(datapath, label, gene, datacond, traceinfo)
    trace = read_tracefiles_grid(datapath, datacond, traceinfo)
    # nframes = traceinfo[3] < 0 ? floor(Int, (720 - traceinfo[2] + traceinfo[1]) / traceinfo[1]) : floor(Int, (traceinfo[3] - traceinfo[2] + traceinfo[1]) / traceinfo[1])
    return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo[1], (trace, Vector[], 0.0, 1))
end

"""
    load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna)

return data structure
"""
function load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna)
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
    elseif datatype == "dwelltime"
        bins, DT = read_dwelltimes(datapath)
        return DwellTimeData(label, gene, bins, DT, dttype)
    elseif occursin("trace", datatype)
        if datatype == "tracejoint"
            load_data_tracejoint(datapath, label, gene, datacond, traceinfo)
        elseif datatype == "tracegrid"
            load_data_tracegrid(datapath, label, gene, datacond, traceinfo)
        else
            load_data_trace(datapath, label, gene, datacond, traceinfo, datatype)
        end
    else
        throw("$datatype not included")
    end
end

"""
    GRSMgridmodel(data, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, grid)

TBW
"""
function GRSMgridmodel(data, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, grid::Int)
    reporter, components = make_reporter_components(data, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
    n = num_rates(transitions, R, S, insertstep)
    raterange = 1:n
    noiserange = n+1:n+length(noisepriors)
    gridrange = n+length(noisepriors)+1:n+length(noisepriors)+1
    priord = prior_distribution(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors)
    propcv < 0 && (propcv = getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root))
    GRSMgridmodel{typeof(r),Nothing,typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, raterange, noiserange, gridrange, grid, nothing, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

"""
    load_model_hierarchical(data::AbstractExperimentalData, r, rm, transitions, G::Int, R::Int, S::Int, insertstep::Int, nalleles, priorcv, splicetype, propcv, fittedparam, fixedeffects, method, components, reporter, noisepriors, hierarchical::Tuple)

hierarchical model
"""
# function load_model_hierarchical(data::AbstractExperimentalData, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, nalleles, priorcv, onstates, decayrate, propcv, splicetype, probfn, noisepriors, method, hierarchical::Tuple)

function GRSMhierarchicalmodel(data::AbstractExperimentalData, r, rm, fittedparam, fixedeffects, transitions, G::Int, R::Int, S::Int, insertstep::Int, splicetype, nalleles, priorcv, propcv, method, noisepriors, hierarchical::Tuple, components, reporter)
    nhypersets = hierarchical[1]
    nrates = num_rates(transitions, R, S, insertstep) + reporter.n
    nindividuals = length(data.trace[1])
    nparams = length(hierarchical[2])
    ratestart = nhypersets * nrates + 1
    paramstart = length(fittedparam) + nhypersets * nparams + 1
    fittedparam, fittedhyper, fittedpriors = make_fitted_hierarchical(fittedparam, hierarchical[1], hierarchical[2], nrates, nindividuals)
    fixedeffects = make_fixed(fixedeffects, hierarchical[3], nrates, nindividuals)
    rprior = rm[1:nhypersets*nrates]
    priord = prior_distribution(rprior, transitions, R, S, insertstep, fittedpriors, priorcv, noisepriors)
    pop = Hierarchy(nhypersets, nrates, nparams, nindividuals, ratestart, paramstart, fittedhyper)
    GRSMhierarchicalmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, pop, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

"""
    GRSMgridhierarchicalmodel(data::AbstractExperimentalData, r, rm, fittedparam, fixedeffects, transitions, G::Int, R::Int, S::Int, insertstep::Int, splicetype, nalleles, priorcv, propcv, method, noisepriors, hierarchical::Tuple, components, reporter)

TBW
"""
function GRSMgridhierarchicalmodel(data::AbstractExperimentalData, r, rm, fittedparam, fixedeffects, transitions, G::Int, R::Int, S::Int, insertstep::Int, splicetype, nalleles, priorcv, propcv, method, noisepriors, hierarchical::Tuple, components, reporter)
    nhypersets = hierarchical[1]
    nrates = num_rates(transitions, R, S, insertstep) + reporter.n
    raterange = 1:nrates
    noiserange = n+1:n+length(noisepriors)
    gridrange = n+length(noisepriors)+1:n+length(noisepriors)+1
    nindividuals = length(data.trace[1])
    nparams = length(hierarchical[2])
    ratestart = nhypersets * nrates + 1
    paramstart = length(fittedparam) + nhypersets * nparams + 1
    fittedparam, fittedhyper, fittedpriors = make_fitted_hierarchical(fittedparam, hierarchical[1], hierarchical[2], nrates, nindividuals)
    fixedeffects = make_fixed(fixedeffects, hierarchical[3], nrates, nindividuals)
    rprior = rm[1:nhypersets*nrates]
    priord = prior_distribution(rprior, transitions, R, S, insertstep, fittedpriors, priorcv, noisepriors)
    hyper = Hierarchy(nhypersets, nrates, nparams, nindividuals, ratestart, paramstart, fittedhyper)
    GRSMgridhierarchicalmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, hyper, raterange, noiserange, gridrange, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

"""
    GRSMcoupledmodel(r::Vector, rm::Vector, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, probfn, noisepriors, method, coupling)

Coupled model
"""
function GRSMcoupledmodel(r::Vector, rm::Vector, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, probfn, noisepriors, method, coupling)
    # function load_model_coupled(r, rm, fittedparam, fixedeffects, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, coupling::Tuple, nalleles, priorcv, propcv, splicetype, probfn, noisepriors, method=Tsit5())
    println(coupling)
    reporter = HMMReporter[]
    !(probfn isa Union{Tuple,Vector}) && (probfn = fill(probfn, length(coupling[1])))
    n_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
    for i in eachindex(G)
        nnoise = length(noisepriors[i])
        weightind = occursin("Mixture", "$(probfn)") ? num_rates(transitions[i], R[i], S[i], insertstep[i]) + nnoise : 0
        push!(reporter, HMMReporter(nnoise, n_per_state[i], probfn[i], weightind, off_states(n_per_state[i])))
    end
    priord = prior_distribution_coupling(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors)
    # components = make_components_TRGCoupled(coupling, transitions, G, R, S, insertstep, "")
    components = TRGCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    GRSMcoupledmodel{typeof(r),Int,typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, coupling[5], transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

function GRSMcoupledhierarchicalmodel(r::Vector, rm::Vector, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, probfn, noisepriors, method, coupling)
    # function load_model_coupled(r, rm, fittedparam, fixedeffects, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, coupling::Tuple, nalleles, priorcv, propcv, splicetype, probfn, noisepriors, method=Tsit5())
    println(coupling)
    reporter = HMMReporter[]
    !(probfn isa Union{Tuple,Vector}) && (probfn = fill(probfn, length(coupling[1])))
    n_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
    for i in eachindex(G)
        nnoise = length(noisepriors[i])
        weightind = occursin("Mixture", "$(probfn)") ? num_rates(transitions[i], R[i], S[i], insertstep[i]) + nnoise : 0
        push!(reporter, HMMReporter(nnoise, n_per_state[i], probfn[i], weightind, off_states(n_per_state[i])))
    end
    priord = prior_distribution_coupling(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors)
    # components = make_components_TRGCoupled(coupling, transitions, G, R, S, insertstep, "")
    components = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    GRSMcoupledhierarchicalmodel{typeof(r),Int,typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, coupling[5], transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

function GRSMcoupledmodel(data::AbstractTraceData, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, probfn, noisepriors, method, coupling)
    reporter, components = make_reporter_components(data, transitions, G, R, S, insertstep, probfn, noisepriors, coupling)
    priord = prior_distribution_coupling(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors)
    GRSMcoupledmodel{typeof(r),Int,typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, coupling[5], transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end

function GRSMcoupledmodel(data::DwellTimeData, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, probfn, noisepriors, method, coupling)
    reporter, components = make_reporter_components(transitions::Tuple, G::Tuple, R, S, insertstep, onstates, data.DTtype, splicetype, coupling)
    priord = prior_distribution_coupling(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors)
    GRSMcoupledmodel{typeof(r),Int,typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, coupling[5], transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
end
"""
    load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical)

TBW
"""
function load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical)
    reporter, components = make_reporter_components(data, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
    if isempty(hierarchical)
        checklength(r, transitions, R, S, insertstep, reporter)
        priord = prior_distribution(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors)
        if propcv < 0
            propcv = getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root)
        end
        if R == 0
            return GMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, nalleles, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        else
            return GRSMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        end
    else
        !isa(method, Tuple) && throw("method not a Tuple")
        GRSMhierarchicalmodel(data, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, method, noisepriors, hierarchical, components, reporter)
        # load_model_hierarchical(data, r, rm, transitions, G, R, S, insertstep, nalleles, priorcv, splicetype, propcv, fittedparam, fixedeffects, method, components, reporter, noisepriors, hierarchical)
    end
end
"""
    load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling)

return model structure
"""
function load_model(data, r, rm, fittedparam::Vector, fixedeffects::Tuple, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid)
    if !isempty(coupling)
        GRSMcoupledmodel(r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, probfn, noisepriors, method, coupling)
    elseif !isnothing(grid)
        if !isempty(hierarchical)
            !isa(method, Tuple) && throw("method not a Tuple")
            GRSMgridhierarchicalmodel(data, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, propcv, method, noisepriors, hierarchical, components, reporter)
        else
            GRSMgridmodel(data, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, prob_Gaussian_grid, noisepriors, method, grid)
        end
    else
        load_model(data, r, rm, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical)
    end
end
"""
    make_reporter_components(data, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)

TBW
"""

function make_reporter_components(data::AbstractRNAData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
    reporter = onstates
    components = MComponents(transitions, G, R, data.nRNA, decayrate, splicetype)
    return reporter, components
end

function make_reporter_components(data::Union{AbstractTraceData,AbstractTraceHistogramData}, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
    nnoise = length(noisepriors)
    weightind = occursin("Mixture", "$(probfn)") ? num_rates(transitions, R, S, insertstep) + nnoise : 0
    reporter = HMMReporter(nnoise, num_reporters_per_state(G, R, S, insertstep), probfn, weightind, off_states(G, R, S, insertstep))
    if typeof(data) <: TraceRNAData
        components = MTComponents(transitions, G, R, S, insertstep, data.nRNA, decayrate, splicetype)
    else
        components = TComponents(transitions, G, R, S, insertstep, splicetype)
    end
    return reporter, components
end

function make_reporter_components(transitions::Tuple, G::Int, R, S, insertstep, onstates, dttype, splicetype)
    sojourn = sojourn_states(onstates, G, R, S, insertstep, dttype)
    components = TDComponents(transitions::Tuple, G, R, S, insertstep, sojourn, dttype, splicetype)
    nonzeros = nonzero_rows(components)
    return (sojourn, nonzeros), components
end

function make_reporter_components(transitions::Tuple, G::Tuple, R, S, insertstep, onstates, dttype, splicetype, coupling)
    sojourn = sojourn_states(onstates, G, R, S, insertstep, dttype)
    components = TDCoupledComponents(coupling, transitions, G, R, S, insertstep, sojourn, dttype, splicetype)
    sojourn = coupled_states(sojourn, coupling, components, G)
    nonzeros = coupled_states(nonzero_rows(components), coupling, components, G)
    return (sojourn, nonzeros), components
end

function make_reporter_components(data::DwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
    make_reporter_components(transitions, G, R, S, insertstep, onstates, Data.DTtypes, splicetype)
end


function make_reporter_components(data::RNADwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
    reporter, tcomponents = make_reporter_components(transitions, G, R, S, insertstep, onstates, data.DTtypes, "")
    components = MTDComponents(MComponents(transitions, G, R, data.nRNA, decayrate, splicetype), tcomponents)
    return reporter, components
end

function make_reporter_components(data::RNAOnOffData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors)
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
            components = MTDComponents(transitions, G, R, S, insertstep, onstates, data.DTtypes, data.nRNA, decayrate, splicetype)
        else
            throw("length of onstates and data.DTtypes not the same")
        end
    else
        components = MTAIComponents(transitions, G, R, S, insertstep, onstates, data.nRNA, decayrate, splicetype)
    end
    return reporter, components
end

function make_reporter_components(data::AbstractTraceData, transitions, G::Tuple, R, S, insertstep, probfn, noisepriors, coupling::Tuple)
    println(coupling)
    reporter = HMMReporter[]
    !(probfn isa Union{Tuple,Vector}) && (probfn = fill(probfn, length(coupling[1])))
    n_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1])
    for i in eachindex(G)
        nnoise = length(noisepriors[i])
        weightind = occursin("Mixture", "$(probfn)") ? num_rates(transitions[i], R[i], S[i], insertstep[i]) + nnoise : 0
        push!(reporter, HMMReporter(nnoise, n_per_state[i], probfn[i], weightind, off_states(n_per_state[i])))
    end
    # components = make_components_TRGCoupled(coupling, transitions, G, R, S, insertstep, "")
    components = TCoupledComponents(coupling, transitions, G, R, S, insertstep, "")
    return reporter, components
end

"""
    checklength(r, transitions, R, S, insertstep, reporter)

TBW
"""
function checklength(r, transitions, R, S, insertstep, reporter)
    n = num_rates(transitions, R, S, insertstep)
    if typeof(reporter) <: HMMReporter
        (length(r) != n + reporter.n) && throw("r has wrong length")
    else
        (length(r) != n) && throw("r has wrong length")
    end
    nothing
end

function prior_hypercv(transitions, R::Int, S, insertstep, noisepriors)
    [fill(1.0, length(transitions)); 1.0; fill(0.1, R - 1); 1.0; fill(1.0, max(0, S - insertstep + 1)); 1.0; fill(0.1, length(noisepriors))]
end

function prior_hypercv(transitions, R::Tuple, S, insertstep, noisepriors)
    rm = Float64[]
    for i in eachindex(R)
        append!(rm, prior_hypercv(transitions[i], R[i], S[i], insertstep[i], noisepriors[i]))
    end
    rm
end


"""
    prior_ratemean(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors::Vector, elongationtime::Float64)

default priors for rates (includes all parameters, fitted or not)
"""
function prior_ratemean(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors::Vector, elongationtime::Float64, initprior::Float64=0.1)
    [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.1, max(0, S - insertstep + 1)); decayrate; noisepriors]
end

"""
    prior_ratemean(transitions, R::Tuple, S::Tuple, insertstep::Tuple, decayrate, noisepriors::Union{Vector,Tuple}, elongationtime::Union{Vector,Tuple}, coupling, initprior=[0.1, 0.1])

TBW
"""
function prior_ratemean(transitions, R::Tuple, S::Tuple, insertstep::Tuple, decayrate, noisepriors::Union{Vector,Tuple}, elongationtime::Union{Vector,Tuple}, coupling, initprior=[0.1, 0.1])
    rm = Float64[]
    for i in eachindex(R)
        append!(rm, prior_ratemean(transitions[i], R[i], S[i], insertstep[i], decayrate, noisepriors[i], elongationtime[i], initprior[i]))
    end
    [rm; fill(0.0, coupling[5])]
end

"""
    prior_ratemean_hierarchical(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, nhypersets, coupling=tuple(), cv::Float64=1.0)

default priors for hierarchical models, arranged into a single vector, shared and hyper parameters come first followed by individual parameters
"""
function prior_ratemean_hierarchical(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, nhypersets, coupling=tuple(), cv::Float64=1.0)
    r = isempty(coupling) ? prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime) : prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling)
    hypercv = prior_hypercv(transitions, R, S, insertstep, noisepriors)
    # hypercv = [fill(1.0, length(transitions)); 1.0; fill(0.1, R - 1); 1.0; fill(1.0, max(0, S - insertstep + 1)); 1.0; fill(.1,length(noisepriors))]
    append!(r, hypercv)
    for i in 3:nhypersets
        append!(r, cv .* rm)
    end
    r
end

"""
    prior_ratemean_grid(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors::Vector, elongationtime::Float64)

TBW
"""
function prior_ratemean_grid(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime)
    [prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime); 0.5]
    # [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.1, max(0, S - insertstep + 1)); decayrate; noisepriors; 0.5]
end

"""
    set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid)

set priormean if empty
"""
function set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid)
    if !isempty(priormean)
        return priormean
    else
        if !isnothing(grid)
            return prior_ratemean_grid(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime)
        end
        if isempty(hierarchical)
            if isempty(coupling)
                return prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime)
            else
                return prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling)
            end
        else
            return prior_ratemean_hierarchical(transitions, R, S, insertstep, decayrate, noisepriors, hierarchical[1], elongationtime, coupling)
        end
    end
end

"""
    prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, priorcv, noisepriors, factor=10)
    prior_distribution(rm, transitions, R::Tuple, S::Tuple, insertstep::Tuple, fittedparam::Vector, priorcv, noisepriors, factor=10)


"""
function prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, priorcv, noisepriors, factor=10)
    if isempty(rm)
        throw("No prior mean")
        # rm = prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors)
    end
    if priorcv isa Number
        n = num_rates(transitions, R, S, insertstep)
        rcv = fill(priorcv, length(rm))
        rcv[n] = 0.1
        rcv[n+1:n+length(noisepriors)] /= factor
    else
        rcv = priorcv
    end
    if length(rcv) == length(rm)
        return distribution_array(log.(rm[fittedparam]), sigmalognormal(rcv[fittedparam]), Normal)
    else
        throw("priorcv not the same length as prior mean")
    end
end
function prior_distribution_coupling(rm, transitions, R::Tuple, S::Tuple, insertstep::Tuple, fittedparam::Vector, priorcv, noisepriors, factor=10)
    if isempty(rm)
        throw("No prior mean")
        # rm = prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors)
    end
    if priorcv isa Number
        rcv = fill(priorcv, length(rm))
        s = 0
        for i in eachindex(R)
            lnp = isempty(noisepriors) ? 0 : length(noisepriors[i])
            n = num_rates(transitions[i], R[i], S[i], insertstep[i])
            rcv[s+n] = 0.1
            rcv[s+n+1:s+n+lnp] ./= factor
            s += n + lnp
        end
    else
        rcv = priorcv
    end
    if length(rcv) == length(rm)
        return distribution_array(transform_array(rm[fittedparam], length(rm), fittedparam, logv, log_shift1), sigmalognormal(rcv[fittedparam]), Normal)
    else
        throw("priorcv not the same length as prior mean")
    end
end

# function set_rinit(infolder, inlabel,  gene, priormean, transitions,G, R, S, insertstep, nalleles, ratetype, hierarchical)
"""
    set_rinit(r, priormean)

set rinit to prior if empty
"""
function set_rinit(r, priormean)
    if isempty(r)
        println("No rate file, set rate to prior")
        r = priormean
    end
    println(r)
    r
end


"""
    set_rinit(rm, transitions, R::Int, S::Int, insertstep, noisepriors, nindividuals)

set rinit for hierarchical models
"""
function set_rinit(r, priormean, transitions, R::Int, S::Int, insertstep, noisepriors, nindividuals)
    if isempty(r)
        println("No rate file, set rate to prior")
        r = copy(priormean)
        nrates = num_rates(transitions, R, S, insertstep) + length(noisepriors)
        for i in 1:nindividuals
            append!(r, priormean[1:nrates])
        end
    end
    r
end



"""
    num_noiseparams(datatype, noisepriors)

TBW
"""
function num_noiseparams(datatype, noisepriors)
    if occursin("trace", lowercase(datatype))
        if eltype(noisepriors) <: Number
            return length(noisepriors)
        else
            return length.(noisepriors)
        end
    else
        return 0
    end
end



"""
    default_fitted(datatype::String, transitions, R::Tuple, S::Tuple, insertstep::Tuple, noiseparams::Tuple, coupling)

create vector of fittedparams that includes all rates except the decay time
"""
function default_fitted(datatype::String, transitions, R::Tuple, S::Tuple, insertstep::Tuple, noiseparams::Tuple, coupling, grid)
    fittedparam = Int[]
    totalrates = 0
    for i in eachindex(R)
        fittedparam = vcat(fittedparam, totalrates .+ default_fitted(datatype, transitions[i], R[i], S[i], insertstep[i], noiseparams[i], coupling, grid))
        totalrates += num_rates(transitions[i], R[i], S[i], insertstep[i]) + noiseparams[i]
    end
    [fittedparam; collect(fittedparam[end]+1:fittedparam[end]+coupling[5])]
end

"""
    default_fitted(datatype::String, transitions, R::Int, S, insertstep, noiseparams, coupling)

TBW
"""
function default_fitted(datatype::String, transitions, R::Int, S, insertstep, noiseparams, coupling, grid)
    n = num_rates(transitions, R, S, insertstep)
    fittedparam = collect(1:n-1)
    if occursin("trace", datatype)
        if isnothing(grid)
            isempty(noiseparams) && throw("noisepriors cannot be empty for trace data")
            fittedparam = vcat(fittedparam, collect(n+1:n+noiseparams))
        else
            fittedparam = vcat(fittedparam, collect(n+1:n+noiseparams+1))
        end
    end
    fittedparam
end

"""
    set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noiseparams, coupling)

TBW
"""
function set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noisepriors, coupling, grid)
    if isempty(fittedparam)
        return default_fitted(datatype, transitions, R, S, insertstep, num_noiseparams(datatype, noisepriors), coupling, grid)
    else
        return fittedparam
    end
end



"""
    make_fitted_hierarchical(fittedparams,N)

make fittedparams vector for hierarchical model

returns 
-`f`: all fitted parameters
-`fhyper`: fitted hyper parameters (needed for hyper distribution)
-`fpriors`: fitted shared and hyper parameters (needed for specifying prior distributions)

Hierarchical models have shared, individual, and hyper parameters.  shared parameters pertain to all individuals, individual parameters are fitted to each individual, and
hyper parameters specify the hyper distribution for the individual parameters. Currently, all hyper parameters are fitted as determined by the fitted individual parameters
shared and hyper parameters have priors, whereas individual parameters are drawn from the hyper distribution
"""
function make_fitted_hierarchical(fittedshared, nhypersets, fittedindividual, nrates, nindividuals)
    f = [fittedshared; fittedindividual] # shared parameters come first followed by hyper parameters and individual parameters
    f = sort(f)
    fhyper = [fittedindividual] # fitted hyper parameters correspond to fitted individual parameters
    for i in 1:nhypersets-1
        append!(f, fittedindividual .+ i * nrates)
        push!(fhyper, fittedindividual .+ i * nrates)
    end
    fpriors = copy(f) # priors only apply to shared and hyper parameters
    for i in 1:nindividuals
        append!(f, fittedindividual .+ (i + nhypersets - 1) * nrates)
    end
    f, fhyper, fpriors
end

##### Under construction - specifying fitted hyper parameters independently of individual parameters
# function make_fitted_hierarchical(fittedshared, fittedhyper, nhypersets, fittedindividual, nrates, nindividuals)
#     fittedall = [fittedshared; fittedhyper; fittedindividual] # shared parameters come first followed by hyper parameters and individual parameters
#     fittedall = sort(fittedall)
#     fhyper = [fittedindividual]
#     for i in 1:nhypersets-1
#         append!(fittedall, fittedindividual .+ i * nrates)
#         push!(fhyper, fittedindividual .+ i * nrates)
#     end
#     fpriors = copy(f)
#     for i in 1:nindividuals
#         append!(fittedall, fittedindividual .+ (i + nhypersets - 1) * nrates)
#     end
#     fittedall, fhyper, fpriors
# end

"""
    make_fixed(fixedpop, fixedindividual, nrates, nindividuals)

make fixed effects tuple for hierarchical model
"""
function make_fixed(fixedshared, fixedindividual, nrates, nindividuals)
    fixed = Vector{Int}[]
    for f in fixedshared
        push!(fixed, f)
    end
    for h in fixedindividual
        push!(fixed, [h + i * nrates for i in 0:nindividuals-1])
    end
    tuple(fixed...)
end

"""
    make_fixedfitted(fixedeffects::String,transitions,R,S,insertstep)

make default fixedeffects tuple and fittedparams vector from fixedeffects String
"""
function make_fixedfitted(datatype, fixedeffects::String, transitions, R, S, insertstep, noiseparams, coupling, grid)
    fittedparam = default_fitted(datatype, transitions, R, S, insertstep, noiseparams, coupling, grid)
    fixed = split(fixedeffects, "-")
    if length(fixed) > 1
        fixed = parse.(Int, fixed)
        deleteat!(fittedparam, fixed[2:end])
        fixed = tuple(fixed)
    else
        fixed = tuple()
    end
    return fixed, fittedparam
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
    M = make_mat_M(MComponents([(2, 1)], 2, R, total, r[end], ""), r)
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
function check_genename(gene::String, p1)
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
    check_genename(gene::Vector, p1)


"""
function check_genename(gene::Vector, p1)
    for i in eachindex(gene)
        gene[i] = check_genename(gene[i], p1)
    end
    return gene
end

"""
    print_ll(param, data, model, message)

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
# function printinfo(gene, G, datacond, datapath, infolder, resultfolder, maxtime)
#     println("Gene: ", gene, " G: ", G, " Treatment:  ", datacond)
#     println("data: ", datapath)
#     println("in: ", infolder, " out: ", resultfolder)
#     println("maxtime: ", maxtime)
# end

function printinfo(gene, G::Int, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime)
    if R == 0
        println("Gene: ", gene, " G: ", G, " Treatment:  ", datacond)
    else
        println("Gene: ", gene, "Label: ", datacond, " G R S insertstep: ", G, R, S, insertstep)
    end
    println("data: ", datapath)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime)
end

function printinfo(gene, G::Tuple, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime)
    println("Gene: ", gene)
    for i in eachindex(G)
        println("Label: ", datacond[i], " G R S insertstep: ", G[i], R[i], S[i], insertstep[i])
    end
    println("data: ", datapath)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime)
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

"""
    getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root, verbose=true)

TBW
"""
function getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root, verbose=true)
    paramfile = path_Gmodel("param-stats", gene, G, nalleles, inlabel, infolder, root)
    if isfile(paramfile)
        cv = read_covlogparam(paramfile)
        cv = float.(cv)
        if !isposdef(cv) || size(cv)[1] != length(fittedparam)
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
    get_decay(a, gene::String)
    get_decay(a::Float64)

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
    if !isnothing(ind)
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
