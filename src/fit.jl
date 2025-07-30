# This file is part of StochasticGene.jl   

# fit.jl
#
# Fit GRS models (generalized telegraph models) to RNA abundance and live cell imaging data
#
"""
    genes_hbec()

Return a vector of HBEC gene names.

# Returns
- `Vector{String}`: Array of HBEC gene names

# Notes
- Contains predefined list of genes used in HBEC cell line studies
- Used for gene-specific parameter lookups
"""
genes_hbec() = ["CANX"; "DNAJC5"; "ERRFI1"; "KPNB1"; "MYH9"; "Rab7a"; "RHOA"; "RPAP3"; "Sec16A"; "SLC2A1"]

"""
    genelength_hbec()

Return a dictionary mapping HBEC gene names to their lengths.

# Returns
- `Dict{String, Int}`: Mapping of gene names to gene lengths in base pairs

# Notes
- Contains gene length information for HBEC genes
- Used for gene-specific calculations and parameter estimation
"""
genelength_hbec() = Dict([("Sec16A", 42960); ("SLC2A1", 33802); ("ERRFI1", 14615); ("RHOA", 52948); ("KPNB1", 33730); ("MYH9", 106741); ("DNAJC5", 40930); ("CANX", 32710); ("Rab7a", 88663); ("RPAP3", 44130); ("RAB7A", 88663); ("SEC16A", 42960)])

"""
    MS2end_hbec()

Return a dictionary mapping HBEC gene names to their MS2 end positions.

# Returns
- `Dict{String, Int}`: Mapping of gene names to MS2 end positions

# Notes
- Contains MS2 binding site end positions for HBEC genes
- Used for MS2-based transcription analysis
"""
MS2end_hbec() = Dict([("Sec16A", 5220); ("SLC2A1", 26001); ("ERRFI1", 5324); ("RHOA", 51109); ("KPNB1", 24000); ("MYH9", 71998); ("DNAJC5", 14857); ("CANX", 4861); ("Rab7a", 83257); ("RPAP3", 38610); ("SEC16A", 5220); ("RAB7A", 83257)])

"""
    halflife_hbec()

Return a dictionary mapping HBEC gene names to their half-lives.

# Returns
- `Dict{String, Float64}`: Mapping of gene names to half-lives in minutes

# Notes
- Contains mRNA half-life information for HBEC genes
- Used for decay rate calculations and model parameterization
"""
halflife_hbec() = Dict([("CANX", 50.0), ("DNAJC5", 5.0), ("ERRFI1", 1.35), ("KPNB1", 9.0), ("MYH9", 10.0), ("Rab7a", 50.0), ("RHOA", 50.0), ("RPAP3", 7.5), ("Sec16A", 8.0), ("SLC2A1", 5.0), ("RAB7A", 50.0), ("SEC16A", 8.0)])

"""
    get_transitions(G::Int, TransitionType)

Return the default transitions tuple for a given number of gene states `G` and a transition type string.

# Arguments
- `G::Int`: Number of gene states.
- `TransitionType`: String describing the transition type (e.g., "KP", "cyclic").

# Returns
- Tuple of allowed transitions for the specified model.
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
        throw(ArgumentError("Transition type unknown for G = $G"))
    end
end


const SUPPORTED_DATATYPES = Set([
    :rna, :rnacount, :rnaonoff, :rnadwelltime,
    :dwelltime, :trace, :tracerna, :tracejoint, :tracegrid
])

"""
    normalize_datatype(datatype::Symbol) -> Symbol

Normalize and validate a user-supplied or internal data type symbol.

# Arguments
- `datatype::Symbol`: Data type symbol (e.g., :rna, :tracegrid).

# Returns
- Validated lowercase `Symbol`.

# Throws
- `ArgumentError` if the type is unsupported.
"""
function normalize_datatype(datatype::Symbol)
    datatype ∈ SUPPORTED_DATATYPES || error(
        "Unsupported datatype '$datatype'. Supported types: $(join(SUPPORTED_DATATYPES, ", "))"
    )
    return datatype
end

"""
    normalize_datatype(datatype::AbstractString) -> Symbol

Normalize and validate a user-supplied or internal data type string.

# Arguments
- `datatype::AbstractString`: Data type string (e.g., "RNA", "tracegrid").

# Returns
- Validated lowercase `Symbol`.

# Throws
- `ArgumentError` if the type is unsupported.
"""
function normalize_datatype(datatype::AbstractString)
    return normalize_datatype(Symbol(lowercase(datatype)))
end

"""
    get_transitions(G::Tuple, TransitionType)

Return a tuple of transitions for each element in `G` (for coupled models).

# Arguments
- `G::Tuple`: Tuple of gene state counts.
- `TransitionType`: String describing the transition type.

# Returns
- Tuple of transitions for each model.
"""
function get_transitions(G::Tuple, TransitionType)
    Tuple([get_transitions(G, TransitionType) for G in G])
end

"""
    make_coupling(source::UnitRange{Int64}=1:3, target::UnitRange{Int64}=1:3)

Construct a default coupling structure for coupled models.

# Arguments
- `source::UnitRange{Int64}`: Source unit indices (default 1:3)
- `target::UnitRange{Int64}`: Target unit indices (default 1:3)

# Returns
- `Array`: Array of coupling tuples

# Notes
- Creates coupling structure for models with multiple interacting units
- Each coupling tuple contains source and target information
- Used for specifying how different model units influence each other
"""
function make_coupling(source::UnitRange{Int64}=1:3, target::UnitRange{Int64}=1:3)
    coupling = []
    for s in source, t in target
        push!(coupling, ((1, 2), (tuple(), tuple(1)), (source, 0), (0, target), 1))
    end
    return coupling
end

"""
    set_decayrate(decayrate::Float64, gene, cell, root)

Return the decay rate, using the provided value or looking it up if negative.

# Arguments
- `decayrate::Float64`: Decay rate (if < 0, will be looked up)
- `gene`: Gene name
- `cell`: Cell type
- `root`: Root directory

# Returns
- `Float64`: Decay rate value

# Notes
- If decayrate < 0, looks up gene-specific decay rate from halflife files
- Otherwise returns the provided decayrate value
- Used for setting mRNA degradation rates in models
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

Return a tuple of decay rates, using the provided values or looking them up if negative.

# Arguments
- `decayrate::Tuple`: Tuple of decay rates (if any < 0, will be looked up)
- `gene`: Gene name(s)
- `cell`: Cell type
- `root`: Root directory

# Returns
- `Tuple`: Tuple of decay rates

# Notes
- Processes each decay rate in the tuple individually
- For coupled models with multiple genes, handles gene-specific lookups
- If any decay rate < 0, looks up from halflife files
- Used for multi-unit models requiring different decay rates
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

Adjust the number of splice sites S if it exceeds R - insertstep + 1.

# Arguments
- `S::Int`: Number of splice sites
- `R::Int`: Number of pre-RNA steps
- `insertstep::Int`: Reporter insertion step

# Returns
- `Int`: Adjusted value of S

# Notes
- Ensures S does not exceed the maximum possible value given R and insertstep
- Maximum S = R - insertstep + 1 (number of RNA steps after insertion)
- Prints warning message if adjustment is made
- Used to prevent invalid model configurations
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

Adjust each element of S if it exceeds R - insertstep + 1 for coupled models.

# Arguments
- `S::Tuple`: Tuple of splice site counts
- `R::Tuple`: Tuple of pre-RNA step counts
- `insertstep::Tuple`: Tuple of reporter insertion steps

# Returns
- `Tuple`: Tuple of adjusted S values

# Notes
- Processes each unit in coupled models individually
- Converts tuple to array for modification, then back to tuple
- Prints warning message for each unit that requires adjustment
- Used for multi-unit models to ensure valid configurations
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

Set nalleles to 1 if coupling is nonempty (for coupled models).

# Arguments
- `nalleles`: Number of alleles
- `coupling`: Coupling structure (tuple)

# Returns
- Adjusted nalleles value

# Notes
- For coupled models, forces nalleles to 1 since coupling handles multiple units
- Prints warning message when adjustment is made
- Used to prevent conflicts between coupling and allele counting
- In coupled models, multiple units are treated as separate "alleles"
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

Fit steady state or transient GM/GRSM model to RNA data for a single gene, write the result (through function finalize), and return fit results and diagnostics.

For coupled transcribing units, arguments transitions, G, R, S, insertstep, and trace become tuples of the single unit type, e.g. If two types of transcription models are desired with G=2 and G=3 then G = (2,3).

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
- `datacol=3`: column of data to use, default is 3 for rna data
- `datatype::String=""`: String that describes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracerna", "tracejoint", "tracegrid"
- `datacond=""`: string or vector of strings describing data, e.g. "WT", "DMSO" or ["DMSO","AUXIN"], ["gene","enhancer"]
- `datapath=""`: path to data file or folder or array of files or folders
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `ejectnumber=1`: number of mRNAs produced per burst, default is 1, if Int then deterministic, if Tuple = (r, p) then stochastic obeying NegativeBinomial(r, p)
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
- `maxtime=60`: maximum wall time for run, default = 60 min
- `method=Tsit5()`: DifferentialEquations.jl numerical method (e.g. Tsit5(), Tsit5(),...); use a tuple for hierarchical models: method = tuple(method, Bool) = (numerical method (currently not used), true if transition rates are shared)
- `nalleles=1`: number of alleles, value in alleles folder will be used if it exists, for coupled models, nalleles is only used when computing steady state RNA histograms and considered uncoupled.  For add coupled alleles as units and set nalleles to 1.
- `nchains::Int=2`: number of MCMC chains = number of processors called by Julia, default = 2
- `noisepriors=[]`: priors of observation noise (use empty set if not fitting traces), superseded if priormean is set
- `onstates=Int[]`: vector of on or sojourn states, e.g. [2], if multiple onstates are desired, use vector of vectors, e.g. [[2,3],Int[]], use empty vector for R states or vector of empty vectors for R states in coupled models, do not use Int[] for R=0 models
- `optimize=false`: use optimizer to compute maximum likelihood value
- `priormean=Float64[]`: mean rates of prior distribution (must set priors for all rates including those that are not fitted)
- `priorcv=10.`: (vector or number) coefficient of variation(s) for the rate prior distributions, default is 10.
- `probfn=prob_Gaussian`: probability function for HMM observation probability (i.e., noise distribution), tuple of functions for each unit, e.g. (prob_Gaussian, prob_Gaussian) for coupled models, use 1 for forced (e.g. one unit drives the other)
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
- `traceinfo=(1.0, 1., -1, 1., 0.5)`: 5 tuple = (frame interval of intensity traces in minutes, starting frame time in minutes, ending frame time (use -1 for last index), fraction of observed active traces, background mean)
    for simultaneous joint traces, the fraction of active traces is a vector of the active fractions for each trace, e.g. (1.0, 1., -1, [.5, .7], [0.5,0.5]) 
    If active fraction is 1.0, then traceinfo can be a 3-tuple, e.g. (1.0, 1., -1) since background correction is not needed
    Note that all traces are scaled by the maximum of the medians of all the traces, the traces are all scaled by the same factor since the signal amplitude should be the same
- `TransitionType=""`: String describing G transition type, e.g. "3state", "KP" (kinetic proofreading), "cyclic", or if hierarchical, coupled
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2-state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3-state kinetic proofreading model, empty for G=1
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `writesamples=false`: write out MH samples if true, default is false
- `zeromedian=false`: if true, subtract the median of each trace from each trace, then scale by the maximum of the medians

# Returns
- `fits`: MCMC fit results (posterior samples, log-likelihoods, etc.)
- `stats`: Summary statistics for parameters
- `measures`: Diagnostic measures (including WAIC and its standard error, which is now for the total WAIC and scaled by sqrt(n_obs))
- `data`, `model`, `options`: The data, model, and options structures used

# Notes
- If `propcv < 0`, proposal covariance is read from previous run if available.
- WAIC standard error is for the total WAIC (not per observation), and is scaled by sqrt(n_obs).
- File and folder conventions may have changed; see README for details.

# Example
If you are in the folder where data/HCT116_testdata is installed, you can fit the mock RNA histogram running 4 MCMC chains with:

```julia
fits = fit(
    G = 2,
    R = 0,
    transitions = ([1,2], [2,1]),
    datatype = "rna",
    datapath = "data/HCT116_testdata/",
    gene = "MYC",
    datacond = "MOCK"
)
```

Trace data fit:

```julia
fits = fit(
    G = 3,
    R = 2,
    S = 2,
    insertstep = 1,
    transitions = ([1,2], [2,1], [2,3], [3,1]),
    datatype = "trace",
    datapath = "data/testtraces",
    cell = "TEST",
    gene = "test",
    datacond = "testtrace",
    traceinfo = (1.0, 1., -1, 1.),
    noisepriors = [40., 20., 200., 10.],
    nchains = 4
)
```
"""
function fit(; rinit=nothing, nchains::Int=2, datatype::String="rna", dttype=String[], datapath="HCT116_testdata/", gene="MYC", cell="HCT116", datacond="MOCK", traceinfo=(1.0, 1, -1, 1.0), infolder::String="HCT116_test", resultfolder::String="HCT116_test", inlabel::String="", label::String="", fittedparam=Int[], fixedeffects=tuple(), transitions=([1, 2], [2, 1]), G=2, R=0, S=0, insertstep=1, coupling=tuple(), TransitionType="nstate", grid=nothing, root=".", elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, maxtime=60, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, TransitionType)
    if isnothing(rinit)
        fit(nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber)
    else
        fit(rinit, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber)
    end
end

"""
    fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene::String, cell::String, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fixedeffects::String, G::String, R::String, S::String, insertstep::String, TransitionType="", grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)

"""
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene::String, cell::String, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fixedeffects::String, G::String, R::String, S::String, insertstep::String, TransitionType="", grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)
    transitions = get_transitions(G, TransitionType)
    fixedeffects, fittedparam = make_fixedfitted(datatype, fixedeffects, transitions, parse(Int, R), parse(Int, S), parse(Int, insertstep), length(noisepriors), coupling, grid)
    println("transitions: ", transitions)
    println("fixedeffects: ", fixedeffects)
    println("fittedparam: ", fittedparam)
    fit(nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects, transitions, parse(Int, G), parse(Int, R), parse(Int, S), parse(Int, insertstep), tuple(), root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype,
        propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber)
end

"""
    fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)


"""
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)
    S = reset_S(S, R, insertstep)
    nalleles = alleles(gene, cell, root, nalleles=nalleles)
    propcv = get_propcv(propcv, folder_path(infolder, root, "results"), inlabel, gene, G, R, S, insertstep, nalleles)
    fit(readrates(folder_path(infolder, root, "results"), inlabel, gene, G, R, S, insertstep, nalleles, ratetype), nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber)
end

"""
    fit(rinit, nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)

"""
function fit(rinit, nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)
    println(now())
    printinfo(gene, G, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime, nalleles, propcv)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    data, model, options = make_structures(rinit, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, method, zeromedian, datacol, ejectnumber)
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
    make_structures(rinit, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)

Create and configure data, model, and options structures for fitting.

# Arguments
- `rinit`: Initial rate parameters
- `datatype::String`: Type of data being analyzed
- `dttype::Vector`: Dwell time types
- `datapath`: Path to data files
- `gene`: Gene name
- `cell`: Cell type
- `datacond`: Data condition
- `traceinfo`: Trace information tuple
- `infolder::String`: Input folder path
- `label::String`: Label for the analysis
- `fittedparam`: Parameters to fit
- `fixedeffects`: Fixed effects structure
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `coupling::Tuple`: Coupling structure (default: empty tuple)
- `grid`: Grid parameter (default: nothing)
- `root`: Root directory (default: ".")
- `maxtime`: Maximum runtime in minutes (default: 60)
- `elongationtime`: RNA elongation time (default: 6.0)
- `priormean`: Prior means for parameters (default: empty array)
- `priorcv`: Prior coefficient of variation (default: 10.0)
- `nalleles`: Number of alleles (default: 1)
- `onstates`: ON states for the model (default: empty array)
- `decayrate`: Decay rate (default: -1.0)
- `splicetype`: Splicing type (default: "")
- `probfn`: Probability function (default: prob_Gaussian)
- `noisepriors`: Noise priors (default: empty array)
- `hierarchical`: Hierarchical structure (default: empty tuple)
- `ratetype`: Rate type (default: "median")
- `propcv`: Proposal coefficient of variation (default: 0.01)
- `samplesteps`: Number of sampling steps (default: 1000000)
- `warmupsteps`: Number of warmup steps (default: 0)
- `annealsteps`: Number of annealing steps (default: 0)
- `temp`: Temperature (default: 1.0)
- `tempanneal`: Annealing temperature (default: 100.0)
- `temprna`: RNA temperature (default: 1.0)
- `method`: Numerical method (default: Tsit5())
- `zeromedian`: Whether to zero median (default: false)
- `datacol`: Data column (default: 3)
- `ejectnumber`: Ejection number (default: 1)

# Returns
- `Tuple`: (data, model, options) structures

# Notes
- Validates and adjusts model parameters
- Loads and processes data according to datatype
- Sets up priors and initial conditions
- Creates appropriate model structure based on parameters
- Handles hierarchical, coupled, and grid models
"""
function make_structures(rinit, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)
    gene = check_genename(gene, "[")
    S = reset_S(S, R, insertstep)
    if G == 1 && !isempty(transitions)
        println("G=1, transitions are ignored")
        transitions = tuple()
    end
    nalleles = reset_nalleles(nalleles, coupling)
    infolder = folder_path(infolder, root, "results")
    datapath = folder_path(datapath, root, "data")
    data = load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna, datacol, zeromedian)
    decayrate = set_decayrate(decayrate, gene, cell, root)
    priormean = set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), coupling, grid)
    fittedparam = set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noisepriors, coupling, grid)
    model = load_model(data, rinit, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid, zeromedian, ejectnumber)
    if samplesteps > 0
        options = MHOptions(samplesteps, warmupsteps, annealsteps, Float64(maxtime), temp, tempanneal)
    else
        throw(ArgumentError("samplesteps must be greater than 0"))
    end
    return data, model, options
end


"""
    set_trace_background(traceinfo)

Extract background information from traceinfo tuple.

# Arguments
- `traceinfo`: Trace information tuple containing background data

# Returns
- Background information as vector or scalar

# Notes
- Extracts background from 5th element of traceinfo if present
- Handles both scalar and vector background values
- Throws ArgumentError if background information is missing
- Used for trace data preprocessing and normalization
"""
function set_trace_background(traceinfo)
    if length(traceinfo) > 4
        if eltype(traceinfo[5]) <: AbstractVector
            background = Vector[]
            for t in traceinfo[5]
                push!(background, t)
            end
        else
            background = traceinfo[5]
        end
    else
        throw(ArgumentError("Must include trace background"))
        # background = Vector[]
    end
    return background
end

"""
    set_trace_weight(traceinfo)

Calculate trace weights from active fraction information.

# Arguments
- `traceinfo`: Trace information tuple containing active fraction data

# Returns
- Weight values as vector or scalar

# Notes
- Extracts active fraction from 4th element of traceinfo
- Calculates weight as (1 - fraction) / fraction
- Handles both scalar and vector active fractions
- Used for trace data weighting in likelihood calculations
"""
function set_trace_weight(traceinfo)
    if traceinfo[4] isa Vector
        weight = Float64[]
        for f in traceinfo[4]
            push!(weight, max((1 - f) / f, 0.0))
        end
    else
        weight = max((1 - traceinfo[4]) / traceinfo[4], 0.0)
    end
    return weight
end

"""
    zero_median(tracer::Vector{T}, zeromedian::Bool) where {T<:AbstractVector}

Zero-center and scale trace data by subtracting median values.

# Arguments
- `tracer::Vector{T}`: Vector of trace vectors
- `zeromedian::Bool`: Whether to zero-center the traces

# Returns
- `Tuple`: (processed traces, scale factor)

# Notes
- Calculates median and MAD for each trace
- If zeromedian=true, subtracts median and scales by maximum median
- If zeromedian=false, returns traces unchanged with scale=1
- Used for trace normalization in fluorescence data analysis
"""
function zero_median(tracer::Vector{T}, zeromedian::Bool) where {T<:AbstractVector}
    medians = [median(t) for t in tracer]
    mads = [mad(t, normalize=false) for t in tracer]
    scale = max(1., maximum(medians))
    madscale = maximum(mads)
    if zeromedian
        trace = similar(tracer)
        for i in eachindex(tracer)
            trace[i] = (tracer[i] .- medians[i]) ./ scale
        end
    else
        trace = deepcopy(tracer)
        scale = 1.
    end
    return trace, madscale ./ scale
end

"""
    zero_median(tracer::Vector{T}, zeromedian::Bool) where {T<:AbstractMatrix}

Zero-center and scale matrix trace data by subtracting median values.

# Arguments
- `tracer::Vector{T}`: Vector of trace matrices
- `zeromedian::Bool`: Whether to zero-center the traces

# Returns
- `Tuple`: (processed traces, scale factors)

# Notes
- Calculates median and MAD for each column of each matrix
- If zeromedian=true, subtracts median and scales by maximum median per column
- If zeromedian=false, returns traces unchanged with scale=1
- Used for multi-channel trace normalization
"""
function zero_median(tracer::Vector{T}, zeromedian::Bool) where {T<:AbstractMatrix}
    medians = [median(t, dims=1) for t in tracer]
    # Calculate MAD for each column of each matrix
    mads = [reshape([mad(t[:, i], normalize=false) for i in 1:size(t, 2)], 1, :) for t in tracer]
    scale = max.(1., maximum(vcat(medians...), dims=1))
    madscale = median(vcat(mads...), dims=1)
    if zeromedian
        trace = similar(tracer)
        for i in eachindex(tracer)
            trace[i] = (tracer[i] .- medians[i]) ./ scale
        end
    else
        trace = deepcopy(tracer)
        scale = 1.
    end
    return trace, madscale ./ scale
end

"""
    zero_median(tracer::Vector{T}, zeromedian::Vector{Bool}) where {T<:AbstractMatrix}

Zero-center and scale matrix trace data with per-column control.

# Arguments
- `tracer::Vector{T}`: Vector of trace matrices
- `zeromedian::Vector{Bool}`: Boolean vector controlling zero-centering per column

# Returns
- `Tuple`: (processed traces, scale factors)

# Notes
- Applies zero-centering only to columns where zeromedian[i] = true
- Calculates median and MAD for each column of each matrix
- Scales by maximum median per column
- Used for selective trace normalization in multi-channel data
"""
function zero_median(tracer::Vector{T}, zeromedian::Vector{Bool}) where {T<:AbstractMatrix}
    medians = [median(t, dims=1) for t in tracer]
    # Calculate MAD for each column of each matrix
    mads = [reshape([mad(t[:, i], normalize=false) for i in 1:size(t, 2)], 1, :) for t in tracer]
    scale = max.(1., maximum(vcat(medians...), dims=1))
    madscale = median(vcat(mads...), dims=1)
    trace = deepcopy(tracer)
    println("zeromedian: ", zeromedian)
    for z in eachindex(zeromedian)
        if zeromedian[z]
            for i in eachindex(tracer)
                trace[i][:, z] = (tracer[i][:, z] .- medians[i][z]) ./ scale[z]
            end
        else
            scale[z] = 1.
        end
    end
    return trace, madscale ./ scale
end

# function robust_zscore(tracer::Vector, zeromedian=false)
#     if zeromedian
#         trace = similar(tracer)
#         medians = [median(t) for t in tracer]
#         mads = [mad(t) for t in tracer]
#         for i in eachindex(tracer)
#             trace[i] = (tracer[i] .- medians[i]) ./ mad(mads)
#         end
#     else
#         trace = tracer
#     end
#     return trace
# end

"""
    load_data_trace(datapath, label, gene, datacond, traceinfo, datatype::Symbol, col=3, zeromedian=false)

Load and process trace data from files.

# Arguments
- `datapath`: Path to trace data files
- `label`: Label for the dataset
- `gene`: Gene name
- `datacond`: Data condition
- `traceinfo`: Trace information tuple
- `datatype::Symbol`: Type of trace data (:trace, :tracejoint, :tracerna)
- `col`: Column to read from trace files (default: 3)
- `zeromedian`: Whether to zero-center traces (default: false)

# Returns
- Trace data structure (TraceData or TraceRNAData)

# Notes
- Reads trace files using read_tracefiles
- Applies zero-median normalization if requested
- Handles background and weight calculations
- Supports single and joint trace data types
- Throws error if no traces are found
"""
function load_data_trace(datapath, label, gene, datacond, traceinfo, datatype::Symbol, col=3, zeromedian=false)
    if typeof(datapath) <: String
        tracer = read_tracefiles(datapath, datacond, traceinfo, col)
    else
        tracer = read_tracefiles(datapath[1], datacond, traceinfo, col)
    end
    (length(tracer) == 0) && throw("No traces")
    trace, tracescale = zero_median(tracer, zeromedian)
    println("number of traces: ", length(trace))
    println("datapath: ", datapath)
    println("datacond: ", datacond)
    println("traceinfo: ", traceinfo)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo) > 3 && traceinfo[4] < 1.0
        weight = set_trace_weight(traceinfo)
        background = set_trace_background(traceinfo)
    else
        weight = 0.0
        background = 0.0
    end
    if datatype == :trace || datatype == :tracejoint
        return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo[1], (trace, background, weight, nframes, tracescale))
    elseif datatype == :tracerna
        len, h = read_rna(gene, datacond, datapath[2])
        return TraceRNAData(label, gene, traceinfo[1], (trace, background, weight, nframes), len, h)
    else
        throw(ArgumentError("Unsupported datatype '$datatype'"))
    end
end


"""
    load_data_grid(datapath, label, gene, datacond, traceinfo)

Load grid data for analysis.

# Arguments
- `datapath`: Path to grid data files
- `label`: Label for the dataset
- `gene`: Gene name
- `datacond`: Data condition
- `traceinfo`: Trace information tuple

# Returns
- GridData structure for analysis
"""
function load_data_grid(datapath, label, gene, datacond, traceinfo)
    # Implementation for grid data loading
    error("load_data_grid not yet implemented")
end

"""
    load_data_tracegrid(datapath, label, gene, datacond, traceinfo)

Load trace data for grid-based analysis.

# Arguments
- `datapath`: Path to trace data files
- `label`: Label for the dataset
- `gene`: Gene name
- `datacond`: Data condition
- `traceinfo`: Trace information tuple

# Returns
- TraceData structure for grid analysis

# Notes
- Reads trace files using read_tracefiles_grid
- Creates TraceData structure with grid-specific parameters
- Used for grid-based model fitting approaches
"""
function load_data_tracegrid(datapath, label, gene, datacond, traceinfo)
    trace = read_tracefiles_grid(datapath, datacond, traceinfo)
    # nframes = traceinfo[3] < 0 ? floor(Int, (720 - traceinfo[2] + traceinfo[1]) / traceinfo[1]) : floor(Int, (traceinfo[3] - traceinfo[2] + traceinfo[1]) / traceinfo[1])
    return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo[1], (trace, Vector[], 0.0, 1))
end

const TRACE_DATATYPES = Set([
    :trace, :tracegrid, :tracejoint, :tracerna
])

"""
    load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna, datacol=3, zeromedian=false)

Load RNA or trace data based on the provided `datatype` string or symbol.

Supports multiple data formats, including steady-state RNA histograms,
dwell time distributions, ON/OFF state durations, and fluorescence traces.

# Arguments
- `datatype`: String or Symbol describing the data type (e.g. "rna", "tracegrid")
- `dttype`: Dwell time type (used only for rnadwelltime and dwelltime)
- `datapath`: Path(s) to the data file(s)
- `label`: Label for the dataset
- `gene`: Gene name
- `datacond`: Experimental condition
- `traceinfo`: Tuple of trace metadata
- `temprna`: Integer divisor for histogram normalization
- `datacol`: Column of trace data to extract (default = 3)
- `zeromedian`: If true, zero-center each trace before fitting (default = false)

# Returns
- A concrete data structure subtype (e.g., `RNAData`, `TraceRNAData`, `DwellTimeData`)

# Throws
- `ArgumentError` if `datatype` is unsupported
"""
function load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna, datacol=3, zeromedian=false)
    dt = normalize_datatype(datatype)

    if dt == :rna
        len, h = read_rna(gene, datacond, datapath)
        return RNAData(label, gene, len, h)

    elseif dt == :rnacount
        countsRNA, yieldfactor, nRNA = read_rnacount(gene, datacond, datapath)
        return RNACountData(label, gene, nRNA, countsRNA, yieldfactor)

    elseif dt == :rnaonoff
        len, h = read_rna(gene, datacond, datapath[1])
        h = div.(h, temprna)
        LC = readfile(gene, datacond, datapath[2])
        return RNAOnOffData(label, gene, len, h, LC[:, 1], LC[:, 2], LC[:, 3])

    elseif dt == :rnadwelltime
        len, h = read_rna(gene, datacond, datapath[1])
        h = div.(h, temprna)
        bins, DT = read_dwelltimes(datapath[2:end])
        return RNADwellTimeData(label, gene, len, h, bins, DT, dttype)

    elseif dt == :dwelltime
        bins, DT = read_dwelltimes(datapath)
        return DwellTimeData(label, gene, bins, DT, dttype)

    elseif dt ∈ TRACE_DATATYPES
        if dt == :tracegrid
            return load_data_tracegrid(datapath, label, gene, datacond, traceinfo)
        else
            return load_data_trace(datapath, label, gene, datacond, traceinfo, dt, datacol, zeromedian)
        end

    else
        throw(ArgumentError("Unsupported datatype '$datatype'"))
    end
end

"""
    make_reporter_components(transitions::Tuple, G::Int, R, S, insertstep, onstates, dttype, splicetype)

Create reporter components for dwell time analysis.

# Arguments
- `transitions::Tuple`: Model transitions
- `G::Int`: Number of gene states
- `R`: Number of RNA steps
- `S`: Number of splice states
- `insertstep`: Insertion step
- `onstates`: ON states for the model
- `dttype`: Dwell time types
- `splicetype`: Splicing type

# Returns
- `Tuple`: (sojourn states, nonzero rows), components

# Notes
- Creates TDComponents for dwell time analysis
- Determines sojourn states based on onstates and model structure
- Identifies nonzero rows in the component matrix
- Used for dwell time distribution modeling
"""
function make_reporter_components(transitions::Tuple, G::Int, R, S, insertstep, onstates, dttype, splicetype)
    make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, onstates, dttype)
end

"""
    make_reporter_components_DT(transitions, G::Int, R::Int, S::Int, insertstep, splicetype, onstates, dttype, coupling=tuple())

Create reporter components for dwell time analysis.

# Arguments
- `transitions`: Model transitions
- `G::Int`: Number of gene states
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep`: Insertion step
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `dttype`: Dwell time types
- `coupling`: Coupling structure (default: empty tuple)

# Returns
- `Tuple`: (sojourn states, nonzero rows), components

# Notes
- Creates TDComponents for dwell time analysis
- Determines sojourn states based on onstates and model structure
- Identifies nonzero rows in the component matrix
- Used for dwell time distribution modeling
"""
function make_reporter_components_DT(transitions, G::Int, R::Int, S::Int, insertstep, splicetype, onstates, dttype, coupling=tuple())
    sojourn = sojourn_states(onstates, G, R, S, insertstep, dttype)
    components = TDComponents(transitions, G, R, S, insertstep, sojourn, dttype, splicetype)
    nonzeros = nonzero_rows(components)
    return (sojourn, nonzeros), components
end

"""
    make_reporter_components_DT(transitions, G::Tuple, R::Tuple, S::Tuple, insertstep, splicetype, onstates, dttype, coupling)

Create reporter components for coupled dwell time analysis.

# Arguments
- `transitions`: Model transitions for each unit
- `G::Tuple`: Number of gene states for each unit
- `R::Tuple`: Number of RNA steps for each unit
- `S::Tuple`: Number of splice states for each unit
- `insertstep`: Insertion steps for each unit
- `splicetype`: Splicing type
- `onstates`: ON states for each unit
- `dttype`: Dwell time types
- `coupling`: Coupling structure

# Returns
- `Tuple`: (coupled sojourn states, coupled nonzero rows), components

# Notes
- Creates TDCoupledComponents for multi-unit dwell time analysis
- Applies coupling to sojourn states and nonzero rows
- Handles interactions between different model units
- Used for coupled dwell time distribution modeling
"""
function make_reporter_components_DT(transitions, G::Tuple, R::Tuple, S::Tuple, insertstep, splicetype, onstates, dttype, coupling)
    sojourn = sojourn_states(onstates, G, R, S, insertstep, dttype)
    components = TDCoupledComponents(coupling, transitions, G, R, S, insertstep, sojourn, dttype, splicetype)
    sojourn = coupled_states(sojourn, coupling, components, G)
    nonzeros = coupled_states(nonzero_rows(components), coupling, components, G)
    return (sojourn, nonzeros), components
end

"""
    make_reporter_components_original(transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, splicetype, onstates, probfn, noisepriors, coupling=tuple())

Create original reporter components for trace analysis.

# Arguments
- `transitions::Tuple`: Model transitions
- `G::Int`: Number of gene states
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep::Int`: Insertion step
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure (default: empty tuple)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Creates HMMReporter for hidden Markov model analysis
- Handles mixture models if probfn contains "Mixture"
- Creates TComponents for trace analysis
- Used for fluorescence trace data modeling
"""
function make_reporter_components_original(transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, splicetype, onstates, probfn, noisepriors, coupling=tuple())
    nnoise = length(noisepriors)
    n = num_rates(transitions, R, S, insertstep)
    weightind = occursin("Mixture", "$(probfn)") ? n + nnoise : 0
    reporter = HMMReporter(nnoise, num_reporters_per_state(G, R, S, insertstep, onstates), probfn, weightind, off_states(G, R, S, insertstep, onstates), collect(n+1:n+nnoise))
    components = TComponents(transitions, G, R, S, insertstep, splicetype)
    return reporter, components
end

"""
    make_reporter_components(transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, splicetype, onstates, probfn, noisepriors, coupling=tuple())

Create reporter components for trace analysis with coupling support.

# Arguments
- `transitions::Tuple`: Model transitions
- `G::Int`: Number of gene states
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep::Int`: Insertion step
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure (default: empty tuple)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Creates HMMReporter for hidden Markov model analysis
- Handles mixture models if probfn contains "Mixture"
- Creates TComponents or TForcedComponents based on coupling
- Used for fluorescence trace data modeling with optional coupling
"""
function make_reporter_components(transitions::Tuple, G::Int, R::Int, S::Int, insertstep::Int, splicetype, onstates, probfn, noisepriors, coupling=tuple())
    nnoise = length(noisepriors)
    n = num_rates(transitions, R, S, insertstep)
    weightind = occursin("Mixture", "$(probfn)") ? n + nnoise : 0
    reporter = HMMReporter(nnoise, num_reporters_per_state(G, R, S, insertstep, onstates), probfn, weightind, off_states(G, R, S, insertstep, onstates), collect(n+1:n+nnoise))
    if isempty(coupling)
        components = TComponents(transitions, G, R, S, insertstep, splicetype)
    else
        components = TForcedComponents(coupling, transitions, G, R, S, insertstep, splicetype)
    end
    return reporter, components
end

"""
    make_reporter_components(transitions::Tuple, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype, onstates, probfn, noisepriors, coupling)

Create reporter components for coupled trace analysis.

# Arguments
- `transitions::Tuple`: Model transitions for each unit
- `G::Tuple`: Number of gene states for each unit
- `R::Tuple`: Number of RNA steps for each unit
- `S::Tuple`: Number of splice states for each unit
- `insertstep::Tuple`: Insertion steps for each unit
- `splicetype`: Splicing type
- `onstates`: ON states for each unit
- `probfn`: Probability function(s)
- `noisepriors`: Noise priors for each unit
- `coupling`: Coupling structure

# Returns
- `Tuple`: (reporters, components)

# Notes
- Creates HMMReporter for each unit in coupled model
- Handles different probability functions per unit
- Creates TCoupledComponents for multi-unit trace analysis
- Used for coupled fluorescence trace data modeling
"""
function make_reporter_components(transitions::Tuple, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype, onstates, probfn, noisepriors, coupling)
    reporter = HMMReporter[]
    !(probfn isa Union{Tuple,Vector}) && (probfn = fill(probfn, length(coupling[1])))
    n_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1], onstates)
    for i in eachindex(G)
        nnoise = length(noisepriors[i])
        n = num_rates(transitions[i], R[i], S[i], insertstep[i])
        weightind = occursin("Mixture", "$(probfn)") ? n + nnoise : 0
        push!(reporter, HMMReporter(nnoise, n_per_state[i], probfn[i], weightind, off_states(n_per_state[i]), collect(n+1:n+nnoise)))
    end
    components = TCoupledComponents(coupling, transitions, G, R, S, insertstep, splicetype)
    return reporter, components
end


"""
    make_reporter_components(data::AbstractRNAData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)

Create reporter components for RNA data analysis.

# Arguments
- `data::AbstractRNAData`: RNA data structure
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Creates MComponents for RNA abundance modeling
- Uses onstates as reporter structure
- Handles RNA decay and splicing processes
- Used for steady-state RNA histogram analysis
"""
function make_reporter_components(data::AbstractRNAData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)
    reporter = onstates
    components = MComponents(transitions, G, R, data.nRNA, decayrate, splicetype, ejectnumber)
    return reporter, components
end

"""
    make_reporter_components(data::RNAOnOffData, transitions, G::Int, R::Int, S::Int, insertstep::Int, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)

Create reporter components for RNA ON/OFF data analysis.

# Arguments
- `data::RNAOnOffData`: RNA ON/OFF data structure
- `transitions`: Model transitions
- `G::Int`: Number of gene states
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep::Int`: Insertion step
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Tuple`: (onstates, components)

# Notes
- Creates MTAIComponents for RNA ON/OFF analysis
- Automatically determines ON states if not provided
- Handles RNA abundance and ON/OFF state transitions
- Used for combined RNA histogram and ON/OFF state analysis
"""
function make_reporter_components(data::RNAOnOffData, transitions, G::Int, R::Int, S::Int, insertstep::Int, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)
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
    return onstates, MTAIComponents(transitions, G, R, S, insertstep, onstates, data.nRNA, decayrate, splicetype, ejectnumber)
end

"""
    make_reporter_components(data::DwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)

Create reporter components for dwell time data analysis.

# Arguments
- `data::DwellTimeData`: Dwell time data structure
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Delegates to make_reporter_components_DT for dwell time analysis
- Uses dwell time types from data structure
- Creates appropriate components for dwell time distribution modeling
- Used for analyzing time spent in different model states
"""
function make_reporter_components(data::DwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)
    make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, onstates, data.DTtypes, coupling)
end

"""
    make_reporter_components(data::RNADwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)

Create reporter components for combined RNA and dwell time data analysis.

# Arguments
- `data::RNADwellTimeData`: Combined RNA and dwell time data structure
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Creates MTComponents combining MComponents and TDComponents
- Handles both RNA abundance and dwell time distributions
- Combines steady-state and kinetic information
- Used for comprehensive transcription analysis
"""
function make_reporter_components(data::RNADwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)
    reporter, tcomponents = make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, onstates, data.DTtypes, coupling)
    mcomponents = MComponents(transitions, G, R, data.nRNA, decayrate, splicetype, ejectnumber)
    # components = MTDComponents(MComponents(transitions, G, R, data.nRNA, decayrate, splicetype, ejectnumber), tcomponents)
    return reporter, MTComponents{typeof(mcomponents),typeof(tcomponents)}(mcomponents, tcomponents)
end

"""
    make_reporter_components(data::AbstractTraceData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)

Create reporter components for trace data analysis.

# Arguments
- `data::AbstractTraceData`: Trace data structure
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Delegates to make_reporter_components for trace analysis
- Creates appropriate components for fluorescence trace modeling
- Handles HMM-based analysis of time series data
- Used for live-cell imaging data analysis
"""
function make_reporter_components(data::AbstractTraceData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)
    make_reporter_components(transitions, G, R, S, insertstep, splicetype, onstates, probfn, noisepriors, coupling)
end

"""
    make_reporter_components(data::AbstractTraceHistogramData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)

Create reporter components for combined trace and histogram data analysis.

# Arguments
- `data::AbstractTraceHistogramData`: Combined trace and histogram data structure
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Tuple`: (reporter, components)

# Notes
- Creates MTComponents combining trace and RNA abundance analysis
- Handles both time series and steady-state data
- Combines fluorescence traces with RNA histograms
- Used for comprehensive transcription analysis with multiple data types
"""
function make_reporter_components(data::AbstractTraceHistogramData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1)
    reporter, tcomponents = make_reporter_components(transitions, G, R, S, insertstep, splicetype, onstates, probfn, noisepriors, coupling)
    mcomponents = MComponents(transitions, G, R, data.nRNA, decayrate, splicetype, ejectnumber)
    return reporter, MTComponents{typeof(mcomponents),typeof(tcomponents)}(mcomponents, tcomponents)
end

"""
    make_grid(transitions, R::Int, S, insertstep, noisepriors, grid)

Create parameter ranges for grid-based analysis.

# Arguments
- `transitions`: Model transitions
- `R::Int`: Number of RNA steps
- `S`: Number of splice states
- `insertstep`: Insertion step
- `noisepriors`: Noise priors
- `grid`: Grid parameter

# Returns
- `Tuple`: (raterange, noiserange, gridrange)

# Notes
- Calculates parameter indices for different parameter types
- raterange: indices for rate parameters
- noiserange: indices for noise parameters
- gridrange: indices for grid parameters
- Used for grid-based parameter search and optimization
"""
function make_grid(transitions, R::Int, S, insertstep, noisepriors, grid)
    n = num_rates(transitions, R, S, insertstep)
    raterange = 1:n
    noiserange = n+1:n+length(noisepriors)
    gridrange = n+length(noisepriors)+1:n+length(noisepriors)+Int(!isnothing(grid))
    raterange, noiserange, gridrange
end

"""
    grid_indices(transitions, R, S, insertstep, reporter, coupling, grid)

Calculate grid parameter indices.

# Arguments
- `transitions`: Model transitions
- `R, S, insertstep`: Model structure parameters
- `reporter`: Reporter structure
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector`: Grid parameter indices

# Notes
- Calculates indices for grid parameters in the full parameter vector
- Used for grid-based model fitting
- Returns single-element array with grid parameter index
"""
function grid_indices(transitions, R, S, insertstep, reporter, coupling, grid)
    [num_all_parameters(transitions, R, S, insertstep, reporter, coupling, grid)]
end

"""
    coupling_indices(transitions, R, S, insertstep, reporter, coupling, grid)

Calculate coupling parameter indices.

# Arguments
- `transitions`: Model transitions
- `R, S, insertstep`: Model structure parameters
- `reporter`: Reporter structure
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector`: Coupling parameter indices

# Notes
- Calculates indices for coupling parameters in the full parameter vector
- Accounts for grid parameters if present
- Returns range of indices for coupling parameters
- Used for coupled model parameter identification
"""
function coupling_indices(transitions, R, S, insertstep, reporter, coupling, grid)
    n = num_all_parameters(transitions, R, S, insertstep, reporter, coupling, grid)
    g = isnothing(grid) ? 0 : 1
    collect(n-g-coupling[5]+1:n-g)
end



"""
    make_fitted_hierarchical(fittedshared, nhypersets, fittedindividual, nallparams, nindividuals)

Create fitted parameter vectors for hierarchical models.

# Arguments
- `fittedshared`: Indices of shared parameters to fit
- `nhypersets`: Number of hyperparameter sets
- `fittedindividual`: Indices of individual parameters to fit
- `nallparams`: Total number of parameters per individual
- `nindividuals`: Number of individuals

# Returns
- `Tuple`: (f, fhyper, fpriors)

# Notes
- Creates comprehensive parameter indexing for hierarchical models
- f: all fitted parameters (shared + hyper + individual)
- fhyper: fitted hyper parameters for each hyperparameter set
- fpriors: fitted parameters that have priors (shared + hyper)
- Used for hierarchical model parameter organization
"""
function make_fitted_hierarchical(fittedshared, nhypersets, fittedindividual, nallparams, nindividuals)
    f = union(fittedshared, fittedindividual) # shared parameters come first followed by hyper parameters and individual parameters
    f = sort(f)
    fhyper = [fittedindividual] # fitted hyper parameters correspond to fitted individual parameters
    for i in 1:nhypersets-1
        append!(f, fittedindividual .+ i * nallparams)
        push!(fhyper, fittedindividual .+ i * nallparams)
    end
    fpriors = sort(union(fittedshared, fhyper...)) # priors only apply to shared and hyper parameters
    for i in 1:nindividuals
        append!(f, fittedindividual .+ (i + nhypersets - 1) * nallparams)
    end
    # fhyper =[findall(x -> x in sub_target, f) for sub_target in fhyper]
    f, fhyper, fpriors
end


"""
    make_hierarchical(data, rmean, fittedparam, fixedeffects, transitions, R, S, insertstep, priorcv, noisepriors, hierarchical::Tuple, reporter, coupling=tuple(), couplingindices=nothing, grid=nothing, factor=10, ratetransforms=nothing, zeromedian=false)

Construct hierarchical model traits and priors.

# Arguments
- `data`: Data structure
- `rmean`: Prior means
- `fittedparam`: Fitted parameter indices
- `fixedeffects`: Fixed effects structure
- `transitions, R, S, insertstep`: Model structure parameters
- `priorcv`: Prior coefficient of variation
- `noisepriors`: Noise priors
- `hierarchical::Tuple`: Hierarchical structure tuple
- `reporter`: Reporter structure
- `coupling`: Coupling structure (default: empty tuple)
- `couplingindices`: Coupling parameter indices (default: nothing)
- `grid`: Grid parameter (default: nothing)
- `factor`: Prior scaling factor (default: 10)
- `ratetransforms`: Rate transformations (default: nothing)
- `zeromedian`: Whether to zero median (default: false)

# Returns
- `Tuple`: (hierarchical trait, fittedparam, fixedeffects, priord, ratetransforms)

# Notes
- Creates HierarchicalTrait for multi-level models
- Organizes parameters into shared, hyper, and individual components
- Sets up appropriate prior distributions
- Handles parameter indexing for hierarchical structure
"""
function make_hierarchical(data, rmean, fittedparam, fixedeffects, transitions, R, S, insertstep, priorcv, noisepriors, hierarchical::Tuple, reporter, coupling=tuple(), couplingindices=nothing, grid=nothing, factor=10, ratetransforms=nothing, zeromedian=false)
    fittedindividual = hierarchical[2]
    fittedshared = setdiff(fittedparam, fittedindividual)
    nhypersets = hierarchical[1]
    n_all_params = num_all_parameters(transitions, R, S, insertstep, reporter, coupling, grid)
    nindividualparams = length(fittedindividual) # number of fitted params per individual
    nindividuals = length(data.trace[1])
    ratestart = nhypersets * n_all_params + 1
    paramstart = length(fittedshared) + nhypersets * nindividualparams + 1
    fittedparam, fittedhyper, fittedpriors = make_fitted_hierarchical(fittedshared, hierarchical[1], hierarchical[2], n_all_params, nindividuals)
    hierarchy = HierarchicalTrait(nhypersets, n_all_params, nindividualparams, nindividuals, ratestart, paramstart, fittedhyper, fittedshared, fittedpriors)
    fixedeffects = make_fixed(fixedeffects, hierarchical[3], n_all_params, nindividuals)
    # rprior = rmean[1:nhypersets*n_all_params]
    # priord = prior_distribution(rprior, transitions, R, S, insertstep, fittedpriors, priorcv, noisepriors, couplingindices, factor, ratetransforms)
    priord = prior_distribution(rmean, priorcv, ratetransforms, fittedpriors)
    return hierarchy, fittedparam, fixedeffects, priord, ratetransforms
end


"""
    rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates::Int, reporter, zeromedian)

Populate transformation functions for rates and noise parameters.

# Arguments
- `ftransforms`: Array to populate with forward transformations
- `invtransforms`: Array to populate with inverse transformations
- `sigmatransforms`: Array to populate with sigma transformations
- `nrates::Int`: Number of rate parameters
- `reporter`: Reporter structure
- `zeromedian`: Whether to zero-center traces

# Returns
- Nothing (modifies arrays in place)

# Notes
- Adds log/exp transformations for rate parameters
- Handles noise parameters for HMMReporter
- Applies identity transformation for zero-median noise parameters
- Used for parameter transformation in MCMC sampling
"""
function rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates::Int, reporter, zeromedian)
    for i in 1:nrates
        push!(ftransforms, log)
        push!(invtransforms, exp)
        push!(sigmatransforms, sigmalognormal)
    end
    if typeof(reporter) <: HMMReporter
        for i in eachindex(reporter.noiseparams)
            # Handle both Bool and Vector{Bool} for zeromedian
            should_zero = zeromedian isa Bool ? zeromedian :
                          (zeromedian isa Vector{Bool} && length(zeromedian) > 0 ? zeromedian[2] : false)
            if isodd(i) && should_zero
                push!(ftransforms, identity)
                push!(invtransforms, identity)
                push!(sigmatransforms, sigmanormal)
            else
                push!(ftransforms, log)
                push!(invtransforms, exp)
                push!(sigmatransforms, sigmalognormal)
            end
        end
    end
end

"""
    rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates::Vector, reporter, zeromedian)

Populate transformation functions for multi-unit models.

# Arguments
- `ftransforms`: Array to populate with forward transformations
- `invtransforms`: Array to populate with inverse transformations
- `sigmatransforms`: Array to populate with sigma transformations
- `nrates::Vector`: Vector of rate counts for each unit
- `reporter`: Reporter structure for each unit
- `zeromedian`: Whether to zero-center traces

# Returns
- Nothing (modifies arrays in place)

# Notes
- Processes each unit in coupled models separately
- Delegates to single-unit rate_transforms! for each unit
- Used for multi-unit model parameter transformations
"""
function rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates::Vector, reporter, zeromedian)
    for n in eachindex(nrates)
        rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates[n], reporter[n], zeromedian)
    end
end


"""
    make_ratetransforms(data, nrates, transitions, G, R, S, insertstep, reporter, coupling, grid, hierarchical, zeromedian)

Create transformation functions for all model parameters.

# Arguments
- `data`: Data structure
- `nrates`: Number of rate parameters
- `transitions, G, R, S, insertstep`: Model structure parameters
- `reporter`: Reporter structure
- `coupling`: Coupling structure
- `grid`: Grid parameter
- `hierarchical`: Hierarchical structure
- `zeromedian`: Whether to zero-center traces

# Returns
- `Transformation`: Struct with forward/inverse/sigma transforms

# Notes
- Creates comprehensive transformation functions for all parameters
- Handles rate parameters, noise parameters, coupling parameters, and grid parameters
- Extends transformations for hierarchical models
- Used for parameter transformation in MCMC sampling
"""
function make_ratetransforms(data, nrates, transitions, G, R, S, insertstep, reporter, coupling, grid, hierarchical, zeromedian)
    ftransforms = Function[]
    invtransforms = Function[]
    sigmatransforms = Function[]

    if isa(G, Int) && !isempty(coupling)
        # for forced model, assume that forced data is zero median
        rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates, reporter, true)
    else
        rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates, reporter, zeromedian)
    end

    # rate_transforms!(ftransforms, invtransforms, sigmatransforms, nrates, reporter, zeromedian)

    if !isempty(coupling)
        couplingindices = coupling_indices(transitions, R, S, insertstep, reporter, coupling, grid)
        for i in eachindex(couplingindices)
            push!(ftransforms, log_shift1)
            push!(invtransforms, invlog_shift1)
            push!(sigmatransforms, sigmalognormal)
        end
    end
    if !isnothing(grid)
        gridindices = grid_indices(transitions, R, S, insertstep, noisepriors, coupling, grid)
        for i in eachindex(gridindices)
            push!(ftransforms, log)
            push!(invtransforms, exp)
            push!(sigmatransforms, sigmalognormal)
        end
    end
    if !isempty(hierarchical)
        fset = ftransforms
        iset = invtransforms
        nindividuals = length(data.trace[1])
        for i in 1:hierarchical[1]-1
            ftransforms = vcat(ftransforms, repeat([log], length(fset)))
            invtransforms = vcat(invtransforms, repeat([exp], length(iset)))
            sigmatransforms = vcat(sigmatransforms, repeat([sigmalognormal], length(fset)))
        end
        ftransforms = vcat(ftransforms, repeat(fset, nindividuals))
        invtransforms = vcat(invtransforms, repeat(iset, nindividuals))
        sigmatransforms = vcat(sigmatransforms, repeat(sigmatransforms, nindividuals))
    end
    Transformation(ftransforms, invtransforms, sigmatransforms)
end


"""
    load_model(data, r, rmean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid, zeromedian=false, ejectnumber=1, factor=10)

Construct and return the appropriate model struct for the given data and options.

# Arguments
- `data`: Data structure
- `r`: Initial rate parameters
- `rmean`: Prior means
- `fittedparam`: Fitted parameter indices
- `fixedeffects`: Fixed effects structure
- `transitions, G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type
- `nalleles`: Number of alleles
- `priorcv`: Prior coefficient of variation
- `onstates`: ON states for the model
- `decayrate`: Decay rate
- `propcv`: Proposal coefficient of variation
- `probfn`: Probability function
- `noisepriors`: Noise priors
- `method`: Numerical method
- `hierarchical`: Hierarchical structure
- `coupling`: Coupling structure
- `grid`: Grid parameter
- `zeromedian`: Whether to zero-center traces (default: false)
- `ejectnumber`: Number of mRNAs per burst (default: 1)
- `factor`: Prior scaling factor (default: 10)

# Returns
- Model struct (e.g., `GMmodel`, `GRSMmodel`)

# Notes
- Creates appropriate model based on parameter combinations
- Handles hierarchical, coupled, and grid models
- Sets up traits and prior distributions
- Returns concrete model type based on model complexity
"""
function load_model(data, r, rmean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid, zeromedian=false, ejectnumber=1, factor=10)
    reporter, components = make_reporter_components(data, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber)

    nrates = num_rates(transitions, R, S, insertstep)
    ratetransforms = make_ratetransforms(data, nrates, transitions, G, R, S, insertstep, reporter, coupling, grid, hierarchical, zeromedian)

    if !isempty(coupling)
        couplingindices = coupling_indices(transitions, R, S, insertstep, reporter, coupling, grid)
        ncoupling = coupling[5]
        couplingtrait = CouplingTrait(ncoupling, couplingindices)
    else
        couplingindices = nothing
    end
    if !isnothing(grid)
        gridindices = grid_indices(transitions, R, S, insertstep, noisepriors, coupling, grid)
        gridtrait = GridTrait(grid, gridindices)
    else
        gridindices = nothing
    end
    if !isempty(hierarchical)
        hierarchicaltrait, fittedparam, fixedeffects, priord = make_hierarchical(data, rmean, fittedparam, fixedeffects, transitions, R, S, insertstep, priorcv, noisepriors, hierarchical, reporter, coupling, couplingindices, grid, factor, ratetransforms, zeromedian)
    else
        # priord = prior_distribution(rmean, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors, couplingindices, factor, ratetransforms)
        priord = prior_distribution(rmean, priorcv, ratetransforms, fittedparam)
    end

    CBool = isempty(coupling)
    GBool = isnothing(grid)
    HBool = isempty(hierarchical)

    if CBool && GBool && HBool
        if R == 0
            if typeof(data) <: AbstractTraceData
                # For trace data with R=0, still use GRSMmodel but with simplified components
                return GRSMmodel{Nothing,typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(nothing, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
            else
                # For non-trace data with R=0, use GMmodel
                return GMmodel{typeof(r),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(r, transitions, G, nalleles, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
            end
        else
            return GRSMmodel{Nothing,typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(nothing, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        end
    elseif CBool && GBool && !HBool
        trait = (hierarchical=hierarchicaltrait,)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    elseif CBool && !GBool && HBool
        trait = (grid=gridtrait,)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    elseif CBool && !GBool && !HBool
        trait = (hierarchical=hierarchicaltrait, grid=gridtrait)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    elseif !CBool && GBool && HBool
        trait = (coupling=couplingtrait,)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    elseif !CBool && GBool && !HBool
        trait = (coupling=couplingtrait, hierarchical=hierarchicaltrait)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
    elseif !CBool && !GBool && HBool
        trait = (coupling=couplingtrait, grid=gridtrait)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        # return GRSMcoupledgridmodel()
    elseif !CBool && !GBool && !HBool
        trait = (coupling=couplingtrait, hierarchical=hierarchicaltrait, grid=gridtrait)
        return GRSMmodel{typeof(trait),typeof(r),typeof(nrates),typeof(G),typeof(priord),typeof(propcv),typeof(fittedparam),typeof(method),typeof(components),typeof(reporter)}(trait, r, ratetransforms, nrates, transitions, G, R, S, insertstep, nalleles, splicetype, priord, propcv, fittedparam, fixedeffects, method, components, reporter)
        # return GRSMcoupledgridhierarchicalmodel()
    end
end

"""
    checklength(r, transitions, R, S, insertstep, reporter)

Check if parameter vector has correct length for the model.

# Arguments
- `r`: Parameter vector
- `transitions`: Model transitions
- `R, S, insertstep`: Model structure parameters
- `reporter`: Reporter structure

# Returns
- Nothing

# Notes
- Validates parameter vector length against expected model parameters
- Includes rate parameters and noise parameters if HMMReporter
- Throws error if parameter vector has wrong length
- Used for parameter validation before model fitting
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

"""
    prior_hypercv(transitions, R::Int, S, insertstep, noisepriors)

Generate default hyperparameter coefficient of variation values.

# Arguments
- `transitions`: Model transitions
- `R::Int`: Number of RNA steps
- `S`: Number of splice states
- `insertstep`: Insertion step
- `noisepriors`: Noise priors

# Returns
- `Vector{Float64}`: Coefficient of variation values for hyperparameters

# Notes
- Creates default CV values for hierarchical model hyperparameters
- Uses 1.0 for most parameters, 0.1 for RNA shifts and noise parameters
- Provides reasonable starting values for hyperparameter priors
- Used for hierarchical model setup
"""
function prior_hypercv(transitions, R::Int, S, insertstep, noisepriors)
    [fill(1.0, length(transitions)); 1.0; fill(0.1, R - 1); 1.0; fill(1.0, max(0, S - insertstep + 1)); 1.0; fill(0.1, length(noisepriors))]
end

"""
    prior_hypercv(transitions, R::Tuple, S, insertstep, noisepriors, coupling)

Generate default hyperparameter CV values for coupled models.

# Arguments
- `transitions`: Model transitions for each unit
- `R::Tuple`: Number of RNA steps for each unit
- `S`: Number of splice states for each unit
- `insertstep`: Insertion steps for each unit
- `noisepriors`: Noise priors for each unit
- `coupling`: Coupling structure

# Returns
- `Vector{Float64}`: Coefficient of variation values for hyperparameters

# Notes
- Processes each unit in coupled model separately
- Appends coupling parameter CV values
- Used for hierarchical coupled model setup
"""
function prior_hypercv(transitions, R::Tuple, S, insertstep, noisepriors, coupling)
    rm = Float64[]
    for i in eachindex(R)
        append!(rm, prior_hypercv(transitions[i], R[i], S[i], insertstep[i], noisepriors[i]))
    end
    [rm; fill(1.0, coupling[5])]
end

"""
    prior_hypercv(transitions, R, S, insertstep, noisepriors, coupling, grid)

Generate default hyperparameter CV values with grid support.

# Arguments
- `transitions`: Model transitions
- `R, S, insertstep`: Model structure parameters
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector{Float64}`: Coefficient of variation values for hyperparameters

# Notes
- Delegates to appropriate prior_hypercv function based on coupling
- Appends grid parameter CV if grid is not nothing
- Used for hierarchical models with optional grid and coupling
"""
function prior_hypercv(transitions, R, S, insertstep, noisepriors, coupling, grid)
    if isempty(coupling)
        pcv = prior_hypercv(transitions, R, S, insertstep, noisepriors)
    else
        pcv = prior_hypercv(transitions, R, S, insertstep, noisepriors, coupling)
    end
    if !isnothing(grid)
        append!(pcv, 1.0)
    end
    pcv
end

# function prior_ratemean_hierarchical(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, nhypersets, coupling=tuple(), cv::Float64=1.0)
#     r = isempty(coupling) ? prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime) : prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling)
#     hypercv = prior_hypercv(transitions, R, S, insertstep, noisepriors, coupling)
#     # hypercv = [fill(1.0, length(transitions)); 1.0; fill(0.1, R - 1); 1.0; fill(1.0, max(0, S - insertstep + 1)); 1.0; fill(.1,length(noisepriors))]
#     append!(r, hypercv)
#     for i in 3:nhypersets
#         append!(r, fill(cv, length(rm)))
#     end
#     r
# end

"""
    prior_ratemean_hierarchical(priormean, hypercv, nhypersets, cv::Float64=1.0)

Create prior means for hierarchical models.

# Arguments
- `priormean`: Base prior means
- `hypercv`: Hyperparameter coefficient of variation values
- `nhypersets`: Number of hyperparameter sets
- `cv::Float64`: Coefficient of variation for additional hyperparameter sets (default: 1.0)

# Returns
- `Vector{Float64}`: Extended prior means for hierarchical model

# Notes
- Extends base prior means with hyperparameter values
- Adds additional hyperparameter sets with specified CV
- Used for hierarchical model prior setup
- Arranges parameters as: shared + hyper + individual
"""
function prior_ratemean_hierarchical(priormean, hypercv, nhypersets, cv::Float64=1.0)
    r = copy(priormean)
    append!(r, hypercv)
    for i in 3:nhypersets
        append!(r, fill(cv, length(priormean)))
    end
    r
end


# function prior_ratemean_grid(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime)
#     [prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime); 0.5]
#     # [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.1, max(0, S - insertstep + 1)); decayrate; noisepriors; 0.5]
# end

"""
    prior_ratemean_grid(priormean)

Add grid parameter to prior means.

# Arguments
- `priormean`: Base prior means

# Returns
- `Vector{Float64}`: Prior means with grid parameter

# Notes
- Appends grid probability parameter (0.5) to prior means
- Used for grid-based model fitting
- Grid parameter represents probability of grid state
"""
function prior_ratemean_grid(priormean)
    [priormean; 0.5]
end

"""
    prior_ratemean(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors::Vector, elongationtime::Float64, initprior::Float64=0.1)

Generate default prior means for rate parameters.

# Arguments
- `transitions`: Model transitions
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep`: Insertion step
- `decayrate`: Decay rate
- `noisepriors::Vector`: Noise priors
- `elongationtime::Float64`: RNA elongation time
- `initprior::Float64`: Initial prior value (default: 0.1)

# Returns
- `Vector{Float64}`: Prior means for all parameters

# Notes
- Creates default prior means for all model parameters
- Uses 0.01 for gene transitions, initprior for initiation
- Calculates RNA processing rates from elongation time
- Uses 0.1 for splicing, provided decay rate, and noise priors
- Used for model initialization when no prior means provided
"""
function prior_ratemean(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors::Vector, elongationtime::Float64, initprior::Float64=0.1)
    [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.1, max(0, S - insertstep + 1)); decayrate; noisepriors]
end


"""
    prior_ratemean(transitions, R::Tuple, S::Tuple, insertstep::Tuple, decayrate, noisepriors::Union{Vector,Tuple}, elongationtime::Union{Vector,Tuple}, coupling, initprior=[0.1, 0.1])

Generate default prior means for coupled models.

# Arguments
- `transitions`: Model transitions for each unit
- `R::Tuple`: Number of RNA steps for each unit
- `S::Tuple`: Number of splice states for each unit
- `insertstep::Tuple`: Insertion steps for each unit
- `decayrate`: Decay rates for each unit
- `noisepriors::Union{Vector,Tuple}`: Noise priors for each unit
- `elongationtime::Union{Vector,Tuple}`: RNA elongation times for each unit
- `coupling`: Coupling structure
- `initprior`: Initial prior values for each unit (default: [0.1, 0.1])

# Returns
- `Vector{Float64}`: Prior means for all parameters including coupling

# Notes
- Processes each unit in coupled model separately
- Appends coupling parameter prior means (0.0)
- Used for coupled model initialization
"""
function prior_ratemean(transitions, R::Tuple, S::Tuple, insertstep::Tuple, decayrate, noisepriors::Union{Vector,Tuple}, elongationtime::Union{Vector,Tuple}, coupling, initprior=[0.1, 0.1])
    rm = Float64[]
    for i in eachindex(R)
        append!(rm, prior_ratemean(transitions[i], R[i], S[i], insertstep[i], decayrate, noisepriors[i], elongationtime[i], initprior[i]))
    end
    [rm; fill(0.0, coupling[5])]
end

"""
    set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid)

Set prior means if empty, using default values.

# Arguments
- `priormean`: Prior means (if empty, will be generated)
- `transitions, R, S, insertstep`: Model structure parameters
- `decayrate`: Decay rate
- `noisepriors`: Noise priors
- `elongationtime`: RNA elongation time
- `hierarchical`: Hierarchical structure
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector{Float64}`: Prior means for all parameters

# Notes
- Generates default prior means if priormean is empty
- Handles coupled, grid, and hierarchical models
- Applies appropriate modifications based on model type
- Used for model initialization and setup
"""
function set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid)
    if !isempty(priormean)
        return priormean
    else
        if !isempty(coupling)
            priormean = prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling)
        else
            priormean = prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime)
        end
        if !isnothing(grid)
            priormean = prior_ratemean_grid(priormean)
        end
        if !isempty(hierarchical)
            priormean = prior_ratemean_hierarchical(priormean, prior_hypercv(transitions, R, S, insertstep, noisepriors, coupling, grid), hierarchical[1])
        end
    end
    priormean
end

"""
    prior_distribution(rprior, priorcv, ratetransforms, fittedparam)

Create prior distribution for fitted parameters.

# Arguments
- `rprior`: Prior means
- `priorcv`: Prior coefficient of variation
- `ratetransforms`: Rate transformations
- `fittedparam`: Indices of fitted parameters

# Returns
- Prior distribution array

# Notes
- Creates prior distributions for fitted parameters only
- Handles scalar and vector priorcv values
- Validates that priorcv length matches rprior length
- Delegates to prior_distribution_array for actual distribution creation
"""
function prior_distribution(rprior, priorcv, ratetransforms, fittedparam)
    if priorcv isa Number
        priorcv = fill(priorcv, length(rprior))
    end
    if length(priorcv) == length(rprior)
        return prior_distribution_array(rprior, priorcv, ratetransforms, fittedparam)
    else
        throw(ArgumentError("priorcv not the same length as prior mean"))
    end
end



# function prior_distribution_coupling(rm, transitions, R::Tuple, S::Tuple, insertstep::Tuple, fittedparam::Vector, priorcv, noisepriors, couplingindices, factor, ratetransforms)
#     if isempty(rm)
#         throw(ArgumentError("No prior mean"))
#         # rm = prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors)
#     end
#     if priorcv isa Number
#         rcv = fill(priorcv, length(rm))
#         s = 0
#         for i in eachindex(R)
#             lnp = isempty(noisepriors) ? 0 : length(noisepriors[i])
#             n = num_rates(transitions[i], R[i], S[i], insertstep[i])
#             rcv[s+n] = 0.1
#             rcv[couplingindices] ./= factor
#             s += n + lnp
#         end
#     else
#         rcv = priorcv
#     end
#     if length(rcv) == length(rm)
#         # return distribution_array(apply_transform(rm[fittedparam], ratetransforms.f[fittedparam]), prior_sigma(rm[fittedparam], rcv[fittedparam], ratetransforms.f_cv[fittedparam]), Normal)
#         return prior_distribution_array(rm, rcv, ratetransforms, fittedparam)
#     else
#         throw(ArgumentError("priorcv not the same length as prior mean"))
#     end
# end

# """
#     prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, priorcv, noisepriors, factor=10)
#     prior_distribution(rm, transitions, R::Tuple, S::Tuple, insertstep::Tuple, fittedparam::Vector, priorcv, noisepriors, factor=10)


# """
# function prior_distribution(rm, transitions, R::Int, S::Int, insertstep, fittedparam::Vector, priorcv, noisepriors, factor, ratetransforms)
#     if isempty(rm)
#         throw(ArgumentError("No prior mean"))
#         # rm = prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors)
#     end
#     if priorcv isa Number
#         n = num_rates(transitions, R, S, insertstep)
#         rcv = fill(priorcv, length(rm))
#         rcv[n] = 0.1
#         rcv[n+1:n+length(noisepriors)] /= factor
#     else
#         rcv = priorcv
#     end
#     if length(rcv) == length(rm)
#         # return distribution_array(apply_transform(rm[fittedparam], ratetransforms.f[fittedparam]), prior_sigma(rm[fittedparam], rcv[fittedparam], ratetransforms.f_cv[fittedparam]), Normal)
#         return prior_distribution_array(rm, rcv, ratetransforms, fittedparam)
#     else
#         throw(ArgumentError("priorcv not the same length as prior mean"))
#     end
# end

# # prior_distribution(rm, transitions, R::Tuple, S::Tuple, insertstep::Tuple, fittedparam::Vector, priorcv, noisepriors, couplingindices, factor=10) = prior_distribution_coupling(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors, couplingindices, factor)

# function prior_distribution(rm, transitions, R, S, insertstep, fittedparam::Vector, priorcv, noisepriors, couplingindices, factor, ratetransforms)
#     if isnothing(couplingindices)
#         prior_distribution(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors, factor, ratetransforms)
#     else
#         prior_distribution_coupling(rm, transitions, R, S, insertstep, fittedparam, priorcv, noisepriors, couplingindices, factor, ratetransforms)
#     end
# end

"""
    prior_sigma(r, cv, sigmatransforms)

Calculate prior standard deviations using transformation functions.

# Arguments
- `r`: Parameter values
- `cv`: Coefficient of variation values
- `sigmatransforms`: Sigma transformation functions

# Returns
- `Vector{Float64}`: Prior standard deviations

# Notes
- Applies sigma transformation functions to calculate standard deviations
- Handles both normal and log-normal parameter transformations
- Used for prior distribution setup in MCMC sampling
"""
function prior_sigma(r, cv, sigmatransforms)
    sigma = Vector{Float64}(undef, length(sigmatransforms))
    for i in eachindex(sigmatransforms)
        f = sigmatransforms[i]
        if f === sigmanormal
            sigma[i] = f(r[i], cv[i])
        else
            sigma[i] = f(cv[i])
        end
    end
    return sigma
end

"""
    prior_distribution_array(rm, rcv, ratetransforms, fittedparam)

Create prior distribution array for fitted parameters.

# Arguments
- `rm`: Prior means
- `rcv`: Prior coefficient of variation values
- `ratetransforms`: Rate transformations
- `fittedparam`: Indices of fitted parameters

# Returns
- Prior distribution array

# Notes
- Applies transformations to prior means and calculates standard deviations
- Creates distribution array for fitted parameters only
- Delegates to prior_distribution_array with transformed parameters
"""
function prior_distribution_array(rm, rcv, ratetransforms, fittedparam)
    means = apply_transform(rm[fittedparam], ratetransforms.f[fittedparam])
    sigmas = prior_sigma(rm[fittedparam], rcv[fittedparam], ratetransforms.f_cv[fittedparam])
    transforms = ratetransforms.f[fittedparam]
    prior_distribution_array(means, sigmas, transforms)
end

"""
    prior_distribution_array(position::Vector, scale::Vector, transforms::Vector{Function}, k=10)

Create array of prior distributions for parameters.

# Arguments
- `position::Vector`: Transformed parameter means
- `scale::Vector`: Parameter standard deviations
- `transforms::Vector{Function}`: Transformation functions
- `k`: Truncation factor for log-normal distributions (default: 10)

# Returns
- `Vector`: Array of prior distributions

# Notes
- Creates appropriate distribution type based on transformation function
- Uses truncated normal for log and log_shift1 transformations
- Uses normal distribution for identity transformations
- Used for MCMC prior setup
"""
function prior_distribution_array(position::Vector, scale::Vector, transforms::Vector{Function}, k=10)
    d = []
    for i in eachindex(transforms)
        if transforms[i] == log || transforms[i] == log_shift1
            push!(d, truncated_normal(position[i], scale[i], k))
        else
            push!(d, Normal(position[i], scale[i]))
        end
    end
    return d
end

# function set_rinit(infolder, inlabel,  gene, priormean, transitions,G, R, S, insertstep, nalleles, ratetype, hierarchical)

"""
    set_rinit(r, priormean, minval=1e-10, maxval=1e10)

Set initial rate parameters to prior if empty or invalid.

# Arguments
- `r`: Initial rate parameters
- `priormean`: Prior means
- `minval`: Minimum valid value (default: 1e-10)
- `maxval`: Maximum valid value (default: 1e10)

# Returns
- `Vector{Float64}`: Valid initial rate parameters

# Notes
- Uses prior means if r is empty
- Checks for NaN and Inf values
- Prints warning messages for invalid parameters
- Used for model initialization
"""
function set_rinit(r, priormean, minval=1e-10, maxval=1e10)
    if isempty(r)
        println("No rate file, set rate to prior")
        r = priormean
    elseif any(isnan.(r)) || any(isinf.(r))
        println("r out of bounds: ", r)
        r = priormean
    end
    println("initial: ", r)
    r
end

"""
    set_rinit(r, priormean, transitions, R, S, insertstep, noisepriors, nindividuals, coupling=tuple(), grid=nothing, minval=1e-7, maxval=300.0)

Set initial parameters for hierarchical models.

# Arguments
- `r`: Initial parameters
- `priormean`: Prior means
- `transitions, R, S, insertstep`: Model structure parameters
- `noisepriors`: Noise priors
- `nindividuals`: Number of individuals
- `coupling`: Coupling structure (default: empty tuple)
- `grid`: Grid parameter (default: nothing)
- `minval`: Minimum valid value (default: 1e-7)
- `maxval`: Maximum valid value (default: 300.0)

# Returns
- `Vector{Float64}`: Valid initial parameters for hierarchical model

# Notes
- Extends parameters for hierarchical models with multiple individuals
- Calculates total parameters including coupling and grid parameters
- Repeats prior means for each individual
- Used for hierarchical model initialization
"""
function set_rinit(r, priormean, transitions, R, S, insertstep, noisepriors, nindividuals, coupling=tuple(), grid=nothing, minval=1e-7, maxval=300.0)
    if isempty(r) || any(isnan.(r)) || any(isinf.(r))
        isempty(r) && println("No rate file, set rate to prior")
        any(isnan.(r)) && println("r contains NaN, set rate to prior")
        any(isinf.(r)) && println("r contains Inf, set rate to prior")
        r = copy(priormean)
        c = isempty(coupling) ? 0 : coupling[5]
        g = isnothing(grid) ? 0 : 1
        # n_all_params = num_rates(transitions, R, S, insertstep) + length(noisepriors)
        n_all_params = num_all_parameters(transitions, R, S, insertstep, noisepriors) + c + g
        for i in 1:nindividuals
            append!(r, priormean[1:n_all_params])
        end
    end
    r
end

"""
    num_noiseparams(datatype, noisepriors)

Calculate the number of noise parameters for a given data type.

# Arguments
- `datatype`: Type of data being analyzed
- `noisepriors`: Noise prior specifications

# Returns
- `Int`: Number of noise parameters

# Notes
- Returns 0 for non-trace data types
- Returns length of noisepriors for trace data types
- Used for parameter counting in model setup
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
    default_fitted(datatype::String, transitions, R::Tuple, S::Tuple, insertstep::Tuple, noiseparams::Tuple, coupling, grid)

Generate default fitted parameter indices for coupled models.

# Arguments
- `datatype::String`: Type of data being analyzed
- `transitions`: Model transitions for each unit
- `R::Tuple`: Number of RNA steps for each unit
- `S::Tuple`: Number of splice states for each unit
- `insertstep::Tuple`: Insertion steps for each unit
- `noiseparams::Tuple`: Number of noise parameters for each unit
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector{Int}`: Default fitted parameter indices

# Notes
- Processes each unit in coupled model separately
- Includes coupling parameters if coupling is not empty
- Includes grid parameter if grid is not nothing
- Used for coupled model parameter setup
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
    default_fitted(datatype::String, transitions, R::Int, S::Int, insertstep::Int, noiseparams, coupling, grid)

Generate default fitted parameter indices for single-unit models.

# Arguments
- `datatype::String`: Type of data being analyzed
- `transitions`: Model transitions
- `R::Int`: Number of RNA steps
- `S::Int`: Number of splice states
- `insertstep::Int`: Insertion step
- `noiseparams`: Number of noise parameters
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector{Int}`: Default fitted parameter indices

# Notes
- Creates indices for all rate parameters except decay rate
- Includes noise parameters for trace data types
- Includes coupling parameters if coupling is not empty
- Includes grid parameter if grid is not nothing
- Used for single-unit model parameter setup
"""
function default_fitted(datatype::String, transitions, R::Int, S::Int, insertstep::Int, noiseparams, coupling, grid)
    n = num_rates(transitions, R, S, insertstep)
    fittedparam = collect(1:n-1)
    if occursin("trace", datatype)
        if isnothing(grid)
            isempty(noiseparams) && throw(ArgumentError("noisepriors cannot be empty for trace data"))
            fittedparam = vcat(fittedparam, collect(n+1:n+noiseparams))
        else
            fittedparam = vcat(fittedparam, collect(n+1:n+noiseparams+1))
        end
    end
    fittedparam
end

"""
    set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noisepriors, coupling, grid)

Set fitted parameter indices if empty, using default values.

# Arguments
- `fittedparam`: Fitted parameter indices (if empty, will be generated)
- `datatype`: Type of data being analyzed
- `transitions, R, S, insertstep`: Model structure parameters
- `noisepriors`: Noise priors
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Vector{Int}`: Fitted parameter indices

# Notes
- Generates default fitted parameter indices if fittedparam is empty
- Handles both single-unit and coupled models
- Delegates to appropriate default_fitted function
- Used for model parameter setup
"""
function set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noisepriors, coupling, grid)
    if isempty(fittedparam)
        return default_fitted(datatype, transitions, R, S, insertstep, num_noiseparams(datatype, noisepriors), coupling, grid)
    else
        return fittedparam
    end
end


"""
    make_fixed(fixedshared, fixedindividual, nallparams, nindividuals)

Create fixed parameter indices for hierarchical models.

# Arguments
- `fixedshared`: Indices of fixed shared parameters
- `fixedindividual`: Indices of fixed individual parameters
- `nallparams`: Total number of parameters per individual
- `nindividuals`: Number of individuals

# Returns
- `Tuple`: (f, fhyper, fpriors)

# Notes
- Creates comprehensive fixed parameter indexing for hierarchical models
- f: all fixed parameters (shared + hyper + individual)
- fhyper: fixed hyper parameters for each hyperparameter set
- fpriors: fixed parameters that have priors (shared + hyper)
- Used for hierarchical model parameter organization
"""
function make_fixed(fixedshared, fixedindividual, nallparams, nindividuals)
    fixed = Vector{Int}[]
    for f in fixedshared
        push!(fixed, f)
    end
    for h in fixedindividual
        push!(fixed, [h + i * nallparams for i in 0:nindividuals-1])
    end
    tuple(fixed...)
end

"""
    make_fixedfitted(datatype, fixedeffects::String, transitions, R, S, insertstep, noiseparams, coupling, grid)

Create fixed and fitted parameter indices based on fixedeffects string.

# Arguments
- `datatype`: Type of data being analyzed
- `fixedeffects::String`: String specifying which parameters to fix
- `transitions, R, S, insertstep`: Model structure parameters
- `noiseparams`: Number of noise parameters
- `coupling`: Coupling structure
- `grid`: Grid parameter

# Returns
- `Tuple`: (fixedeffects, fittedparam)

# Notes
- Parses fixedeffects string to determine which parameters to fix
- Creates complementary fitted parameter indices
- Handles both single-unit and coupled models
- Used for model parameter setup with fixed effects
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
lossfn(x, data, model) = -loglikelihood(x, data, model)[1]

"""
    burstsize(fits, model::AbstractGMmodel)

Calculate burst size from fitted parameters for gene-only models.

# Arguments
- `fits`: Fitted parameters and statistics
- `model::AbstractGMmodel`: Gene-only model

# Returns
- `Float64`: Burst size

# Notes
- Calculates burst size as ratio of initiation rate to gene inactivation rate
- Uses maximum likelihood parameters from fits
- Used for gene-only model analysis
"""
function burstsize(fits, model::AbstractGMmodel)
    if model.G > 1
        b = Float64[]
        for p in eachcol(fits.param)
            r = get_rates(p, model)
            push!(b, r[2*model.G-1] / r[2*model.G-2])
        end
        return BurstMeasures(mean(b), std(b), median(b), mad(b, normalize=false), quantile(b, [0.025; 0.5; 0.975]))
    else
        return 0
    end
end
"""
    burstsize(fits::Fit, model::AbstractGRSMmodel)

Calculate burst size from fitted parameters for GRS models.

# Arguments
- `fits::Fit`: Fitted parameters and statistics
- `model::AbstractGRSMmodel`: GRS model

# Returns
- `Float64`: Burst size

# Notes
- Calculates burst size using model-specific burstsize function
- Uses maximum likelihood parameters from fits
- Delegates to burstsize(r, transitions, G, R, S, insertstep, splicetype, ejectnumber)
- Used for GRS model analysis
"""
function burstsize(fits::Fit, model::AbstractGRSMmodel)
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

burstsize(r, model::AbstractGRSMmodel, ejectnumber=1) = burstsize(r, model.Gtransitions, model.G, model.R, model.S, model.insertstep, model.splicetype, ejectnumber)

"""
    burstsize(r, transitions, G, R, S, insertstep, splicetype="", ejectnumber=1)

Calculate burst size from rate parameters for GRS models.

# Arguments
- `r`: Rate parameters
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type (default: "")
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Float64`: Burst size

# Notes
- Calculates burst size as ratio of initiation rate to gene inactivation rate
- Handles different gene state configurations
- Accounts for splicing and mRNA ejection processes
- Used for GRS model burst size analysis
"""
function burstsize(r, transitions, G, R, S, insertstep, splicetype="", ejectnumber=1)
    ntransitions = num_rates(transitions, R, S, insertstep)
    total = min(Int(div(r[ntransitions+1], r[ntransitions])) * 2, 400)
    components = MComponents(transitions, G, R, total, r[num_rates(transitions, R, S, insertstep)], splicetype, ejectnumber)
    M = make_mat_M(components, r)
    nT = G * 2^R
    L = nT * total
    S0 = zeros(L)
    S0[2] = 1.0
    s = time_evolve_diff([1, 10 / minimum(r[ntransitions:ntransitions+R+1])], M, S0)
    mean_histogram(s[2, collect(1:nT:L)])
end

"""
    check_genename(gene::String, p1)

Check if gene name matches expected pattern.

# Arguments
- `gene::String`: Gene name to check
- `p1`: Expected pattern or reference

# Returns
- `String`: Validated gene name

# Notes
- Validates gene name against expected pattern
- Returns gene name if valid, otherwise returns default
- Used for gene name validation in data processing
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

Check if gene names in vector match expected pattern.

# Arguments
- `gene::Vector`: Vector of gene names to check
- `p1`: Expected pattern or reference

# Returns
- `Vector{String}`: Validated gene names

# Notes
- Validates each gene name in vector against expected pattern
- Returns vector of validated gene names
- Used for gene name validation in coupled model data processing
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
    println("number of rates: ", length(model.rates))
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

function printinfo(gene, G::Int, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime, nalleles, propcv)
    if R == 0
        println("Gene: ", gene, " G: ", G, " Treatment:  ", datacond, ", nalleles: ", nalleles)
    else
        println("Gene: ", gene, ", datacond: ", datacond, ", GRSinsertstep: ", G, R, S, insertstep, ", nalleles: ", nalleles)
    end
    println("data: ", datapath)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime, ", propcv: ", propcv)
end

function printinfo(gene, G::Tuple, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime, nalleles, propcv)
    println("Gene: ", gene)
    for i in eachindex(G)
        println("datacond: ", datacond[i], ", GRSinsertstep: ", G[i], R[i], S[i], insertstep[i], ", nalleles: ", nalleles)
    end
    println("data: ", datapath)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime, ", propcv: ", propcv)
end



"""
    finalize(data,model,fits,stats,measures,temp,resultfolder,optimized,burst,writesamples,root)

write out run results and print out final loglikelihood and deviance
"""
function finalize(data, model, fits, stats, measures, temp, writefolder, optimized, burst, writesamples)
    println("final max ll: ", fits.llml)
    print_ll(transform_rates(vec(stats.medparam), model), data, model, "median ll: ")
    println("Median fitted rates: ", stats.medparam[:, 1])
    println("ML rates: ", inverse_transform_params(fits.parml, model))
    println("Acceptance: ", fits.accept, "/", fits.total, " (", fits.accept / fits.total * 100, "% )    ")
    if is_histogram_compatible(data)
        println("Deviance: ", deviance(fits, data, model))
    end
    println("rhat: ", maximum(measures.rhat))
    println("ess: ", minimum(measures.ess))
    println("geweke: ", maximum(measures.geweke))
    println("mcse: ", maximum(measures.mcse))
    println("waic: ", measures.waic)
    println("aic: ", aic(fits))
    # println("aic_onstates: ", aic_onstates(fits.parml, data, model))
    if optimized != 0
        println("Optimized ML: ", Optim.minimum(optimized))
        println("Optimized rates: ", exp.(Optim.minimizer(optimized)))
    end
    writeall(writefolder, fits, stats, measures, data, temp, model, optimized=optimized, burst=burst, writesamples=writesamples)
end

"""
    get_propcv(propcv, infolder, label, gene, G, R, S, insertstep, nalleles)

Get proposal coefficient of variation from file or use default.

# Arguments
- `propcv`: Proposal coefficient of variation (if negative, will be read from file)
- `infolder`: Input folder path
- `label`: Label for the dataset
- `gene`: Gene name
- `G, R, S, insertstep`: Model structure parameters
- `nalleles`: Number of alleles

# Returns
- Proposal coefficient of variation (scalar or matrix)

# Notes
- If propcv < 0, reads covariance matrix from param-stats file
- Scales covariance by 2.38^2 / n for optimal MCMC acceptance rate
- Returns absolute value of propcv if file doesn't exist or matrix is not positive definite
- Used for MCMC proposal distribution setup
"""
function get_propcv(propcv, infolder, label, gene, G, R, S, insertstep, nalleles)
    if propcv < 0.0
        file = get_resultfile("param-stats", infolder, label, gene, G, R, S, insertstep, nalleles)
        if !isfile(file)
            println(file, " does not exist")
            return abs(propcv)
        end
        cv = read_bottom_float_block(file)
        cv = 2.38^2 * cv / size(cv, 1)
        if isposdef(cv)
            return (cv, abs(propcv))
        else
            return abs(propcv)
        end
    else
        return propcv
    end
end

"""
    getcv(gene, G, nalleles, fittedparam, inlabel, infolder, root, verbose=true)

Obsolete, needs to be fixed
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
        return nalleles
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




