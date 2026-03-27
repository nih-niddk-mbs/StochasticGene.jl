# This file is part of StochasticGene.jl   

# fit.jl
#
# Fit GRS models (generalized telegraph models) to RNA abundance and live cell imaging data
#

"""
    fit(; <keyword arguments> )

Fit steady state or transient GM/GRSM model to RNA data for a single gene, write the result (through function finalize), and return fit results and diagnostics.

For coupled transcribing units, arguments transitions, G, R, S, insertstep, and trace become tuples of the single unit type, e.g. If two types of transcription models are desired with G=2 and G=3 then G = (2,3).

#Arguments
- `annealsteps=0`: number of annealing steps (during annealing temperature is dropped from tempanneal to temp)
- `burst=false`: if true then compute burst frequency
- `cell::String=""`: cell type for halflives and allele numbers
- `coupling=tuple()`: if nonempty, a 2- or 3-tuple: `(unit_model, connections::Vector{ConnectionSpec}[, sign_modes])`. Each connection is `(β, s, α, t) = (source unit, source state, target unit, target transition)`. Use `make_coupling("31", G, R)` or `make_coupling_reciprocal("3131", G, R)` in io.jl to build from a coupling field string. For **three units** with a **hidden latent** third unit (G=3, fully connected G-only dynamics) modulating observed units 1–2 from hidden states **1** and **3**, use `make_coupling_hidden_latent(t1, t2)` or `make_coupling_hidden_latent("H3#t1-t2")` and `transitions_hidden_g3_all_pairs()` for unit 3; see docs *Units and models — hidden latent unit*. Empty `connections` is valid (uncoupled T is still built). Optional third element `sign_modes` constrains the sign of each coupling parameter γ: use a single `Symbol` for all connections, or a vector/tuple of one per connection. Canonical symbols: `:free` (γ ∈ (-1, ∞)), `:activate` (γ ∈ (0, ∞)), `:inhibit` (γ ∈ (-1, 0)). Aliases `:positive` and `:negative` are normalized to `:activate` and `:inhibit`. See `coupling_ranges`.
- `datacol=3`: column of data to use, default is 3 for rna data
- `datatype::String=""`: String that describes data type, choices are "rna", "rnaonoff", "rnadwelltime", "trace", "tracerna", "tracejoint", "tracegrid"
- `datacond=""`: string or vector of strings describing data, e.g. "WT", "DMSO" or ["DMSO","AUXIN"], ["gene","enhancer"]
- `datapath=""`: path to data file or folder or array of files or folders
- `decayrate=1.0`: decay rate of mRNA, if set to -1, value in halflives folder will be used if it exists
- `ejectnumber=1`: number of mRNAs produced per burst, default is 1, if Int then deterministic, if Tuple = (r, p) then stochastic obeying NegativeBinomial(r, p)
- `dttype=String[]`: dwelltime types, choices are "OFF", "ON", for R states and "OFFG", "ONG" for G states
- `elongationtime=6.0`: average time for elongation, vector of times for coupled model
- `fittedparam::Vector=Int[]`: vector of rate indices to be fit, e.g. [1,2,3,5,6,7], if empty then all rates (except decay) are fit (applies to shared rates for hierarchical models, fitted hyper parameters are specified by individual fittedparams)
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
- `insertstep=1`: R step where reporter is inserted. Must be >= 1; when R = 0, insertstep is ignored (no RNA steps).
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
- `resultfolder::String=test`: folder for results of MCMC run. Resolved with `root` as: if `joinpath(root, resultfolder)` exists that path is used, else `joinpath(root, "results", resultfolder)` is used (and created if missing). So results go under `root` or `root/results/`
- `R=0`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `root="."`: name of root directory for project, e.g. "scRNA"
- `samplesteps::Int=1000000`: number of MCMC sampling steps
- `S=0`: number of splice sites (set to 0 for classic telegraph models and R - insertstep + 1 for GRS models)
- `splicetype=""`: RNA pathway for GRS models, (e.g., "offeject" = spliced intron is not viable)
- `temp=1.0`: MCMC temperature
- `tempanneal=100.`: annealing temperature
- `temprna=1.`: reduce RNA counts by temprna compared to dwell times
- `trace_specs=[]`: container of trace specs; each spec is a NamedTuple with at least `unit`, `interval`, `start`, `t_end`, `zeromedian` (and optionally `active_fraction`, `background`). **Coupled `tracejoint`:** when `trace_specs` is empty, `make_structures` fills defaults via `default_trace_specs_for_coupled` so `data.units` lists observed units (required for correct HMM emission masking with hidden units / Rany).
- `dwell_specs=[]`: container of dwell-time specs per unit (e.g. onstates, bins, dttype). When non-empty, used for dwell-time data with multiple units or observation mapping. Legacy dwell data uses single-unit defaults when empty.
- `traceinfo=(1.0, 1., -1, 1., 0.5)`: 5 tuple = (frame interval of intensity traces in minutes, starting frame time in minutes, ending frame time (use -1 for last index), fraction of observed active traces, background mean)
    for simultaneous joint traces, the fraction of active traces is a vector of the active fractions for each trace, e.g. (1.0, 1., -1, [.5, .7], [0.5,0.5]) 
    If active fraction is 1.0, then traceinfo can be a 3-tuple, e.g. (1.0, 1., -1) since background correction is not needed
    Note that all traces are scaled by the maximum of the medians of all the traces, the traces are all scaled by the same factor since the signal amplitude should be the same
- `TransitionType=""`: String describing G transition type, e.g. "3state", "KP" (kinetic proofreading), "cyclic", or if hierarchical, coupled
- `transitions::Tuple=([1,2],[2,1])`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2-state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3-state kinetic proofreading model, empty for G=1
- `warmupsteps=0`: number of MCMC warmup steps to find proposal distribution covariance
- `writesamples=false`: write out MH samples if true, default is false
- `zeromedian=false`: if true, subtract the median of each trace from each trace, then scale by the maximum of the medians
- `key=nothing`: when nothing, fit uses the keyword arguments you pass (and defaults). When a string (e.g. `key=\"33il\"`), fit looks for `info_<key>.toml` in the results folder; if found, loads that spec and overrides with any kwargs you pass (kwargs take precedence). If not found, uses your kwargs and defaults. Results are always written to `info_<stem>.toml`; with a key, that file is also read on the next run when present.

# Returns
- `fits`: MCMC fit results (posterior samples, log-likelihoods, etc.)
- `stats`: Summary statistics for parameters
- `measures`: Diagnostic measures (including WAIC and its standard error, which is now for the total WAIC and scaled by sqrt(n_obs))
- `data`, `model`, `options`: The data, model, and options structures used

# Notes
- If `propcv < 0`, proposal covariance is read from previous run if available.
- WAIC standard error is for the total WAIC (not per observation), and is scaled by sqrt(n_obs).
- File and folder conventions: see the package manual (*Package overview*, *Cluster and batch workflows*) and the [GitHub README](https://github.com/nih-niddk-mbs/StochasticGene.jl#readme).

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
# Set by top-level fit(; key=..., ...) so finalize/writeall get key-based naming and info TOML without kwargs in the call chain.
const _current_run_spec = Ref{Any}(nothing)
const _current_name_override = Ref{Any}(nothing)

const _FIT_DEFAULTS = (
    rinit=nothing,
    nchains=2,
    datatype="rna",
    dttype=String[],
    datapath="HCT116_testdata/",
    gene="MYC",
    cell="HCT116",
    datacond="MOCK",
    traceinfo=(1.0, 1, -1, 1.0),
    infolder="HCT116_test",
    resultfolder="HCT116_test",
    inlabel="",
    label="",
    fittedparam=Int[],
    fixedeffects=tuple(),
    transitions=([1, 2], [2, 1]),
    G=2,
    R=0,
    S=0,
    insertstep=1,
    coupling=tuple(),
    TransitionType="nstate",
    grid=nothing,
    root=".",
    elongationtime=6.0,
    priormean=Float64[],
    priorcv=10.0,
    nalleles=1,
    onstates=Int[],
    decayrate=-1.0,
    splicetype="",
    probfn=nothing,  # resolved to prob_Gaussian at runtime (defined in hmm.jl, included after fit.jl)
    noisepriors=[],
    hierarchical=tuple(),
    ratetype="median",
    propcv=0.01,
    maxtime=60.0,
    samplesteps=1000000,
    warmupsteps=0,
    annealsteps=0,
    temp=1.0,
    tempanneal=100.0,
    temprna=1.0,
    burst=false,
    optimize=false,
    writesamples=false,
    method=Tsit5(),
    zeromedian=false,
    datacol=3,
    ejectnumber=1,
    yieldfactor=1.0,
    trace_specs=[],
    dwell_specs=[],
)

"""
    fit_default_spec() -> Dict{Symbol,Any}

Return a copy of default `fit` keyword arguments as a `Dict`, with `probfn` resolved to `prob_Gaussian`
when the default would be `nothing`. Used by batch utilities (e.g. `makeswarm_models`) to build
run specs consumed by `write_run_spec_preset` / `fit(; key=...)`.
"""
function fit_default_spec()
    d = Dict{Symbol, Any}(pairs(_FIT_DEFAULTS))
    if d[:probfn] === nothing
        d[:probfn] = prob_Gaussian
    end
    return d
end

"""
    fit_coupled_default_spec() -> Dict{Symbol,Any}

Like [`fit_default_spec`](@ref), but baseline keywords for **coupled** / trace-joint batch jobs
(`makeswarmfiles` with CSV keys, `base_keys`, or H3 grids): `datatype=\"tracejoint\"`, more chains,
fewer MCMC steps than the single-gene RNA default, etc. Per-model structure (`G`, `R`, `coupling`,
priors, …) should come from merging an existing `info_<key>.jld2` (see [`makeswarmfiles`](@ref)) or
explicit `kwargs`.
"""
function fit_coupled_default_spec()
    d = fit_default_spec()
    d[:datatype] = "tracejoint"
    d[:nchains] = 16
    d[:samplesteps] = 100_000
    d[:annealsteps] = 0
    d[:warmupsteps] = 0
    d[:writesamples] = false
    return d
end

function fit(; key=nothing, kwargs...)
    defaults = Dict{Symbol, Any}(pairs(_FIT_DEFAULTS))
    if defaults[:probfn] === nothing
        defaults[:probfn] = prob_Gaussian
    end
    kw = Dict{Symbol, Any}(pairs(kwargs))
    merged = copy(defaults)
    merge!(merged, kw)
    if key !== nothing
        resultfolder = merged[:resultfolder]
        root = merged[:root]
        spec_path = joinpath(folder_path(resultfolder, root, "results"), "info_" * key * ".toml")
        if isfile(spec_path)
            try
                spec = read_run_spec(spec_path)
                merged = merge(defaults, spec, kw)
            catch
                # Info file exists but failed to parse (e.g. malformed TOML); use kwargs only
            end
        end
    end
    run_spec = Dict{Symbol, Any}(merged)
    if key !== nothing
        run_spec[:key] = key
    end
    name_override = key !== nothing ? filename(key) : nothing
    rinit = merged[:rinit]
    nchains = merged[:nchains]
    datatype = merged[:datatype]
    dttype = merged[:dttype]
    datapath = merged[:datapath]
    gene = merged[:gene]
    cell = merged[:cell]
    datacond = merged[:datacond]
    traceinfo = merged[:traceinfo]
    infolder = merged[:infolder]
    resultfolder = merged[:resultfolder]
    inlabel = merged[:inlabel]
    label = merged[:label]
    fittedparam = merged[:fittedparam]
    fixedeffects = merged[:fixedeffects]
    transitions = merged[:transitions]
    G = merged[:G]
    R = merged[:R]
    S = merged[:S]
    insertstep = merged[:insertstep]
    coupling = merged[:coupling]
    grid = merged[:grid]
    root = merged[:root]
    elongationtime = merged[:elongationtime]
    priormean = merged[:priormean]
    priorcv = merged[:priorcv]
    nalleles = merged[:nalleles]
    onstates = merged[:onstates]
    decayrate = merged[:decayrate]
    splicetype = merged[:splicetype]
    probfn = merged[:probfn]
    noisepriors = merged[:noisepriors]
    hierarchical = merged[:hierarchical]
    ratetype = merged[:ratetype]
    propcv = merged[:propcv]
    maxtime = merged[:maxtime]
    samplesteps = merged[:samplesteps]
    warmupsteps = merged[:warmupsteps]
    annealsteps = merged[:annealsteps]
    temp = merged[:temp]
    tempanneal = merged[:tempanneal]
    temprna = merged[:temprna]
    burst = merged[:burst]
    optimize = merged[:optimize]
    writesamples = merged[:writesamples]
    method = merged[:method]
    zeromedian = merged[:zeromedian]
    datacol = merged[:datacol]
    ejectnumber = merged[:ejectnumber]
    yieldfactor = merged[:yieldfactor]
    trace_specs = merged[:trace_specs]
    dwell_specs = merged[:dwell_specs]
    label, inlabel = create_label(label, inlabel, datatype, datacond, cell, merged[:TransitionType])
    run_spec[:label] = label
    run_spec[:inlabel] = inlabel
    if rinit === nothing && key !== nothing && key != ""
        rr = folder_path(infolder, root, "results")
        rates_path = joinpath(rr, "rates_" * key * ".txt")
        if isfile(rates_path)
            rinit = readrates(rates_path, get_row(ratetype))
        end
    end
    _current_run_spec[] = run_spec
    _current_name_override[] = name_override
    try
        if rinit === nothing
            fit(nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, inlabel, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)
        else
            fit(rinit, nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)
        end
    finally
        _current_run_spec[] = nothing
        _current_name_override[] = nothing
    end
end

"""
    fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)


"""
function fit(nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, inlabel::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[])
    S = reset_S(S, R, insertstep)
    nalleles = alleles(gene, cell, root, nalleles=nalleles)
    propcv = get_propcv(propcv, folder_path(infolder, root, "results"), inlabel, gene, G, R, S, insertstep, nalleles)
    fit(readrates(folder_path(infolder, root, "results"), inlabel, gene, G, R, S, insertstep, nalleles, ratetype), nchains, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, resultfolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, burst, optimize, writesamples, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)
end

"""
    fit(rinit, nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1)

"""
function fit(rinit, nchains::Int, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, resultfolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, burst=false, optimize=false, writesamples=false, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[])
    println(now())
    printinfo(gene, G, R, S, insertstep, datacond, datapath, infolder, resultfolder, maxtime, nalleles, propcv)
    resultfolder = folder_path(resultfolder, root, "results", make=true)
    data, model, options = make_structures(rinit, datatype, dttype, datapath, gene, cell, datacond, traceinfo, infolder, label, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling, grid, root, maxtime, elongationtime, priormean, priorcv, nalleles, onstates, decayrate, splicetype, probfn, noisepriors, hierarchical, ratetype, propcv, samplesteps, warmupsteps, annealsteps, temp, tempanneal, temprna, method, zeromedian, datacol, ejectnumber, yieldfactor, trace_specs, dwell_specs)
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
    [((1, 2), [(1, Int(s), 2, Int(t))]) for s in source for t in target]
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
    normalize_insertstep(R, insertstep)

Enforce insertstep >= 1 everywhere. When R = 0, insertstep is ignored (no RNA steps) and 1 is used.
When R > 0, insertstep is clamped to at least 1; a warning is issued if it was < 1.
"""
function normalize_insertstep(R::Int, insertstep::Int)
    if R == 0
        return 1
    end
    if insertstep < 1
        @warn "insertstep was $insertstep ( < 1 ); set to 1"
        return 1
    end
    return insertstep
end
function normalize_insertstep(R::Tuple, insertstep::Tuple)
    out = collect(insertstep)
    for i in eachindex(R)
        if R[i] == 0
            out[i] = 1
        elseif insertstep[i] < 1
            @warn "insertstep[$i] was $(insertstep[i]) ( < 1 ); set to 1"
            out[i] = 1
        else
            out[i] = insertstep[i]
        end
    end
    return tuple(out...)
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
    ncoupling(coupling)

Return the number of coupling strength parameters.

Coupling format is `(unit_model, connections::Vector{ConnectionSpec})` or
`(unit_model, connections::Vector{ConnectionSpec}, ranges)`, where `ranges`
encodes per-connection sign modes (see `coupling_ranges`).
"""
ncoupling(coupling) = isempty(coupling) ? 0 : length(coupling[2])

"""
    coupling_ranges(coupling)

Return per-connection coupling sign modes.

- If `coupling` has only two elements `(unit_model, connections)`, all modes
  default to `:free` (no sign constraint, γ ∈ (-1, ∞)).
- If `coupling` has a third element, it is interpreted as:
    * a single `Symbol` (`:free`, `:activate`, or `:inhibit`) applied to all
      connections, or
    * a tuple/vector of such symbols, one per connection.

Valid modes and their semantics:

- `:activate`: γ ∈ (0, ∞)  (log / exp transform)
- `:inhibit`:  γ ∈ (-1, 0) (logit-style transform via `coupling_inhibitory_fwd` / `coupling_inhibitory_inv`)
- `:free`:     γ ∈ (-1, ∞) (shifted-log transform via `log_shift1` / `invlog_shift1`)

The returned value is a `Vector{Symbol}` of length `ncoupling(coupling)`.
"""
@inline function _normalize_coupling_mode(mode::Symbol)
    # Accept historical aliases but normalize to the canonical names.
    mode === :positive && return :activate
    mode === :negative && return :inhibit
    mode
end

function coupling_ranges(coupling)
    n = ncoupling(coupling)
    n == 0 && return Symbol[]
    if length(coupling) < 3
        return fill(:free, n)
    end
    ranges = coupling[3]
    if ranges isa Symbol
        m = _normalize_coupling_mode(ranges)
        return fill(m, n)
    elseif ranges isa Tuple || ranges isa AbstractVector
        length(ranges) == n || throw(ArgumentError("coupling_ranges length ($(length(ranges))) must match number of connections ($n)"))
        return [_normalize_coupling_mode(Symbol(r)) for r in ranges]
    else
        throw(ArgumentError("Unsupported coupling_ranges type $(typeof(ranges)); expected Symbol or tuple/vector of Symbols"))
    end
end

"""
    default_coupling_prior_mean(mode::Symbol) -> Float64

Return a sensible prior mean for a single coupling parameter based on its sign mode.

- `:activate` → `0.1`  (weak positive activation)
- `:inhibit`  → `-0.5` (weak inhibition; value is in the logit-space parameter)
- `:free`     → `0.0`  (no preferred direction)
"""
function default_coupling_prior_mean(mode::Symbol)
    m = _normalize_coupling_mode(mode)
    m === :activate && return 0.1
    m === :inhibit  && return -0.5
    return 0.0
end

"""
    default_coupling_prior_means(coupling) -> Vector{Float64}

Return a vector of default prior means (one per coupling connection) based on the
sign modes embedded in `coupling` (via `coupling_ranges`).
"""
function default_coupling_prior_means(coupling)
    modes = coupling_ranges(coupling)
    return [default_coupling_prior_mean(m) for m in modes]
end


"""
    make_structures(rinit, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, ...; trace_specs=[], dwell_specs=[])

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
- `trace_specs`: Per-unit observation windows for `tracejoint` (default: empty). For **coupled** `tracejoint`, if left empty and `dwell_specs` is empty, defaults are built from `traceinfo` and `zeromedian` (see `default_trace_specs_for_coupled`) so observed units are set on `TraceData`.
- `dwell_specs`: Container of dwell-time specs per unit (default: empty); used for multi-unit or observation-mapped dwell data when non-empty. Legacy single-unit when empty.

# Returns
- `Tuple`: (data, model, options) structures

# Notes
- Validates and adjusts model parameters
- Loads and processes data according to datatype
- Sets up priors and initial conditions
- Creates appropriate model structure based on parameters
- Handles hierarchical, coupled, and grid models
"""
function make_structures(rinit, datatype::String, dttype::Vector, datapath, gene, cell, datacond, traceinfo, infolder::String, label::String, fittedparam, fixedeffects, transitions, G, R, S, insertstep, coupling::Tuple=tuple(), grid=nothing, root=".", maxtime=60, elongationtime=6.0, priormean=Float64[], priorcv=10.0, nalleles=1, onstates=Int[], decayrate=-1.0, splicetype="", probfn=prob_Gaussian, noisepriors=[], hierarchical=tuple(), ratetype="median", propcv=0.01, samplesteps::Int=1000000, warmupsteps=0, annealsteps=0, temp=1.0, tempanneal=100.0, temprna=1.0, method=Tsit5(), zeromedian=false, datacol=3, ejectnumber=1, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[])
    gene = check_genename(gene, "[")
    insertstep = normalize_insertstep(R, insertstep)
    S = reset_S(S, R, insertstep)
    if G == 1 && !isempty(transitions)
        println("G=1, transitions are ignored")
        transitions = tuple()
    end
    nalleles = reset_nalleles(nalleles, coupling)
    infolder = folder_path(infolder, root, "results")
    datapath = folder_path(datapath, root, "data")
    data = load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna, datacol, zeromedian, yieldfactor, trace_specs, dwell_specs)
    decayrate = set_decayrate(decayrate, gene, cell, root)
    priormean, priorcv = set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid, datatype; priorcv=priorcv)
    rinit = isempty(hierarchical) ? set_rinit(rinit, priormean) : set_rinit(rinit, priormean, transitions, R, S, insertstep, noisepriors, length(data.trace[1]), coupling, grid; nhypersets=hierarchical[1])
    fittedparam = set_fittedparam(fittedparam, datatype, transitions, R, S, insertstep, noisepriors, coupling, grid)
    model = load_model(data, rinit, priormean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid, zeromedian, ejectnumber, 10, dwell_specs)
    if samplesteps > 0
        options = MHOptions(samplesteps, warmupsteps, annealsteps, Float64(maxtime), temp, tempanneal)
    else
        throw(ArgumentError("samplesteps must be greater than 0"))
    end
    return data, model, options
end

"""
    n_observed_trace_units(coupling)

Number of units with observed fluorescence traces, for building `trace_specs` / `data.units`.
For a 3-unit hidden-latent layout (`unit_model` = `(1,2,3)`), returns `2` (enhancer + gene only).
"""
function n_observed_trace_units(coupling::Tuple)
    length(coupling) < 1 && return 0
    um = coupling[1]
    if length(um) == 3 && um == (1, 2, 3)
        return 2
    end
    return length(um)
end

"""
    default_trace_specs_for_coupled(traceinfo, zeromedian, n_units::Int)

Construct `trace_specs` for coupled `tracejoint` when the caller omits them: one NamedTuple per
observed unit (`unit`, `interval`, `start`, `t_end`, `zeromedian`), using `traceinfo[1]` as the
sampling interval and `traceinfo[3]` as the trace end time (or a large upper bound if `< 0`).
"""
function default_trace_specs_for_coupled(traceinfo, zeromedian, n_units::Int)
    n_units >= 1 || throw(ArgumentError("n_units must be >= 1"))
    interval = Float64(traceinfo[1])
    tracetime = length(traceinfo) >= 3 ? Float64(traceinfo[3]) : -1.0
    t_end = tracetime < 0 ? 1.0e30 : tracetime
    zm_vec = zeromedian isa AbstractVector ? zeromedian : fill(zeromedian, n_units)
    length(zm_vec) >= n_units || throw(ArgumentError("zeromedian vector length must be >= n_units ($n_units)"))
    return [NamedTuple{(:unit, :interval, :start, :t_end, :zeromedian)}((u, interval, 0.0, t_end, zm_vec[u])) for u in 1:n_units]
end

"""
    default_trace_specs_for_coupled(traceinfo, zeromedian, observed_units::Vector{Int})

Same as `default_trace_specs_for_coupled(traceinfo, zeromedian, length(observed_units))`, but with
`unit` set to `observed_units[i]` (for models where observed units are not contiguous `1:n`, e.g. units 1 and 3).
"""
function default_trace_specs_for_coupled(traceinfo, zeromedian, observed_units::Vector{Int})
    n = length(observed_units)
    n >= 1 || throw(ArgumentError("observed_units must be non-empty"))
    interval = Float64(traceinfo[1])
    tracetime = length(traceinfo) >= 3 ? Float64(traceinfo[3]) : -1.0
    t_end = tracetime < 0 ? 1.0e30 : tracetime
    zm_vec = zeromedian isa AbstractVector ? zeromedian : fill(zeromedian, n)
    length(zm_vec) >= n || throw(ArgumentError("zeromedian vector length must be >= length(observed_units) ($n)"))
    return [NamedTuple{(:unit, :interval, :start, :t_end, :zeromedian)}((observed_units[i], interval, 0.0, t_end, zm_vec[i])) for i in 1:n]
end

"""
    get_units(dwell_specs, trace_specs)

Return the vector of observed unit indices (1-based) from dwell_specs or trace_specs.
Each spec in either container is a NamedTuple with a `unit` field.
When both are empty, returns Int[] (legacy: no unit subset).
When only one is non-empty, returns [spec.unit for spec in that container].
When both non-empty, returns units from dwell_specs (used e.g. for DwellTimeData).
"""
function get_units(dwell_specs, trace_specs)
    if !isempty(dwell_specs)
        return [spec.unit for spec in dwell_specs]
    elseif !isempty(trace_specs)
        return [spec.unit for spec in trace_specs]
    else
        return Int[]
    end
end

"""
    observed_units_from_dwell_specs(dwell_specs)

Return the vector of observed unit indices from dwell_specs (same order as dwell_specs).
Convenience alias for get_units(dwell_specs, []).
"""
observed_units_from_dwell_specs(dwell_specs) = get_units(dwell_specs, Int[])

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
    if length(traceinfo) < 4
        return 0.0
    else
        if traceinfo[4] isa Vector
            weight = Float64[]
            for f in traceinfo[4]
                push!(weight, max((1. - f) / f, 0.0))
            end
        else
            weight = max((1. - traceinfo[4]) / traceinfo[4], 0.0)
        end
        return weight
    end
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
function load_data_trace(datapath, label, gene, datacond, traceinfo, datatype::Symbol, col=3, zeromedian=false, yieldfactor::Float64=1.0, trace_specs=[])
    # When trace_specs is provided, override interval, zeromedian, and units from specs.
    # When empty, use legacy traceinfo / zeromedian parameters unchanged.
    if !isempty(trace_specs)
        interval_eff = trace_specs[1].interval
        zm_eff = [spec.zeromedian for spec in trace_specs]
        length(zm_eff) == 1 && (zm_eff = zm_eff[1])  # scalar for single-unit
        traceinfo_eff = (interval_eff, traceinfo[2:end]...)
        units_out = [spec.unit for spec in trace_specs]
    else
        traceinfo_eff = traceinfo
        zm_eff = zeromedian
        units_out = Int[]
    end
    if typeof(datapath) <: String
        tracer = read_tracefiles(datapath, datacond, traceinfo_eff, col)
    else
        tracer = read_tracefiles(datapath[1], datacond, traceinfo_eff, col)
    end
    (length(tracer) == 0) && throw("No traces")
    trace, tracescale = zero_median(tracer, zm_eff)
    println("number of traces: ", length(trace))
    println("datapath: ", datapath)
    println("datacond: ", datacond)
    println("traceinfo: ", traceinfo_eff)
    nframes = round(Int, mean(size.(trace, 1)))  #mean number of frames of all traces
    if length(traceinfo_eff) < 4
        weight = 0.0
        background = 0.0
    else
        weight = set_trace_weight(traceinfo_eff)
        background = set_trace_background(traceinfo_eff)
    end
    if datatype == :trace || datatype == :tracejoint
        return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo_eff[1], (trace, background, weight, nframes, tracescale), units_out)
    elseif datatype == :tracerna
        len, h = read_rna(gene, datacond, datapath[2])
        yield = yieldfactor < 1.0 ? (yieldfactor, nhist_loss(len, yieldfactor)) : yieldfactor
        return TraceRNAData{typeof((trace, background, weight, nframes)),typeof(h)}(label, gene, traceinfo_eff[1], (trace, background, weight, nframes), len, h, yield, units_out)
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
    return TraceData{typeof(label),typeof(gene),Tuple}(label, gene, traceinfo[1], (trace, Vector[], 0.0, 1), Int[])
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
function load_data(datatype, dttype, datapath, label, gene, datacond, traceinfo, temprna, datacol=3, zeromedian=false, yieldfactor::Float64=1.0, trace_specs=[], dwell_specs=[])
    dt = normalize_datatype(datatype)
    # Units for non-trace data types come from dwell_specs; for trace types they come from trace_specs.
    units = get_units(dwell_specs, trace_specs)

    if dt == :rna
        len, h = read_rna(gene, datacond, datapath)
        yield = yieldfactor < 1.0 ? (yieldfactor, nhist_loss(len, yieldfactor)) : yieldfactor
        return RNAData{typeof(len),typeof(h)}(label, gene, len, h, yield, units)

    elseif dt == :rnacount
        countsRNA, yield_vec, nRNA = read_rnacount(gene, datacond, datapath)
        min_yield = minimum(yield_vec)
        nRNA_computed = min_yield < 1.0 ? nhist_loss(nRNA, min_yield) : nRNA
        return RNACountData(label, gene, nRNA_computed, countsRNA, yield_vec, units)

    elseif dt == :rnaonoff
        len, h = read_rna(gene, datacond, datapath[1])
        h = div.(h, temprna)
        LC = readfile(gene, datacond, datapath[2])
        yield = yieldfactor < 1.0 ? (yieldfactor, nhist_loss(len, yieldfactor)) : yieldfactor
        return RNAOnOffData(label, gene, len, h, LC[:, 1], LC[:, 2], LC[:, 3], yield, units)

    elseif dt == :rnadwelltime
        len, h = read_rna(gene, datacond, datapath[1])
        h = div.(h, temprna)
        bins, DT = read_dwelltimes(datapath[2:end])
        yield = yieldfactor < 1.0 ? (yieldfactor, nhist_loss(len, yieldfactor)) : yieldfactor
        return RNADwellTimeData(label, gene, len, h, bins, DT, dttype, yield, units)

    elseif dt == :dwelltime
        bins, DT = read_dwelltimes(datapath)
        du = isempty(dwell_specs) ? units : observed_units_from_dwell_specs(dwell_specs)
        return DwellTimeData(label, gene, bins, DT, dttype, du)

    elseif dt ∈ TRACE_DATATYPES
        if dt == :tracegrid
            return load_data_tracegrid(datapath, label, gene, datacond, traceinfo)
        else
            return load_data_trace(datapath, label, gene, datacond, traceinfo, dt, datacol, zeromedian, yieldfactor, trace_specs)
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
- Creates `TDCoupledFullComponents` for multi-unit dwell time analysis (default coupled stack)
- Applies coupling to sojourn states and nonzero rows
- Handles interactions between different model units
- Used for coupled dwell time distribution modeling
"""
function make_reporter_components_DT(transitions, G::Tuple, R::Tuple, S::Tuple, insertstep, splicetype, onstates, dttype, coupling)
    sojourn = sojourn_states(onstates, G, R, S, insertstep, dttype)
    components = TDCoupledFullComponents(coupling, transitions, G, R, S, insertstep, sojourn, dttype)
    unit_model = coupling[1]
    nT_vec = collect(T_dimension(G, R, S, unit_model))
    sojourn_c = [
        [begin
             soj_unit = sojourn[unit_model[k]][i]
             if occursin("G", dttype[unit_model[k]][i])
                 α = unit_model[k]
                 soj_unit = g_sojourn_to_T_sojourn(soj_unit, G[α], R[α], S[α])
             end
             full_state_indices_for_unit_sojourn(k, soj_unit, nT_vec)
         end
         for i in eachindex(sojourn[unit_model[k]])]
        for k in eachindex(unit_model)]
    return (sojourn_c, sojourn_c), components
end

"""
    make_reporter_components_DT(transitions, G::Tuple, R::Tuple, S::Tuple, insertstep, splicetype, dwell_specs, coupling)

Coupled dwell time reporter components from per-unit dwell_specs.

dwell_specs lists only observed units (hidden units may be omitted). Builds full
per-unit onstates and dttype in unit_model order; for units not in dwell_specs uses
a valid placeholder (R=0 needs explicit G-state onstates).
"""
function make_reporter_components_DT(transitions, G::Tuple, R::Tuple, S::Tuple, insertstep, splicetype, dwell_specs, coupling)
    unit_model = coupling[1]
    n_units = length(unit_model)
    placeholder = dwell_specs[1]
    onstates_full = Vector{typeof(placeholder.onstates)}(undef, n_units)
    dttype_full = Vector{typeof(placeholder.dttype)}(undef, n_units)
    for k in 1:n_units
        j = findfirst(s -> s.unit == unit_model[k], dwell_specs)
        if j !== nothing
            onstates_full[k] = dwell_specs[j].onstates
            dttype_full[k] = dwell_specs[j].dttype
        else
            dttype_full[k] = placeholder.dttype
            onstates_full[k] = R[k] == 0 ? [[G[k]] for _ in 1:length(placeholder.dttype)] : placeholder.onstates
        end
    end
    make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, onstates_full, dttype_full, coupling)
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
- Creates `TCoupledFullComponents` for multi-unit trace analysis (default coupled stack)
- Used for coupled fluorescence trace data modeling
"""
function make_reporter_components(transitions::Tuple, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype, onstates, probfn, noisepriors, coupling)
    reporter = HMMReporter[]
    nunits = length(G)
    if probfn isa Union{Tuple,Vector}
        # Pad to nunits if shorter (hidden units at the end have nnoise=0 so probfn is unused)
        if length(probfn) < nunits
            probfn = [probfn..., fill(probfn[end], nunits - length(probfn))...]
        end
    else
        probfn = fill(probfn, nunits)
    end
    n_per_state = num_reporters_per_state(G, R, S, insertstep, coupling[1], onstates)
    for i in eachindex(G)
        nnoise = length(noisepriors[i])
        n = num_rates(transitions[i], R[i], S[i], insertstep[i])
        weightind = occursin("Mixture", "$(probfn)") ? n + nnoise : 0
        push!(reporter, HMMReporter(nnoise, n_per_state[i], probfn[i], weightind, off_states(n_per_state[i]), collect(n+1:n+nnoise)))
    end
    components = TCoupledFullComponents(coupling, transitions, G, R, S, insertstep, splicetype)
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
    # Use nRNA_true from data (already computed in load_data)
    # For RNACountData: data.nRNA is already nRNA_true (computed using nhist_loss)
    # For other types: extract nRNA_true from yield tuple (computed in load_data)
    if typeof(data) <: RNACountData
        nRNA_size = data.nRNA  # Already nRNA_true from load_data
    else
        nRNA_size = get_nRNA_true(data.yield, data.nRNA)  # Extract nRNA_true from yield tuple
    end
    components = MComponents(transitions, G, R, nRNA_size, decayrate, splicetype, ejectnumber)
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
    # Use nRNA_true if available (from yield tuple), otherwise use observed nRNA
    nRNA_size = get_nRNA_true(data.yield, data.nRNA)
    return onstates, MTAIComponents(transitions, G, R, S, insertstep, onstates, nRNA_size, decayrate, splicetype, ejectnumber)
end

"""
    full_onstates_dttype_from_dwell_specs(dwell_specs, n_units)

Build full-length onstates and dttype (length n_units) for model construction when some units are hidden.
Observed units get their spec; hidden units get a placeholder (copy of first spec) so coupled TD matrices are built for all units.
"""
function full_onstates_dttype_from_dwell_specs(dwell_specs, n_units)
    observed = observed_units_from_dwell_specs(dwell_specs)
    full_onstates = Vector{typeof(dwell_specs[1].onstates)}(undef, n_units)
    full_dttype = Vector{typeof(dwell_specs[1].dttype)}(undef, n_units)
    placeholder_o = dwell_specs[1].onstates
    placeholder_d = dwell_specs[1].dttype
    for k in 1:n_units
        j = findfirst(==(k), observed)
        if j !== nothing
            full_onstates[k] = dwell_specs[j].onstates
            full_dttype[k] = dwell_specs[j].dttype
        else
            full_onstates[k] = placeholder_o
            full_dttype[k] = placeholder_d
        end
    end
    full_onstates, full_dttype
end

"""
    make_reporter_components(data::DwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1; dttype_full=nothing)

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
- `dttype_full`: Optional full-length dwell type vector for all units (used when some units are hidden).

# Returns
- `Tuple`: (reporter, components)

# Notes
- Delegates to make_reporter_components_DT for dwell time analysis
- Uses dwell time types from data structure (or dttype_full when provided)
- Creates appropriate components for dwell time distribution modeling
- Used for analyzing time spent in different model states
"""
function make_reporter_components(data::DwellTimeData, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber=1; dttype_full=nothing, dwell_specs=[])
    if !isempty(dwell_specs)
        make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, dwell_specs, coupling)
    else
        dtype = (dttype_full !== nothing) ? dttype_full : data.DTtypes
        make_reporter_components_DT(transitions, G, R, S, insertstep, splicetype, onstates, dtype, coupling)
    end
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
    # Use nRNA_true if available (from yield tuple), otherwise use observed nRNA
    nRNA_size = get_nRNA_true(data.yield, data.nRNA)
    mcomponents = MComponents(transitions, G, R, nRNA_size, decayrate, splicetype, ejectnumber)
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
    # Use nRNA_true if available (from yield tuple), otherwise use observed nRNA
    nRNA_size = get_nRNA_true(data.yield, data.nRNA)
    mcomponents = MComponents(transitions, G, R, nRNA_size, decayrate, splicetype, ejectnumber)
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
    c = ncoupling(coupling)
    collect(n-g-c+1:n-g)
end

function coupling_indices_full(transitions, R, S, insertstep, reporter, coupling, grid)
    indices = coupling_indices(transitions, R, S, insertstep, reporter, coupling, grid)
    targets = Vector{Tuple{Int,Int}}(undef, length(coupling[2]))
    for (i, c) in enumerate(coupling[2])
        targets[i] = (coupling[1][c[3]], c[4])
    end
    indices, targets
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
        modes = coupling_ranges(coupling)
        for i in eachindex(couplingindices)
            mode = i <= length(modes) ? modes[i] : :free
            if mode === :activate
                # γ ∈ (0, ∞)
                push!(ftransforms, log)
                push!(invtransforms, exp)
                push!(sigmatransforms, sigmalognormal)
            elseif mode === :inhibit
                # γ ∈ (-1, 0)
                push!(ftransforms, coupling_inhibitory_fwd)
                push!(invtransforms, coupling_inhibitory_inv)
                push!(sigmatransforms, sigmanormal)
            else
                # :free — γ ∈ (-1, ∞)
                push!(ftransforms, log_shift1)
                push!(invtransforms, invlog_shift1)
                push!(sigmatransforms, sigmalognormal)
            end
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
        fset = copy(ftransforms)
        iset = copy(invtransforms)
        sset = copy(sigmatransforms)
        nindividuals = length(data.trace[1])
        for i in 1:hierarchical[1]-1
            ftransforms = vcat(ftransforms, fset)
            invtransforms = vcat(invtransforms, iset)
            sigmatransforms = vcat(sigmatransforms, sset)
        end
        ftransforms = vcat(ftransforms, repeat(fset, nindividuals))
        invtransforms = vcat(invtransforms, repeat(iset, nindividuals))
        sigmatransforms = vcat(sigmatransforms, repeat(sset, nindividuals))
    end
    Transformation(ftransforms, invtransforms, sigmatransforms)
end


"""
    load_model(data, r, rmean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid, zeromedian=false, ejectnumber=1, factor=10, dwell_dttype_full=nothing)

Construct and return the appropriate model struct for the given data and options.
When `dwell_dttype_full` is provided and data is DwellTimeData, it is used (full-length dttype for all units) so that coupled TD matrices are built correctly when some units are hidden.

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
function load_model(data, r, rmean, fittedparam, fixedeffects, transitions, G, R, S, insertstep, splicetype, nalleles, priorcv, onstates, decayrate, propcv, probfn, noisepriors, method, hierarchical, coupling, grid, zeromedian=false, ejectnumber=1, factor=10, dwell_specs=[])
    dwell_specs = (dwell_specs === nothing) ? [] : dwell_specs
    # For coupled dwell models, extract full onstates from dwell_specs here.
    if !isempty(dwell_specs) && !isempty(coupling) && G isa Tuple
        n_units = length(coupling[1])
        full_onstates, _ = full_onstates_dttype_from_dwell_specs(dwell_specs, n_units)
        onstates = full_onstates
    end
    insertstep = normalize_insertstep(R, insertstep)
    reporter, components = if data isa DwellTimeData
        make_reporter_components(data, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber; dwell_specs=dwell_specs)
    else
        make_reporter_components(data, transitions, G, R, S, insertstep, splicetype, onstates, decayrate, probfn, noisepriors, coupling, ejectnumber)
    end

    nrates = num_rates(transitions, R, S, insertstep)
    ratetransforms = make_ratetransforms(data, nrates, transitions, G, R, S, insertstep, reporter, coupling, grid, hierarchical, zeromedian)

    if !isempty(coupling)
        couplingindices = coupling_indices(transitions, R, S, insertstep, reporter, coupling, grid)
        ncp = ncoupling(coupling)
        couplingtrait = CouplingTrait(ncp, couplingindices)
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
    [rm; fill(1.0, ncoupling(coupling))]
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
    AbstractPriorContext

Supertype for which **default prior mean recipe** applies given data (`datatype`) and model layout
(`coupling`, `R` scalar vs tuple, `grid`, `hierarchical`, `transitions`). Used by [`prior_context`](@ref)
and [`set_priormean_empty`](@ref). You may add subtypes and methods to `set_priormean_empty` for custom
workflows.
"""
abstract type AbstractPriorContext end

"""Coupled multi-unit model: use `prior_ratemean(..., coupling)`."""
struct PriorContextCoupled <: AbstractPriorContext end

"""Trace-like single-unit GRSM (`datatype` contains `\"trace\"`), no grid/hierarchical, ≥2 transition pairs: use `prior_ratemean_trace`."""
struct PriorContextTraceSingleUnit <: AbstractPriorContext end

"""Default single-unit GRSM prior (`prior_ratemean`): RNA, hierarchical trace, grid, short transition list, etc."""
struct PriorContextGenericSingle <: AbstractPriorContext end

"""
    prior_context(datatype, coupling, R, grid, hierarchical, transitions) -> AbstractPriorContext

Classify data + model for default **base** prior means (before grid append and hierarchical hyper
rows). **Order:** (1) coupled; (2) trace single-unit recipe; (3) generic single-unit.

See also: [`set_priormean`](@ref), [`set_priormean_empty`](@ref).
"""
function prior_context(datatype::String, coupling, R, grid, hierarchical, transitions)
    if !isempty(coupling)
        return PriorContextCoupled()
    end
    if occursin("trace", datatype) && R isa Int && isnothing(grid) && isempty(hierarchical) && length(transitions) >= 2
        return PriorContextTraceSingleUnit()
    end
    return PriorContextGenericSingle()
end

"""
    set_priormean_empty(ctx, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling, grid)

Build the **base** prior mean vector for empty `priormean`; `set_priormean` then appends grid and
hierarchical blocks. Extend by adding methods on new `AbstractPriorContext` subtypes.
"""
function set_priormean_empty(::PriorContextTraceSingleUnit, transitions, R::Int, S::Int, insertstep, decayrate, noisepriors, elongationtime, coupling, grid)
    prior_ratemean_trace(transitions, R, S, insertstep, Float64(decayrate), noisepriors, Float64(elongationtime))
end

function set_priormean_empty(::PriorContextCoupled, transitions, R::Tuple, S::Tuple, insertstep::Tuple, decayrate, noisepriors, elongationtime, coupling, grid)
    prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling)
end

function set_priormean_empty(::PriorContextGenericSingle, transitions, R::Int, S::Int, insertstep, decayrate, noisepriors, elongationtime, coupling, grid)
    prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime)
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
- For **trace** single-unit fits with empty `priormean`, `set_priormean` uses
  `prior_ratemean_trace` instead (see `datatype` in `make_structures`).
"""
function prior_ratemean(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors::Vector, elongationtime::Float64, initprior::Float64=0.1)
    [fill(0.01, length(transitions)); initprior; fill(R / elongationtime, R); fill(0.1, max(0, S - insertstep + 1)); decayrate; noisepriors]
end

"""
    prior_ratemean_trace(transitions, R::Int, S::Int, insertstep, decayrate, noisepriors, elongationtime, initprior=0.1)

Default **prior means** for single-unit **trace** GRSM fits (MS2-style live traces): first two G
transition edges at `0.001`, remaining edges at `0.01`, then initiation, elongation, splicing, decay,
and noise blocks matching the legacy `set_variables` / `trace_prior_variables` recipe.

This is used automatically by `set_priormean` when `priormean` is empty, `datatype` contains
`"trace"`, and the model is a non-coupled, non-hierarchical, non-grid single unit with at least two
transition pairs. Otherwise `prior_ratemean` (generic GRSM) is used.

Coupled (`tracejoint`), hierarchical trace, or grid models keep the previous generic `prior_ratemean`
path unless you pass explicit `priormean`.
"""
function prior_ratemean_trace(transitions, R::Int, S::Int, insertstep, decayrate::Float64, noisepriors::Vector, elongationtime::Float64, initprior::Float64=0.1)
    n_t = length(transitions)
    if n_t < 2
        return prior_ratemean(transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, initprior)
    end
    [fill(0.001, 2); fill(0.01, n_t - 2); initprior; fill(R / elongationtime, R); fill(0.05, max(0, S - insertstep + 1)); decayrate; noisepriors]
end

"""
    priorcv_trace_grsm(transitions, R, S, insertstep, noisepriors)

Structured **prior CV** vector matching `prior_ratemean_trace` (same block layout). Applied inside
[`set_priormean`](@ref) when `priorcv` is still a scalar and `priormean` was auto-filled for a trace
single-unit fit (`priorcv` keyword passed).
"""
function priorcv_trace_grsm(transitions, R::Int, S::Int, insertstep, noisepriors::Vector)
    nn = length(noisepriors)
    noise_cv = nn == 4 ? [0.2, 0.2, 0.1, 0.1] : fill(0.2, nn)
    [fill(1.0, length(transitions)); 0.2; fill(0.2, R); fill(0.2, max(0, S - insertstep + 1)); 0.1; noise_cv]
end

"""
    trace_prior_variables(transitions, R, S, insertstep, hierarchical, fixed, noisepriors, elongationtime, propcv;
        initprior=0.1, ejectprior=0.05, nfitted=[1, 3], coupling=false, fitcoupling=false)

Construct prior means, prior CVs, fitted-parameter indices, fixed-effects tuple, hierarchical spec, ODE
`method` (`Tsit5()` or `(Tsit5(), true)`), and proposal CV for **trace-oriented GRSM** single-unit models (common for live-cell
traces and single-gene batches). For **non-hierarchical, non-coupling** cases, `priormean` and `priorcv`
match the same `prior_ratemean_trace` / `priorcv_trace_grsm` construction used automatically by
`fit(...)` when `priormean` is empty (see `set_priormean`). For hierarchical or fixed-rates recipes,
this helper still fills `fittedparam`, `fixedeffects`, and hyperparameter blocks.

**Prefer** calling `fit` with empty `priormean` for standard trace fits so defaults stay in one place.

# Returns
`(priormean, priorcv, fittedparam, fixedeffects, hierarchical, method, propcv)` where `hierarchical` is
`tuple()` or a 3-tuple compatible with the `fit` keyword `hierarchical`.

# Notes
- `ejectprior` is reserved for compatibility with older scripts; the current prior vector uses the
  fixed decay constant `0.03165055618995184` in the same position as generic `prior_ratemean` trace priors.
- Does not set `samplesteps`; choose `samplesteps` in your run spec (e.g. `1_000_000` non-hierarchical,
  `100_000` hierarchical) to match your cluster budget.
"""
function trace_prior_variables(transitions, R::Int, S::Int, insertstep, hierarchical::Bool, fixed::Bool, noisepriors, elongationtime, propcv, initprior=0.1, ejectprior=0.05, nfitted=[1, 3], coupling=false, fitcoupling=false)
    n = num_rates(transitions, R, S, insertstep)
    if hierarchical
        method = (Tsit5(), true)
        propcv == 0.0 && (propcv = 0.001)
    else
        method = Tsit5()
        propcv == 0.0 && (propcv = 0.005)
    end
    if fixed
        if hierarchical
            fitted = [collect(1:length(transitions) + 2); collect(length(transitions) + R + 1:n - 1 - max(0, S - 1))]
        else
            fitted = [collect(1:length(transitions) + 2); collect(length(transitions) + R + 1:n - 1 - max(0, S - 1)); n .+ nfitted]
        end
        f = R > 2 ? (collect(length(transitions) + 2:length(transitions) + R),) : tuple()
    else
        if hierarchical
            fitted = collect(1:n - 1 - max(0, S - 1))
        else
            fitted = [collect(1:n - 1 - max(0, S - 1)); n .+ nfitted]
        end
        f = tuple()
    end
    if fitcoupling
        fitted = vcat(fitted, [n + 5])
    end
    priormean = prior_ratemean_trace(transitions, R, S, insertstep, 0.03165055618995184, noisepriors, elongationtime, initprior)
    priorcv = priorcv_trace_grsm(transitions, R, S, insertstep, noisepriors)
    if coupling
        priormean = vcat(priormean, [0.0])
        priorcv = vcat(priorcv, [2.0])
    end
    if hierarchical
        hypercv = [fill(1.0, length(transitions)); 1.0; fill(0.25, R - 1); 1.0; fill(1.0, max(0, S - insertstep + 1)); 1.0; fill(0.25, 4)]
        coupling && (hypercv = vcat(hypercv, [1.0]))
        priormean = [priormean; hypercv]
        priorcv = [priorcv; fill(2.0, length(priorcv))]
        h = (2, n .+ nfitted, tuple())
    else
        h = tuple()
    end
    return priormean, priorcv, fitted, f, h, method, propcv
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
- Appends coupling parameter prior means using coupling sign modes:
  * `:free`     ⇒ 0.0
  * `:activate` ⇒ +0.5
  * `:inhibit`  ⇒ -0.5
- Used for coupled model initialization
"""
function prior_ratemean(transitions, R::Tuple, S::Tuple, insertstep::Tuple, decayrate, noisepriors::Union{Vector,Tuple}, elongationtime::Union{Vector,Tuple}, coupling, initprior=[0.1, 0.1])
    length(initprior) < length(R) && (initprior = vcat(initprior, fill(0.1, length(R) - length(initprior))))
    rm = Float64[]
    for i in eachindex(R)
        append!(rm, prior_ratemean(transitions[i], R[i], S[i], insertstep[i], decayrate, noisepriors[i], elongationtime[i], initprior[i]))
    end
    ranges = coupling_ranges(coupling)
    coupling_means = Float64[
        mode === :activate ? 0.5 :
        mode === :inhibit  ? -0.5 :
        0.0
        for mode in ranges
    ]
    [rm; coupling_means]
end

"""
    set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid, datatype=""; ctx=nothing, priorcv=nothing)

Set default prior means (and optionally prior CVs) when the user left them at defaults.

# Arguments
- `priormean`: Prior means (if empty, filled via [`prior_context`](@ref) and [`set_priormean_empty`](@ref))
- `transitions, R, S, insertstep`, `decayrate`, `noisepriors`, `elongationtime`, `hierarchical`, `coupling`, `grid`, `datatype`: Same role as in [`make_structures`](@ref).
- `ctx`: Optional [`AbstractPriorContext`](@ref); when `nothing`, `prior_context(datatype, coupling, R, grid, hierarchical, transitions)` is used.
- `priorcv`: When `nothing` (default), only prior means are considered and the return value is just `priormean`.
  When passed (e.g. scalar default from `make_structures`), returns `(priormean, priorcv)`; if `priormean` was
  empty and `priorcv` is still a scalar and the context is [`PriorContextTraceSingleUnit`](@ref), `priorcv`
  is replaced by [`priorcv_trace_grsm`](@ref) to match the filled `priormean`.

# Returns
- If `priorcv === nothing`: `Vector{Float64}` prior means only.
- Otherwise: `(priormean, priorcv)` with the same conventions as above.

# Notes
- All branching on empty `priormean`, model/data context, and trace structured CV lives here, not in callers.
"""
function set_priormean(priormean, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, hierarchical, coupling, grid, datatype::String=""; ctx=nothing, priorcv=nothing)
    adjust_priorcv = priorcv !== nothing
    priormean_was_empty = isempty(priormean)
    ctx = ctx === nothing ? prior_context(datatype, coupling, R, grid, hierarchical, transitions) : ctx
    if !priormean_was_empty
        return adjust_priorcv ? (priormean, priorcv) : priormean
    end
    priormean = set_priormean_empty(ctx, transitions, R, S, insertstep, decayrate, noisepriors, elongationtime, coupling, grid)
    if !isnothing(grid)
        priormean = prior_ratemean_grid(priormean)
    end
    if !isempty(hierarchical)
        priormean = prior_ratemean_hierarchical(priormean, prior_hypercv(transitions, R, S, insertstep, noisepriors, coupling, grid), hierarchical[1])
    end
    if adjust_priorcv && !(priorcv isa AbstractVector) && ctx isa PriorContextTraceSingleUnit
        priorcv = priorcv_trace_grsm(transitions, R, S, insertstep, noisepriors)
    end
    return adjust_priorcv ? (priormean, priorcv) : priormean
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
    Vector{Float64}(r)
end

"""
    set_rinit(r, priormean, transitions, R, S, insertstep, noisepriors, nindividuals, coupling=tuple(), grid=nothing; nhypersets=1, ...)

Set initial parameters for hierarchical models.

# Arguments
- `r`: Initial parameters
- `priormean`: Prior means
- `transitions, R, S, insertstep`: Model structure parameters
- `noisepriors`: Noise priors
- `nindividuals`: Number of individuals (trace trials)
- `coupling`: Coupling structure (default: empty tuple)
- `grid`: Grid parameter (default: nothing)
- `nhypersets`: Number of hyperparameter blocks at the start of `r` (must match `hierarchical[1]` in `load_model`; default `1`)
- `minval`: Minimum valid value (default: 1e-7)
- `maxval`: Maximum valid value (default: 300.0)

# Returns
- `Vector{Float64}`: Valid initial parameters for hierarchical model

# Notes
- Builds `r` with `(nhypersets + nindividuals)` blocks of `n_all_params` each, matching `make_hierarchical`.
- If `length(priormean) ≥ nhypersets * n_all_params`, hyper blocks are taken from successive slices of `priormean`; otherwise each hyper block is filled from the first block of `priormean` / loaded rates.
- Used for hierarchical model initialization
"""
function set_rinit(r, priormean, transitions, R, S, insertstep, noisepriors, nindividuals, coupling=tuple(), grid=nothing; nhypersets::Int=1, minval=1e-7, maxval=300.0)
    c = ncoupling(coupling)
    g = isnothing(grid) ? 0 : 1
    n_all_params = num_all_parameters(transitions, R, S, insertstep, noisepriors) + c + g
    # Stacked layout matches `make_hierarchical`: nhypersets full blocks + one block per individual
    total_expected = (nhypersets + nindividuals) * n_all_params
    if isempty(r) || any(isnan.(r)) || any(isinf.(r)) || length(r) < total_expected
        isempty(r) && println("No rate file, set rate to prior")
        any(isnan.(r)) && println("r contains NaN, set rate to prior")
        any(isinf.(r)) && println("r contains Inf, set rate to prior")
        length(r) < total_expected && !isempty(r) && println("Rate file too short ($(length(r)) < $total_expected), expanding from prior")
        # Seed individual inits from loaded rates when available, else fall back to priormean
        seed = Vector{Float64}(
            (length(r) >= n_all_params && !any(isnan.(Float64.(r[1:n_all_params])))) ?
            r[1:n_all_params] : priormean[1:n_all_params]
        )
        # Validate coupling parameters match their specified sign mode
        if c > 0
            modes = coupling_ranges(coupling)
            coupling_start = n_all_params - g - c + 1
            for k in 1:c
                idx = coupling_start + k - 1
                mode = k <= length(modes) ? modes[k] : :free
                γ = seed[idx]
                # Check if initial coupling parameter sign matches the specified mode
                if mode === :activate && γ <= 0.0
                    throw(ArgumentError("Coupling parameter $k (value $γ) is invalid for mode $mode; must be > 0"))
                elseif mode === :inhibit && (γ >= 0.0 || γ <= -1.0)
                    throw(ArgumentError("Coupling parameter $k (value $γ) is invalid for mode $mode; must be in (-1, 0)"))
                elseif mode === :free && γ <= -1.0
                    throw(ArgumentError("Coupling parameter $k (value $γ) is invalid for mode $mode; must be > -1"))
                end
            end
        end
        r = Vector{Float64}()
        if length(priormean) >= nhypersets * n_all_params
            for b in 1:nhypersets
                lo = (b - 1) * n_all_params + 1
                hi = b * n_all_params
                append!(r, priormean[lo:hi])
            end
        else
            for _ in 1:nhypersets
                append!(r, seed)
            end
        end
        for _ in 1:nindividuals
            append!(r, seed)
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
    [fittedparam; collect(fittedparam[end]+1:fittedparam[end]+ncoupling(coupling))]
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

Burst size = mean emission rate (while in state G) / k_off (rate of leaving state G)

# Arguments
- `r`: Rate parameters
- `transitions`: Model transitions
- `G, R, S, insertstep`: Model structure parameters
- `splicetype`: Splicing type (default: "")
- `ejectnumber`: Number of mRNAs per burst (default: 1)

# Returns
- `Float64`: Mean burst size (mean number of mRNA produced per burst)

# Notes
- For classic telegraph models (GM, R=0), burst size = eject_rate / off_rate
- For GRSM models (R>0), this function:
  1. Computes k_off: sum of all transition rates from gene state G to other gene states
  2. Computes steady-state distribution over R states while in gene state G
  3. Computes mean emission rate = sum(P(R config) * emission_rate(R config))
  4. Returns mean emission rate / k_off

The lifetime of state G is exponential with rate k_off, so the mean number of mRNA
produced during a burst is the mean emission rate divided by k_off.
"""
function burstsize(r, transitions, G, R, S, insertstep, splicetype="", ejectnumber=1)
    ntransitions = num_rates(transitions, R, S, insertstep)
    
    # Step 1: Compute k_off (rate of leaving gene state G)
    # Find all transitions that start from gene state G
    k_off = 0.0
    nG_transitions = length(transitions)
    for (i, trans) in enumerate(transitions)
        if trans[1] == G  # Transition starts from state G
            k_off += r[i]
        end
    end
    
    if k_off == 0.0
        # If state G never transitions away, burst size is infinite
        return Inf
    end
    
    # Step 2: Compute mean emission rate while in state G
    # For R=0: emission rate is just the ejection/initiation rate
    if R == 0
        # For R=0, ejection rate is at index ntransitions (before decay)
        # Rate order: G transitions, then ejection, then decay
        emission_rate = r[ntransitions] * ejectnumber
        return emission_rate / k_off
    end
    
    # For R>0: Need to compute steady-state distribution over R states while in gene state G
    # Set up reduced master equation: only R state transitions, no gene state transitions, no mRNA
    # State space: only R configurations with gene state = G
    nR_configs = 2^R  # Number of R configurations
    nT = G * 2^R
    
    # Identify R configurations that correspond to gene state G
    # For gene state G: config indices are G, 2*G, 3*G, ..., nR_configs*G
    # But we need to map these to a reduced state space (just R configs)
    R_configs_in_G = Int[]
    for z in 1:nR_configs
        config_idx = state_index(G, G, z)  # Gene state G, R config z
        push!(R_configs_in_G, config_idx)
    end
    
    # Build reduced transition matrix for R states only (within gene state G)
    # This includes: R step transitions, S transitions, and ejection
    # Need to extract relevant elements from the full TComponents
    tcomponents = TComponents(transitions, G, R, S, insertstep, splicetype)
    T_full = make_mat_T(tcomponents, r)
    
    # Extract submatrix for R states in gene state G
    # Only include transitions within gene state G (R state transitions)
    # Exclude transitions from gene state G to other gene states
    nR = length(R_configs_in_G)
    T_reduced = zeros(nR, nR)
    for i in 1:nR
        for j in 1:nR
            if i != j
                # Off-diagonal: transitions between R configs within gene state G
                T_reduced[i, j] = T_full[R_configs_in_G[i], R_configs_in_G[j]]
            end
        end
        # Diagonal: negative sum of all outgoing transitions within gene state G only
        # (excludes transitions to other gene states, which we don't want in reduced matrix)
        T_reduced[i, i] = -sum(T_reduced[i, :])
    end
    
    # Compute steady-state distribution over R states
    T_reduced_sparse = sparse(T_reduced)
    pss_R = normalized_nullspace(T_reduced_sparse)
    
    # Step 3: Compute mean emission rate
    # Emission occurs from R states where the final R step is occupied
    # For each R config z, check if final R step is occupied and get emission rate
    mean_emission_rate = 0.0
    nG_rates = length(transitions)
    # Rate order: G transitions (1 to nG_rates), initiation (nG_rates+1), 
    # R transitions (nG_rates+2 to nG_rates+1+R), S transitions, ejection (ntransitions), decay (ntransitions+1)
    # Ejection rate is at index ntransitions (before decay)
    eject_rate_idx = ntransitions
    eject_rate = r[eject_rate_idx] * ejectnumber
    
    for (i, config_idx) in enumerate(R_configs_in_G)
        # Extract R config z from config_idx
        z = div(config_idx - 1, G) + 1
        
        # Check if final R step is occupied in config z
        # R config z is represented as binary: digits(z-1, base=2, pad=R)
        R_binary = digits(z - 1, base=2, pad=R)
        if R > 0 && R_binary[R] == 1  # Final R step is occupied
            mean_emission_rate += pss_R[i] * eject_rate
        end
    end
    
    # Step 4: Burst size = mean emission rate / k_off
    mean_emission_rate / k_off
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
        cond = i <= length(datacond) ? datacond[i] : "(hidden)"
        println("datacond: ", cond, ", GRSinsertstep: ", G[i], R[i], S[i], insertstep[i], ", nalleles: ", nalleles)
    end
    println("data: ", datapath)
    println("in: ", infolder, " out: ", resultfolder)
    println("maxtime: ", maxtime, ", propcv: ", propcv)
end



"""
    finalize(data,model,fits,stats,measures,temp,resultfolder,optimized,burst,writesamples)

Write out run results and print final loglikelihood and deviance. When invoked from the top-level `fit(; key=..., ...)`, writes info_<stem>.toml for reproducibility (run_spec/name_override are taken from internal refs).
"""
function finalize(data, model, fits, stats, measures, temp, writefolder, optimized, burst, writesamples)
    run_spec = _current_run_spec[]
    name_override = _current_name_override[]
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
    writeall(writefolder, fits, stats, measures, data, temp, model, optimized=optimized, burst=burst, writesamples=writesamples, name_override=name_override, run_spec=run_spec)
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
            println(gene, " has no decay time, set to 1.0")
            return 1.0
        end
    else
        path = get_file(root, "data/halflives", cell, "csv")
        if isnothing(path)
            println(gene, " has no decay time, set to 1.0")
            return 1.0
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




