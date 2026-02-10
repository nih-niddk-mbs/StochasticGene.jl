# This file is part of StochasticGene.jl

# simulator.jl
# Functions to simulate Markov gene transcription models
# Uses hybrid first and next reaction method

"""
    Reaction

Structure representing a reaction in the GRSM model.

# Fields
- `action::Int`: Type of reaction (1=activateG!, 2=deactivateG!, 3=transitionG!, 4=initiate!, 5=transitionR!, 6=eject!, 7=splice!, 8=decay!)
- `index::Int`: Rate index for the reaction
- `disabled::Vector{Int64}`: Indices of reactions that are disabled by this reaction
- `enabled::Vector{Int64}`: Indices of reactions that are enabled by this reaction
- `initial::Int`: Initial gene/RNA state
- `final::Int`: Final gene/RNA state
"""
struct Reaction
    action::Int
    index::Int
    disabled::Vector{Int64}
    enabled::Vector{Int64}
    initial::Int
    final::Int
end
"""
    ReactionIndices

Structure containing rate indices for different types of reactions.

# Fields
- `grange::UnitRange{Int64}`: Range of indices for gene state transitions
- `irange::UnitRange{Int64}`: Range of indices for initiation reactions
- `rrange::UnitRange{Int64}`: Range of indices for RNA elongation reactions
- `erange::UnitRange{Int64}`: Range of indices for ejection reactions
- `srange::UnitRange{Int64}`: Range of indices for splicing reactions
- `decay::Int`: Index for mRNA decay reaction
"""
struct ReactionIndices
    grange::UnitRange{Int64}
    irange::UnitRange{Int64}
    rrange::UnitRange{Int64}
    erange::UnitRange{Int64}
    srange::UnitRange{Int64}
    decay::Int
end
"""
    set_actions()

Create dictionary mapping reaction names to their action codes.

# Returns
- `Dict{String, Int}`: Dictionary with keys as reaction names and values as action codes:
  - "activateG!" => 1: Gene activation
  - "deactivateG!" => 2: Gene deactivation  
  - "transitionG!" => 3: Gene state transition
  - "initiate!" => 4: Transcription initiation
  - "transitionR!" => 5: RNA elongation
  - "eject!" => 6: RNA ejection
  - "splice!" => 7: RNA splicing
  - "decay!" => 8: mRNA decay
"""
set_actions() = Dict("activateG!" => 1, "deactivateG!" => 2, "transitionG!" => 3, "initiate!" => 4, "transitionR!" => 5, "eject!" => 6, "splice!" => 7, "decay!" => 8)
# invert_dict(D) = Dict(D[k] => k for k in keys(D)) put into utilities

"""
    simulator(r, transitions, G, R, S, insertstep; coupling=tuple(), nalleles=1, nhist=20, onstates=Int[], bins=Float64[], traceinterval::Float64=0.0, probfn=prob_Gaussian, noiseparams=4, totalsteps::Int=10000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", verbose::Bool=false)

Simulate any GRSM model. Returns steady state mRNA histogram. If bins not a null vector will return a vector of the mRNA histogram and ON and OFF time histograms. If traceinterval > 0, it will return a vector containing the mRNA histogram and the traces

#Arguments
- `r`: vector of rates
- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
- `G`: number of gene states
- `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
- `S`: number of splice sites (set to 0 for G (classic telegraph) and GR models and R for GRS models)
- `insertstep`: reporter insertion step

#Named arguments
- `bins::Vector=Float64[]`: vector of time bin vectors for each set of ON and OFF histograms or vector of vectors of time bins (one time bin vector for each onstate)
- `coupling=tuple()`: if nonempty, a 4-tuple where elements are
    1. tuple of model indices corresponding to each unit, e.g. (1, 1, 2) means that unit 1 and 2 use model 1 and unit 3 uses model 2
    2. tuple of vectors indicating source units for each unit, e.g. ([2,3], [1], Int[]) means unit 1 is influenced by source units 2 and 3, unit 2 is influenced by unit 1 and unit 3 is uninfluenced.
    3. source states, e.g. (3,0) means that model 1 influences other units whenever it is in G state 3, while model 2 does not influence any other unit
    4. target transitions, e.g. (0, 4) means that model 1 is not influenced by any source while model 2 is influenced by sources at G transition 4.
    5. Int indicating number of coupling parameters
- `nalleles`: Number of alleles
- `nhist::Int`: Size of mRNA histogram
- `onstates::Vector`: a vector of vector of ON states (use empty set for any R step is ON), ON and OFF time distributions are computed for each ON state set
- `probfn`=prob_Gaussian: reporter distribution
- `reporterfn`=sum: how individual reporters are combined
- `splicetype`::String: splice action
- `tol`::Float64=1e-6: convergence error tolerance for mRNA histogram (not used when simulating traces are made)
- `totalsteps`::Int=10000000: maximum number of simulation steps (not usred when simulating traces)
- `totaltime`::Float64=0.0: total time of simulation
- `traceinterval`: Interval in minutes between frames for intensity traces.  If 0, traces are not made.
- `verbose::Bool=false`: flag for printing state information

#Example:

julia> h=simulator([.1, .1, .1, .1, .1, .1, .1, .1, .1, .01],([1,2],[2,1],[2,3],[3,2]),3,2,2,1,nhist=20,bins=[collect(2.:2.:200),collect(.2:.2:20)],onstates=[Int[],[3]],nalleles=2, totalsteps = 200000)
5-element Vector{Vector}:
 [7823.967526508377, 33289.19787176562, 69902.6774554014, 92942.59127561412, 91816.91189325438, 70259.88319069796, 43895.28637479579, 22426.725922619895, 9005.190755732247, 3043.2417332890695, 1005.6412773072143, 203.11430396725336, 6.420639815427421, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
 [2911, 2568, 2228, 1694, 1354, 1088, 819, 661, 514, 401  …  0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 [622, 643, 592, 596, 553, 524, 520, 489, 448, 437  …  19, 13, 12, 12, 14, 10, 10, 17, 8, 8]
 [550, 584, 569, 555, 498, 510, 497, 495, 489, 487  …  89, 96, 107, 99, 89, 103, 86, 97, 87, 77]
 [593, 519, 560, 512, 492, 475, 453, 468, 383, 429  …  84, 73, 85, 92, 73, 81, 85, 101, 79, 78]

"""
function simulator(rin, transitions, G, R, S, insertstep; warmupsteps=0, coupling=tuple(), nalleles=1, nhist=20, onstates=Int[], bins=Float64[], traceinterval::Float64=0.0, probfn=prob_Gaussian, noiseparams=4, totalsteps::Int=100000000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", a_grid=nothing, verbose::Bool=false, ejectnumber=1)

    r = copy(rin)
    if !isempty(coupling)
        coupling, nalleles, noiseparams, r = prepare_coupled(r, coupling, transitions, G, R, S, insertstep, nalleles, noiseparams)
        nhist = 0
    end

    S = reset_S(S, R, insertstep)

    # if length(r) < num_rates(transitions, R, S, insertstep) + noiseparams * (traceinterval > 0)
    #     throw("r has too few elements")
    # end
    # if insertstep > R > 0
    #     throw("insertstep>R")
    # end

    nhist, mhist, mhist0, m, steps, t, ts, t0, tsample, err = initialize_sim(r, nhist, tol)
    reactions = set_reactions(transitions, G, R, S, insertstep)
    tau, state = initialize(r, G, R, reactions, nalleles)

    if length(bins) < 1
        onoff = false
    else
        onoff = true
        onstates, bins, before, after, ndt, dt, histofftdd, histontdd, tIA, tAI = set_onoff(onstates, bins, nalleles, coupling)
    end

    if traceinterval > 0
        par = set_par(r, noiseparams)
        tracelog = initialize_tracelog(t, state)
        if isempty(onstates) && !isempty(coupling)
            onstates = fill(Int64[], length(state))
        end
    end

    if verbose
        invactions = invert_dict(set_actions())
    end
    if totaltime > 0.0
        err = 0.0
        totalsteps = 0
    end

    if warmupsteps > 0
        for steps in 1:warmupsteps
            t, index, allele = findmin_tau(tau)
            initial, final, disabled, enabled, action = set_arguments(reactions, index)
            m = update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling)
        end # for
        # After warmup, tau values are absolute times (proposed next event times)
        # Shift all tau values back by the warmup time t so they remain absolute times
        # but now relative to the new time 0 (start of main simulation)
        # Note: tau[e] = -log(rand()) / r[e] + t during update!, so tau[e] > t
        # Therefore after shifting, tau[e] - t > 0 (proposed next absolute time is positive)
        # Note: for taui in tau iterates over each Matrix in Vector{Matrix}, and
        # taui .-= t modifies each matrix in place (taui is a reference, not a copy)
        for taui in tau
            taui .-= t
        end
        # Reset simulation time to match t0 initialization from initialize_sim
        t = 0.0
        t0 = 0.0
    end # if

    while (err > tol && steps < totalsteps) || t < totaltime
        steps += 1
        t, index, allele = findmin_tau(tau)   # reaction index and allele for least time transition
        initial, final, disabled, enabled, action = set_arguments(reactions, index)
        dth = t - t0
        t0 = t
        nhist > 0 && update_mhist!(mhist, m, dth, nhist)

        if t - ts > tsample && traceinterval == 0 && nhist > 0
            err, mhist0 = update_error(mhist, mhist0)
            ts = t
        end

        if onoff
            before = set_before(before, onstates, state, allele, G, R, insertstep)
        end

        if verbose
            println("---")
            println(steps)
            println("m: ", m)
            println("state: ", state)
            println("tau: ", tau)
            println("t: ", t)
            println("index: ", index, ", allele: ", allele)
            println(invactions[action])
            println("enabled: ", enabled, ", disabled: ", disabled)
            println(initial, "->", final)
        end

        m = update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling, ejectnumber, verbose)

        if onoff
            set_after!(histofftdd, histontdd, tAI, tIA, dt, ndt, before, after, t, onstates, state, allele, G, R, insertstep, verbose)
        end
        if traceinterval > 0
            set_tracelog!(tracelog, t, state)

        end
    end  # while
    verbose && println(steps)

    results = []
    nhist > 0 && push!(results, prune_mhist(mhist, nhist))
    if onoff
        if !isempty(coupling)
            for i in eachindex(histontdd)
                push!(results, vcat([[histontdd[i][j], histofftdd[i][j]] for j in eachindex(histontdd[i])]...))
            end
        else
            for i in eachindex(histontdd)
                push!(results, histontdd[i])
                push!(results, histofftdd[i])
            end
        end
    end
    if traceinterval > 0.0
        if isnothing(a_grid)
            push!(results, make_trace(tracelog, G, R, S, insertstep, onstates, traceinterval, par, probfn, reporterfn))
        else
            push!(results, make_trace(tracelog, G, R, S, insertstep, onstates, traceinterval, par, probfn, reporterfn, a_grid))
        end
        if verbose
            push!(results, tracelog)
        end
    end
    return results
end



"""
    simulator_ss(r, transitions, G, R, S, insertstep; warmupsteps=0, coupling=tuple(), nalleles=1, nhist=20, onstates=Int[], bins=Float64[], traceinterval::Float64=0.0, probfn=prob_Gaussian, noiseparams=4, totalsteps::Int=100000000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", a_grid=nothing, verbose::Bool=false)

TBW
"""
function simulator_ss(rin, transitions, G, R, S, insertstep; warmupsteps=0, coupling=tuple(), nalleles=1, nhist=20, onstates=Int[], bins=Float64[], traceinterval::Float64=0.0, probfn=prob_Gaussian, noiseparams=4, totalsteps::Int=100000000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", a_grid=nothing, verbose::Bool=false)
    r = copy(rin)
    if !isempty(coupling)
        coupling, nalleles, noiseparams, r = prepare_coupled(r, coupling, transitions, G, R, S, insertstep, nalleles, noiseparams)
        nhist = 0
    end

    S = reset_S(S, R, insertstep)

    sshist = initialize_ss(G, R, S)

    _, _, _, m, steps, t, _, t0, _, _ = initialize_sim(r, nhist, tol)
    reactions = set_reactions(transitions, G, R, S, insertstep)
    tau, state = initialize(r, G, R, reactions, nalleles)

    if verbose
        invactions = invert_dict(set_actions())
    end
    if totaltime > 0.0
        err = 0.0
        totalsteps = 0
    end

    while steps < totalsteps
        steps += 1
        t, index, allele = findmin_tau(tau)   # reaction index and allele for least time transition
        initial, final, disabled, enabled, action = set_arguments(reactions, index)
        dth = t - t0
        t0 = t

        if verbose
            println("---")
            println(steps)
            println("state: ", state, ": ", state_index(state, G, R, S))
            println("state indices: ", findall(!iszero, vec(state[index[1], 1])), ", coupling: ", coupling[3])
            println("tau: ", tau)
            println("t: ", t)
            println("index: ", index, " ", allele)
            println(invactions[action])
            println("enabled: ", enabled, ", disabled: ", disabled)
            println(initial, "->", final)
        end

        update_sshist!(sshist, state, dth, G, R, S)

        update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling, ejectnumber, verbose)

    end  # while
    # verbose && println(steps)

    return sshist
end

function initialize_ss(G, R, S)
    zeros(prod(T_dimension(G, R, S)))
end

function update_sshist!(sshist, state, dt, G, R, S)
    # println(state_index(state, G, R, S))
    sshist[state_index(state, G, R, S)] += dt
end
"""
    initialize(r, G::Int, R, reactions, nalleles, initstate=1, initreaction=1)

Initialize simulation state and reaction times for a single gene model.

# Arguments
- `r::Vector{Float64}`: Vector of reaction rates
- `G::Int`: Number of gene states
- `R`: Number of pre-RNA steps
- `reactions::Vector{Reaction}`: Vector of reaction definitions
- `nalleles::Int`: Number of alleles to simulate
- `initstate::Int=1`: Initial gene state
- `initreaction::Int=1`: Initial reaction to schedule

# Returns
- `Tuple{Matrix{Float64}, Matrix{Int}}`: (tau, states) where tau contains next reaction times and states contains current system state
"""
function initialize(r, G::Int, R, reactions, nalleles, initstate=1, initreaction=1)
    nreactions = length(reactions)
    tau = fill(Inf, nreactions, nalleles)
    states = zeros(Int, G + R, nalleles)
    for n in 1:nalleles
        tau[initreaction, n] = -log(rand()) / r[1]
        states[initstate, n] = 1
    end
    return tau, states
end
"""
    initialize(r, G::Tuple, R, reactions, nalleles, initstate=1, initreaction=1)

Initialize simulation state and reaction times for coupled gene models.

# Arguments
- `r::Vector{Vector{Float64}}`: Vector of rate vectors for each gene
- `G::Tuple`: Tuple of gene state counts for each gene
- `R::Tuple`: Tuple of pre-RNA step counts for each gene
- `reactions::Vector{Vector{Reaction}}`: Vector of reaction definitions for each gene
- `nalleles::Int`: Number of alleles to simulate
- `initstate::Int=1`: Initial gene state
- `initreaction::Int=1`: Initial reaction to schedule

# Returns
- `Tuple{Vector{Matrix{Float64}}, Vector{Matrix{Int}}}`: (tau, states) where tau contains next reaction times and states contains current system state for each gene
"""
function initialize(r, G::Tuple, R, reactions, nalleles, initstate=1, initreaction=1)
    tau = Matrix[]
    states = Matrix[]
    for i in eachindex(G)
        t, s = initialize(r[i], G[i], R[i], reactions[i], 1, initstate, initreaction)
        push!(tau, t)
        push!(states, s)
    end
    return tau, states
end


"""
    initialize_sim(r::Vector{Vector}, nhist, tol, samplefactor=20.0, errfactor=10.0)

Initialize simulation parameters for coupled models.

# Arguments
- `r::Vector{Vector{Float64}}`: Vector of rate vectors for each gene
- `nhist`: Histogram size(s) for mRNA counting
- `tol::Float64`: Convergence tolerance
- `samplefactor::Float64=20.0`: Factor for sampling interval calculation
- `errfactor::Float64=10.0`: Factor for error tolerance calculation

# Returns
- `Tuple`: (nhist, mhist, mhist0, m, steps, t, ts, t0, tsample, err) containing simulation state variables
"""
function initialize_sim(r::Vector{Vector}, nhist, tol, samplefactor=20.0, errfactor=10.0)
    if nhist isa Number && nhist > 0
        nhist = fill(nhist, length(r) - 1)
    end
    mhist = Vector[]
    mhist0 = similar(mhist)
    rmin = Inf
    for i in eachindex(nhist)
        rmin = min(rmin, minimum(r[i]))
        push!(mhist, zeros(nhist[i] + 1))
        push!(mhist0, ones(nhist[1] + 1))
    end
    return nhist, mhist, mhist0, zeros(Int, length(r) - 1), 0, 0.0, 0.0, 0.0, samplefactor / rmin, errfactor * tol
end

"""
    initialize_sim(r::Vector{Float64}, nhist, tol, samplefactor=20.0, errfactor=10.0)

Initialize simulation parameters for single gene models.

# Arguments
- `r::Vector{Float64}`: Vector of reaction rates
- `nhist::Int`: Histogram size for mRNA counting
- `tol::Float64`: Convergence tolerance
- `samplefactor::Float64=20.0`: Factor for sampling interval calculation
- `errfactor::Float64=10.0`: Factor for error tolerance calculation

# Returns
- `Tuple`: (nhist, mhist, mhist0, m, steps, t, ts, t0, tsample, err) containing simulation state variables
"""
initialize_sim(r::Vector{Float64}, nhist, tol, samplefactor=20.0, errfactor=10.0) = nhist, zeros(nhist + 1), ones(nhist + 1), 0, 0, 0.0, 0.0, 0.0, samplefactor / minimum(r), errfactor * tol

"""
    set_par(r, noiseparams::Int)
    set_par(r, noiseparams::Vector)

Extract noise parameters from rate vector.

# Arguments
- `r::Vector{Float64}`: Vector of rates including noise parameters
- `noiseparams::Int`: Number of noise parameters at the end of the rate vector

# Returns
- `Vector{Float64}`: Vector containing only the noise parameters

# Arguments (Vector version)
- `r::Vector{Vector{Float64}}`: Vector of rate vectors
- `noiseparams::Vector{Int}`: Number of noise parameters for each rate vector

# Returns (Vector version)
- `Vector{Vector{Float64}}`: Vector of noise parameter vectors
"""
set_par(r, noiseparams::Int) = r[end-noiseparams+1:end]

function set_par(r, noiseparams::Vector)
    par = Vector[]
    for i in eachindex(noiseparams)
        push!(par, set_par(r[i], noiseparams[i]))
    end
    par
end


"""
    targets(coupling)

Find coupled targets for each unit in a coupled model.

# Arguments
- `coupling::Tuple`: Coupling tuple containing (models, sources, source_states, target_transitions, n_coupling_params)

# Returns
- `Vector{Vector{Int}}`: Vector where targets[i] contains indices of units that are influenced by unit i
"""
function targets(coupling)
    targets = Vector{Int}[]
    models = coupling[1]
    sources = coupling[2]
    for m in models
        t = Int[]
        for i in eachindex(sources)
            if m ∈ sources[i]
                push!(t, i)
            end
        end
        push!(targets, t)
    end
    targets
end

"""
    prepare_coupled(r, coupling, transitions, G, R, S, insertstep, nalleles, noiseparams)

Prepare parameters for coupled model simulation.

# Arguments
- `r::Vector{Float64}`: Vector of rates including coupling parameters
- `coupling::Tuple`: Coupling configuration tuple
- `transitions::Tuple`: Transition definitions for each gene
- `G::Tuple`: Number of gene states for each gene
- `R::Tuple`: Number of pre-RNA steps for each gene
- `S::Tuple`: Number of splice sites for each gene
- `insertstep::Tuple`: Reporter insertion steps for each gene
- `nalleles::Int`: Number of alleles
- `noiseparams::Union{Int, Vector{Int}}`: Number of noise parameters

# Returns
- `Tuple`: (coupling, nalleles, noiseparams, r_prepared) where r_prepared contains separated rate vectors for each gene
"""
function prepare_coupled(r, coupling, transitions, G, R, S, insertstep, nalleles, noiseparams)
    ncoupling = coupling[5]
    if length(r) >= ncoupling && any(r[end-ncoupling+1:end] .< -1.0)
        throw(ArgumentError("all coupling strengths must be > -1.0"))
    end
    if noiseparams isa Number
        noiseparams = fill(noiseparams, length(G))
    end
    coupling = (coupling..., targets(coupling))
    return coupling, nalleles, noiseparams, prepare_rates_sim(r, coupling, transitions, R, S, insertstep, noiseparams)
end


"""

    prepare_rates_sim(rates, coupling, transitions, R, S, insertstep, n_noise)

Prepare rate vectors for coupled model simulation by separating rates for each gene.

# Arguments
- `rates::Vector{Float64}`: Combined vector of rates for all genes
- `coupling::Tuple`: Coupling configuration tuple
- `transitions::Tuple`: Transition definitions for each gene
- `R::Tuple`: Number of pre-RNA steps for each gene
- `S::Tuple`: Number of splice sites for each gene
- `insertstep::Tuple`: Reporter insertion steps for each gene
- `n_noise::Vector{Int}`: Number of noise parameters for each gene

# Returns
- `Vector{Vector{Float64}}`: Vector of rate vectors, one for each gene, plus coupling parameters
"""
function prepare_rates_sim(rates, coupling, transitions, R, S, insertstep, n_noise)
    r = Vector[]
    couplingStrength = Float64[]
    j = 1
    for i in eachindex(R)
        n = num_rates(transitions[i], R[i], S[i], insertstep[i]) + n_noise[i]
        push!(r, rates[j:j+n-1])
        j += n
    end
    # Canonical order: to_connections(coupling), one γ per connection
    conns = to_connections(coupling)
    for k in 1:length(conns)
        push!(couplingStrength, rates[j])
        j += 1
    end
    push!(r, couplingStrength)
    r
end

"""
    simulate_trace(r, transitions, G, R, S, insertstep, interval, onstates; totaltime=1000.0, reporterfn=sum)

Simulate a single intensity trace for a GRSM model.

# Arguments
- `r::Vector{Float64}`: Vector of reaction rates
- `transitions`: Transition definitions
- `G`: Number of gene states
- `R`: Number of pre-RNA steps
- `S`: Number of splice sites
- `insertstep`: Reporter insertion step
- `interval::Float64`: Time interval between frames
- `onstates::Vector{Int}`: Vector of ON states
- `totaltime::Float64=1000.0`: Total simulation time
- `reporterfn::Function=sum`: Function to combine reporter signals

# Returns
- `Matrix{Float64}`: Trace matrix with columns [time, intensity, reporters, state]
"""
simulate_trace(r, transitions, G, R, S, insertstep, interval, onstates; totaltime=1000.0, reporterfn=sum) = simulator(r, transitions, G, R, S, insertstep, nalleles=1, nhist=0, onstates=onstates, traceinterval=interval, reporterfn=reporterfn, totaltime=totaltime, warmupsteps=100)[1][1:end-1, :]

"""
    simulate_trace_data(datafolder::String; ntrials::Int=10, r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231, 30, 20, 200, 100, 0.8], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=Int[], interval=1.0, totaltime=1000.0)

Create simulated trace files in the specified data folder.

# Arguments
- `datafolder::String`: Path to folder where trace files will be saved
- `ntrials::Int=10`: Number of simulation trials
- `r::Vector{Float64}`: Vector of reaction rates
- `transitions`: Transition definitions
- `G::Int=3`: Number of gene states
- `R::Int=2`: Number of pre-RNA steps
- `S::Int=2`: Number of splice sites
- `insertstep::Int=1`: Reporter insertion step
- `onstates::Vector{Int}=Int[]`: Vector of ON states
- `interval::Float64=1.0`: Time interval between frames
- `totaltime::Float64=1000.0`: Total simulation time
- `reporterfn::Function=sum`: Function to combine reporter signals
- `a_grid`: Optional grid transition matrix

# Returns
- Creates `.trk` files in the specified folder containing simulated trace data
"""
function simulate_trace_data(datafolder::String; a_grid=nothing, ntrials::Int=10, r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231, 40, 20, 200, 10], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=Int[], interval=1.0, totaltime=1000.0, reporterfn=sum)
    isdir(datafolder) || mkdir(datafolder)
    trace = simulate_trace_vector(r, transitions, G, R, S, insertstep, interval, totaltime, ntrials, onstates=onstates, reporterfn=reporterfn, a_grid=a_grid)
    for i in eachindex(trace)
        l = length(trace[i])
        if isnothing(a_grid)
            writedlm(joinpath(datafolder, "testtrace$i.trk"), [zeros(l) zeros(l) trace[i] collect(1:l)])
        else
            writedlm(joinpath(datafolder, "testtrace$i.trk"), trace[i][:, 1:end-1])
        end
    end
end
"""
    simulate_trace_vector(r, transitions, G, R, S, interval, totaltime, ntrials; insertstep=1, onstates=Int[], reporterfn=sum)

- `hierarchical`: tuple of (background mean index, vector of background means)
"""
function simulate_trace_vector(rin, transitions, G::Int, R, S, insertstep, interval, totaltime, ntrials; onstates=Int[], reporterfn=sum, a_grid=nothing, hierarchical=tuple(), warmupsteps=100, col=2)
    trace = Vector{Vector{Float64}}(undef, ntrials)
    r = copy(rin)
    for i in eachindex(trace)
        if !isempty(hierarchical)
            r[hierarchical[1]] = hierarchical[2][i]
        end
        if isnothing(a_grid)
            trace[i] = simulator(r, transitions, G, R, S, insertstep, onstates=onstates, traceinterval=interval, totaltime=totaltime, nhist=0, reporterfn=reporterfn, a_grid=a_grid, warmupsteps=warmupsteps)[1][1:end-1, col]
        else
            trace[i] = simulator(r, transitions, G, R, S, insertstep, onstates=onstates, traceinterval=interval, totaltime=totaltime, nhist=0, reporterfn=reporterfn, a_grid=a_grid, warmupsteps=warmupsteps)[1]
        end
    end
    trace
end


"""
    simulate_trace_vector(r, transitions, G::Tuple, R, S, insertstep, coupling::Tuple, interval, totaltime, ntrials; onstates=Int[], reporterfn=sum)

TBW
"""
# function simulate_trace_vector(r, transitions, G::Tuple, R, S, insertstep, coupling::Tuple, interval, totaltime, ntrials; onstates=Int[], reporterfn=sum, a_grid=nothing, hierarchical=tuple(), col=2)
function simulate_trace_vector(rin, transitions, G::Tuple, R, S, insertstep, coupling::Tuple, interval, totaltime, ntrials; onstates=Int[], reporterfn=sum, a_grid=nothing, hierarchical=tuple(), col=2, totalsteps=100000, verbose=false, warmupsteps::Int=100)
    trace = Array{Array{Float64}}(undef, ntrials)
    r = copy(rin)
    if verbose
        totaltime = 0.0
        totalsteps = 10
    end
    # Convert warmup_time to warmup steps (events). For coupled models, we want sufficient warmup.
    # Use conservative estimate: if average event rate is ~1 per time unit, use warmup_time steps.
    # But ensure minimum of 1000 steps for stability.
    # warmup_steps = warmup_time > 0.0 ? max(Int(round(warmup_time)), 1000) : 0
    for i in eachindex(trace)
        if !isempty(hierarchical)
            r[hierarchical[1]] = hierarchical[2][i]
        end
        t = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, onstates=onstates, traceinterval=interval, totaltime=totaltime, totalsteps=totalsteps, nhist=0, reporterfn=reporterfn, warmupsteps=warmupsteps, verbose=verbose)[1]
        if col isa Vector
            # When col is a vector, extract multiple columns and combine them
            tr = Vector[]
            for t in t
                extracted = t[1:end-1, col]  # This is a Matrix when col is a vector
                # Convert each column to a vector and push them
                for j in 1:size(extracted, 2)
                    push!(tr, extracted[:, j])
                end
            end
            trace[i] = hcat(tr...)
        else
            # When col is a scalar, extract single column
            tr = Vector[]
            for t in t
                push!(tr, t[1:end-1, col])
            end
            trace[i] = hcat(tr...)
        end
    end
    trace
end


"""
    make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Int}, interval, par, probfn, reporterfn)

TBW
"""
function make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Int}, interval, par, probfn, reporterfn)
    n = length(tracelog)
    trace = Matrix{Float64}(undef, 0, 4)
    state = tracelog[1][2]
    frame = interval
    reporters = num_reporters_per_state(G, R, S, insertstep, onstates, reporterfn)
    i = 2
    base = S > 0 ? 3 : 2
    d = probfn(par, reporters, G * base^R)
    while i < n
        while tracelog[i][1] <= frame && i < n
            state = tracelog[i][2]
            i += 1
        end
        trace = vcat(trace, [frame intensity(state, G, R, S, d) reporters[state_index(state, G, R, S)] state_index(state, G, R, S)])
        frame += interval
    end
    return trace
end

"""
    make_trace(tracelog, G::Tuple, R, S, insertstep, onstates, interval, par, probfn, reporterfn=sum)

Return array of frame times and intensities

- `tracelog`: Vector if Tuples of (time,state of allele 1)
- `interval`: Number of minutes between frames
- `onstates`: Vector of G on states, empty for GRS models
- `G` and `R` as defined in simulator

"""
function make_trace(tracelog, G::Tuple, R, S, insertstep, onstates, interval, par, probfn, reporterfn=sum)
    trace = Matrix[]
    for i in eachindex(tracelog)
        push!(trace, make_trace(tracelog[i], G[i], R[i], S[i], insertstep[i], onstates[i], interval, par[i], probfn, reporterfn))
    end
    trace
end
"""
    make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Vector}, interval, par, probfn, reporterfn=sum)

TBW
"""
function make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Vector{Int}}, interval, par, probfn, reporterfn)
    traces = Matrix[]
    for o in onstates
        push!(traces, make_trace(tracelog, G, R, S, insertstep, o, interval, par, probfn, reporterfn))
    end
    traces
end



function make_trace_grid(trace::Vector{T}, a_grid, d_background) where {T<:Array}
    gridtrace = Matrix[]
    for t in trace
        push!(gridtrace, make_trace_grid(t, a_grid, d_background))
    end
    gridtrace
end

function make_trace_grid(trace::Matrix, a_grid, d_background)
    make_trace_grid(trace[:, 2], a_grid, d_background)
end

function make_trace_grid(trace::Vector{Float64}, a_grid, d_background)
    Nstate = size(a_grid, 1)
    tracematrix = Matrix{Float64}(undef, Nstate, length(trace))
    cdf = cumsum(a_grid, dims=2)
    position = rand(1:Nstate)
    for t in eachindex(trace)
        position = searchsortedfirst(cdf[position, :], rand())
        tracematrix[position, t] = trace[t]
        for p in 1:Nstate
            if p != position
                tracematrix[p, t] = max(rand(d_background), 0.0)
            end
        end
    end
    tracematrix
end

function make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Int}, interval, par, probfn, reporterfn, a_grid)
    trace = make_trace(tracelog, G, R, S, insertstep, onstates::Vector{Int}, interval, par, probfn, reporterfn)
    make_trace_grid(trace, a_grid, probfn(par, 0))
end
"""
    intensity(state, G, R, S, d)

Returns the trace intensity given the state of a system

For R = 0, the intensity is occupancy of any onstates
For R > 0, intensity is the number of reporters in the nascent mRNA

"""
function intensity(state, G, R, S, d)
    stateindex = state_index(state, G, R, S)
    # max(rand(d[stateindex]), 0)
    rand(d[stateindex])
end

"""
    gstate(G, state, allele)

TBW
"""
gstate(G, state, allele) = argmax(state[1:G, allele])

"""
    state_index(state::Array, G, allele)
    state_index(state::Array, G, R, S,allele=1)

returns state index given state vector
"""
state_index(state::Array, G, allele) = argmax(state[1:G, allele])

"""
    state_index(state::Array, G, R, S, allele=1)

TBW
"""
function state_index(state::Array, G::Int, R, S, allele=1)
    Gstate = gstate(G, state, allele)
    if R == 0
        return Gstate
    else
        if S > 0
            base = 3
            Rstate = state[G+1:end, allele]
        else
            base = 2
            Rstate = Int.(state[G+1:end, allele] .> 0)
        end
        return Gstate + G * decimal(Rstate, base)
    end
end

function state_index(state, G::Tuple, R::Tuple, S::Tuple, allele=1)
    si = Vector{Int}(undef, 2)
    for i in eachindex(G)
        si[i] = state_index(state[i, 1], G[i], R[i], S[i], allele)
    end
    (si[1] - 1) * T_dimension(G[2], R[2], S[2]) + si[2]
end

function coupled_state_index(jointstate::Vector, G, R, S)
    (jointstate[1] - 1) * T_dimension(G[2], R[2], S[2]) + jointstate[2]
end

"""
    initialize_tracelog(t, state::Vector{Matrix})


TBW
"""
function initialize_tracelog(t, state::Vector{Matrix})
    tracelog = Vector{Vector}(undef, length(state))
    for i in eachindex(state)
        tracelog[i] = initialize_tracelog(t, state[i])
    end
    tracelog
end
"""
    initialize_tracelog(t, state::Matrix)

TBW
"""
initialize_tracelog(t, state::Matrix) = [(t, state[:, 1])]


"""
    set_tracelog!(tracelog, t, state::Matrix)

TBW
"""
function set_tracelog!(tracelog, t, state::Matrix)
    push!(tracelog, (t, state[:, 1]))
end
"""
    set_tracelog(tracelog, t, state::Vector{Matrix})

TBW
"""
function set_tracelog!(tracelog, t, state::Vector{Matrix})
    for i in eachindex(state)
        set_tracelog!(tracelog[i], t, state[i])
    end
end





"""
    findmin_tau(tau::Vector{Matrix})

TBW
"""
function findmin_tau(tau::Vector{Matrix})
    τ = Inf
    index = 1
    allele = 1
    for i in eachindex(tau)
        t, ind, a = findmin_tau(tau[i])
        if t < τ
            τ = t
            index = (i, ind)
            # allele = a
        end
    end
    τ, index, allele
end

function findmin_tau(tau::Matrix)
    t, ind = findmin(tau)
    return t, ind[1], ind[2]
end




"""
    set_arguments(reaction, index::Tuple)

return values of fields of Reaction structure
"""
set_arguments(reaction, index::Tuple) = set_arguments(reaction[index[1]], index[2])

"""
    set_arguments(reaction, index::Int)

TBW
"""
set_arguments(reaction, index::Int) = (reaction[index].initial, reaction[index].final, reaction[index].disabled, reaction[index].enabled, reaction[index].action)


"""
    set_reactionindices(Gtransitions, R::Int, S)

return structure of ranges for each type of transition
"""
function set_reactionindices(Gtransitions, R::Int, S)
    nG = length(Gtransitions)
    g = 1:nG
    i = nG+1:nG+Int(R > 0)
    r = nG+2:nG+R
    e = nG+R+1:nG+R+1
    s = nG+1+R+1:nG+1+R+S
    d = nG + 1 + R + S + 1
    ReactionIndices(g, i, r, e, s, d)
end
"""
    set_reactions(Gtransitions, G::Tuple, R, S, insertstep)

return vector of Reaction structures for each transition in GRS model

reporter first appears at insertstep
"""
function set_reactions(Gtransitions, G::Tuple, R, S, insertstep)
    reactions = Vector{Reaction}[]
    for i in eachindex(G)
        push!(reactions, set_reactions(Gtransitions[i], G[i], R[i], S[i], insertstep[i]))
    end
    reactions
end
"""
    set_reactions(Gtransitions, G::Int, R, S, insertstep)

TBW
"""
function set_reactions(Gtransitions, G::Int, R, S, insertstep)
    actions = set_actions()
    indices = set_reactionindices(Gtransitions, R, S)
    reactions = Reaction[]
    nG = length(Gtransitions)
    Sstride = R - insertstep + 1
    for g in eachindex(Gtransitions)
        u = Int[] #disabled (upstream)
        d = Int[] #enabled (downstream)
        ginitial = Gtransitions[g][1]
        gfinal = Gtransitions[g][2]
        for s in eachindex(Gtransitions)
            # if ginitial == Gtransitions[s][1] && gfinal != Gtransitions[s][2]
            if ginitial == Gtransitions[s][1]
                push!(u, s)
            end
            if gfinal == Gtransitions[s][1]
                push!(d, s)
            end
        end
        if gfinal == G
            push!(d, nG + 1)
            push!(reactions, Reaction(actions["activateG!"], g, u, d, ginitial, gfinal))
        elseif ginitial == G
            push!(u, nG + 1)
            push!(reactions, Reaction(actions["deactivateG!"], g, u, d, ginitial, gfinal))
        else
            push!(reactions, Reaction(actions["transitionG!"], g, u, d, ginitial, gfinal))
        end
    end
    for i in indices.irange
        if S > 0 && insertstep == 1
            push!(reactions, Reaction(actions["initiate!"], i, [], [nG + 2; nG + 2 + R], G, G + 1))
        else
            push!(reactions, Reaction(actions["initiate!"], i, [], [nG + 2], G, G + 1))
        end
    end
    i = G
    for r in indices.rrange  # Corrected for allowing S < R + insertstep - 1
        i += 1
        if S == 0 || i > G + S
            push!(reactions, Reaction(actions["transitionR!"], r, Int[r], [r - 1; r + 1], i, i + 1))
        else
            if i < G + insertstep - 1
                push!(reactions, Reaction(actions["transitionR!"], r, Int[r], [r - 1; r + 1], i, i + 1))
            elseif i == G + insertstep - 1
                push!(reactions, Reaction(actions["transitionR!"], r, Int[r], [r - 1; r + 1; r + 1 + Sstride], i, i + 1))
            elseif i > G + insertstep - 1
                push!(reactions, Reaction(actions["transitionR!"], r, Int[r; r + Sstride], [r - 1; r + 1; r + 1 + Sstride], i, i + 1))
            end
        end
    end
    for e in indices.erange
        if S > 0 && e + Sstride <= nG + R + S + 1 # Corrected for allowing S < R + insertstep - 1
            push!(reactions, Reaction(actions["eject!"], e, Int[e, e+Sstride], Int[e-1, indices.decay], G + R, 0))
        elseif R > 0
            push!(reactions, Reaction(actions["eject!"], e, Int[e], Int[e-1, indices.decay], G + R, 0))
        else
            push!(reactions, Reaction(actions["eject!"], e, Int[], Int[e, indices.decay], G + 1, 0))
        end
    end
    j = G + insertstep - 1
    for s in indices.srange
        j += 1
        push!(reactions, Reaction(actions["splice!"], s, Int[], Int[], j, j))
    end
    push!(reactions, Reaction(actions["decay!"], indices.decay, Int[], Int[indices.decay], 0, 0))
    return reactions
end
"""
    update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling, ejectnumber=1, verbose=false)

updates proposed next reaction time and state given the selected action and returns updated number of mRNA

(uses if-then statements because that executes faster than an element of an array of functions)

Arguments are same as defined in simulator

"""
function update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling, ejectnumber=1, verbose=false)
    !isempty(coupling) && (initialstate = copy(state[index[1], 1]))
    if action < 5
        if action < 3
            if action == 1
                activateG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
            else
                deactivateG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
            end
        else
            if action == 3
                transitionG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
            else
                initiate!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
            end
        end
    else
        if action < 7
            if action == 5
                transitionR!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
            else
                m = eject!(tau, state, index, t, m, r, allele, G, R, S, disabled, enabled, initial, ejectnumber)
            end
        else
            if action == 7
                splice!(tau, state, index, t, m, r, allele, G, R, initial)
            else
                m = decay!(tau, index, t, m, r)
            end
        end
    end
    # println("taup: ",tau)
    !isempty(coupling) && update_coupling!(tau, state, index[1], t, r, enabled, initialstate, coupling, verbose)
    return m
end

"""
    coupling!(tau, state, unit::Int, t, r, disabled, enabled, initial, final, coupling)

TBW
"""

function _state_in_s(state, s)
    s isa Int && return s in state
    return !isdisjoint(state, s)
end

function update_coupling!(tau, state, unit::Int, t, r, enabled, initialstate, coupling, verbose=false)
    conns = to_connections(coupling)
    coupling_strength = r[end]  # Vector of length ncoupling, same order as conns
    oldstate = findall(!iszero, vec(initialstate))
    newstate = findall(!iszero, vec(state[unit, 1]))

    verbose && println("unit: ", unit, ", oldstate: ", oldstate, ", newstate: ", newstate, ", conns: ", conns)
    verbose && println("tau1: ", tau)

    # unit as source (β): for each connection (unit, α, s, t_trans), update tau[α][t_trans]
    for k in 1:length(conns)
        β, α, s, t_trans = conns[k]
        β != unit && continue
        if isfinite(tau[α][t_trans, 1])
            γc = coupling_strength[k]
            s_set = s isa Int ? [s] : s
            if isdisjoint(oldstate, s_set) && _state_in_s(newstate, s)
                tau[α][t_trans, 1] = 1 / (1 + γc) * (tau[α][t_trans, 1] - t) + t
            elseif _state_in_s(oldstate, s) && isdisjoint(newstate, s_set)
                tau[α][t_trans, 1] = (1 + γc) * (tau[α][t_trans, 1] - t) + t
            end
        end
    end

    # unit as target (α): for each connection (β, unit, s, t_trans), update tau[unit][t_trans] when source β in state s
    for k in 1:length(conns)
        β, α, s, t_trans = conns[k]
        α != unit && continue
        t_trans ∈ enabled || continue
        source_state_vec = findall(!iszero, vec(state[β, 1]))
        _state_in_s(source_state_vec, s) || continue
        γc = coupling_strength[k]
        tau[unit][t_trans, 1] = 1 / (1 + γc) * (tau[unit][t_trans, 1] - t) + t
    end
    verbose && println("tau2: ", tau)
end

"""
    transitionG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)

update tau and state for G transition

"""
function transitionG!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, disabled, enabled, initial, final, coupling)
    transitionG!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], 1, G[index[1]], R[index[1]], disabled, enabled, initial, final, coupling)
    # update_coupling!(tau, state, index[1], t, r, disabled, enabled, initial, final, coupling)
end
"""
    transitionG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())

TBW
"""
function transitionG!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, disabled, enabled, initial, final, coupling=tuple())
    # println(enabled)
    for e in enabled
        tau[e, allele] = -log(rand()) / r[e] + t
    end
    for d in disabled
        tau[d, allele] = Inf
    end
    state[final, allele] = 1
    state[initial, allele] = 0
    # println("tG: ",tau)
end
"""
    activateG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)

"""
function activateG!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, disabled, enabled, initial, final, coupling)

    # println("before mutation:\n", tau)
    # taui = tau[index[1]]
    activateG!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], disabled, enabled, initial, final, coupling)
    # println("after mutation:\n", tau)
    # tau[index[1]] .= taui
    # update_coupling!(tau, state, index[1], t, r, disabled, enabled, initial, final, coupling)
end
"""
    activateG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())

TBW
"""
function activateG!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, disabled, enabled, initial, final, coupling=tuple())
    transitionG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
    if R > 0 && state[G+1, allele] > 0
        tau[enabled[end], allele] = Inf
    end
end
"""
    deactivateG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)

"""
function deactivateG!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, disabled, enabled, initial, final, coupling)
    deactivateG!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], 1, G[index[1]], R[index[1]], disabled, enabled, initial, final, coupling)
    # update_coupling!(tau, state, index[1], t, r, disabled, enabled, initial, final, coupling)
end
"""
    deactivateG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())

TBW
"""
function deactivateG!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, disabled, enabled, initial, final, coupling=tuple())
    transitionG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
end
"""
    initiate!(tau, state, index::Tuple, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
    initiate!(tau, state, index::Int, t, m, r, allele, G, R, S, disabled, enabled, initial, final,nsertstep)

"""
function initiate!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, S, insertstep, disabled, enabled, initial, final)
    initiate!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], 1, G[index[1]], R[index[1]], S[index[1]], insertstep[index[1]], disabled, enabled, initial, final)
end
"""
    initiate!(tau, state, index::Int, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)

TBW
"""
function initiate!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, S, insertstep, disabled, enabled, initial, final)
    if final + 1 > G + R || state[final+1, allele] == 0
        tau[enabled[1], allele] = -log(rand()) / (r[enabled[1]]) + t
    end
    if insertstep == 1
        state[final, allele] = 2
        if S > 0
            tau[enabled[2], allele] = -log(rand()) / (r[enabled[2]]) + t
        end
    else
        state[final, allele] = 1
    end
    tau[index, allele] = Inf
end
"""
    transitionR!(tau, state, index::Tuple, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)

"""
function transitionR!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, S, insertstep, disabled, enabled, initial, final)
    transitionR!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], S[index[1]], insertstep[index[1]], disabled, enabled, initial, final)
end
"""
    transitionR!(tau, state, index::Int, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)

TBW
"""
function transitionR!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, S, insertstep, disabled, enabled, initial, final)
    if state[initial-1, allele] > 0
        tau[enabled[1], allele] = -log(rand()) / r[enabled[1]] + t
    end
    if final + 1 > G + R || state[final+1, allele] == 0
        tau[enabled[2], allele] = -log(rand()) / r[enabled[2]] + t
    end
    if S > 0 && final >= G + insertstep && final <= G + S # Corrected for allowing S < R + insertstep - 1
        tau[enabled[3], allele] = -log(rand()) / r[enabled[3]] + t
        # if final == insertstep + G
        #     tau[enabled[3], allele] = -log(rand()) / r[enabled[3]] + t
        # elseif state[initial, allele] > 1
        #     tau[enabled[3], allele] = r[enabled[3]-1] / r[enabled[3]] * (tau[enabled[3]-1, allele] - t) + t
        # end
    end
    for d in disabled
        tau[d, allele] = Inf
    end
    if final == insertstep + G
        state[final, allele] = 2
    else
        state[final, allele] = state[initial, allele]
    end
    state[initial, allele] = 0
end

"""
    eject!(tau, state, index::Tuple, t, m, r, allele, G, R, S, disabled, enabled, initial)


"""
function eject!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, S, disabled, enabled, initial, ejectnumber=1)
    m[index[1]] = eject!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], S[index[1]], disabled, enabled, initial, ejectnumber)
    m
end
"""
    eject!(tau, state, index::Int, t, m, r, allele, G, R, S, disabled, enabled, initial)

TBW
"""
function eject!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, S, disabled, enabled, initial, ejectnumber=1)
    if state[initial-1, allele] > 0
        tau[enabled[1], allele] = -log(rand()) / (r[enabled[1]]) + t
    end
    for d in disabled
        tau[d, allele] = Inf
    end
    if R > 0
        state[initial, allele] = 0
    end
    set_decay!(tau, enabled[end], t, m, r, ejectnumber)
end
"""
    splice!(tau, state, index::Tuple, t, m, r, allele, G, R, initial)

"""
function splice!(tau::Vector, state, index::Tuple, t, m, r, allele, G::Tuple, R::Tuple, initial)
    splice!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], initial)
end
"""
    splice!(tau, state, index::Int, t, m, r, allele, G, R, initial)

TBW
"""
function splice!(tau::Matrix, state, index::Int, t, m, r, allele, G::Int, R::Int, initial)
    state[initial, allele] = 1
    tau[index, allele] = Inf
end
"""
    decay!(tau, index::Tuple, t, m, r)
    decay!(tau, index::Int, t, m, r)

"""
function decay!(tau::Vector, index::Tuple, t, m, r)
    m[index[1]] = decay!(tau[index[1]], index[2], t, m[index[1]], r[index[1]])
    m
end
"""
    decay!(tau, index::Int, t, m, r)

TBW
"""
function decay!(tau::Matrix, index::Int, t, m, r)
    m -= 1
    if m == 0
        tau[index, 1] = Inf
    else
        tau[index, 1] = -log(rand()) / (m * r[index]) + t
    end
    m
end

"""
    set_decay!(tau, index::Tuple, t, m, r)
    set_decay!(tau, index::Int, t, m, r)

update tau matrix for decay rate

"""
function set_decay!(tau::Vector, index::Tuple, t, m, r, ejectnumber=1)
    m[index[1]] = set_decay!(tau[index[1]], index[2], t, m[index[1]], r[index[1]], ejectnumber)
    m
end
"""
    set_decay!(tau, index::Int, t, m, r)

TBW
"""
function set_decay!(tau::Matrix, index::Int, t, m, r, ejectnumber=1)
    m += ejectnumber
    tau[index, 1] = -log(rand()) / (m * r[index]) + t
    # if m == 1
    #     tau[index, 1] = -log(rand()) / r[index] + t
    # else
    #     tau[index, 1] = (m - ejectnumber) / m * (tau[index, 1] - t) + t
    # end
    m
end


"""
    update_error(mhist::Vector{Vector}, mhist0)

TBW
"""
function update_error(mhist::Vector{Vector}, mhist0)
    n = 0
    for i in eachindex(mhist)
        n = max(n, norm(mhist[i] / sum(mhist[i]) - mhist0[i] / sum(mhist0[i]), Inf))
    end
    return n, copy(mhist)
end


update_error(mhist, mhist0) = (norm(mhist / sum(mhist) - mhist0 / sum(mhist0), Inf), copy(mhist))
"""
    update_mhist!(mhist, m::Vector, dt, nhist)
    update_mhist!(mhist, m::Int, dt, nhist)

"""
function update_mhist!(mhist, m::Vector, dt, nhist)
    for i in eachindex(m)
        update_mhist!(mhist[i], m[i], dt, nhist[i])
    end
end

"""
    update_mhist!(mhist, m::Int, dt, nhist)

TBW
"""
function update_mhist!(mhist, m::Int, dt, nhist)
    if m + 1 <= nhist
        mhist[m+1] += dt
    else
        mhist[nhist+1] += dt
    end
end


function prune_mhist(mhist, nhist)
    if eltype(mhist) <: Vector
        for i in eachindex(mhist)
            mhist[i] = mhist[i][1:nhist[i]]
        end
    else
        mhist = mhist[1:nhist]
    end
    return mhist
end

"""
    set_onoff(onstates, bins, nalleles, coupling)

TBW
"""
function set_onoff(onstates, bins, nalleles, coupling)
    if isempty(coupling)
        if ~(eltype(onstates) <: Vector)
            onstates = [onstates]
        end
        if ~(eltype(bins) <: Vector)
            bins = [bins]
        end
        if length(bins) != length(onstates)
            throw(ArgumentError("Number of time bin vectors ($(length(bins))) does not match number of onstates ($(length(onstates)))"))
        end
        before, after, ndt, dt, histofftdd, histontdd, tIA, tAI = set_onoff(onstates, bins, nalleles)
        return onstates, bins, before, after, ndt, dt, histofftdd, histontdd, tIA, tAI
    else
        if !(eltype(onstates) <: Vector{Vector{T}} where {T})
            throw(ArgumentError("onstates must be a Vector{Vector{Vector}}, got $(typeof(onstates))"))
        end
        if ~(eltype(bins) <: Vector && eltype(eltype(bins)) <: Vector)
            throw(ArgumentError("bins must be a Vector{Vector{Vector}}, got $(typeof(bins))"))
        end

        tIA = Vector[]
        tAI = Vector[]
        before = Vector[]
        after = Vector[]
        ndt = Vector[]
        dt = Vector[]
        histofftdd = Vector[]
        histontdd = Vector[]

        for i in eachindex(onstates)
            b, a, nndt, ddt, off, on, IA, AI = set_onoff(onstates[i], bins[i], 1)
            push!(before, b)
            push!(after, a)
            push!(ndt, nndt)
            push!(dt, ddt)
            push!(histofftdd, off)
            push!(histontdd, on)
            push!(tIA, IA)
            push!(tAI, AI)
        end
        return onstates, bins, before, after, ndt, dt, histofftdd, histontdd, tIA, tAI
    end
end
"""
    set_onoff(onstates, bins, nalleles)

TBW
"""
function set_onoff(onstates, bins, nalleles)
    nn = length(onstates)
    tIA = Vector{Float64}[]
    tAI = Vector{Float64}[]
    before = Vector{Int}(undef, nn)
    after = Vector{Int}(undef, nn)
    ndt = Int[]
    dt = Float64[]
    histofftdd = Vector{Int}[]
    histontdd = Vector{Int}[]
    for i in eachindex(onstates)
        push!(ndt, length(bins[i]))
        push!(dt, bins[i][2] - bins[i][1])
        push!(histofftdd, zeros(Int, ndt[i]))
        push!(histontdd, zeros(Int, ndt[i]))
        push!(tIA, zeros(nalleles))
        push!(tAI, zeros(nalleles))
    end
    return before, after, ndt, dt, histofftdd, histontdd, tIA, tAI
end


"""
    set_before(onstates, state, allele, G, R, inserstep)

find before and after states for the same allele to define dwell time histograms
"""
function set_before(before, onstates, state, allele, G::Int, R, insertstep)
    for i in eachindex(onstates)
        before[i] = isempty(onstates[i]) ? num_reporters(state, allele, G, R, insertstep) : Int(gstate(G, state, allele) ∈ onstates[i])
    end
    before
end

function set_before(before, onstates, state, allele, G::Tuple, R, insertstep)
    for i in eachindex(onstates)
        for j in eachindex(onstates[i])
            before[i][j] = isempty(onstates[i][j]) ? num_reporters(state[i], 1, G[i], R[i], insertstep[i]) : Int(gstate(G[i], state[i], 1) ∈ onstates[i][j])
        end
    end
    before
end

"""
    set_after!(histofftdd, histontdd, tAI, tIA, dt, ndt, before, after, t, onstates, state, allele, G::Int, R, insertstep, verbose)

TBW
"""
function set_after!(histofftdd, histontdd, tAI, tIA, dt, ndt, before, after, t, onstates, state, allele, G::Int, R, insertstep, verbose)
    for i in eachindex(onstates)
        after[i] = isempty(onstates[i]) ? num_reporters(state, allele, G, R, insertstep) : Int(gstate(G, state, allele) ∈ onstates[i])
        firstpassagetime!(histofftdd[i], histontdd[i], tAI[i], tIA[i], t, dt[i], ndt[i], allele, before[i], after[i], verbose)
    end
    verbose && println(tAI)
    verbose && println(tIA)
end

"""
    set_after!(histofftdd, histontdd, tAI, tIA, dt, ndt, before, after, t, onstates, state, allele, G::Tuple, R, insertstep, verbose)

TBW
"""
function set_after!(histofftdd, histontdd, tAI, tIA, dt, ndt, before, after, t, onstates, state, allele, G::Tuple, R, insertstep, verbose)
    for i in eachindex(onstates)
        for j in eachindex(onstates[i])
            after[i][j] = isempty(onstates[i][j]) ? num_reporters(state[i], 1, G[i], R[i], insertstep[i]) : Int(gstate(G[i], state[i], 1) ∈ onstates[i][j])
            firstpassagetime!(histofftdd[i][j], histontdd[i][j], tAI[i][j], tIA[i][j], t, dt[i][j], ndt[i][j], 1, before[i][j], after[i][j], verbose)
        end
    end
    verbose && println(tAI)
    verbose && println(tIA)
end



"""
    num_reporters(state, allele, G::Int, R, insertstep)

return number of states with R steps > 1
"""
num_reporters(state, allele, G::Int, R, insertstep) = sum(state[G+insertstep:G+R, allele] .> 1)


"""
    firstpassagetime!(histofftdd, histontdd, tAI, tIA, t, dt, ndt, allele, before, after, verbose)

decide if transition exits or enters sojourn states then in place update appropriate histogram
"""
function firstpassagetime!(histofftdd, histontdd, tAI, tIA, t, dt, ndt, allele, before, after, verbose)
    if before == 1 && after == 0  # turn off
        firstpassagetime!(histontdd, tAI, tIA, t, dt, ndt, allele)
        if verbose
            println("off:", allele)
        end
    elseif before == 0 && after == 1 # turn on
        firstpassagetime!(histofftdd, tIA, tAI, t, dt, ndt, allele)
        if verbose
            println("on:", allele)
        end
    end
end

"""
    firstpassagetime!(hist, t1, t2, t, dt, ndt, allele)

in place update of first passage time histograms
"""
function firstpassagetime!(hist, t1, t2, t, dt, ndt, allele)
    t1[allele] = t
    t12 = (t - t2[allele]) / dt
    if t12 <= ndt && t12 > 0 && t2[allele] > 0
        hist[ceil(Int, t12)] += 1
    end
end





"""
    compute_transition_matrix(trace, G, R, S, insertstep; coupling=tuple())

Compute the state-to-state transition probability matrix from frame-to-frame using trace data.

# Arguments
- `trace`: Output from make_trace function
- `G`: Number of gene states (Int for uncoupled, Tuple for coupled models)
- `R`: Number of pre-RNA steps
- `S`: Number of splice sites
- `insertstep`: Reporter insertion step
- `coupling`: Optional coupling tuple for coupled models

# Returns
- For uncoupled models: A matrix P where P[i,j] is probability of transitioning from state i to j
- For coupled models: A matrix P where P[i,j] is probability of transitioning from combined state i to j
"""

# jointstates = [StochasticGene.coupled_state_index(collect(r),G,R,S) for r in eachrow(states)]


function compute_transition_matrix(states, nstates::Int)
    states = Int.(states)
    a = zeros(nstates, nstates)
    p0 = zeros(nstates)
    for i in 1:length(states)-1
        a[states[i], states[i+1]] += 1.0
        p0[states[i]] += 1.0
    end
    a ./ max.(sum(a, dims=2), eps()), p0 / sum(p0)
end


function compute_transition_matrix(trace, G::Tuple, R, S)
    states = [trace[1][1][:, 4] trace[1][2][:, 4]]
    jointstates = [StochasticGene.coupled_state_index(collect(r), G, R, S) for r in eachrow(states)]
    nstates = prod(T_dimension.(G, R, S))
    compute_transition_matrix(jointstates, nstates)
end

# function compute_transition_matrixa(trace, G, R, S, insertstep; coupling=tuple())
#     if isempty(coupling)
#         # Single unit case
#         n_states = T_dimension(G, R, S)
#         P = zeros(n_states, n_states)

#         # Get state indices from trace data
#         gstates = trace[:, "Gstate1"]
#         rstates = trace[:, "Rstate1"]

#         # Count transitions
#         for i in 1:(length(gstates)-1)
#             # Create state vectors for current and next states
#             current_state = zeros(Int, G + R, 1)
#             next_state = zeros(Int, G + R, 1)

#             # Set G states
#             current_state[gstates[i], 1] = 1
#             next_state[gstates[i+1], 1] = 1

#             # Set R states if they exist
#             if R > 0
#                 current_state[G+1:end, 1] = rstates[i]
#                 next_state[G+1:end, 1] = rstates[i+1]
#             end

#             # Convert to indices
#             from_state = state_index(current_state, G, R, S)
#             to_state = state_index(next_state, G, R, S)

#             P[from_state, to_state] += 1
#         end

#         # Normalize to get probabilities
#         for i in 1:n_states
#             row_sum = sum(P[i,:])
#             if row_sum > 0
#                 P[i,:] ./= row_sum
#             end
#         end

#         return P
#     else
#         # Coupled model case
#         if !(G isa Tuple)
#             error("G must be a tuple for coupled models")
#         end

#         # Calculate total number of combined states
#         n_states = prod(T_dimension.(G, R, S))
#         P = zeros(n_states, n_states)

#         # For coupled models, process each unit's states
#         for i in 1:length(trace)-1
#             # Create state matrices for both units
#             current_states = Vector{Matrix{Int}}()
#             next_states = Vector{Matrix{Int}}()

#             for unit in 1:length(G)
#                 # Get state data for current unit
#                 gstate_col = Symbol("Gstate1_$unit")
#                 rstate_col = Symbol("Rstate1_$unit")

#                 # Create current state matrix
#                 current_state = zeros(Int, G[unit] + R[unit], 1)
#                 current_state[trace[i, gstate_col], 1] = 1
#                 if R[unit] > 0
#                     current_state[G[unit]+1:end, 1] = trace[i, rstate_col]
#                 end
#                 push!(current_states, current_state)

#                 # Create next state matrix
#                 next_state = zeros(Int, G[unit] + R[unit], 1)
#                 next_state[trace[i+1, gstate_col], 1] = 1
#                 if R[unit] > 0
#                     next_state[G[unit]+1:end, 1] = trace[i+1, rstate_col]
#                 end
#                 push!(next_states, current_state)
#             end

#             # Convert to indices using the coupled state_index
#             from_state = state_index(current_states, G, R, S)
#             to_state = state_index(next_states, G, R, S)

#             P[from_state, to_state] += 1
#         end

#         # Normalize to get probabilities
#         for i in 1:n_states
#             row_sum = sum(P[i,:])
#             if row_sum > 0
#                 P[i,:] ./= row_sum
#             end
#         end

#         return P
#     end
# end
