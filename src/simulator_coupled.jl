# This file is part of StochasticGene.jl

# simulator.jl
# Functions to simulate Markov gene transcription models
# Uses hybrid first and next reaction method

"""
	Reaction

structure for reaction type

action: type of reaction
index: rate index for reaction
disabled: reactions that are disabled by current reaction
enabled: reactions that are enabled by current reaction
initial: initial GR state
final: final GR state
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

structure for rate indices of reaction types
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

create dictionary for all the possible transitions
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
- `nalleles`: Number of alleles, set to 1 if coupling nonempty
- `nhist::Int`: Size of mRNA histogram
- `onstates::Vector`: a vector of vector of ON states (use empty set for any R step is ON)
- `probfn`=prob_GaussianMixture: reporter distribution
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
function simulator(r, transitions, G, R, S, insertstep; coupling=tuple(), nalleles=1, nhist=20, onstates=Int[], bins=Float64[], traceinterval::Float64=0.0, probfn=prob_Gaussian, noiseparams=4, totalsteps::Int=100000000, totaltime::Float64=0.0, tol::Float64=1e-6, reporterfn=sum, splicetype="", verbose::Bool=false)

    if !isempty(coupling)
        coupling, nalleles, noiseparams, r = prepare_coupled(r, coupling, transitions, G, R, S, insertstep, nalleles, noiseparams)
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

    while (err > tol && steps < totalsteps) || t < totaltime
        steps += 1
        t, index, allele = findmin_tau(tau)   # reaction index and allele for least time transition
        initial, final, disabled, enabled, action = set_arguments(reactions, index)
        dth = t - t0
        t0 = t
        update_mhist!(mhist, m, dth, nhist)

        if t - ts > tsample && traceinterval == 0
            err, mhist0 = update_error(mhist, mhist0)
            ts = t
        end

        if onoff
            before = set_before(before, onstates, state, allele, G, R, insertstep)
        end

        if verbose
            println("---")
            println("m:", m)
            println(state)
            println(tau)
            println("t:", t)
            println(index, " ", allele)
            println(invactions[action], " ", allele)
            println(initial, "->", final)
            println(steps)
        end

        m = update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling)

        if onoff
            set_after!(histofftdd, histontdd, tAI, tIA, dt, ndt, before, after, t, onstates, state, allele, G, R, insertstep, verbose)
        end
        if traceinterval > 0
            set_tracelog!(tracelog, t, state)

        end
    end  # while
    verbose && println(steps)
    # counts = max(sum(mhist), 1)
    # mhist /= counts

    results = []
    push!(results, prune_mhist(mhist, nhist))
    if onoff
        for i in eachindex(histontdd)
            push!(results, histontdd[i])
            push!(results, histofftdd[i])
        end
    end
    if traceinterval > 0.0
        push!(results, make_trace(tracelog, G, R, S, insertstep, onstates, traceinterval, par, probfn, reporterfn))
    end
    if verbose
        push!(results, tracelog)
    end
    return results

    # if onoff && traceinterval == 0.
    #     dwelltimes = Vector[]
    #     push!(dwelltimes, prune_mhist(mhist, nhist))
    #     for i in eachindex(histontdd)
    #         push!(dwelltimes, histontdd[i])
    #         push!(dwelltimes, histofftdd[i])
    #     end
    #     return dwelltimes
    # elseif traceinterval > 0.0 && !onoff
    #     return [prune_mhist(mhist, nhist), make_trace(tracelog, G, R, S, insertstep, onstates, traceinterval, par, probfn, reporterfn), tracelog]
    # elseif onoff && traceinterval > 0
    #     return [prune_mhist(mhist, nhist), histontdd, histofftdd, make_trace(tracelog, G, R, S, insertstep, onstates, traceinterval, par, probfn, reporterfn)]
    # else
    #     return prune_mhist(mhist, nhist)
    # end
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
            throw("number of time bin vectors do not match number of onstates")
        end
        before, after, ndt, dt, histofftdd, histontdd, tIA, tAI = set_onoff(onstates, bins, nalleles)
        return onstates, bins, before, after, ndt, dt, histofftdd, histontdd, tIA, tAI
    else
        if ~(eltype(onstates) <: Vector && eltype(eltype(onstates)) <: Vector)
            throw("onstates not a vector of vectors of vectors")
        end
        if ~(eltype(bins) <: Vector && eltype(eltype(bins)) <: Vector)
            throw("bins not a vector of vectors of vectors")
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
    set_tracelog(tracelog, t, state::Vector{Matrix})

TBW
"""
function set_tracelog!(tracelog, t, state::Vector{Matrix})
    for i in eachindex(state)
        set_tracelog!(tracelog[i], t, state[i])
    end
end

"""
    set_tracelog!(tracelog, t, state::Matrix)

TBW
"""
function set_tracelog!(tracelog, t, state::Matrix)
    push!(tracelog, (t, state[:, 1]))
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
    set_par(r, noiseparams::Vector)
    set_par(r, noiseparams::Int)

TBW
"""
function set_par(r, noiseparams::Vector)
    par = Vector[]
    for i in eachindex(noiseparams)
        push!(par, set_par(r[i], noiseparams[i]))
    end
    par
end
set_par(r, noiseparams::Int) = r[end-noiseparams+1:end]

"""
    targets(coupling)

Find coupled targets for each unit
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

TBW
"""
function prepare_coupled(r, coupling, transitions, G, R, S, insertstep, nalleles, noiseparams)
    # if nalleles isa Number
    #     nalleles = fill(1, length(G))
    # end
    if noiseparams isa Number
        noiseparams = fill(noiseparams, length(G))
    end
    coupling = (coupling..., targets(coupling))
    return coupling, nalleles, noiseparams, prepare_rates(r, coupling, transitions, R, S, insertstep, noiseparams)
end


"""
    prepare_rates(r, coupling, transitions, R, S, insertstep, noiseparams)

TBW
"""
function prepare_rates(r, coupling, transitions, R, S, insertstep, noiseparams)
    rv = Vector[]
    n = 0
    for i in eachindex(R)
        num = num_rates(transitions[i], R[i], S[i], insertstep[i]) + noiseparams[i]
        push!(rv, r[n+1:n+num])
        n += num
    end
    push!(rv, r[end-coupling[5]+1:end])
    rv
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
    simulate_trace_data(datafolder::String; ntrials::Int=10, r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231, 30, 20, 200, 100, 0.8], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=Int[], interval=1.0, totaltime=1000.0)

create simulated trace files in datafolder
"""
function simulate_trace_data(datafolder::String; ntrials::Int=10, r=[0.038, 1.0, 0.23, 0.02, 0.25, 0.17, 0.02, 0.06, 0.02, 0.000231, 30, 20, 200, 100, 0.8], transitions=([1, 2], [2, 1], [2, 3], [3, 1]), G=3, R=2, S=2, insertstep=1, onstates=Int[], interval=1.0, totaltime=1000.0)
    ~ispath(datafolder) && mkpath(datafolder)
    for i in 1:ntrials
        trace = simulator(r, transitions, G, R, S, insertstep, onstates=onstates, traceinterval=interval, totaltime=totaltime)[2][1:end-1, 2]
        l = length(trace)
        datapath = joinpath(datafolder, "testtrace$i.trk")
        writedlm(datapath, [zeros(l) zeros(l) trace collect(1:l)])
    end
end
"""
    simulate_trace_vector(r, transitions, G, R, S, interval, totaltime, ntrials; insertstep=1, onstates=Int[], reporterfn=sum)

TBW
"""
function simulate_trace_vector(r, transitions, G, R, S, interval, totaltime, ntrials; insertstep=1, onstates=Int[], reporterfn=sum)
    trace = Array{Array{Float64}}(undef, ntrials)
    for i in eachindex(trace)
        trace[i] = simulator(r, transitions, G, R, S, insertstep, onstates=onstates, traceinterval=interval, totaltime=totaltime)[2][1:end-1, 2]
    end
    trace
end

function simulate_trace_vector(r, transitions, G::Tuple, R, S, insertstep, coupling::Tuple, interval, totaltime, ntrials; onstates=Int[], reporterfn=sum)
    trace = Array{Array{Float64}}(undef, ntrials)
    for i in eachindex(trace)
        t = simulator(r, transitions, G, R, S, insertstep, coupling=coupling, onstates=onstates, traceinterval=interval, totaltime=totaltime)[2]
        tr = Vector[]
        for t in t
            tr = push!(tr, t[1:end-1,2])
        end
        trace[i] = hcat(tr...)
        end
    trace
end

"""
    simulate_trace(r, transitions, G, R, S, insertstep, interval, onstates; totaltime=1000.0, reporterfn=sum)

simulate a trace
"""
simulate_trace(r, transitions, G, R, S, insertstep, interval, onstates; totaltime=1000.0, reporterfn=sum) = simulator(r, transitions, G, R, S, insertstep, nalleles=1, nhist=2, onstates=onstates, traceinterval=interval, reporterfn=reporterfn, totaltime=totaltime)[2][1:end-1, :]


"""
    initialize(r, G::Tuple, R, reactions, nalleles, initstate=1, initreaction=1)

return initial proposed next reaction times and states
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
    initialize(r, G::Int, R, reactions, nalleles, initstate=1, initreaction=1)

TBW
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
    initialize_sim(r::Vector{Vector}, nhist, tol, samplefactor=20.0, errfactor=10.0)

"""
function initialize_sim(r::Vector{Vector}, nhist, tol, samplefactor=20.0, errfactor=10.0)
    if nhist isa Number
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

TBW
"""
initialize_sim(r::Vector{Float64}, nhist, tol, samplefactor=20.0, errfactor=10.0) = nhist, zeros(nhist + 1), ones(nhist + 1), 0, 0, 0.0, 0.0, 0.0, samplefactor / minimum(r), errfactor * tol

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
function make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Vector{Int}}, interval, par, probfn, reporterfn=sum)
    traces = []
    for o in onstates
        push!(traces, make_trace(tracelog, G, R, S, insertstep, o, interval, par, probfn, reporterfn))
    end
    traces
end

"""
    make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Int}, interval, par, probfn, reporterfn=sum)

TBW
"""
function make_trace(tracelog, G::Int, R, S, insertstep, onstates::Vector{Int}, interval, par, probfn, reporterfn=sum)
    n = length(tracelog)
    trace = Matrix(undef, 0, 4)
    state = tracelog[1][2]
    frame = interval
    if isempty(onstates)
        reporters = num_reporters_per_state(G, R, S, insertstep, reporterfn)
    else
        reporters = num_reporters_per_state(G, onstates)
    end
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
    intensity(state, G, R, S, d)

Returns the trace intensity given the state of a system

For R = 0, the intensity is occupancy of any onstates
For R > 0, intensity is the number of reporters in the nascent mRNA

"""
function intensity(state, G, R, S, d)
    stateindex = state_index(state, G, R, S)
    max(rand(d[stateindex]), 0)
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
function state_index(state::Array, G, R, S, allele=1)
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
    set_reactionindices(Gtransitions, R::Int, S, insertstep)

return structure of ranges for each type of transition
"""
function set_reactionindices(Gtransitions, R::Int, S, insertstep)
    if S > 0
        S = R - insertstep + 1
    end
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
    indices = set_reactionindices(Gtransitions, R, S, insertstep)
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
            push!(reactions, Reaction(actions["initiate!"], i, [], [nG + 2; nG + 2 + S], G, G + 1))
        else
            push!(reactions, Reaction(actions["initiate!"], i, [], [nG + 2], G, G + 1))
        end
    end
    i = G
    for r in indices.rrange
        i += 1
        if S == 0
            push!(reactions, Reaction(actions["transitionR!"], r, Int[r], [r - 1; r + 1], i, i + 1))
        end
        if S > 0 && insertstep == 1
            push!(reactions, Reaction(actions["transitionR!"], r, Int[r; r + Sstride], [r - 1; r + 1; r + 1 + Sstride], i, i + 1))
        end
        if S > 0 && insertstep > 1 && i > G + insertstep - 1
            push!(reactions, Reaction(actions["transitionR!"], r, Int[r; r + Sstride], [r - 1; r + 1; r + 1 + Sstride], i, i + 1))
        end
        if S > 0 && insertstep > 1 && i == G + insertstep - 1
            push!(reactions, Reaction(actions["transitionR!"], r, Int[r], [r - 1; r + 1; r + 1 + Sstride], i, i + 1))
        end
        if S > 0 && insertstep > 1 && i < G + insertstep - 1
            push!(reactions, Reaction(actions["transitionR!"], r, Int[r], [r - 1; r + 1], i, i + 1))
        end
    end
    for e in indices.erange
        if S > 0
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
        push!(reactions, Reaction(actions["splice!"], s, Int[], Int[], j, 0))
    end
    push!(reactions, Reaction(actions["decay!"], indices.decay, Int[], Int[indices.decay], 0, 0))
    return reactions
end
"""
    update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling=0)

updates proposed next reaction time and state given the selected action and returns updated number of mRNA

(uses if-then statements because that executes faster than an element of an array of functions)

Arguments are same as defined in simulator

"""
function update!(tau, state, index, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final, action, coupling=0)
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
                m = eject!(tau, state, index, t, m, r, allele, G, R, S, disabled, enabled, initial)
            end
        else
            if action == 7
                splice!(tau, state, index, t, m, r, allele, G, R, initial)
            else
                m = decay!(tau, index, t, m, r)
            end
        end
    end
    return m
end

"""
    couplingG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)

TBW
"""
function couplingG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
    models = coupling[1]
    unit = models[index]
    sources = coupling[2][unit]
    sstate = coupling[3]
    ttrans = coupling[4]
    targets = coupling[6][unit]

    # unit just changed state

    # unit as source
    # new state moves into a primed source state

    for target in targets
        if isfinite(tau[target][ttrans[target], 1])
            if initial != sstate[unit] && final == sstate[unit]
                tau[target][ttrans[target], 1] = 1 / (1 + r[end][models[unit]]) * (tau[target][ttrans[target], 1] - t) + t
            elseif initial == sstate[unit] && final != sstate[unit]
                tau[target][ttrans[target], 1] = (1 + r[end][unit]) * (tau[target][ttrans[target], 1] - t) + t
            end
        end
       
    end
    # unit as target
    # new state moves to a target transition
    for source in sources
        if ttrans[index] ∈ enabled
            if state[source][sstate[source], 1] > 0
                tau[unit][ttrans[unit], 1] = 1 / (1 + r[end][source]) * (tau[unit][ttrans[unit], 1] - t) + t
            end
        end
    end

end


"""
    transitionG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)

update tau and state for G transition

"""
function transitionG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
    transitionG!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], 1, G[index[1]], R[index[1]], disabled, enabled, initial, final, coupling)
    couplingG!(tau, state, index[1], t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
end
"""
    transitionG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())

TBW
"""
function transitionG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())
    for e in enabled
        tau[e, allele] = -log(rand()) / r[e] + t
    end
    for d in disabled
        tau[d, allele] = Inf
    end
    state[final, allele] = 1
    state[initial, allele] = 0
end
"""
    activateG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)

"""
function activateG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
    activateG!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], disabled, enabled, initial, final, coupling)
    couplingG!(tau, state, index[1], t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
end
"""
    activateG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())

TBW
"""
function activateG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())
    transitionG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
    if R > 0 && state[G+1, allele] > 0
        tau[enabled[end], allele] = Inf
    end
end
"""
    deactivateG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
 
"""
function deactivateG!(tau, state, index::Tuple, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
    deactivateG!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], 1, G[index[1]], R[index[1]], disabled, enabled, initial, final, coupling)
    couplingG!(tau, state, index[1], t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
end
"""
    deactivateG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())

TBW
"""
function deactivateG!(tau, state, index::Int, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling=tuple())
    transitionG!(tau, state, index, t, m, r, allele, G, R, disabled, enabled, initial, final, coupling)
end
"""
    initiate!(tau, state, index::Tuple, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
    initiate!(tau, state, index::Int, t, m, r, allele, G, R, S, disabled, enabled, initial, final,nsertstep)

"""
function initiate!(tau, state, index::Tuple, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
    initiate!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], 1, G[index[1]], R[index[1]], S[index[1]], insertstep[index[1]], disabled, enabled, initial, final)
end
"""
    initiate!(tau, state, index::Int, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)

TBW
"""
function initiate!(tau, state, index::Int, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
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
function transitionR!(tau, state, index::Tuple, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
    transitionR!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], S[index[1]], insertstep[index[1]], disabled, enabled, initial, final)
end
"""
    transitionR!(tau, state, index::Int, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)

TBW
"""
function transitionR!(tau, state, index::Int, t, m, r, allele, G, R, S, insertstep, disabled, enabled, initial, final)
    if state[initial-1, allele] > 0
        tau[enabled[1], allele] = -log(rand()) / r[enabled[1]] + t
    end
    if final + 1 > G + R || state[final+1, allele] == 0
        tau[enabled[2], allele] = -log(rand()) / r[enabled[2]] + t
    end
    if S > 0 && final >= G + insertstep
        # tau[enabled[3], allele] = -log(rand()) / r[enabled[3]] + t
        if final == insertstep + G
            tau[enabled[3], allele] = -log(rand()) / r[enabled[3]] + t
        elseif state[initial, allele] > 1
            tau[enabled[3], allele] = r[enabled[3]-1] / r[enabled[3]] * (tau[enabled[3]-1, allele] - t) + t
        end
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
function eject!(tau, state, index::Tuple, t, m, r, allele, G, R, S, disabled, enabled, initial)
    m[index[1]] = eject!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], S[index[1]], disabled, enabled, initial)
    m
end
"""
    eject!(tau, state, index::Int, t, m, r, allele, G, R, S, disabled, enabled, initial)

TBW
"""
function eject!(tau, state, index::Int, t, m, r, allele, G, R, S, disabled, enabled, initial)
    if state[initial-1, allele] > 0
        tau[enabled[1], allele] = -log(rand()) / (r[enabled[1]]) + t
    end
    for d in disabled
        tau[d, allele] = Inf
    end
    if R > 0
        state[initial, allele] = 0
    end
    set_decay!(tau, enabled[end], t, m, r)
end
"""
    splice!(tau, state, index::Tuple, t, m, r, allele, G, R, initial)
   
"""
function splice!(tau, state, index::Tuple, t, m, r, allele, G, R, initial)
    splice!(tau[index[1]], state[index[1]], index[2], t, m[index[1]], r[index[1]], allele, G[index[1]], R[index[1]], initial)
end
"""
    splice!(tau, state, index::Int, t, m, r, allele, G, R, initial)

TBW
"""
function splice!(tau, state, index::Int, t, m, r, allele, G, R, initial)
    state[initial, allele] = 1
    tau[index, allele] = Inf
end
"""
    decay!(tau, index::Tuple, t, m, r)
    decay!(tau, index::Int, t, m, r)

"""
function decay!(tau, index::Tuple, t, m, r)
    m[index[1]] = decay!(tau[index[1]], index[2], t, m[index[1]], r[index[1]])
    m
end
"""
    decay!(tau, index::Int, t, m, r)

TBW
"""
function decay!(tau, index::Int, t, m, r)
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
function set_decay!(tau, index::Tuple, t, m, r)
    m[index[1]] = set_decay!(tau[index[1]], index[2], t, m[index[1]], r[index[1]])
    m
end
"""
    set_decay!(tau, index::Int, t, m, r)

TBW
"""
function set_decay!(tau, index::Int, t, m, r)
    m += 1
    if m == 1
        tau[index, 1] = -log(rand()) / r[index] + t
    else
        tau[index, 1] = (m - 1) / m * (tau[index, 1] - t) + t
    end
    m
end
