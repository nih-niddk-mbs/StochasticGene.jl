# simulator.jl
# Functions to simulate Markov gene transcription models
# Uses simplified next reaction method

"""
	Reaction

	structure for reaction type
"""
struct Reaction
    action::Int
    index::Int
    upstream::Vector{Int64}
    downstream::Vector{Int64}
    initial::Int
    final::Int
end


"""
	ReactionIndices

	structure for rate indices of reaction types
"""
struct ReactionIndices
    grange::UnitRange{Int64}
    rrange::UnitRange{Int64}
    srange::UnitRange{Int64}
    decay::Int
end

"""
	set_actions()

	create dictionary for all the possible transitions
"""
set_actions() = Dict("activateG!" => 1, "deactivateG!" => 2, "transitionG!" => 3, "initiate!" => 4, "transitionR!" => 5, "eject!" => 6, "splice!" => 7, "decay!" => 8)
invert_dict(D) = Dict(D[k] => k for k in keys(D))

"""
	simulator(r::Vector{Float64},transitions,G::Int,R::Int,S::Int,nhist::Int,nalleles::Int;onstates::Vector{Int}=[G],range::Vector{Float64}=Float64[],total::Int=10000000,tol::Float64=1e-6, par=[30, 14, 200, 75], verbose::Bool=false)

	Simulate any GRSM model. Returns steady state mRNA histogram and if range not a null vector will return ON and OFF time histograms.
    If trace is set to true, it returns a nascent mRNA trace

	Arguments
	- `r`: vector of rates
	- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
	- `G`: number of gene states
    - `R`: number of pre-RNA steps (set to 0 for classic telegraph models)
    - `S`: number of splice sites (set to 0 for G (classic telegraph) and GR models and R for GRS models)
	- `nhist::Int`: Size of mRNA histogram
	- `nalleles`: Number of alleles

	Named arguments
    - `onstates::Vector`: a vector of ON G states
	- `range::Vector{Float64}=Float64[]`: vector of time bins for ON and OFF histograms
	- `totalsteps::Int=10000000`: maximum number of simulation steps
	- `tol::Float64=1e-6`: convergence error tolerance for mRNA histogram (not used when traces are made)
    - `traceinterval`: Interval in minutes between frames for intensity traces.  If 0, traces are not made.
	- `par=[30, 14, 200, 75, 0.2]`: Vector of 5 parameters for noise model [background mean, background std, signal mean, signal std, weight of background <= 1]
    - `verbose::Bool=false`: flag for printing state information


    Examples:

    julia> trace = simulator([.1,.02,.1,.05,.01,.01],([1,2],[2,1],[2,3],[3,1]),3,0,0,100,1,onstates=[2,3],traceinterval=100.,totalsteps = 1000)

 	julia> hoff,hon,mhist = simulator([.1,.02,.1,.05,.01,.01],([1,2],[2,1],[2,3],[3,1]),3,0,0,20,1,onstates=[2,3],range=collect(1.:100.))

    
"""
function simulator(r::Vector{Float64}, transitions::Tuple, G::Int, R::Int, S::Int, nhist::Int, nalleles::Int; onstates::Vector{Int}=[G], range::Vector{Float64}=Float64[], totalsteps::Int=10000000, tol::Float64=1e-6, traceinterval::Float64=0.0, par=[50, 20, 250, 75], verbose::Bool=false)
    mhist, mhist0, m, steps, t, ts, t0, tsample, err = initialize_sim(r, nhist, tol)
    reactions = set_reactions(transitions, G, R, S)
    tau, state = initialize(r, G, R, length(reactions), nalleles)
    tIA = zeros(Float64, nalleles)
    tAI = zeros(Float64, nalleles)
    if length(range) < 1
        onoff = false
    else
        onoff = true
        ndt = length(range)
        dt = range[2] - range[1]
        histofftdd = zeros(Int, ndt)
        histontdd = zeros(Int, ndt)
    end
    if traceinterval > 0
        tracelog = [(t, state[:, 1])]
    end
    if verbose
        invactions = invert_dict(set_actions())
    end
    while err > tol && steps < totalsteps
        steps += 1
        t, rindex = findmin(tau)
        index = rindex[1]
        allele = rindex[2]
        initial, final, upstream, downstream, action = set_arguments(reactions[index])
        dth = t - t0
        t0 = t
        update_mhist!(mhist, m, dth, nhist)
        if t - ts > tsample && traceinterval == 0
            err, mhist0 = update_error(mhist, mhist0)
            ts = t
        end
        if verbose
            println("---")
            println(state)
            if R > 0
                println(num_reporters(state, allele, G, R))
            end
            println(tau)
            println(rindex)
            println(invactions[action])
            println(initial, "->", final)
        end
        if onoff
            if R == 0
                if initial ∈ onstates && final ∉ onstates && final > 0  # turn off
                    firstpassagetime!(histontdd, tAI, tIA, t, dt, ndt, allele)
                elseif initial ∉ onstates && final ∈ onstates && final > 0 # turn on
                    firstpassagetime!(histofftdd, tIA, tAI, t, dt, ndt, allele)
                end
            else
                if num_reporters(state, allele, G, R) == 1 && ((action == 6 && state[G+R] == 2) || action == 7 )  # turn off
                    firstpassagetime!(histontdd, tAI, tIA, t, dt, ndt, allele); if verbose; println("off"); end;
                elseif action == 4 && num_reporters(state, allele, G, R) == 0   # turn on
                    firstpassagetime!(histofftdd, tIA, tAI, t, dt, ndt, allele); if verbose; println("on"); end;
                end
            end
        end
		m = update!(tau, state, index, t, m, r, allele, G, R, S, upstream, downstream, initial, final, action)
        if traceinterval > 0
            push!(tracelog, (t, state[:, 1]))
        end
    end  # while
    counts = max(sum(mhist), 1)
    mhist /= counts
    if onoff
        return histofftdd / max(sum(histofftdd), 1), histontdd / max(sum(histontdd), 1), mhist[1:nhist]
    elseif traceinterval > 0.0
        make_trace(tracelog, G, R, onstates, traceinterval, par)
    else
        return mhist[1:nhist]
    end
end

"""
update!(tau, state, index, t, m, r, allele, G, R, S, upstream, downstream, initial, final, action)

updates proposed next reaction time and state given the selected action and returns updated number of mRNA

(uses if-then statements because that executes faster than an element of an array of functions)

Arguments are same as defined in simulator

"""
function update!(tau, state, index, t, m, r, allele, G, R, S, upstream, downstream, initial, final, action)
    if action < 5
        if action < 3
            if action == 1
                activateG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
            else
                deactivateG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
            end
        else
            if action == 3
                transitionG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
            else
                initiate!(tau, state, index, t, m, r, allele, G, R, S, downstream)
            end
        end
    else
        if action < 7
            if action == 5
                transitionR!(tau, state, index, t, m, r, allele, G, R, S, upstream, downstream, initial, final)
            else
                m = eject!(tau, state, index, t, m, r, allele, G, R, S, upstream, downstream)
            end
        else
            if action == 7
                splice!(tau, state, index, t, m, r, allele, G, R, initial)
            else
                m = decay!(tau, state, index, t, m, r)
            end
        end
    end
    return m
end

"""
make_trace(tracelog, G, R, onstates, interval=100.0)

Return array of frame times and intensities

- `tracelog`: Vector if Tuples of (time,state of allele 1)
- `interval`: Number of minutes between frames
- `onstates`: Vector of G on states
- `G` and `R` as defined in simulator

"""
function make_trace(tracelog, G, R, onstates, interval=100.0, par=[30, 14, 200, 75])
    n = length(tracelog)
    trace = Matrix(undef, 0, 2)
    state = tracelog[1][2]
    frame = interval
    if R > 0
        onstates = collect(G+1:G+R)
    end
    i = 2
    d = prob_Gaussian(par, onstates, 2)
    while i < n
        while tracelog[i][1] <= frame && i < n
            state = tracelog[i][2]
            i += 1
        end
        trace = vcat(trace, [frame intensity(state, onstates, G, R, d, par)])
        frame += interval
    end
    return trace
end

"""
intensity(state,onstates,G,R)

Returns the trace intensity given the state of a system

For R = 0, the intensity is occupancy of any onstates
For R > 0, intensity is the number of reporters in the nascent mRNA

"""
function intensity(state, onstates, G, R, d, par)
    if R == 0
        return max(rand(d[sum(state[onstates] .== 1)+1]),0)
    else
        s = sum(state[G+1:G+R] .> 1)
        return rand(d[Int(s > 0)+1])
    end
end

function num_reporters(state, allele, G, R)
    d = 0
    for i in G+1:G+max(R, 1)
        d = d + Int(state[i, allele] > 1)
    end
    d
end

function firstpassagetime!(hist, t1, t2, t, dt, ndt, allele)
    t1[allele] = t
    t12 = (t - t2[allele]) / dt
    if t12 <= ndt && t12 > 0 && t2[allele] > 0
        hist[ceil(Int, t12)] += 1
    end
end

set_arguments(reaction) = (reaction.initial, reaction.final, reaction.upstream, reaction.downstream, reaction.action)

function set_reactionindices(Gtransitions, R, S)
    g = 1:length(Gtransitions)
    r = length(Gtransitions)+1:length(Gtransitions)+1+R
    s = length(Gtransitions)+1+R+1:length(Gtransitions)+1+R+S
    d = length(Gtransitions) + 1 + R + S + 1
    ReactionIndices(g, r, s, d)
end
"""
set_reactions(Gtransitions,G,R,S,indices,actions)

create a vector of Reaction structures all the possible reactions
"""
function set_reactions(Gtransitions, G, R, S)
    actions = set_actions()
    indices = set_reactionindices(Gtransitions, R, S)
    reactions = Reaction[]
    nG = length(Gtransitions)
    for g in eachindex(Gtransitions)
        u = Int[]
        d = Int[]
        ginitial = Gtransitions[g][1]
        gfinal = Gtransitions[g][2]
        for s in eachindex(Gtransitions)
            if ginitial == Gtransitions[s][1] && gfinal != Gtransitions[s][2]
                push!(u, s)
            end
            if gfinal == Gtransitions[s][1]
                push!(d, s)
            end
        end
        if gfinal == G
            push!(d, length(Gtransitions) + 1)
            push!(reactions, Reaction(actions["activateG!"], g, u, d, ginitial, gfinal))
        elseif ginitial == G
            push!(u, length(Gtransitions) + 1)
            push!(reactions, Reaction(actions["deactivateG!"], g, u, d, ginitial, gfinal))
        else
            push!(reactions, Reaction(actions["transitionG!"], g, u, d, ginitial, gfinal))
        end
    end
    if R > 0
        # set downstream to splice reaction
        push!(reactions, Reaction(actions["initiate!"], indices.rrange[1], Int[], [nG + 2 + S], G, G + 1))
    end
    i = G
    for r in indices.rrange
        if r < length(Gtransitions) + R
            i += 1
            push!(reactions, Reaction(actions["transitionR!"], r + 1, [r], [r + 2], i, i + 1))
        end
    end
    push!(reactions, Reaction(actions["eject!"], indices.rrange[end], Int[nG+R], Int[nG+R+S+2], G + R, 0))
    j = G
    for s in indices.srange
        j += 1
        push!(reactions, Reaction(actions["splice!"], s, Int[], Int[], j, 0))
    end
    push!(reactions, Reaction(actions["decay!"], indices.decay, Int[], Int[], 0, 0))
    return reactions
end

"""
initialize_sim(r, nhist, tol, samplefactor=20.0, errfactor=10.0)

"""
initialize_sim(r, nhist, tol, samplefactor=20.0, errfactor=10.0) = zeros(nhist + 1), ones(nhist + 1), 0, 0, 0.0, 0.0, 0.0, samplefactor / minimum(r), errfactor * tol

"""
    update_error(mhist, mhist0)

TBW
"""
update_error(mhist, mhist0) = (norm(mhist / sum(mhist) - mhist0 / sum(mhist0), Inf), copy(mhist))
"""
update_mhist!(mhist,m,dt,nhist)

"""
function update_mhist!(mhist, m, dt, nhist)
    if m + 1 <= nhist
        mhist[m+1] += dt
    else
        mhist[nhist+1] += dt
    end
end
"""
initialize(r, G, R, nreactions, nalleles, initstate=1, initreaction=1)

return initial proposed next reaction times and states

"""
function initialize(r, G, R, nreactions, nalleles, initstate=1, initreaction=1)
    tau = fill(Inf, nreactions, nalleles)
    states = zeros(Int, G + max(R, 1), nalleles)
    for n in 1:nalleles
        tau[initreaction, n] = -log(rand()) / r[1]
        states[initstate, n] = 1
    end
    return tau, states
end
"""
	transitionG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)

	update tau and state for G transition

"""
function transitionG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
    tau[index, allele] = Inf
    state[final, allele] = 1
    state[initial, allele] = 0
    for d in downstream
        tau[d, allele] = -log(rand()) / r[d] + t
    end
    for u in upstream
        tau[u, allele] = Inf
    end
end
"""
	activateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)

"""
function activateG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
    transitionG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
    if R == 0
        state[G+1, allele] = 2
    elseif state[G+1, allele] > 0
        tau[downstream[end], allele] = Inf
    end
end
"""
	deactivateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)

"""
function deactivateG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
    transitionG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
    if R == 0
        state[G+1, allele] = 0
    end
end

"""
	initiate!(tau,state,index,t,m,r,allele,G,R,S,downstream)

"""
function initiate!(tau, state, index, t, m, r, allele, G, R, S, downstream)
    tau[index, allele] = Inf
    state[G+1, allele] = 2
    if R < 2 || state[G+2, allele] < 1
        tau[index+1, allele] = -log(rand()) / (r[index+1]) + t
    end
    if S > 0
        tau[downstream[1], allele] = -log(rand()) / (r[downstream[1]]) + t
    end
end
"""
	transitionR!(tau,state,index,t,m,r,allele,G,R,S,upstream,downstream,initial,final)

"""
function transitionR!(tau, state, index, t, m, r, allele, G, R, S, u, d, initial, final)
    tau[index, allele] = Inf
    if S > 0 && isfinite(tau[index+S, allele])
        tau[index+S+1, allele] = -log(rand()) / (r[index+S+1]) + t
        tau[index+S, allele] = Inf
    end
    if state[initial-1, allele] > 0
        tau[u[1], allele] = -log(rand()) / r[u[1]] + t
    else
        tau[u[1], allele] = Inf
    end
    if final == G + R || state[final+1, allele] < 1
        tau[d[1], allele] = -log(rand()) / r[d[1]] + t
    else
        tau[d[1], allele] = Inf
    end
    state[final, allele] = state[initial, allele]
    state[initial, allele] = 0
end
"""
eject!

"""
function eject!(tau, state, index, t, m, r, allele, G, R, S, upstream, downstream)
    m += 1
    set_decay!(tau, downstream[end], t, m, r)
    if S > 0 && isfinite(tau[index+S, allele])
        tau[index+S, allele] = Inf
    end
    if R > 0
        tau[index, allele] = Inf
        state[G+R, allele] = 0
        if state[G+R-1, allele] > 0
            tau[upstream[1], allele] = -log(rand()) / (r[upstream[1]]) + t
        end
    else
        tau[index, allele] = -log(rand()) / (r[index]) + t
    end
    m
end
"""
splice!

"""
function splice!(tau, state, index, t, m, r, allele, G, R, initial)
    state[initial, allele] = 1
    tau[index, allele] = Inf
end
"""
decay!

"""
function decay!(tau, state, index, t, m, r)
    m -= 1
    tau[index, 1] = -log(rand()) / (m * r[index]) + t
    m
end

"""
set_decay!(tau,reaction,t,m)
update tau matrix for decay rate

"""
function set_decay!(tau, index, t, m, r)
    if m == 1
        tau[index, 1] = -log(rand()) / r[index] + t
    else
        tau[index, 1] = (m - 1) / m * (tau[index, 1] - t) + t
    end
end
