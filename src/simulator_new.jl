# simulator.jl
# Functions to simulate Markov gene transcription models
# Use Gibson Bruck Next Reaction method

"""
	ReactionIndices

	Structure that specifies the range of indices for each type of reaction in the rate vector

	- `grange`: range for G reactions
	- `rrange`: range for R reactions
	- `srange`: range for S reactions
	- `decay`: index for mRNA decay rate

	(eject rate follows grange for GM models, or is last index of rrange)


"""
struct ReactionIndices
    grange::UnitRange{Int64}
    rrange::UnitRange{Int64}
    srange::UnitRange{Int64}
    decay::Int
end
"""
	Reaction

	Structure that specifies a reaction

	- `action`: type of reaction (e.g. G state transition, mRNA ejection, mRNA decay, etc.)
	- `index`: index of reaction (used for selecting reaction in next reaction scheme)
	- `upstream`: vector of upstream reactions that are affected when reaction occurs
	- `downstream`: vector of downstream reactions
	- `initial`: initial state of reaction
	- `final`: final state of reaction
    - `f`: function to execute action
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
min_hepify!(a,tau,i)

- `a`: binary heap of indices of tau
- `tau`: matrix of proposed next reaction times
- `i`: index of a to be percolated down (based on tau[a])


"""
function min_heapify!(heap,node,tau,i)
    left = 2*i
    right = 2*i + 1
    smallest = i

    if left <= length(heap) && tau[heap[left]] < tau[heap[smallest]]
        smallest = left
    end

    if right <= length(heap) && tau[heap[right]] < tau[heap[smallest]]
        smallest = right
    end

    if smallest != i
        swap!(heap,node,i,smallest)
        min_heapify!(heap,tau,smallest)
    end
end

"""
build_heap(tau)

Returns an array containing a binary minimum heap of a priority index of tau

- `tau`: matrix of proposed times for next reaction

returned array is structured such that
parent is located at index div(i,2)
left child is located at index 2*i
right child is located at index 2*i+1

"""
function build_heap(tau)
    n = length(tau)
    heap = vec(collect(1:n))
    node = similar(heap)
    for i in div(n,2):-1:1
	    min_heapify!(heap,node,tau,i)
    end
    for (j,h) in enumerate(heap)
        node[h] = j
    end
    return heap,node
end
"""
update_heap!(heap,node,tau,i,i1)


"""
function update_heap!(heap,node,tau,index,taunew)
    tau[heap[index]] = taunew
    update_heap!(heap,node,tau,index)
end


function update_heap!(heap,node,tau,i)
    ip = div(i,2)
    ic = min_child(a,tau,i)
    if ip > 0 && tau[heap[i]] < tau[heap[ip]]
        swap(a,i,ip)
        update_heap!(heap,node,tau,ip)
    elseif tau[heap[i]] > tau[heap[ic]]
        swap!(heap,node,i,ic)
        update_heap!(heap,node,tau,ic)
    else
        nothing
    end
end

function min_child(heap,tau,i)
    l = 2*i
    r = 2*i + 1
    tau[heap[l]] < tau[heap[r]] ? l : r
end

function swap!(heap,node,i,ic)
    t = heap[i]
    heap[i] = heap[ic]
    heap[ic] = t
    node[heap[i]] = ic
    node[heap[ic]] = i
end

linear_index(i,j,nreactions) = (nreactions-1)*(j-1) + i


"""
initialize(r, G, R, nreactions, nalleles, initstate=1, initreaction=1)

Returns tau,states
- `tau` is an nreactions X nalleles matrix of proposed initial next reaction times
- `states` is an G + max(R,1) + 1 X nalleles matrix of occupancy of all states in the system

In all models, first G rows are for G states, and the last row is the number of free mRNA
In GM models, the rows after the G state is an auxiliary state, in GRM models, the next R rows are R states.
In GRM models these states are set to 0 for no pre-mRNA and  1 for pre-mRNA but in GRSM models they are 0,1,or 2, with 2 indicating the presence of the intron reporter and 1
indicating the pre-mRNA without the intron.
- 

Arguments
	- `r`: vector of rates
	- `G`: number of gene states
    - `R`: number of pre-RNA steps (set to 0 for GM models)
    - `nreactions`:
	- `nalleles`: Number of alleles
    - `initstate`: Index of initial occupied G state
	- `initreaction`: Index of initial reaction


"""
function initialize(r, G, R, nreactions, nalleles, initstate=1, initreaction=1,minitial=0)
    tau = fill(Inf, (nreactions-1)*nalleles + 1)
    states = zeros(Int, (G + max(R, 1))*nalleles + 1)
    for n in 1:nalleles
        tau[initreaction * n] = -log(rand()) / r[1]
        states[initstate * n] = 1
    end
    states[end] = minitial
    return build_heap(tau), tau, states
end


"""
	set_actions()

	dictionary for index of each type of transition function

    - `transitionG!`: move forward or backwards one G state 
    - `eject!`: emit an mRNA, either from active G state in GM models or last R state in GR models
    - `decay!`: removal of one mRNA
    - `initiate!`: Load first R (1 for GRM models) or RS state (2 for GRSM models)
    - `transitionR!`: move to next R state (carry S state in GRSM model)
    - `splice!`: eject splice state (only in GRS model, state goes from 2 to 1)
"""
set_actions() = Dict("transitionG!" => 1, "eject!" => 2, "decay!" => 3, "initiate!" => 4, "transitionR!" => 5, "splice!" => 6)
invert_dict(D) = Dict(D[k] => k for k in keys(D))

get_fields(reaction) = (reaction.action, reaction.initial, reaction.final, reaction.upstream, reaction.downstream)

"""
mhist, mhist0, steps, t, ts, t0, tsample, err = initialize_sim(r, nhist, tol)

- `nhist`: length of mRNA histogram
- `tol`: tolerance for histogram

- `mhist`: mRNA histogram, initially set to zeros

- `steps`: maximum number of simulation steps
- `t`: initial time set to zero
- `ts`: set to zero
- `t0`: set to zero
- `tsample`:
- `err`:
"""
initialize_sim(r, nhist, tol, samplefactor=20.0, errfactor=10.0) = zeros(nhist + 1), ones(nhist + 1), 0, 0.0, 0.0, 0.0, samplefactor / minimum(r), errfactor * tol


function simulator(model)


end

"""
	simulator(r::Vector{Float64},transitions,G::Int,R::Int,S::Int,nhist::Int,nalleles, onstates::Vector;::Int;range::Vector{Float64}=Float64[],total::Int=10000000,tol::Float64=1e-6,count=false,verbose=false)

	Simulate any GRSM model. Returns steady state mRNA histogram and if range not a null vector will return ON and OFF time histograms.

	Arguments
	- `r`: vector of rates
	- `transitions`: tuple of vectors that specify state transitions for G states, e.g. ([1,2],[2,1]) for classic 2 state telegraph model and ([1,2],[2,1],[2,3],[3,1]) for 3 state kinetic proof reading model
	- `G`: number of gene states
    - `R`: number of pre-RNA steps (set to 0 for GM models)
    - `S`: number of splice sites (set to 0 for GM models and R for GRS models)
	- `nhist::Int`: Size of mRNA histogram
	- `nalleles`: Number of alleles
	- `onstates`: vector of active states (for ON and OFF time distributions)

	Named arguments
	- `range::Vector{Float64}=Float64[]`: vector of time bins for ON and OFF histograms
	- `total::Int=10000000`: maximum number of simulation steps
	- `tol::Float64=1e-6`: convergence error tolerance for mRNA histogram
	- `verbose::Bool=false`: flag for printing state information

"""


function simulator(r::Vector{Float64}, transitions::Tuple, G::Int, R::Int, S::Int, nhist::Int, nalleles::Int; onstates::Vector=[G],range::Vector{Float64}=Float64[], total::Int=10000000, tol::Float64=1e-6, verbose::Bool=false)
    mhist, mhist0, steps, t, ts, t0, tsample, err = initialize_sim(r, nhist, tol)
    reactions = set_reactions(transitions, G, R, S)
    nreactions = length(reactions)
    heap, tau, state = initialize(r, G, R, nreactions, nalleles)
    tIA = zeros(Float64, nalleles)
    tAI = zeros(Float64, nalleles)
    if length(range) < 1
        count = false
    else
        count = true
        ndt = length(range)
        dt = range[2] - range[1]
        histofftdd = zeros(Int, ndt)
        histontdd = zeros(Int, ndt)
    end
    if verbose
        invactions = invert_dict(set_actions())
    end
    while err > tol && steps < total
        steps += 1
        # t, rindex = findmin(tau)
        tauindex = heap[1]
        t = tau[tauindex]
        index = mod(tauindex,nalleles) + 1
        allele = div(tauindex,nreactions) + 1
        action, initial, final, upstream, downstream  = get_fields(reactions[index])
        update!(heap,tau,state,index,nreactions,t,r,allele,G,R,upstream,downstream,initial,final,action)
        dth = t - t0
        t0 = t
        update_mhist!(mhist, state[end], dth, nhist)
        update_histogram()
        if t - ts > tsample
            err, mhist0 = update_error(mhist, mhist0)
            ts = t
        end
        if verbose
            println(state)
            println(num_introns(state, allele, G, R))
            println(tau)
            println(rindex)
            println(invactions[action])
        end
    end  # while
    counts = max(sum(mhist), 1)
    mhist /= counts
    if count
        return histofftdd / sum(histofftdd), histontdd / sum(histontdd), mhist[1:nhist]
    else
        return mhist[1:nhist]
    end
end


function update!(heap, tau, state, index, nreactions, t, r, allele, G, R, upstream, downstream, initial, final, action)
    if action < 4
        if action == 1
            transitionG!(heap,tau, state, index, nreactions, t, r, allele, G, R, upstream, downstream, initial, final)
        else
            if action == 2
                eject!(heap,tau, state, index, nreactions, t, r, allele, G, R, S, upstream, downstream)
            else
                decay!(heap,tau, state, index, nreactions, t, r)
            end
        end
    else
        if action == 4
            initiate!(heap,tau, state, index, nreactions, t, r, allele, G, R, S, downstream)
        else
            if action == 5
                transitionR!(heap,tau, state, index, nreactions, t, r, allele, G, R, S, upstream, downstream, initial, final)
            else
                splice!(heap,tau, state, index, nreactions, t, r, allele, G, R, initial)
            end
        end
    end
end

function num_introns(state, allele, G, R)
    d = 0
    for i in G+1:G+max(R, 1)
        d = d + Int(state[i, allele] > 1)
    end
    d
end

function ontime!(histon, tIA, tAI, t, dt, ndt, state, allele, G, R)
    if num_introns(state, allele, G, R) == 1
        firstpassagetime!(histon, tAI, tIA, t, dt, ndt, allele)
    end
end

function offtime!(histoff, tIA, tAI, t, dt, ndt, state, allele, G, R)
    if num_introns(state, allele, G, R) == 0
        firstpassagetime!(histoff, tIA, tAI, t, dt, ndt, allele)
    end
end

function firstpassagetime!(hist, t1, t2, t, dt, ndt, allele)
    t1[allele] = t
    t12 = (t - t2[allele]) / dt
    if t12 <= ndt && t12 > 0 && t2[allele] > 0
        hist[ceil(Int, t12)] += 1
    end
end

function set_reactionindices(Gtransitions, R, S)
    g = 1:length(Gtransitions)
    r = length(Gtransitions)+1:length(Gtransitions)+1+R
    s = length(Gtransitions)+1+R+1:length(Gtransitions)+1+R+S
    d = length(Gtransitions) + 1 + R + S + 1
    ReactionIndices(g, r, s, d)
end
"""
set_reactions(Gtransitions,G,R,S,indices,actions)

create a vector of Reactions
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
		push!(reactions, Reaction(actions["transitionG!"], g, u, d, ginitial, gfinal))
    end
    if R > 0
        # set downstream to splice reaction
        push!(reactions, Reaction(actions["initiate!"], indices.rrange[1], Int[], [nG + 2 + S], G, G + 1))
    end
    i = G
    for r in indices.rrange
        d = Int[]
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


function update_tau!(heap,tau,node,nreactions,index,allele,taunew)
    ind = linear_index(index, allele,nreactions)
    tau[ind] = taunew
    update_heap!(heap,node,tau,ind,taunew)
end

function update_tau!(heap,tau,node,nreactions,index,taunew)
    tau[end] = taunew
    update_heap!(heap,node,tau,end,taunew)
end


"""
	transitionG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)

"""
function transitionG!(heap, tau, state, index, nreactions, t, r, allele, G, R, upstream, downstream, initial, final)
    update_tau!(heap,tau,node,nreactions,index,allele,Inf)
    # ind = linear_index(index, allele,nreactions)
    # tau[ind] = Inf
    # update_heap!(heap,node,tau,ind,Inf)
    state[linear_index(final, allele,nreactions)] = 1
    state[linear_index(initial, allele,nreactions)] = 0
    for d in downstream
        update_tau!(heap,tau,node,nreactions,d,allele,-log(rand()) / r[d] + t)
        # taunew = -log(rand()) / r[d] + t
        # ind = linear_index(d, allele,nreactions)
        # tau[ind] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
    end
    for u in upstream
        update_tau!(heap,tau,node,nreactions,u,allele,Inf)
        # ind = linear_index(u, allele,nreactions)
        # tau[ind] = Inf
        # update!(heap,node,tau,ind,Inf)
    end
end
# """
# 	activateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)

# """
# function activateG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
#     transitionG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
#     if R == 0
#         state[G+1, allele] = 2
#     elseif state[G+1, allele] > 0
#         tau[downstream[end], allele] = Inf
#     end
# end
# """
# 	deactivateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)

# """
# function deactivateG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
#     transitionG!(tau, state, index, t, m, r, allele, G, R, upstream, downstream, initial, final)
#     if R == 0
#         state[G+1, allele] = 0
#     end
# end

"""
	initiate!(tau,state,index,t,m,r,allele,G,R,S,downstream)

"""
function initiate!(tau, state, index, nreactions, t, r, allele, G, R, S, downstream)
    update_tau!(heap,tau,node,nreactions,index,allele,Inf)
    # ind = linear_index(index, allele,nreactions)
    # tau[ind] = Inf
    # update_heap!(heap,node,tau,ind,Inf)
    state[linear_index(G+1, allele,nreactions)] = 2
    if R < 2 || state[linear_index(G+2, allele,nreactions)] < 1
        update_tau!(heap,tau,node,nreactions,index+1,allele,-log(rand()) / (r[index+1]) + t)
        # taunew = -log(rand()) / (r[index+1]) + t
        # ind = linear_index(index+1, allele,nreactions)
        # tau[ind] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
    end
    if S > 0
        update_tau!(heap,tau,node,nreactions,donwstream[1],allele,-log(rand()) / (r[downstream[1]]) + t)
        # taunew = -log(rand()) / (r[downstream[1]]) + t
        # ind = linear_index(downstream[1], allele,nreactions)
        # tau[linear_index(downstream[1], allele,nreactions)] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
    end
end
"""
	transitionR!(tau,state,index,t,m,r,allele,G,R,S,upstream,downstream,initial,final)

"""
function transitionR!(tau, state, index, nreactions, t, r, allele, G, R, S, u, d, initial, final)
    update_tau!(heap,tau,node,nreactions,index,allele,Inf)
    # ind = linear_index(index, allele,nreactions)
    # tau[ind] = Inf
    # update_heap!(heap,node,tau,ind,Inf)
    if S > 0 && isfinite(tau[index+S, allele])
        update_tau!(heap,tau,node,nreactions,index+S+1,allele,tau[linear_index(index+S, allele,nreactions)])
        update_tau!(heap,tau,node,nreactions,index+S,allele,Inf)

        # ind = linear_index(index+S+1, allele,nreactions)
        # indS = linear_index(index+S, allele,nreactions)
        # taunew = tau[indS]
        # tau[ind] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
        # tau[indS] = Inf
        # update_heap!(heap,node,tau,indS,Inf)
    end
    if state[linear_index(initial-1, allele,nreactions)] > 0
        update_tau!(heap,tau,node,nreactions,u[1],allele,-log(rand()) / r[u[1]] + t)
        # taunew = -log(rand()) / r[u[1]] + t
        # ind = linear_index(u[1], allele,nreactions)
        # tau[ind] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
    else
        update_tau!(heap,tau,node,nreactions,u[1],allele,Inf)
        # ind = linear_index(u[1],allele,nreactions)
        # tau[ind] = Inf
        # update_heap!(heap,node,tau,ind,Inf)
    end
    if final == G + R || state[final+1, allele] < 1
        update_tau!(heap,tau,node,nreactions,d[1],allele,-log(rand()) / r[d[1]] + t)
        # taunew = -log(rand()) / r[d[1]] + t
        # ind = linear_index(d[1], allele,nreactions)
        # tau[ind] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
    else
        update_tau!(heap,tau,node,nreactions,d[1],allele,Inf)
        # ind = linear_index(d[1], allele,nreactions)
        # tau[ind] = Inf
        # update_heap!(heap,node,tau,ind,Inf)
    end
    state[final, allele] = state[initial, allele]
    state[initial, allele] = 0
end
"""
eject!

"""
function eject!(tau, state, index, t, r, allele, G, R, S, upstream, downstream)
    m = state[end]
    m += 1
    set_decay!(tau, downstream[end], t, m, r)
    if S > 0 && isfinite(tau[index+R, allele])
        update_tau!(heap,tau,node,nreactions,index+R,allele,Inf)
        # ind = linear_index(index+R, allele,nreactions)
        # tau[inde] = Inf
        # update_heap!(heap,node,tau,ind,Inf)
    end
    if R > 0
        update_tau!(heap,tau,node,nreactions,index,allele,Inf)
        # ind = linear_index(index, allele,nreactions)
        # tau[index, allele] = Inf
        # update_heap!(heap,node,tau,ind,Inf)
        state[G+R, allele] = 0
        if state[G+R-1, allele] > 0
            taunew = -log(rand()) / (r[upstream[1]]) + t
            ind = linear_index(upstream[1], allele,nreactions)
            tau[ind] = -log(rand()) / (r[upstream[1]]) + t
            update_heap!(heap,node,tau,ind,taunew)
        end
    else
        update_tau!(heap,tau,node,nreactions,index,allele, -log(rand()) / (r[index]) + t)
        # taunew = -log(rand()) / (r[index]) + t
        # ind = linear_index(index, allele,nreactions)
        # tau[ind] = taunew
        # update_heap!(heap,node,tau,ind,taunew)
    end
    state[end] = m
end
"""
splice!

"""
function splice!(tau, state, index, t, m, r, allele, G, R, initial)
    state[initial, allele] -= 1
    update_tau!(heap,tau,node,nreactions,index,allele,Inf)
    # tau[index, allele] = Inf
    # update_heap!(heap,node,tau,index,Inf)
end
"""
decay!

"""
function decay!(tau, state, index, t, r)
    m = state[end] -= 1
    update_tau!(heap,tau,node,nreactions, -log(rand()) / (state[end] * r[index]) + t)
    # taunew = -log(rand()) / (m * r[index]) + t
    # tau[index, 1] = taunew
    # update_heap!(heap,node,tau,index,taunew)
end

"""
set_decay!(tau,reaction,t,m)
update tau matrix for decay rate

"""
function set_decay!(tau, index, t, m, r)
    if m == 1
        update_tau!(heap,tau,node,nreactions, -log(rand()) / r[index] + t)
        # taunew = -log(rand()) / r[index] + t
        # tau[index, 1] = taunew
        # update_heap!(heap,node,tau,index,taunew)
    else
        update_tau!(heap,tau,node,nreactions, (m - 1) / m * (tau[end] - t) + t)
        # taunew = = (m - 1) / m * (tau[index, 1] - t) + t
        # tau[index, 1] = taunew
        # update_heap!(heap,node,tau,index,taunew)
    end
end
