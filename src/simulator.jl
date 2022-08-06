# simulator.jl
# Functions to simulate Markov gene transcription models


struct ReactionIndices
	grange::UnitRange{Int64}
	rrange::UnitRange{Int64}
	srange::UnitRange{Int64}
	decay::Int
end

struct Reaction
	action::Int
	index::Int
	upstream::Vector{Int64}
	downstream::Vector{Int64}
	initial::Int
	final::Int
end

actions = Dict("activateG!" => 1, "deactivateG!" => 2, "transitionG!" => 3, "initiate!" => 4, "transitionR!" => 5, "eject!" => 6, "splice!" => 7, "decay!" => 8)

function simulator(r::Vector,transitions,G::Int,R::Int,S::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
	indices = set_indices(transitions,R,S)
	reactions = set_reactions(transitions,G,R,S,indices)
	tau,state = initialize(r,G,R,length(reactions),nalleles)
	while err > tol && steps < total
		steps += 1
		t,rindex = findmin(tau)
		reaction = reactions[rindex[1]]
		dt = t-t0
		t0 = t
		update_mhist!(mhist,m,dt,nhist)
		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
		end
		m = reaction.action(tau,state,reaction,t,m,r,rindex[2],G,R)
	end  # while
	counts = max(sum(mhist),1)
	mhist /= counts
	if count
		return mhist[1:nhist],counts,steps,err
	else
		return mhist[1:nhist]
	end
end

function simulator2(r::Vector{Float64},transitions,G::Int,R::Int,S::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
	indices = set_indices(transitions,R,S)
	reactions = set_reactions(transitions,G,R,S,indices,actions)
	tau,state = initialize(r,G,R,length(reactions),nalleles)
	while err > tol && steps < total
		steps += 1
		t,rindex = findmin(tau)
		index = rindex[1]
		allele = rindex[2]
		initial,final,upstream,downstream,action = set_arguments(reactions[index])
		dt = t-t0
		t0 = t
		update_mhist!(mhist,m,dt,nhist)
		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
		end
		# println(tau)
		# println(rindex)
		# println(action)
		if action < 5
			if action < 3
				if action == 1
					activateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
				else
					deactivateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
				end
			else
				if action == 3
					transitionG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
				else
					initiate!(tau,state,index,t,m,r,allele,G,R,S,downstream)
				end
			end
		else
			if action < 7
				if action == 5
					transitionR!(tau,state,index,t,m,r,allele,G,R,S,upstream,downstream,initial,final)
				else
					m = eject!(tau,state,index,t,m,r,allele,G,R,S,upstream,downstream,initial,final)
				end
			else
				if action == 7
					splice!(tau,state,index,t,m,r,allele,G,R,initial)
				else
					m = decay!(tau,state,index,t,m,r)
				end
			end
		end
	end  # while
	counts = max(sum(mhist),1)
	mhist /= counts
	if count
		return mhist[1:nhist],counts,steps,err
	else
		return mhist[1:nhist]
	end
end

set_arguments(reaction) = (reaction.initial,reaction.final,reaction.upstream,reaction.downstream,reaction.action)

function set_indices(Gtransitions,R,S)
	g = 1:length(Gtransitions)
	r = length(Gtransitions) + 1 : length(Gtransitions) + 1 + R
	s = length(Gtransitions) + 1 + R + 1 : length(Gtransitions) + 1 + R + S
	d = length(Gtransitions) + 1 + R + S + 1
	ReactionIndices(g,r,s,d)
end

function set_reactions(Gtransitions,G,R,S,indices,actions)
	reactions = Reaction[]
	nG = length(Gtransitions)
	for g in eachindex(Gtransitions)
		u = Int[]
		d = Int[]
		ginitial = Gtransitions[g][1]
		gfinal = Gtransitions[g][2]
		for s in eachindex(Gtransitions)
			if ginitial == Gtransitions[s][2] && gfinal != Gtransitions[s][1]
				push!(u,s)
			end
			if gfinal == Gtransitions[s][1]
				push!(d,s)
			end
		end
		if gfinal == G
			push!(d,length(Gtransitions)+1)
			push!(reactions,Reaction(actions["activateG!"],g,u,d,ginitial,gfinal))
		elseif ginitial == G
			push!(u,length(Gtransitions)+1)
			push!(reactions,Reaction(actions["deactivateG!"],g,u,d,ginitial,gfinal))
		else
			push!(reactions,Reaction(actions["transitionG!"],g,u,d,ginitial,gfinal))
		end
	end
	if R > 0
		push!(reactions,Reaction(actions["initiate!"],indices.rrange[1],Int[],[nG+1],G,G+1))
	end
	i = G
	for r in indices.rrange
		d = Int[]
		if r  < indices.grange[end] + R
			i += 1
			push!(reactions,Reaction(actions["transitionR!"],r+1,[r],[r+2],i,i+1))
		end
	end
	push!(reactions,Reaction(actions["eject!"],indices.rrange[end],Int[nG+R],Int[nG+R+S+2],G+R+1,0))
	for s in indices.srange
		push!(reactions,Reaction(actions["splice!"],s,Int[],Int[],s-R-1,0))
	end
	push!(reactions,Reaction(actions["decay!"],indices.decay,Int[],Int[],0,0))
	return reactions
end

initialize_sim(r,nhist,tol,samplefactor=20.,errfactor=10.) = zeros(nhist+1),ones(nhist+1),0,0,0.,0.,0.,samplefactor/minimum(r),errfactor*tol

update_error(mhist,mhist0) = (norm(mhist/sum(mhist)-mhist0/sum(mhist0),Inf),copy(mhist))

function update_mhist!(mhist,m,dt,nhist)
	if m + 1 <= nhist
		mhist[m+1] += dt
	else
		mhist[nhist+1] += dt
	end
end
"""
initialize_times(r,n,nalleles)

All alleles are initialized to state 0
"""
function initialize_times(r,n,nalleles)
	m = 0
	tau = fill(Inf,2*n+2,nalleles)
	for n in 1:nalleles
		tau[1,n] = -log(rand())/r[1]
		tau[end,n] = -log(rand())/(m*r[end])
	end
	return tau
end
function initialize_times(r,n,nr,nalleles)
	tau = fill(Inf,2*n+nr+2,nalleles)
	for n in 1:nalleles
		tau[1,n] = -log(rand())/r[1]
	end
	return tau,zeros(Int,max(nr,1),nalleles)
end
function initialize(r,G,R,nreactions,nalleles,initstate=1,initreaction=1)
	tau = fill(Inf,nreactions,nalleles)
	for n in 1:nalleles
		tau[initreaction,n] = -log(rand())/r[1]
	end
	states = zeros(Int,G+max(R,1),nalleles)
	states[initstate] = 1
	return tau,states
end

function activateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
	transitionG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
	if R == 0
		state[G+1,allele] = 2
	elseif state[G+1,allele] > 1
		tau[downstream[end],allele] = Inf
	end
end

function deactivateG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
	transitionG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
	if R == 0
		state[G+1,allele] = 0
	end
end

function transitionG!(tau,state,index,t,m,r,allele,G,R,upstream,downstream,initial,final)
	tau[index,allele] = Inf
	state[final,allele] = 1
	state[initial,allele] = 0
	for d in downstream
		tau[d,allele] = -log(rand())/r[d] + t
	end
	for u in upstream
		tau[u,allele] = Inf
	end
end

function initiate!(tau,state,index,t,m,r,allele,G,R,S,downstream)
	tau[index,allele] = Inf
	state[G+1,allele] = 2
	if R > 1
		tau[downstream[1],allele] = -log(rand())/(r[downstream[1]])+ t
	end
	if R < 2 || state[G+2,allele] < 1
		tau[index+1,allele] =  -log(rand())/(r[index+1])+ t
	end
	if S > 0


	end
end
function transitionR!(tau,state,index,t,m,r,allele,G,R,S,upstream,downstream,initial,final)
	tau[index,allele] = Inf
	if S > 0 && isfinite(tau[index+R+1])
		tau[index+R+2] = -log(rand())/r[index+R+2] + t
		tau[index+R+1] = Inf
	end
	for d in downstream
		if state[d,allele] > 0
			tau[d,allele] = Inf
		else
			tau[d,allele] = -log(rand())/r[index+1] + t
		end
	end
	for u in upstream
		if state[u,allele] > 0
			tau[u,allele] = -log(rand())/r[index+1] + t
		else
			tau[u,allele] = Inf
		end
	end
	state[final,allele] = state[initial,allele]
	state[initial,allele] = 0
end
function eject!(tau,state,index,t,m,r,allele,G,R,S,upstream,downstream,initial,final)
	m += 1
	set_decay!(tau,downstream[end],t,m,r)
	if S > 0
		# tau[splice] = Inf
	end
	if R > 0
		tau[index,allele] = Inf
		state[G+R,allele] = 0
		if state[G+R-1,allele] > 0
			tau[upstream[1],allele] = -log(rand())/(r[upstream[1]])+ t
		end
	else
		tau[index,allele] = -log(rand())/(r[index]) + t
	end
	m
end
# function eject!(tau,state,reaction,t,m,r,n,nr,nallele)
# 				m += 1
# 				tau[reaction,nallele] = Inf
# 				state[nr,nallele] = 0
# 				if nr == 1
# 					set_initiate!(tau,t,r,n,nallele)
# 				elseif state[nr-1,nallele]
# 					tau[reaction-1,nallele] = -log(rand())/(r[reaction-1])+ t
# 				end
# 				set_decay!(tau,reaction,t,m,r)
# end

function splice!(tau,state,index,t,m,r,allele,G,R,initial)
	state[initial,allele] = 1
	tau[index,allele] = Inf
end
function decay!(tau,state,index,t,m,r)
	m -= 1
	tau[index,1] = -log(rand())/(m*r[index]) + t
	m
end
# function decay!(tau,state,reaction::Reaction,t,m,r,allele,G,R)
# 	m -= 1
# 	tau[reaction.index,1] = -log(rand())/(m*r[reaction.index]) + t
# 	m
# end
"""
set_initiate!(tau,t,r,n,nallele)

update tau matrix for initiation rate

"""
function set_initiate!(tau,t,r,n,nallele)
	if isfinite(tau[2*n,nallele])
		tau[2*n+1,nallele] = -log(rand())/(r[2*n+1])+ t
	end
end
"""
set_decay!(tau,reaction,t,m)
update tau matrix for decay rate

"""
function set_decay!(tau,index,t,m,r)
	if m == 1
		tau[index,1] = -log(rand())/r[index]+ t
	else
		tau[index,1] = (m-1)/m*(tau[index,1] - t) + t
	end
end
function set_decay(rdecay,t,m)
	if m == 1
		return -log(rand())*(taudecay)+ t
	else
		return (m-1)/m*(taudecay - t) + t
	end
end
"""
simulatorGM(r::Vector,n::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)

Modified Gibson and Bruch next reaction algorithm to simulate G state telegraph model (classic telegraph is 2 state)

r = vector of rates in order of state transitions (forward, backward alternating), creation rate, decay rate
n  = number of states - 1, labeled 0,1,...,n
nhist = number of mRNA bins in histogram
total = maximum number of simulator steps
tolerance in histogram change

tau = matrix that keeps track of the putative reaction times of each transition, rows are transitions and columns are alleles
t = time of next reaction and reaction[1] is which reaction and reaction[2] is which allele

"""
function simulatorGM(r::Vector,n::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)

	tau = initialize_times(r,n,nalleles)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)

	while err > tol && steps < total
		steps += 1
		t,reaction = findmin(tau)
		dt = t-t0
		t0 = t

		update_mhist!(mhist,m,dt,nhist)

		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
		end
		if reaction[1] <= 2*n
			if isodd(reaction[1])
				gforward!(tau,reaction[1],t,r,reaction[2])
			else
				greverse!(tau,reaction[1],t,r,reaction[2])
			end
		else
			if reaction[1] == 2*n + 1
				m = eject!(tau,reaction[1],t,m,r,reaction[2])
			else
				m = decay!(tau,reaction[1],t,m,r)
			end
		end
	end  # while

	counts = max(sum(mhist),1)
	mhist /= counts

	if count
		return mhist[1:nhist],counts,steps,err
	else
		return mhist[1:nhist]
	end

end


function simulatorGRM(r::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)
	tau,state = initialize_times(r,n,nr,nalleles)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
	while err > tol && steps < total
		steps += 1
		t,reaction = findmin(tau)
		dt = t-t0
		t0 = t
		update_mhist!(mhist,m,dt,nhist)
		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
		end
		if reaction[1] <= 2*n
			if isodd(reaction[1])
				gforward!(tau,state,reaction[1],t,r,n,reaction[2])
			else
				greverse!(tau,reaction[1],t,r,reaction[2])
			end
		else
			if reaction[1] <= 2*n + nr
				if reaction[1] == 2*n + 1
					initiate!(tau,state,reaction[1],t,r,nr,reaction[2])
				else
					rstep!(tau,state,reaction[1],t,r,n,nr,reaction[2])
				end
			else
				if reaction[1] == 2*n + nr + 1
					m = eject!(tau,state,reaction[1],t,m,r,n,nr,reaction[2])
				else
					m = decay!(tau,reaction[1],t,m,r)
				end
			end
		end
	end  # while
	counts = max(sum(mhist),1)
	mhist /= counts
	if count
		return mhist[1:nhist],counts,steps,err
	else
		return mhist[1:nhist]
	end
end


"""
gforward!(tau,reaction,t,r,nallele)

update tau matrix for forward transition
"""
function gforward!(tau,reaction,t,r,nallele)
		tau[reaction,nallele] = Inf
		tau[reaction+1,nallele] = -log(rand())/r[reaction+1] + t
		if reaction > 1
			tau[reaction-1,nallele] = Inf
		end
		tau[reaction+2,nallele] = -log(rand())/r[reaction+2] + t
		nothing
end

function gforward!(tau,state,reaction,t,r,n,nallele)
		tau[reaction,nallele] = Inf
		tau[reaction+1,nallele] = -log(rand())/r[reaction+1] + t
		if reaction > 1
			tau[reaction-1,nallele] = Inf
		end
		if reaction == 2*n-1 && state[1,nallele]
			tau[2*n+1,nallele] = Inf
		else
			tau[reaction+2,nallele] = -log(rand())/r[reaction+2] + t
		end
		nothing
end

"""
	greverse!(tau,reaction,t,r,nallele)

	update tau matrix for reverse state transition
"""
function greverse!(tau,reaction,t,r,nallele)
			tau[reaction,nallele] = Inf
			tau[reaction+1,nallele] = Inf
			tau[reaction-1,nallele]= -log(rand())/r[reaction-1] + t
			if reaction > 2
				tau[reaction-2,nallele] = -log(rand())/r[reaction-2] + t
			end
			nothing
end

"""
		eject!(tau,reaction,t,m,r,nallele)

		update tau matrix for mRNA ejection

"""
function eject!(tau,reaction,t,m,r,nallele)
				m += 1
				tau[reaction,nallele] = -log(rand())/r[reaction] + t
				set_decay!(tau,reaction,t,m,r)
end

function eject!(tau,state,reaction,t,m,r,n,nr,nallele)
				m += 1
				tau[reaction,nallele] = Inf
				state[nr,nallele] = 0
				if nr == 1
					set_initiate!(tau,t,r,n,nallele)
				elseif state[nr-1,nallele]
					tau[reaction-1,nallele] = -log(rand())/(r[reaction-1])+ t
				end
				set_decay!(tau,reaction,t,m,r)
end



"""
decay!(tau,reaction,t,m,r)

update tau matrix for decay transition

"""
function decay!(tau,reaction,t,m,r)
	m -= 1
	tau[reaction,1] = -log(rand())/(m*r[reaction]) + t
	m
end

function decay(t,m,r)
	m -= 1
	# tau[reaction,1] = -log(rand())/(m*r[reaction]) + t
	m,-log(rand())/(m*r[end]) + t
end

"""
initiate!(tau,state,reaction,t,r,nr,nallele)

update tau matrix for initiation reaction from G to R

"""
function initiate!(tau,state,reaction,t,r,nr,nallele)
	tau[reaction,nallele] = Inf
	state[1,nallele] = 1
	if nr < 2 || ~state[2,nallele]
		tau[reaction+1,nallele] =  -log(rand())/(r[reaction+1])+ t
	end
end
"""
rstep!(tau,state,reaction,t,r,n,nr,nallele)

update tau matrix for R step transition

"""
# function rstep!(tau,state,reaction::Int,t,r,n,nr,nallele)
# 	step = reaction - 2*n - 1   #current Rstep
# 	tau[reaction,nallele] = Inf
# 	state[step,nallele] = 0
# 	state[step+1,nallele] = 1
# 	if  step == nr-1 || ~state[step+2,nallele]
# 		tau[reaction+1,nallele] =  -log(rand())/(r[reaction+1])+ t
# 	end
# 	if step > 1 && state[step-1,nallele]
# 		tau[reaction-1,nallele] =  -log(rand())/(r[reaction-1])+ t
# 	elseif step == 1
# 		set_initiate!(tau,t,r,n,nallele)
# 	end
# end
"""
splice!(tau,state,reaction,t,r,n,nr,nallele)

update tau matrix for splice transition

"""
# function splice!(tau,state,reaction,t,r,n,nr,nallele)
# 	tau[reaction,nallele] = Inf
# end
