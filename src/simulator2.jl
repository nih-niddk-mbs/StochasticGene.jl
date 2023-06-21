# simulator.jl
# Functions to simulate Markov gene transcription models

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

Inputs:  G,R,S, Gtransitions

Initialize
1. index all transitions reactions
2. get initial and final state of each reaction
3. specify state Gstate = int, Rstates = vector of length R with entries 0, 1, 2
4. find upstream and downstream reactions for each reaction
5. Specify initial state
6. generate tau for each reaction


Simulation

Find reaction with minimum tau
find initial and final states
update state
update propensity
Set tau of upstream reactions
Set tau of downstream reactions
determine if start or end of ON state
update rna histogram


Reaction:  a -> b

upstream   x -> a
downstream b -> y

if x is occupied then x -> a has finite tau
if y is occupied then b -> y has infinite tau


need to go from init and final states to reaction index

initial[reaction]  : gives initial state of reaction
final[reaction]:
reactions_initial[state] : gives list of reactions with initial state state
reactions_final[state]: list of reactions with final state state


if pick reaction a -> b

then move state from a to b

find upstream and downstream reactions
find initial state of upstream reactions
find final state of downstream reactions



reactants:  G1,G2,G3...,R1,R2,R3...,S1,S2,S3...,M
reactions: [G1,G2],[G2,G3],[G3,R1],[R1,R2],[S1,L],...,[R3,M],[M,M-]



if propensity[reaction] == 0
	tau = Inf
else
	tau = -log(rand())/r[reaction]
end

tau[reaction]

swap
update
build
reactants
products
init
final
dependson[tau[index]]
affects[index]

Rstate = Array{Array{Int,1},1}(undef,nalleles)


function affects(reactions)
	affected = Int[]
	for r in reactions
		for s in eachindex(reactions)
			if r[1] == reaction[s][1] || r[1] == reactions[s][2] ||  r[2] == reaction[s][1] || r[2] == reactions[s][2]
				push!(affected,s)
			end
		end
	return affected
end

function status(reactions,states)
	for r in eachindex(reactions)
		for s in states
		if occupied(reactions[r][1])
			initial[r] = 1
		end
		if ~occupied(reactions[r][2])
			final[r] = 1
		end
	end
end

function update_occupancy!(occupancy,reaction)
	occupancy[reaction[1]] -= 1
	occupancy[reaction[2]] += 1
	nothing
end

function set_reactions(Gtransitions,G,R,S)
	reaction = Vector{Vector}(undef,0)
	for g in Gtransitions
		push!(reaction,"G" .* string.(g))
	end
	if R == 0
		push!(reaction,["G" * string(G),"m"])
	else
		push!(reaction,["G" * string(G),"R1"])
		for r in 2:R
			push!(reaction,["R" * string(r-1),"R" * string(r)])
		end
		push!(reaction,["R" * string(R),"M"])
		for s in 1:S
			push!(reaction,["S" * string(s),"L"])
		end
	end
	push!(reaction,["M","M-"])
	return reaction
end

function make_reactants(G,R,S)
	reactants = String[]
	for r in 1:G
		push!(reactants,"G" * string(r))
	end
	for r in 1:R
		push!(reactants,"R" * string(r))
	end
	for s in 1:S
		push!(reactants,"S" * string(s))
	end
	push!(reactants,"M")
	push!(reactants,"M-")
	reactants
end

function initial_occupancy(reactants)
	o = Dict(i => 0 for i in reactants)
	o["G1"] = 1
	o
end

struct Reaction
	type
	initstate
	finalstate
end


function updateRstate!(state,a,b,nallele)
	state[b,nallele] = state[a,nallele]
	state[a,nallele] = 0
end

function splice!(state,a,nallele)
	state[a,nallele] = 1
end


function update!(tau,propensity,t,reaction,upstream,downstream)
	for i in df
		tau[i,nallele] = -log(rand())/(r[i]*initial[i]*final[i]) + t
	end
	for i in di
		tau[i,nallele] = Inf
	end
	for i in ui
		tau[i,nallele] = -log(rand())/(r[i]*propensity[i]) + t
	end
	for i in uf
		tau[i,nallele] = Inf
	end
	updatepropensity
end

function update_propensity!(propensity,oldstate,newstate,downstream,upstream)
	if occursin("m",newstate)
		if newstate == "m"
			propensity[decay reaction,nallele] += 1
		elseif newstate == "m-1"
			protensity[decay reaction,nallele ] -= 1
		end
	else
		for i in ui
			propensity[i,nallele] = 1
		end
		for i in uf
			propensity[i,nallele] = 0
		end
	end
	propensity[reactions ending with old state ] = 0

	end
end

function update!(tau,state,t,reaction,nallele,action::Rreaction)
	for d in action.downstream
		if state[action.next] == 0
			tau[d,nallele] = -log(rand())/r[reaction+1] + t
		end
	end
	for u in action.upstream
		tau[u,nallele] = Inf
	end
	update_occpancy!(occupancy,a,b,nallele)
end



function update!(tau,state,t,reaction,nallele,action::Sreaction)
	tau[reaction,nallele] = Inf
	state[reaction.a,nallele] = 1
end

function update!(tau,state,t,reaction,nallele,action::Eject)
	for d in action.downstream
		if state[action.next] == 0
			tau[d,nallele] = -log(rand())/r[reaction+1] + t
		end
	end
	update

end

function update!(tau,state,t,reaction,nallele,action::Decay)
	for d in action.downstream
		if state[action.next] == 0
			tau[d,nallele] = -log(rand())/r[reaction+1] + t
		end
	end
	if eject(reaction)
		update_decay
	end
	for u in action.upstream
		tau[u,nallele] = Inf
	end
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




end

G = 2
R = 3
S = 3

Gtransitions = ([1,2],[2,3],[3,2],[3,1])
Rtransitions = ([0,1],[1,2],[2,3],[3,0])
Stransitions = ([1,0],[2,0],[3,0])
Mtransitions = ([m,m+1],[m,m-1])


ntransitions = length(Gtransitions) + length(Rtransitions) + length(Stransitions) + length(Mtransitions)

struct Transition
	index
	a
	b
	parents
	children
end

Gstate = 1
Rstate = [1,0,1]
State = [0,0,1]



function simulator(state,Gtransitions,G,R,S,nalleles,total,tol)
function simulatorGM(r::Vector,n::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)


		tau = initialize_times(r,n,nalleles)
		mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
		mhist,mhist0,state,action

		while err > tol && steps < total
			steps += 1
			t,reaction = findmin(tau)
			dt = t-t0
			t0 = t

			burs()
			update_mhist!(mhist,m,dt,nhist)
			update!(tau,reaction,state,upstream[reaction],downstream[reaction])

			if t-ts > tsample
				err,mhist0 = update_error(mhist,mhist0)
				ts = t
			end

		end
			# if reaction[1] <= 2*n
			# 	if isodd(reaction[1])
			# 		gforward!(tau,reaction[1],t,r,reaction[2])
			# 	else
			# 		greverse!(tau,reaction[1],t,r,reaction[2])
			# 	end
			# else
			# 	if reaction[1] == 2*n + 1
			# 		m = eject!(tau,reaction[1],t,m,r,reaction[2])
			# 	else
			# 		m = decay!(tau,reaction[1],t,m,r)
			# 	end
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

function update_tau!(tau,state,reaction,t,r)
	tau[reaction] = Inf
	tau[reaction[1]+1,reaction[2]] = -log(rand())/r[index(reaction[1])] + t

	affected_reactions(reaction)




end

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
	tau,Rstep = initialize_times(r,n,nr,nalleles)
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
				gforward!(tau,Rstep,reaction[1],t,r,n,reaction[2])
			else
				greverse!(tau,reaction[1],t,r,reaction[2])
			end
		else
			if reaction[1] <= 2*n + nr
				if reaction[1] == 2*n + 1
					initiate!(tau,Rstep,reaction[1],t,r,nr,reaction[2])
				else
					rstep!(tau,Rstep,reaction[1],t,r,n,nr,reaction[2])
				end
			else
				if reaction[1] == 2*n + nr + 1
					m = eject!(tau,Rstep,reaction[1],t,m,r,n,nr,reaction[2])
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

function greaction!(tau,Rstep,reaction,t,r,n,allele)
	if isodd(reaction)
		gforward!(tau,Rstep,reaction,t,r,n,allele)
	else
		greverse!(tau,reaction,t,r,allele)
	end
end

function rreaction!(tau,Rstep,reaction,t,r,n,nr,allele)
	if reaction == 2*n + 1
		initiate!(tau,Rstep,reaction,t,r,nr,allele)
	else
		rstep!(tau,Rstep,reaction,t,r,n,nr,allele)
	end
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
	return tau,zeros(Bool,nr,nalleles)
end
function initialize(r,n,nr,nalleles)
tau = fill(Inf,2*n+nr+1,nalleles)
for n in 1:nalleles
	tau[1,n] = -log(rand())/r[1]
end
tau,tau[1,:],ones(Int,nalleles),zeros(Bool,nr,nalleles),Inf
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
# function eject1!(tau,state,reaction,taudecay,t,m,r,n,nr,nallele)
# 	m += 1
# 	tau[reaction,nallele] = Inf
# 	state[nr,nallele] = 0
# 	if nr == 1
# 		set_initiate!(tau,t,r,n,nallele)
# 	elseif state[nr-1,nallele]
# 		tau[reaction-1,nallele] = -log(rand())/(r[reaction-1])+ t
# 	end
# 	m,set_decay(taudecay,t,m)
# end

"""
set_decay!(tau,reaction,t,m)

update tau matrix for decay rate

"""

function set_decay!(tau,reaction,t,m,r)
	if m == 1
		tau[reaction+1,1] = -log(rand())/(r[reaction+1])+ t
	else
		tau[reaction+1,1] = (m-1)/m*(tau[reaction+1] - t) + t
	end
	m
end
function set_decay(taudecay,t,m)
	if m == 1
		return -log(rand())/(taudecay)+ t
	else
		return (m-1)/m*(taudecay - t) + t
	end
end
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
function rstep!(tau,state,reaction,t,r,n,nr,nallele)
	step = reaction - 2*n - 1   #current Rstep
	tau[reaction,nallele] = Inf
	state[step,nallele] = 0
	state[step+1,nallele] = 1
	if  step == nr-1 || ~state[step+2,nallele]
		tau[reaction+1,nallele] =  -log(rand())/(r[reaction+1])+ t
	end
	if step > 1 && state[step-1,nallele]
		tau[reaction-1,nallele] =  -log(rand())/(r[reaction-1])+ t
	elseif step == 1
		set_initiate!(tau,t,r,n,nallele)
	end
end
"""
splice!(tau,state,reaction,t,r,n,nr,nallele)

update tau matrix for splice transition

"""
function splice!(tau,state,reaction,t,r,n,nr,nallele)
	tau[reaction,nallele] = Inf
end
