# simulator.jl
# Functions to simulate Markov gene transcription models


struct Reaction
	initial::String
	final::String
	action::Function
end

function set_reactions(Gtransitions,G,R,S)
	reaction = Vector{Vector}(undef,0)
	for g in Gtransitions
		if g[2] == G
			push!(reaction,Reaction(g[1],g[2],activateG!)
		elseif g[1] == G
			push!(reaction,Reaction(g[1],g[2],deactivateG!)
		else
			push!(reaction,Reaction(g[1],g[2],transitionG!)
		end
	end
	for r in 1:R-1
		push!(reaction,Reaction(r,r+1,transitionR!)
	end
	push!(reaction,Reaction(R,0,transitionR!)
	for s in 1:S
		push!(reaction,Reaction(s,0,transitionR!)
	end
	push!(reaction,push!(reaction,Reaction(0,0,decay!)
	return reaction
end


function activateG!(tau,state,reaction,upstream,downstream,t,m,r,G,R,allele)
	transitionG!(tau,state,reaction,upstream,downstream,t,m,r,allele)
	if R == 0
		state[1,allele] = 2
		tau[eject,allele] = -log(rand())/r[eject] + t
	elseif state[1,allele] == 0
		tau[initiation,allele] = -log(rand())/r[initiate] + t
	end
end
function deactivateG!(tau,state,reaction,upstream,downstream,t,m,r,G,R,allele)
	transitionG!(tau,state,reaction,upstream,downstream,t,m,r,allele)
	tau[reaction+1,allele] = Inf
	if R == 1
		state[,allele] = 0
	end
end

function transitionG!(tau,state,reaction,upstream,downstream,t,m,r,G,R,allele)
	tau[reaction,allele] = Inf
	for u in upstream
		tau[u,allele] = Inf
	end
	for d in downstream
		tau[d,allele] = -log(rand())/r[d] + t
	end
end

function transitionR!(tau,state,reaction,upstream,downstream,t,m,r,G,R,allele)
	tau[reaction,nallele] = Inf
	step = length(upstream) + 1
	state[step,nallele] = 0
	state[step+1,nallele] = 2
	for u in upstream
		if step > 1 && state[step-1,nallele]
			tau[reaction-1,nallele] =  -log(rand())/(r[reaction-1])+ t
	end
	if isempty(downstream)
		tau[reaction+1,nallele] =  -log(rand())/(r[reaction+1])+ t
	else
		for d in downstream
		if  step == R-1 || ~state[step+2,nallele]
			tau[reaction+1,nallele] =  -log(rand())/(r[reaction+1])+ t
		end
	end
end
function eject!(tau,state,reaction,upstream,downstream,t,m,r,G,R,allele)
	m += 1
	tau[reaction,nallele] = Inf
	if R = 0
		state[1,nallele] = 0
	else
		state[R,nallele] = 0
	end
	if R == 1
		set_initiate!(tau,t,r,n,nallele)
	elseif state[R-1,nallele]
		tau[reaction-1,nallele] = -log(rand())/(r[reaction-1])+ t
	end
	set_decay!(tau,reaction,t,m,r)
end
function splice!(occupancy,reaction,allele)
	occupancy[allele][reaction[1]] -= 1
end

function decay!(tau,reaction,t,m,r)
	if m == 1
		tau[reaction+1,1] = -log(rand())/(r[reaction+1])+ t
	else
		tau[reaction+1,1] = (m-1)/m*(tau[reaction+1] - t) + t
	end
	m
end



function set_reactions(Gtransitions,G,R,S)
	reaction = Vector{Vector}(undef,0)
	for g in Gtransitions
		push!(reaction,"G" .* string.(g))
	end
	if R == 0
		push!(reaction,["G" * string(G),"M"])
	else
		push!(reaction,["G" * string(G),"R1"])
		for r in 2:R
			push!(reaction,["R" * string(r-1),"R" * string(r)])
		end
		push!(reaction,["R" * string(R),"M"])
	end
	for s in 1:S
		push!(reaction,["R" * string(s),"L"])
	end
	push!(reaction,["M","M-"])
	return reaction
end

function set_reactants(G,R)
	reactants = String[]
	for r in 1:G
		push!(reactants,"G" * string(r))
	end
	push!(reactants,"R" * "1")
	for r in 2:R
		push!(reactants,"R" * string(r))
	end
	push!(reactants,"M")
	reactants
end

function affects(reactions)
	affected = Vector{Vector}(undef,0)
	for r in reactions
		v = Int[]
		# if r[1] == "M"
		# 	push!(v,length(reactions))
		# else
		if r[1] != "M"
			for s in eachindex(reactions)
				if r != reactions[s] && (r[1] == reactions[s][1] || r[1] == reactions[s][2] ||  r[2] == reactions[s][1] || r[2] == reactions[s][2])
					push!(v,s)
				end
			end
		end
		# end
		push!(affected,v)
	end
	return affected
end

function set_rules(reactions)
	rules = Vector{Function}(undef,0)
	for r in reactions
		if occursin("G",r[2])
			push!(rules,slide!)
		elseif occursin("R1",r[2])
			push!(rules,initiate!)
		elseif occursin("R",r[1]) && occursin("R",r[2])
			push!(rules,slide!)
		elseif r[2] == "M"
			push!(rules,eject!)
		elseif occursin("L",r[2])
			push!(rules,splice!)
		elseif r[2] == "M-"
			push!(rules,decay!)
		end
	end
	rules
end

function set_rulesG(reactions,active)
	rules = Vector{Function}(undef,0)
	for r in reactions
		if occursin("G",r[2]) && r[1] != "G" * string(active) && r[2] != "G" * string(active)
			push!(rules,slide!)
		elseif occursin("G" * string(active),r[1]) && occursin("G",r[2])
				push!(rules,deactivateG!)
		elseif occursin("G" * string(active),r[2]) && occursin("G",r[1])
				push!(rules,activateG!)
		elseif r[2] == "M"
			push!(rules,birth!)
		elseif occursin("S",r[1])
			push!(rules,splice!)
		elseif r[2] == "M-"
			push!(rules,decay!)
		end
	end
	rules
end



function slide!(occupancy,reaction,allele)
	occupancy[allele][reaction[2]] = occupancy[allele][reaction[1]]
	occupancy[allele][reaction[1]] = 0
end

function eject!(occupancy,reaction,allele)
	occupancy[allele][reaction[1]] = 0
	birth!(occupancy,reaction,allele)
end

function birth!(occupancy,reaction,allele)
	occupancy[1]["M"] += 1
end

function initiate!(occupancy,reaction,allele)
	occupancy[allele][reaction[2]] = 2
end

function activateG!(occupancy,reaction,allele)
	slide!(occupancy,reaction,allele)
	occupancy[allele]["R1"] = 2
end

function deactivateG!(occupancy,reaction,allele)
	slide!(occupancy,reaction,allele)
	occupancy[allele]["R1"] = 0
end

function decay!(occupancy,reaction,allele)
	occupancy[1]["M"] -= 1
end

function splice!(occupancy,reaction,allele)
	occupancy[allele][reaction[1]] -= 1
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


function initial_occupancy(reactants,nalleles,start)
	o = Dict[]
	for i in 1:nalleles
		d = Dict(i => 0 for i in reactants)
		d["G1"] = start
		push!(o,d)
	end
	o
end

function initial_state(r,reactions,decay,reactants,nalleles,start = 1)
	occupancy = initial_occupancy(reactants,nalleles,start)
	tau = fill(Inf,length(reactions),nalleles)
	for a in 1:length(reactions), allele in 1:nalleles
		set_tau!(tau,occupancy,r,0,a,allele,reactions,decay)
	end
	return tau,occupancy
end
function set_decay(taudecay,decayrate,t,m)
	if m == 0
		return Inf
	elseif m == 1
		return -log(rand())/decayrate+ t
	else
		return (m-1)/m*(taudecay - t) + t
	end
end

function set_tau!(tau,m::Int,r,t,reaction,reactionindex,allele,decay)
	if reaction[1] == "M"
		tau[decay,1] = -log(rand())/(r[decay]*m) + t
	elseif reaction[2] == "M"
		tau[reactionindex,allele] = -log(rand())/r[reactionindex] + t
		tau[decay,1] = set_decay(tau[decay,1],r[decay],t,m)
	elseif reaction[2] != "L"
		tau[reactionindex,allele] = Inf
	end
end

function set_tau!(tau,occupancy::Vector,r,t,a,allele,reactions,decay)
	m = occupancy[1]["M"]
	if a == decay
		tau[a,1] = -log(rand())/(r[a]*m) + t
	elseif occupancy[allele][reactions[a][1]] > 0 && (occupancy[allele][reactions[a][2]] < 1 || reactions[a][2] == "L" || reactions[a][2] == "M")
		tau[a,allele] = -log(rand())/(r[a]) + t
	else
		tau[a,allele] = Inf
	end
end

function update!(tau,occupancy,r,t,reaction,allele,reactions,decay,affect,rule!)
	rule!(occupancy,reactions[reaction],allele)
	set_tau!(tau,occupancy[1]["M"],r,t,reactions[reaction],reaction,allele,decay)
	for a in affect
		set_tau!(tau,occupancy,r,t,a,allele,reactions,decay)
	end
end

initialize_sim(r,nhist,tol,samplefactor=20.,errfactor=10.) = zeros(nhist+1),ones(nhist+1),0,0,0.,0.,0.,samplefactor/minimum(r),errfactor*tol

update_error(mhist,mhist0) = (norm(mhist/sum(mhist)-mhist0/sum(mhist0),Inf),copy(mhist))

function update_mhist!(mhist,m,dt,nhist)
	if 0 < m + 1 <= nhist
		mhist[m+1] += dt
	else
		mhist[nhist+1] += dt
	end
end

function simulator(r::Vector,Gtransitions::Tuple,G,R,S,nalleles,nhist::Int,total::Int=10000000,tol::Float64=1e-6,count=false)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
	reactants = set_reactants(G,R)
	reactions = set_reactions(Gtransitions,G,R,S)
	println(reactants)
	println(reactions)
	decayind = length(reactions)
	if R == 0
		rules = set_rulesG(reactions,G)
	else
		rules = set_rules(reactions)
	end
	affect = affects(reactions)
	tau,occupancy = initial_state(r,reactions,decayind,reactants,nalleles)
	println(rules)
	println(affect)
	println(occupancy)
	println(tau)
	while err > tol && steps < total
		steps += 1
		t,reaction = findmin(tau)
		dt = t-t0
		t0 = t
		update_mhist!(mhist,occupancy[1]["M"],dt,nhist)
		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
		end
		update!(tau,occupancy,r,t,reaction[1],reaction[2],reactions,decayind,affect[reaction[1]],rules[reaction[1]])
		# println(reaction)
		# println(occupancy)
		# println(tau," ",t)
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

# function decay(t,m,r)
# 	m -= 1
# 	# tau[reaction,1] = -log(rand())/(m*r[reaction]) + t
# 	m,-log(rand())/(m*r[end]) + t
# end

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
