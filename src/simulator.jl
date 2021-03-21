"""
simulator(r::Vector,n::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)

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

function simulator1(r::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)
	tau,taumin,reactmin,Rstep,taudecay = initialize(r,n,nr,nalleles)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
	while err > tol && steps < total
		steps += 1
		t,allele = findmin(taumin)
		reaction = reactmin[allele]
		if taudecay < t
			t = taudecay
			dt = t-t0
			update_mhist!(mhist,m,dt,nhist)
			m,taudecay = decay(t,m,r)
		else
			dt = t-t0
			update_mhist!(mhist,m,dt,nhist)
			if reaction <= 2*n
				if isodd(reaction)
					gforward!(tau,Rstep,reaction,t,r,n,allele)
				else
					greverse!(tau,reaction,t,r,allele)
				end
			else
				if reaction <= 2*n + nr
					if reaction == 2*n + 1
						initiate!(tau,Rstep,reaction,t,r,nr,allele)
					else
						rstep!(tau,Rstep,reaction,t,r,n,nr,allele)
					end
				else
					m,taudecay = eject1!(tau,Rstep,reaction,taudecay,t,m,r,n,nr,allele)
				end
			end
			taumin[allele],reactmin[allele] = findmin(tau[:,allele])
		end
		t0 = t
		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
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

function simulator2(r::Vector,n::Int,nr::Int,nhist::Int,nalleles::Int,total::Int=10000000,tol::Float64=1e-6,count=false)
	tau,taumin,reactmin,Rstep,taudecay = initialize(r,n,nr,nalleles)
	mhist,mhist0,m,steps,t,ts,t0,tsample,err = initialize_sim(r,nhist,tol)
	while err > tol && steps < total
		steps += 1
		t,allele = findmin(taumin)
		reaction = reactmin[allele]
		if taudecay < t
			t = taudecay
			dt = t-t0
			update_mhist!(mhist,m,dt,nhist)
			m,taudecay = decay(t,m,r)
		else
			dt = t-t0
			update_mhist!(mhist,m,dt,nhist)
			if reaction <= 2*n
				greaction!(tau,Rstep,reaction,t,r,n,allele)
			else
				if reaction <= 2*n + nr
					rreaction!(tau,Rstep,reaction,t,r,n,nr,allele)
				else
					m,taudecay = eject1!(tau,Rstep,reaction,taudecay,t,m,r,n,nr,allele)
				end
			end
			taumin[allele],reactmin[allele] = findmin(tau[:,allele])
		end
		t0 = t
		if t-ts > tsample
			err,mhist0 = update_error(mhist,mhist0)
			ts = t
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
technical_loss(m,yieldfactor)

Accounts for loss by sampling from a Binomial distribution with probability yieldfactor
"""
function technical_loss(m,yield)
	d = Binomial(m,clamp(yield,0.,1.))
	rand(d)
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
	tau = fill(Inf,2*n+2,nalleles)
	for n in 1:nalleles
		tau[1,n] = -log(rand())/r[1]
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
	set_decay!(tau,reaction,t,m)
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
	set_decay!(tau,reaction,t,m)
end
function eject1!(tau,state,reaction,taudecay,t,m,r,n,nr,nallele)
	m += 1
	tau[reaction,nallele] = Inf
	state[nr,nallele] = 0
	if nr == 1
		set_initiate!(tau,t,r,n,nallele)
	elseif state[nr-1,nallele]
		tau[reaction-1,nallele] = -log(rand())/(r[reaction-1])+ t
	end
	m,set_decay(taudecay,t,m)
end

"""
set_decay!(tau,reaction,t,m)

update tau matrix for decay rate

"""

function set_decay!(tau,reaction,t,m)
	if m == 1
		tau[reaction+1,1] = -log(rand())/(r[reaction+1])+ t
	else
		tau[reaction+1,1] = (m-1)/m*(tau[reaction+1] - t) + t
	end
	m
end
function set_decay(taudecay,t,m)
	if m == 1
		return -log(rand())/(r[end])+ t
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
