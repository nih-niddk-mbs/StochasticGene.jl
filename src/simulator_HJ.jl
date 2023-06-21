using Distributions
using Random
using DelimitedFiles

"""
simulator(r::Vector,G::Int,nalleles::Int,yield,final,minitial::Vector)

Modified Gibson and Bruch next reaction algorithm to simulate G state telegraph model (classic telegraph is 2 state)

r = vector of rates in order of state transitions (forward, backward alternating), creation rate, decay rate
G = number of states, labeled 0,1,...,G-1
yield = probability of capturing mRNA in scRNA experiment
final = output time
minitial = initial mRNA number

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
				m = create!(tau,reaction[1],t,m,r,reaction[2])
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

# rng = MersenneTwister(31415926);
function simulator(r::Vector,G::Int,nalleles::Int,yield,final,minitial::Vector)
    mlist = [0.0] # array for mRNA count
    tlist = [0.0] # array for time points
    Nlist = [0.0] # array for the cell number

    i = 1
    N = length(minitial)
    while i <= N
        m = floor(minitial[i])
        #m = 0.
        t = 0.
        push!(tlist, t)
        push!(mlist, m)
        push!(Nlist, i)

        n = G-1
        tau = initialize_times(r,n,nalleles)
        while t < final
            t,reaction = findmin(tau)
            if reaction[1] <= 2*n
                if isodd(reaction[1])
                    gforward!(tau,reaction[1],t,r,reaction[2])
                else
                    greverse!(tau,reaction[1],t,r,reaction[2])
                end
            else
                if reaction[1] == 2*n + 1
                    m = create!(tau,reaction[1],t,m,r,reaction[2])
                else
                    m = decay!(tau,reaction[1],t,m,r)
                end
            end
            push!(tlist, t)
            push!(mlist, m)
            push!(Nlist, i)
        end # while
        technical_loss(m,yield)
        i += 1
    end # while
    mtlist = [mlist Nlist tlist]
    return mtlist
end # function

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
		tau[1,n] = -log(rand(rng))/r[1]
	end
	return tau
end

"""
initialize_times(r,n,nalleles)

All alleles are initialized to equilibrium distribution
"""
function initialize_times(r,n,nalleles)
	tau = fill(Inf,2*n+2,nalleles)
	for n in 1:nalleles
		d = Categorical(p_G(r,n))
		state = rand(d)
		forward = 2*state - 1
		reverse = 2*state - 2
		tau[forward,n] = -log(rand(rng))/r[forward]
		if state > 1
			tau[reverse,n] = -log(rand(rng))/r[reverse]
		end
	end
	return tau
end
"""
p_G(r::Vector,n::Int)

Compute state probability for equilibrium

"""
function p_G(r::Vector,n::Int)
    Gss = Array{Float64,2}(undef,1,n+1)
    Gss[1,1] = 1.
    for k in 1:n
        Gss[1,k+1] = Gss[1,k]*r[2*k-1]/r[2*k]
    end
    Gss ./= sum(Gss)
end

"""
gforward!(tau,reaction,t,r,nallele)

update tau matrix for forward transition

"""
function gforward!(tau,reaction,t,r,nallele)
	tau[reaction,nallele] = Inf
	tau[reaction+1,nallele] = -log(rand(rng))/r[reaction+1] + t
	tau[reaction+2,nallele] = -log(rand(rng))/r[reaction+2] + t
	if reaction > 1
		tau[reaction-1,nallele] = Inf
	end
	nothing
end

"""
gbackward!(tau,reaction,t,r,nallele)

update tau matrix for reverse state transition

"""
function greverse!(tau,reaction,t,r,nallele)
	tau[reaction,nallele] = Inf
	tau[reaction+1,nallele] = Inf
	tau[reaction-1,nallele]= -log(rand(rng))/r[reaction-1] + t
	if reaction > 2
		tau[reaction-2,nallele] = -log(rand(rng))/r[reaction-2] + t
	end
	nothing
end

"""
create!(tau,reaction,t,m,r,nallele)

update tau matrix for mRNA creation transition

"""
function create!(tau,reaction,t,m,r,nallele)
	m += 1
	tau[reaction,nallele] = -log(rand(rng))/r[reaction] + t
	# tau[reaction+1,1] = -log(rand(rng))/(m*r[reaction+1])+ t
	if m == 1
		tau[reaction+1,1] = -log(rand(rng))/(r[reaction+1])+ t
	else
		tau[reaction+1,1] = (m-1)/m*(tau[reaction+1,1] - t) + t
	end
	m
end

"""
decay!(tau,reaction,t,m,r)

update tau matrix for decay transition

"""
function decay!(tau,reaction,t,m,r)
	m -= 1
	tau[reaction,1] = -log(rand(rng))/(m*r[reaction]) + t
	m
end
