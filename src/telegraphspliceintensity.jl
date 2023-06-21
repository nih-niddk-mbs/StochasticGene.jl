
function prep_param_telegraph(r,n,zeta)
	# Transition rates
	gammap = zeros(n+1)
	gamman = zeros(n+1)
	nu = zeros(zeta+1)
	eps = zeros(zeta)
	for i = 1:n
		gammap[i] = r[2*(i-1)+1]
		gamman[i+1] = r[2*i]
	end
	for i = 1:zeta+1
		nu[i] = r[2*n+i]
	end
    # Splice rate is non decreasing, rates in r represent increase in splice rate from previous step
	eps[1] = r[2*n + 1 + zeta + 1]
	for i = 2:zeta
		eps[i] = eps[i-1] + r[2*n + 1 + zeta + i]
	end
	return gammap,gamman,nu,eps
end


#Compute dwelltime and steady state mRNA distributions
function telegraphsplice0(range::Array,nhist::Int,n::Int,zeta::Int,rin::Vector,total::Int,tol::Float64,nalleles::Int)

	# decay reaction index
	reactiondecay = 2*n + 2*zeta + 2
	# dwell time arrays
	ndt = length(range)
	dt = range[2]-range[1]
	histofftdd = zeros(Int,ndt)
	histontdd  = zeros(Int,ndt)
	# ss mRNA array
	mhist = zeros(nhist+1)
	# tIA and TAI are times when gene state transitions between active and inactive states
	# they are used to compute dwell times
	tIA = zeros(Float64,nalleles)
	tAI = zeros(Float64,nalleles)
	# active start site
	mu = n
	# Get rates
	gammaf,gammab,nu,eta = prep_param_telegraph(rin,n,zeta)
	# assign decay rate
	delta = rin[reactiondecay]
	# assign correlation rate
	if length(rin) > reactiondecay
		rcorr = rin[reactiondecay+1]
	else
		rcorr = 0.
	end
	#initialize mRNA number to 0
	m = 0
	mproduced = 0
	# Initialize gene states
	state = ones(Int,nalleles)
	# pre-mRNA states
	z = Array{Array{Int,1},1}(undef,nalleles)
	# Intron states
	intron = Array{Array{Int,1},1}(undef,nalleles)
	# initialize propensities
	af = Array{Float64,1}(undef,nalleles)
	ab = Array{Float64,1}(undef,nalleles)
	apre1 = Array{Float64,1}(undef,nalleles)
	az = Array{Array{Float64,1},1}(undef,nalleles)
	aintron = Array{Array{Float64,1},1}(undef,nalleles)
	apremid = zeros(nalleles)     # sum(az[2:zeta]) middle pre-mRNAstep
	afree = zeros(nalleles)       # az[zeta+1]  create free mRNA
	for i = 1:nalleles
		z[i] = zeros(Int,zeta)  # Intialize pre-mRNA state to all 0
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
		az[i] = zeros(zeta+1)
		# apre1[i] = float(1-z[i][1])*nu[1]*float(state[i]==n+1)   #first pre-mRNA step
		apre1[i] = float(state[i] > mu)*float(1-z[i][1])*nu[1]   #first pre-mRNA step
		intron[i] = zeros(Int,zeta)
		aintron[i] = zeros(zeta)
	end
	astate = af + ab
	adecay = 0.       #m*delta
	amrna = apre1 + apremid + afree
	t=0.  # time
	ts=0. # time from last sample
	tsample = 20/minimum(rin)  # sample every 100 decay times
	steps = 0
	err = 10*tol
	hoff0 = ones(ndt)/ndt
	hon0  = ones(ndt)/ndt
	mhist0 = ones(nhist+1)
	intensity = Matrix{Float64}(undef,0,2)
	# Begin simulation
	while err > tol && steps < total
		steps += 1
		intensity = vcat(intensity,[t sum(intron[1])])
		asum = sum(amrna) + sum(astate) + adecay
		if asum > 0.
			# choose which state to update and time of update
			r = asum*rand()
			tau = -log(rand())/asum
			# update time
			t += tau
			#print(tau,' ',asum,' ',state,'\n')
			if m + 1 <= nhist
				mhist[m+1] += tau
			else
				mhist[nhist+1] += tau
			end
			if r > sum(amrna + astate)
				# mRNA decay
				m -= 1
				adecay = m*delta
			else
				agenecs = amrna[1] + astate[1]
				agene0 = 0.
				l = 1
				while r > agenecs
					agene0 = agenecs
					agenecs += amrna[1+l] + astate[l+1]
					l += 1
				end
				if l > 1
					r -= agene0
				end
				if r > amrna[l]
					# state update
					if r > ab[l] + amrna[l]
						# update forward reaction
						state[l] += 1
					else
						# update backward reaction
						state[l] -= 1
					end
					# update propensities
					af[l] = gammaf[state[l]]
					ab[l] = gammab[state[l]]
					astate[l] = af[l] + ab[l]
					# Activate coupling if any alleles are in state mu + 1
					if rcorr > 0.
						ba = false
						for k = 1:nalleles
							ba |= state[k] > mu
						end
						if ba
							for k = 1:nalleles
								if state[k] == mu
									af[k] = gammaf[state[k]] + rcorr
									astate[k] = af[k] + ab[k]
								end
							end
						else
							for k = 1:nalleles
								if state[k] == mu
									af[k] = gammaf[state[k]]
									astate[k] = af[k] + ab[k]
								end
							end
						end
					end
					if state[l] > mu
						apre1[l] = float(1-z[l][1])*nu[1]
						# apre1 = az[l][1]
					else
						apre1[l] = 0.
					end
					amrna[l] = apremid[l] + apre1[l] + afree[l] + sum(aintron[l])
				else
					# premRNA update
					if r > apremid[l] + apre1[l] + sum(aintron[l])
						# increase free mRNA
						m += 1
						adecay = m*delta
						mproduced += 1
						z[l][zeta] = 0
						#az[l][zeta+1] = 0
						afree[l] = 0.
						# ON state ends, OFF state begins
						if sum(intron[l]) == 1 && intron[l][zeta] == 1
							tAI[l] = t
							tactive = (tAI[l] - tIA[l])/dt
							if tactive <= ndt && tactive > 0
								tactive = ceil(Int,tactive)
								histontdd[tactive] += 1
							end
						end
						intron[l][zeta] = 0
						aintron[l][zeta] = 0.
						if zeta == 1
							if state[l] > mu
								#az[l][1] = nu[1]
								apre1[l] = nu[1]
							else
								#az[l][1] = 0
								apre1[l] = 0.
							end
						else
							az[l][zeta] = float(z[l][zeta-1])*nu[zeta]
							apremid[l] = sum(az[l][2:zeta])
						end
						amrna[l] = apre1[l] + apremid[l] + afree[l] + sum(aintron[l])
					elseif r > apremid[l] + apre1[l]
						# eject intron
						aintronsum = apremid[l] + apre1[l]
						i = 0
						while r > aintronsum
							i += 1
							aintronsum += aintron[l][i]
						end
						intron[l][i] = 0
						aintron[l][i] = 0
						# ON state ends, OFF state begins
						if sum(intron[l]) == 0 #&& nhist == 0
							tAI[l] = t
							tactive = (tAI[l] - tIA[l])/dt
							if tactive <= ndt && tactive > 0
								tactive = ceil(Int,tactive)
								histontdd[tactive] += 1
							end
						end
						amrna[l] = apremid[l] + apre1[l] + afree[l] + sum(aintron[l])
						# println("*",apre1)
					elseif r > apremid[l]
						# increase first pre-mRNA state
						# OFF state ends, ON state begins
						if sum(intron[l]) == 0 #&& nhist == 0
							tIA[l] = t
							tinactive = (tIA[l] - tAI[l])/dt
							if tinactive <= ndt && tinactive > 0
								tinactive = ceil(Int,tinactive)
								histofftdd[tinactive] += 1
							end
						end
						z[l][1] = 1
						#az[l][1] = 0
						apre1[l] = 0.
						intron[l][1] = 1
						aintron[l][1] = eta[1]
						if zeta == 1
							#az[l][2] = nu[2]
							afree[l] = nu[2]
						else
							az[l][2] = float(1-z[l][2])*nu[2]
							apremid[l] = sum(az[l][2:zeta])
						end
						amrna[l] = apre1[l] + apremid[l] + afree[l] + sum(aintron[l])
					else
						# update mid pre-mRNA state
						azsum = az[l][2]
						k = 1
						while r > azsum
							azsum += az[l][2+k]
							k += 1
						end
						i = k + 1
						z[l][i] = 1
						z[l][i-1] = 0
						az[l][i] = 0.
						if intron[l][i-1] == 1
							intron[l][i] = 1
							intron[l][i-1] = 0
							aintron[l][i] = eta[i]
							aintron[l][i-1] = 0.
						else
							intron[l][i] = 0
							aintron[l][i] = 0.
						end
						if i == zeta
							#az[l][i+1] = nu[i+1]
							afree[l] = nu[i+1]
						else
							az[l][i+1] = float(1-z[l][i+1])*nu[i+1]
						end
						if i == 2
							if state[l] > mu
								#az[l][1] = nu[1]
								apre1[l] = nu[1]
							else
								az[l][1] = 0.
								apre1[l] = 0.
							end
						else
							az[l][i-1] = float(z[l][i-2])*nu[i-1]
						end
						apremid[l] = sum(az[l][2:zeta])
						amrna[l] = apre1[l] + apremid[l] + afree[l] + sum(aintron[l])
					end # if r > apremid[l] + apre1[l]
				end # if r > amrna[l]
			end #if r > sum(amrna + adecay)
		end  #if asum > 0
		# Compute change in histograms from last sample
		ts += tau
		if ts > tsample
			hoff = histofftdd/sum(histofftdd)
			hon = histontdd/sum(histontdd)
			err = max(norm(hoff-hoff0,Inf),norm(hon-hon0,Inf))
			hoff0 = copy(hoff)
			hon0 = copy(hon)
			errss = norm(mhist/sum(mhist) - mhist0/sum(mhist0),Inf)
			err = max(err,errss)
			mhist0 = copy(mhist)
			ts = 0.
		end
	end  # while bursts < nbursts && steps < total
	sumoff = sum(histofftdd)
	sumon = sum(histontdd)
	counts = max(sum(mhist),1)
	mhist /= counts
	# print(mproduced/t)
	return histofftdd/max(sumoff,1), histontdd/max(sumon,1), mhist[1:nhist], intensity
end
