#
# Gillespie Algorithm functions to compute:
# Live Cell ON/OFF distributions
# smFISH steady state histograms
# G state and R step Occupancy probabilities
# Mean dwell time for being ON
#

# Compute ON OFF time histograms
function telegraphprefast(range::Array,n::Int,zeta::Int,rin::Vector,total::Int,nbursts::Int,nalleles::Int)

	# Generalized n-zeta telegraph model
	# Compute dwell time histograms
	# n+1 total gene states
	# gene states are labeled from 1 to n+1
	# zeta pre-mRNA steps
	# There are 2n gene reactions, zeta pre-mRNA reactions, and 1 mRNA decay reaction
	# rin = reaction rates
	# forward gene reactions are labeled by odd numbers <2n, backward reactions are even numbers <= 2n
	# reactions [2n+1, 2n+zeta] are pre-mRNA recations,organized by starting spot, reaction 2n(zeta+1)+1 is mRNA decay reaction

	# Use Gillespie direct method with speed up from Gibson and Bruck

	# copy rates
	# r = copy(rin)
	# nbursts = 100000
	# nbursts = 1000

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

	# dwell time arrays
	ndt = length(range)
	dt = range[2]-range[1]
	histofftdd = zeros(Int,ndt)
	histontdd  = zeros(Int,ndt)

	# time
	t=0.

	# tIA and TAI are times when gene state transitions between active and inactive states
	# they are used to compute dwell times
	tIA = zeros(Float64,nalleles)
	tAI = zeros(Float64,nalleles)

	#Initialize output arrays
	# mRNA = Array{Int,1}(total)
	# times = Array{Float64,1}(total)

	# active start site
	mu = n
	#print(n," ",mu," ",zeta,'\n')

	# Initialize rates
	gammaf = zeros(n+1)
	gammab = zeros(n+1)
	nu = Array{Float64,1}(undef,zeta+1)

	# assign gene rates
	for i = 1:n
		gammaf[i] = rin[2*(i-1)+1]
		gammab[i+1] = rin[2*i]
	end

	# assign pre-mRNA rates
	for i = 1:zeta+1
		nu[i] = rin[2*n + i]
	end

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

	# Initialize gene states
	state = ones(Int,nalleles)
	#state[1] = mu + 1
	# state[1] = 2
	# nindex = n + 1

	# pre-mRNA states
	z = Array{Array{Int,1},1}(undef,nalleles)

	# initialize propensities
	af = Array{Float64,1}(undef,nalleles)
	ab = Array{Float64,1}(undef,nalleles)
	apre1 = Array{Float64,1}(undef,nalleles)
	az = Array{Array{Float64,1},1}(undef,nalleles)

	apremid = zeros(nalleles)     # sum(az[2:zeta]) middle pre-mRNAstep
	afree = zeros(nalleles)       # az[zeta+1]  create free mRNA

	for i = 1:nalleles
		z[i] = zeros(Int,zeta)  # Intialize pre-mRNA state to all 0
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
		az[i] = zeros(zeta+1)
		# apre1[i] = float(1-z[i][1])*nu[1]*float(state[i]==n+1)   #first pre-mRNA step
		apre1[i] = 0.   #first pre-mRNA step
	end

	astate = af + ab

	adecay = 0.       #m*delta

	amrna = apre1 + apremid + afree
	bursts = 0
	steps = 0

	# for steps = 1:total
	while bursts < nbursts && steps < total

		steps += 1

		asum = sum(amrna) + sum(astate) + adecay

		if asum > 0.
			# choose which state to update and time of update
			r = asum*rand()
			tau = -log(rand())/asum

			# update time
			t += tau
			#print(tau,' ',asum,' ',state,'\n')

			# times[steps] = tau
			# mRNA[steps] = m

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

					amrna[l] = apremid[l] + apre1[l] + afree[l]

				else

					# premRNA update
					if r > apremid[l] + apre1[l]
						# increase free mRNA
						m += 1
						adecay = m*delta

						z[l][zeta] = 0
						#az[l][zeta+1] = 0
						afree[l] = 0.

						if sum(z[l]) == 0 #&& nhist == 0
							tAI[l] = t
							tactive = ceil(Int,(tAI[l] - tIA[l])/dt)
							if tactive <= ndt && tactive > 0
								histontdd[tactive] += 1
							end
						end

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

						amrna[l] = apre1[l] + apremid[l] + afree[l]

					elseif r > apremid[l]
						# increase first pre-mRNA state

						if sum(z[l]) == 0 #&& nhist == 0
							tIA[l] = t
							tinactive = ceil(Int,(tIA[l] - tAI[l])/dt)
							if tinactive <= ndt && tinactive > 0
								histofftdd[tinactive] += 1
								bursts += 1
							end
						end

						z[l][1] = 1
						#az[l][1] = 0
						apre1[l] = 0.

						if zeta == 1
							#az[l][2] = nu[2]
							afree[l] = nu[2]
						else
							az[l][2] = float(1-z[l][2])*nu[2]
							apremid[l] = sum(az[l][2:zeta])
						end

						amrna[l] = apre1[l] + apremid[l] + afree[l]

					else
						# update mid pre-mRNA state
						azsum = az[l][2]
						k = 1
						while r > azsum
							azsum += az[l][2+k]
							k += 1
						end

						i = k + 1

						if i <= zeta

							z[l][i] = 1
							z[l][i-1] = 0
							az[l][i] = 0.

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
						end
						apremid[l] = sum(az[l][2:zeta])
						amrna[l] = apre1[l] + apremid[l] + afree[l]

					end # if r > apremid[l] + apre1[l]
				end # if r > amrna[l]
			end #if r > sum(amrna + adecay)

		end  #if asum > 0

	end  # for steps = 1:total

	sumoff = sum(histofftdd)
	sumon = sum(histontdd)

	if sumon > 0 && sumoff > 0
		#return [histofftdd[2]/sumoff; histontdd[2]/sumon], mhistn
		return cumsum(histofftdd/sumoff), cumsum(histontdd/sumon)
	else
		return zeros(length(histofftdd)),zeros(length(histontdd))
	end

end



# Compute smFISH histogram
function telegraphprefast(n::Int,zeta::Int,rin::Vector,total::Int,tmax::Float64,nhist::Int,nalleles::Int,count=false)

	# Generalized n, zeta, telegraph model,
	# Compute steady state mRNA histogram
	# n+1 total gene states
	# gene states are labeled from 1 to n+1
	# zeta pre-mRNA steps
	# There are 2n gene reactions, zeta pre-mRNA reactions, and 1 mRNA decay reaction
	# rin = reaction rates
	# forward gene reactions are labeled by odd numbers <2n, backward reactions are even numbers <= 2n
	# reactions [2n+1, 2n+zeta] are pre-mRNA reactions,organized by starting spot, reaction 2n(zeta+1)+1 is mRNA decay reaction

	# Use Gillespie direct method with speed up from Gibson and Bruck

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

	# time
	t=0.

	# tIA and TAI are times when gene state transitions between active and inactive states
	# they are used to compute dwell times
	# tIA = zeros(Float64,nalleles)
	# tAI = zeros(Float64,nalleles)

	#Initialize m histogram
	mhist=zeros(nhist+1)
	# mRNA = Array{Int,1}(total)
	# times = Array{Float64,1}(total)

	# active start site
	mu = n
	#print(n," ",mu," ",zeta,'\n')

	# Initialize rates
	gammaf = zeros(n+1)
	gammab = zeros(n+1)
	nu = Array{Float64,1}(undef,zeta+1)

	# assign gene rates
	for i = 1:n
		gammaf[i] = rin[2*(i-1)+1]
		gammab[i+1] = rin[2*i]
	end

	# assign pre-mRNA rates
	for i = 1:zeta+1
		nu[i] = rin[2*n + i]
	end

	# assign decay rate
	delta = rin[reactiondecay]

	# assign correlation rate
	if length(rin) > reactiondecay
		rcorr = rin[reactiondecay+1]
	else
		rcorr = 0.
	end
	rcorr = 0.

	#initialize mRNA number to 0
	m = 0

	# Initialize gene states
	state = ones(Int,nalleles)
	#state[1] = mu + 1
	# state[1] = 2
	# nindex = n + 1

	# pre-mRNA states
	z = Array{Array{Int,1},1}(undef,nalleles)

	# initialize propensities
	af = Array{Float64,1}(undef,nalleles)
	ab = Array{Float64,1}(undef,nalleles)
	apre1 = Array{Float64,1}(undef,nalleles)
	az = Array{Array{Float64,1},1}(undef,nalleles)

	apremid = zeros(nalleles)     # sum(az[2:zeta]) middle pre-mRNAstep
	afree = zeros(nalleles)       # az[zeta+1]  create free mRNA

	for i = 1:nalleles
		z[i] = zeros(Int,zeta)  # Intialize pre-mRNA state to all 0
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
		az[i] = zeros(zeta+1)
		apre1[i] = 0.   #first pre-mRNA step
	end

	astate = af + ab

	adecay = 0.       #m*delta

	amrna = apre1 + apremid + afree
	bursts = 0
	steps = 0

	# Begin simulation
	while t < tmax && steps < total

		steps += 1

		asum = sum(amrna) + sum(astate) + adecay

		if asum > 0.
			# choose which state to update and time of update
			r = asum*rand()
			tau = -log(rand())/asum

			# update time
			t += tau
			#print(tau,' ',asum,' ',state,'\n')

			# times[steps] = tau
			# mRNA[steps] = m

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

					amrna[l] = apremid[l] + apre1[l] + afree[l]

				else

					# premRNA update
					if r > apremid[l] + apre1[l]
						# increase free mRNA
						m += 1
						adecay = m*delta

						z[l][zeta] = 0
						#az[l][zeta+1] = 0
						afree[l] = 0.

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

						amrna[l] = apre1[l] + apremid[l] + afree[l]

					elseif r > apremid[l]
						# increase first pre-mRNA state

						z[l][1] = 1
						#az[l][1] = 0
						apre1[l] = 0.

						if zeta == 1
							#az[l][2] = nu[2]
							afree[l] = nu[2]
						else
							az[l][2] = float(1-z[l][2])*nu[2]
							apremid[l] = sum(az[l][2:zeta])
						end

						amrna[l] = apre1[l] + apremid[l] + afree[l]

					else
						# update mid pre-mRNA state
						azsum = az[l][2]
						k = 1
						while r > azsum
							azsum += az[l][2+k]
							k += 1
						end

						i = k + 1

						if i <= zeta

							z[l][i] = 1
							z[l][i-1] = 0
							az[l][i] = 0.

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
						end
						apremid[l] = sum(az[l][2:zeta])
						amrna[l] = apre1[l] + apremid[l] + afree[l]

					end # if r > apremid[l] + apre1[l]
				end # if r > amrna[l]
			end #if r > sum(amrna + adecay)

		end  #if asum > 0

	end  # for steps = 1:total

	counts = max(sum(mhist),1)
	mhist /= counts

	if count
		return mhist[1:nhist],counts
	else
		return mhist[1:nhist]
	end

end


# Compute steady state occupancy probabilities for all G states and R steps
function telegraphprefastSS(n::Int,zeta::Int,rin::Vector,total::Int,nalleles::Int)

	#Steady state arrays
	Gss = zeros(n+1)
	Rss = 0.
	Rssc = 0.

	# copy rates
	r = copy(rin)

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

	# time
	t=0.

	# active start site
	mu = n

	# Initialize rates
	gammaf = zeros(n+1)
	gammab = zeros(n+1)
	nu = Array{Float64}(undef,zeta+1)

	# assign gene rates
	for i = 1:n
		gammaf[i] = r[2*(i-1)+1]
		gammab[i+1] = r[2*i]
	end

	# assign pre-mRNA rates
	for i = 1:zeta+1
		nu[i] = r[2*n + i]
	end

	# assign decay rate
	delta = r[reactiondecay]

	# assign correlation rate
	if length(r) > reactiondecay
		rcorr = r[reactiondecay+1]
	else
		rcorr = 0.
	end

	#initialize mRNA number to 0
	m=0

	# Initialize to gene state 1, all other states are zero
	state = ones(Int,nalleles)
	#state[1] = mu + 1
	#state[2] = mu
	nindex = n + 1

	# pre-mRNA states
	z = Array{Array{Int,1}}(undef,nalleles)

	# initialize propensities
	af = Array{Float64}(undef,nalleles)
	ab = similar(af)
	apre1 = similar(af)
	az = Array{Array{Float64,1}}(undef,nalleles)

	apremid = zeros(nalleles)     # sum(az[2:zeta]) middle pre-mRNAstep
	afree = zeros(nalleles)       # az[zeta+1]  create free mRNA

	for i = 1:nalleles
		z[i] = zeros(Int8,zeta)  # Intialize pre-mRNA state to all 0
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
		az[i] = zeros(zeta+1)
		apre1[i] = float(1-z[i][1])*nu[1]*float(state[i]==n+1)   #first pre-mRNA step
	end
	astate = af + ab

	adecay = 0.       #m*delta

	amrna = apre1 + apremid + afree

	for steps = 1:total

		asum = sum(amrna) + sum(astate) + adecay

		if asum > 0.
			# choose which state to update and time of update
			r = asum*rand()
			tau = -log(rand())/asum

			# update time
			t += tau
			#print(tau,' ',asum,' ',state,'\n')

			#times[steps]=tau
			#mRNA[steps] = m
			for j = 1:nalleles
				Gss[state[j]] += tau
				if state[j]> mu && z[j][1] == 1 #&& sum(z[j]) == 1
					Rss += tau
				end
			end
			if reduce(|,state) > 1 && reduce(|,z[:][1]) == 1
				Rssc += tau
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
					# if ba
					# 	for k = 1:nalleles
					# 		af[k] = gammaf[state[k]] + rcorr*float(state[k]==mu)
					# 		astate[k] = af[k] + ab[k]
					# 	end
					# else
					# 	for k = 1:nalleles
					# 		af[k] = gammaf[state[k]]
					# 		astate[k] = af[k] + ab[k]
					# 	end
					# end

					if state[l] > mu
						apre1[l] = float(1-z[l][1])*nu[1]
						# apre1 = az[l][1]
					else
						apre1[l] = 0.
					end
					amrna[l] = apremid[l] + apre1[l] + afree[l]
				else

					# premRNA update
					if r > apremid[l] + apre1[l]
						# increase free mRNA
						m += 1
						adecay = m*delta

						z[l][zeta] = 0
						#az[l][zeta+1] = 0
						afree[l] = 0.

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

						amrna[l] = apre1[l] + apremid[l] + afree[l]

					elseif r > apremid[l]
						# increase first pre-mRNA state
						z[l][1] = 1
						#az[l][1] = 0
						apre1[l] = 0.

						if zeta == 1
							#az[l][2] = nu[2]
							afree[l] = nu[2]
						else
							az[l][2] = float(1-z[l][2])*nu[2]
							apremid[l] = sum(az[l][2:zeta])
						end

						amrna[l] = apre1[l] + apremid[l] + afree[l]

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
						amrna[l] = apre1[l] + apremid[l] + afree[l]

					end # if r > apremid[l] + apre1[l]
				end # if r > amrna[l]
			end #if r > sum(amrna + adecay)

		end  #if asum > 0

	end  # for steps = 1:total

	println(Gss)
	println(Rss,",",t)

	Gss /= t*nalleles
	Rss /= t*nalleles
	Rssc /= t*nalleles

	eff = Rss/Gss[n+1]
	effc = Rssc/Gss[n+1]

	println(sum(Gss))

	return  Gss, eff, effc

end

# Compute mean dwell times
function telegraphmdt(range::Array,n::Int,zeta::Int,rin::Vector,total::Int,nbursts::Int,nalleles::Int)

	# decay reaction index
	reactiondecay = 2*n + zeta + 2

	# dwell time arrays
	ndt = length(range)
	dt = range[2]-range[1]
	histofftdd = zeros(Int,ndt)
	histontdd  = zeros(Int,ndt)

	# time
	t=0.

	# tIA and TAI are times when gene state transitions between active and inactive states
	# they are used to compute dwell times
	tIA = zeros(Float64,nalleles)
	tAI = zeros(Float64,nalleles)

	ONtimes = Array{Float64}(undef,0)

	#Initialize output arrays
	# mRNA = Array{Int,1}(total)
	# times = Array{Float64,1}(total)

	# active start site
	mu = n
	#print(n," ",mu," ",zeta,'\n')

	# Initialize rates
	gammaf = zeros(n+1)
	gammab = zeros(n+1)
	nu = Array{Float64,1}(undef,zeta+1)

	# assign gene rates
	for i = 1:n
		gammaf[i] = rin[2*(i-1)+1]
		gammab[i+1] = rin[2*i]
	end

	# assign pre-mRNA rates
	for i = 1:zeta+1
		nu[i] = rin[2*n + i]
	end

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

	# Initialize gene states
	state = ones(Int,nalleles)
	#state[1] = mu + 1
	state[1] = 2
	# nindex = n + 1

	# pre-mRNA states
	z = Array{Array{Int,1},1}(undef,nalleles)

	# initialize propensities
	af = Array{Float64,1}(undef,nalleles)
	ab = Array{Float64,1}(undef,nalleles)
	apre1 = Array{Float64,1}(undef,nalleles)
	az = Array{Array{Float64,1},1}(undef,nalleles)

	apremid = zeros(nalleles)     # sum(az[2:zeta]) middle pre-mRNAstep
	afree = zeros(nalleles)       # az[zeta+1]  create free mRNA

	for i = 1:nalleles
		z[i] = zeros(Int,zeta)  # Intialize pre-mRNA state to all 0
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
		az[i] = zeros(zeta+1)
		# apre1[i] = float(1-z[i][1])*nu[1]*float(state[i]==n+1)   #first pre-mRNA step
		apre1[i] = 0.   #first pre-mRNA step
	end

	astate = af + ab

	adecay = 0.       #m*delta

	amrna = apre1 + apremid + afree
	bursts = 0
	steps = 0

	# for steps = 1:total
	while bursts < nbursts && steps < total

		steps += 1

		asum = sum(amrna) + sum(astate) + adecay

		if asum > 0.
			# choose which state to update and time of update
			r = asum*rand()
			tau = -log(rand())/asum

			# update time
			t += tau
			#print(tau,' ',asum,' ',state,'\n')

			# times[steps] = tau
			# mRNA[steps] = m

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

					amrna[l] = apremid[l] + apre1[l] + afree[l]

				else

					# premRNA update
					if r > apremid[l] + apre1[l]
						# increase free mRNA
						m += 1
						adecay = m*delta

						z[l][zeta] = 0
						#az[l][zeta+1] = 0
						afree[l] = 0.

						if sum(z[l]) == 0 #&& nhist == 0
							tAI[l] = t
							tactive = ceil(Int,(tAI[l] - tIA[l])/dt)
							push!(ONtimes,tAI[l] - tIA[l])
							if tactive <= ndt && tactive > 0
								histontdd[tactive] += 1
							end
						end

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

						amrna[l] = apre1[l] + apremid[l] + afree[l]

					elseif r > apremid[l]
						# increase first pre-mRNA state

						if sum(z[l]) == 0 #&& nhist == 0
							tIA[l] = t
							tinactive = ceil(Int,(tIA[l] - tAI[l])/dt)
							if tinactive <= ndt && tinactive > 0
								histofftdd[tinactive] += 1
								bursts += 1
							end
						end

						z[l][1] = 1
						#az[l][1] = 0
						apre1[l] = 0.

						if zeta == 1
							#az[l][2] = nu[2]
							afree[l] = nu[2]
						else
							az[l][2] = float(1-z[l][2])*nu[2]
							apremid[l] = sum(az[l][2:zeta])
						end

						amrna[l] = apre1[l] + apremid[l] + afree[l]

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
						amrna[l] = apre1[l] + apremid[l] + afree[l]

					end # if r > apremid[l] + apre1[l]
				end # if r > amrna[l]
			end #if r > sum(amrna + adecay)

		end  #if asum > 0

	end  # for steps = 1:total

	return mean(ONtimes), std(ONtimes)/sqrt(length(ONtimes)), ONtimes

end
