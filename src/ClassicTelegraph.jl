

# Compute smFISH histogram - classic telegraph model for arbtirary numbers of gene states
function telegraph(n::Int,rin::Vector,total::Int,tol::Float64,nhist::Int,nalleles::Int,count=false)

	# n-state telegraph model for nalleles alleles
	# Compute steady state mRNA histogram
	# n+1 total gene states
	# gene states are labeled from 1 to n+1
	# There are 2n gene reactions, 1 mRNA creation reaction, and 1 mRNA decay reaction
	# rin = reaction rates
	# forward gene reactions are labeled by odd numbers <2n, backward reactions are even numbers <= 2n
	# reactions 2n+1 mRNA creation reaction, reaction 2n+2 is mRNA decay reaction

	# Gillespie direct method with speed up from Gibson and Bruck

	#Initialize m histogram
	mhist = zeros(nhist+1)

	# active start site is n+1
	mu = n
	#print(n," ",mu," ",zeta,'\n')

	# Initialize rates
	gammaf = zeros(n+1)
	gammab = zeros(n+1)

	# assign gene rates
	for i = 1:n
		gammaf[i] = rin[2*(i-1)+1]
		gammab[i+1] = rin[2*i]
	end

	# mRNA creation rate
	nu = rin[2*n + 1]

	# assign decay rate
	# decay reaction index
	reactiondecay = 2*n + 2
	delta = rin[reactiondecay]

	# allele coupling strength
	if length(rin) > reactiondecay
		rcorr = rin[reactiondecay+1]
	else
		rcorr = 0.
	end

	#initialize mRNA number to 0
	m = 0

	# Initialize gene state to most inactive state
	state = ones(Int,nalleles)

	# initialize propensities
	af = Array{Float64,1}(undef,nalleles)  # forward
	ab = Array{Float64,1}(undef,nalleles)  # backward
	afree = nu*ones(nalleles).*(state .== 1)    # mRNA creation
	adecay = 0.       # mRNA decay

	for i = 1:nalleles
		af[i] = gammaf[state[i]]
		ab[i] = gammab[state[i]]
	end
	astate = af + ab

	steps = 0

	# Begin simulation
	t=0.
	tsample = 20/minimum(rin)
	sample = 0
	err = 10*tol
	mhist0 = ones(nhist+1)

	while err > tol && steps < total

		steps += 1

		# sum of all propensities
		asum = sum(afree) + sum(astate) + adecay

		if asum > 0.
			# choose which state to update and time of update
			r = asum*rand()
			tau = -log(rand())/asum

			# update sample time
			t += tau
			#print(tau,' ',asum,' ',state,'\n')

			# update m histogram
			if m + 1 <= nhist
				mhist[m+1] += tau
			else
				mhist[nhist+1] += tau
			end

			# Compute change in histogram from last sample

			if t > tsample
				err = norm(mhist/sum(mhist)-mhist0/sum(mhist0),Inf)
				mhist0 = copy(mhist)
				t = 0.
			end

			if r > sum(afree + astate)
				# mRNA decay reaction
				m -= 1
				adecay = m*delta
			else
				# gene state reaction
				# find gene that updated
				agenecs = afree[1] + astate[1]
				agene0 = 0.
				l = 1
				while r > agenecs
					agene0 = agenecs
					agenecs += afree[1+l] + astate[l+1]
					l += 1
				end
				if l > 1
					r -= agene0
				end

				if r > afree[l]
					# state update
					if r > ab[l] + afree[l]
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

					# gene state is active if > mu
					if state[l] > mu
						afree[l] = nu
					else
						afree[l] = 0.
					end

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
				else
					# increase free mRNA
					m += 1
					adecay = m*delta
				end # if r > free[l]
			end #if r > sum(afree + astate)

		end  #if asum > 0

	end  # while

	counts = max(sum(mhist),1)
	mhist /= counts

	if count
		return mhist[1:nhist],counts,steps,err
	else
		return mhist[1:nhist]
	end

end


function prep_params(n,r)
	gammap = zeros(n+2)
	gamman = zeros(n+2)
	for i = 1:n
		gammap[i+1] = r[2*(i-1)+1]
		gamman[i+1] = r[2*i]
	end
	return gammap, gamman
end

function transitiontelegraph(n,rin,nhist)
	nT = n+1
	total = nhist + 2
	r = rin/rin[2*n + 2]
	gammap,gamman = prep_params(n,r)  # State Transition rates
	nu = r[2*n+1] # mRNA creation rate
	T = zeros(nT,nT)
	B = zeros(nT,nT)
	S = zeros(total,total)
	Sminus = similar(S)
	Splus = similar(S)
	# Generate T = A + B matrix and B matrix
	for ip=1:nT, i=1:nT
		T[i,ip] =  gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip)
		B[i,ip] =  nu*float(i==ip)*float(i==n+1)
	end
	# Generate matrices for m transitions
	for m=1:total,mp=1:total
		S[m,mp] = -(mp)*float(m==mp)
		Sminus[m,mp] = float(m==mp+1)
		Splus[m,mp] = (mp-1)*float(m==mp-1)
	end
	# Construct m truncated Master equation rate matrix
	M = kron(S,Matrix{Float64}(I,nT,nT)) + kron(Matrix{Float64}(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,Matrix{Float64}(I,nT,nT))
	M[end,end]+=nu  # boundary condition to ensure probability is conserved
	return M
end

function transitiontelegraph(r,nhist)
	nu = r[1]/r[2]
	total = nhist + 2
	S = zeros(total,total)
	for m=1:total,mp=1:total
		S[m,mp] = nu *(m==mp+1) - nu*(m==mp) + (mp-1)*(m==mp-1) - (mp-1)*(m==mp)
	end
	S[end,end] += nu
	return S
end

# function marginalize(p,n,nhist)
# 	mhist = zeros(nhist)
# 	nT = n+1
# 	# Marginalize over G states
# 	for m in 1:nhist
# 		i = (m-1)*nT
# 		mhist[m] = sum(p[i+1:i+nT])
# 	end
# 	return mhist
# end

function transientMaster(t::Float64,n::Int,P0,Mvals::Vector,Mvects::Matrix)
	Lambda = Diagonal(exp.(Mvals*t))
	Q = V/P0
	P = V*Lambda*Q
end

# function sstel(n,r,nhist,nalleles)
# 	M = transitiontelegraph(n,r,nhist)
# 	p = normalizednullspace(M)
# 	mhist = marginalize(p,n,nhist)
# 	alleleConvolve(mhist,nalleles)
# end
#
# function sstel(r,nhist,nalleles)
# 	S = transitiontelegraph(r,nhist)
# 	mhist = normalizednullspace(S)
# 	alleleConvolve(mhist,nalleles)
# end

function sstel2(n::Int,rin::Vector{Float64},nhist::Int,nalleles::Int)
	# Closed form steady state mRNA distribution of generalized n-state classic telegraph model
	nT = n+1
	r = rin/rin[2*n + 2]
	total = nhist + 2
	# Transition rates
	gammap,gamman = prep_params(n,r)
	# mRNA creation rate
	nu = r[2*n+1]
	# Declare Transition matrices
	T = zeros(nT,nT)
	B = zeros(nT,nT)
	S = zeros(total,total)
	Sminus = similar(S)
	Splus = similar(S)
	# Generate T = A + B matrix and B matrix
	for i=1:n+1,  ip=1:n+1
		T[i,ip] =  gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip)
		B[i,ip] =  nu*float(i==ip)*float(i==n+1)
	end
	# Generate matrices for m transitions
	for m=1:total,mp=1:total
		S[m,mp] = -(mp-1)*float(m==mp)
		Sminus[m,mp] = float(m==mp+1)
		Splus[m,mp] = (mp-1)*float(m==mp-1)
	end
	# Construct m truncated Master equation rate matrix
	M = kron(S,Matrix{Float64}(I,nT,nT)) + kron(Matrix{Float64}(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,Matrix{Float64}(I,nT,nT))
	# boundary condition to ensure probability is conserved
	M[end,end]+=nu

	# Compute steady state solution of Master equation
	# MS = try svd(M)
	# p = nullspace(M)[:,1]
	# p /= sum(p)
	p = normalizednullspace(M)
	# Marginalize over G states
	mhist = marginalize(p,n,nhist)
	alleleConvolve(mhist,nalleles)
end


function sstelFISH(n::Int,rin::Vector{Float64},nhist::Int,nalleles::Int)
	# Closed form steady state mRNA distribution of generalized n-state classic telegraph model
	r = copy(rin)
	nT = n+1
	reactiondecay = 2*n + 2
	r /= rin[reactiondecay]
	total = nhist + 10
	mhist = zeros(nhist)
	# Transition rates
	gammap = zeros(n+2)
	gamman = zeros(n+2)
	for i = 1:n
		gammap[i+1] = r[2*(i-1)+1]
		gamman[i+1] = r[2*i]
	end
	# mRNA creation rate
	nu = r[2*n+1]
	if n > 0
		# Declare Transition matrices
		T = zeros(nT,nT)
		B = zeros(nT,nT)
		S = zeros(total,total)
		Sminus = similar(S)
		Splus = similar(S)
		# Generate T = A + B matrix and B matrix
		for i=1:n+1,  ip=1:n+1
			T[i,ip] =  gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip)
			B[i,ip] =  nu*float(i==ip)*float(i==n+1)
		end
		# Generate matrices for m transitions
		for m=1:total,mp=1:total
			S[m,mp] = -(mp-1)*float(m==mp)
			Sminus[m,mp] = float(m==mp+1)
			Splus[m,mp] = (mp-1)*float(m==mp-1)
		end
		# Construct m truncated Master equation rate matrix
		M = kron(S,Matrix{Float64}(I,nT,nT)) + kron(Matrix{Float64}(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,Matrix{Float64}(I,nT,nT))
		# boundary condition to ensure probability is conserved
		M[end,end]+=nu

		# Compute steady state solution of Master equation
		# MS = try svd(M)
		p = nullspace(M)[:,1]
		# catch
		# 	tol = 1e-8
		# 	return telegraph(n,r,100000000,tol,nhist,nalleles)
		# end
		# p = MS.Vt[end,:]
		p /= sum(p)
		# Marginalize over G states
		for n in 1:nhist
			i = (n-1)*nT
			mhist[n] = sum(p[i+1:i+nT])
		end
	else
		# 0-state telegraph model
		# total = nhist + 10
		S = zeros(total,total)
		for m=1:total,mp=1:total
			S[m,mp] = nu *(m==mp+1) - nu*(m==mp) + (mp-1)*(m==mp-1) - (mp-1)*(m==mp)
		end
		S[end,end] += nu
		mhist = nullspace(S)[:,1]
		# catch
		# 	tol = 1e-5
		# 	mhist = telegraph(n,r,10000000,tol,nhist,2)
		# end
		mhist /= sum(mhist)
	end
	# Convolve to compute distribution for contributions from multiple alleles
	mhists = Array{Array{Float64,1}}(undef,nalleles)
	mhists[1] = convert.(Float64,mhist)
	for i = 2:nalleles
		mhists[i] = zeros(nhist)
		for m = 0:nhist-1
			for m2 = 0:min(nhist-1,m)
				mhists[i][m+1] += mhists[i-1][m-m2+1]*mhist[m2+1]
			end
		end
	end
	return mhists[nalleles]
end


# function transientFISH(t::Float64,n::Int,P0,Mvals,Mvects,rin::Vector{Float64},nhist::Int,nalleles::Int)
# 	# Closed form steady state mRNA distribution of generalized n-state classic telegraph model
# 	P = transientMaster(t,n,Mvals,Mvects)
# 	mhist= marginalize(P,nhist,n)
# 	alleleConvolve(mhist,nalleles)
# end
