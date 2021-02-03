"""
Transition rate matrices for G-R-S transcription model
n+1 G states, nr R steps, and nr S sites
"""

"""
transition_rate_mat(nu::Float64,nhist::Int)
Full master equation transition rate matrix for (G=1)M (n=0) transcription model
(i.e. Poisson model)
returns S
"""
function transition_rate_mat(nu::Float64,nhist::Int)
	total = nhist + 2
	S = zeros(total,total)
	for m=1:total,mp=1:total
		S[m,mp] = nu *(m==mp+1) - nu*(m==mp) + (mp-1)*(m==mp-1) - (mp-1)*(m==mp)
	end
	S[end,end] += nu
	return S
end

"""
transition_rate_mat(n,gammap,gamman,nu,mu,nhist)
Transition rate matrices for G transcription model (i.e. classic telegraph model with n+1 states)
mRNA decay rate = mu
returns M
"""
function transition_rate_mat(n::Int,gammap::Vector,gamman::Vector,nu::Float64,mu::Float64,nhist::Int)
	T,B = transition_rate_mat(n,gammap,gamman,nu)
	transition_rate_mat(n,T,B,nu,mu,nhist)
end
"""
transition_rate_mat(n,gammap,gamman,nu,nhist)
Full master equation transition rate matrix for GM transcription model (i.e. classic telegraph model with n+1 states)
with mRNA decay time set to 1
returns M
"""
function transition_rate_mat(n::Int,gammap::Vector,gamman::Vector,nu::Float64,nhist::Int)
	T,B = transition_rate_mat(n,gammap,gamman,nu)
	transition_rate_mat(n,T,B,nu,1.,nhist)
end
"""
transition_rate_mat(n::Int,gammap::Vector,gamman::Vector,nu::Float64)

transition rate matrices for G model
returns T and B
"""
function transition_rate_mat(n::Int,gammap::Vector,gamman::Vector,nu::Float64)
    nT = n+1
	T = zeros(nT,nT)
	B = zeros(nT,nT)
	T[1,1] = -(gamman[1]+gammap[2])
	T[1,2] = gamman[2]
	for i=2:nT-1
		# T[i,ip] =  gammap[i]*(ip==i-1) - (gamman[i]+gammap[i+1])*(ip==i) + gamman[i+1]*(ip==i+1)
		T[i,i-1] = gammap[i]
		T[i,i] = -(gamman[i]+gammap[i+1])
		T[i,i+1] = gamman[i+1]
	end
	T[nT,nT-1] = gammap[nT]
	T[nT,nT] = -(gamman[nT]+gammap[nT+1])
	B[nT,nT] = nu
	return T, B
end
"""
transition_rate_mat(n::Int,T::Matrix,B::Matrix,nu::Float64,mu::Float64,nhist::Int)

Construct full GM transition rate matrix using outer product (kron) of T and B matrices
returns M
"""
function transition_rate_mat(n::Int,T::Matrix,B::Matrix,nu::Float64,mu::Float64,nhist::Int)
	nT = n+1
	total = nhist + 2
	S = zeros(total,total)
	Sminus = zeros(total,total)
	Splus = zeros(total,total)
	# Generate matrices for m transitions
	Splus[1,2] = mu
	for m=2:total-1
		S[m,m] = -mu*(m-1)
		Sminus[m,m-1] = 1
		Splus[m,m+1] = mu*m
	end
	S[total,total] = -mu*(total-1)
	Sminus[total,total-1] = 1
	# Construct m truncated Master equation rate matrix
	M = kron(S,Matrix{Float64}(I,nT,nT)) + kron(Matrix{Float64}(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,Matrix{Float64}(I,nT,nT))
	M[end,end]+=nu  # boundary condition to ensure probability is conserved
	return M
end



"""
transition_rate_mat(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector)
transtion matrix for GR (marginalized over m) master equation
returns T and B
"""
function transition_rate_mat(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector)
	# standard GR intron reporter model
	nT = (n+1)*2^nr
	# Declare Transition matrices
	T = transition_rate_mat_G(n,nr,gammap,gamman,2)
	B = zeros(nT,nT)
	# Generate T = A + B transition matrix and B transition matrix
	# R states are a Base 2 number, e.g. 0100, means presence of polymerase at step 2, other steps are empty
	for w=1:2^nr, z=1:2^nr, i=1:n+1
		a = i + (n+1)*(z-1)
		b = i + (n+1)*(w-1)
		zdigits = digits(z-1,base=2,pad=nr)
		wdigits = digits(w-1,base=2,pad=nr)
		z1 = zdigits[1]
		w1 = wdigits[1]
		zr = zdigits[nr]
		wr = wdigits[nr]
		zbar1 = zdigits[2:nr]
		wbar1 = wdigits[2:nr]
		zbarr = zdigits[1:nr-1]
		wbarr = wdigits[1:nr-1]
		T[a,b] +=  nu[1]*((i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0))
		T[a,b] += nu[nr+1]*((zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1))
		B[a,b] =  nu[nr+1]*((zbarr==wbarr)*(zr==0)*(wr==1))
		for j = 1:nr-1
			zbarj = zdigits[[1:j-1;j+2:nr]]
			wbarj = wdigits[[1:j-1;j+2:nr]]
			zj = zdigits[j]
			zj1 = zdigits[j+1]
			wj = wdigits[j]
			wj1 = wdigits[j+1]
			T[a,b] += nu[j+1]*((zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
		end
	end
	return T, B
end
"""
transition_rate_mat_G(n::Int,nr::Int,gammap::Vector,gamman::Vector,base=3)

transition rate matrix for GR(S) model, base = 2 for GR and base = 3 for GRS
returns T
"""
function transition_rate_mat_G(n::Int,nr::Int,gammap::Vector,gamman::Vector,base=3)
	# G state transition matrix for GR and GRS models
	nA = (n+1)*base^nr
	T = zeros(nA,nA)
	# R states are a base 3 number, e.g. 01200, which means pol2 no intron at step 2, pol2 and intron at step 3, and empty elsehwere
	for a = 0 : n+1 :nA-1
		# a = (n+1)*w
		T[a+1,a+1] = -(gamman[1]+gammap[2])
		T[a+1,a+2] = gamman[2]
		for i = 2:n
			T[a+i,a+i-1] = gammap[i]
			T[a+i,a+i] = -(gamman[i]+gammap[i+1])
			T[a+i,a+i+1] = gamman[i+1]
		end
		T[a+n+1,a+n] = gammap[n+1]
		T[a+n+1,a+n+1] = -(gamman[n+1]+gammap[n+2])
	end
	return T
end

"""
transition_rate_mat(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)

transtion matrices for GRS Master equation
(introns are spliced)
returns T, TA, and TI
"""
function transition_rate_mat(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)
	T = transition_rate_mat_T(n,nr,gammap,gamman,nu,eta)
	transition_rate_mat(T,n,nr)
end

"""
transition_rate_mat(T,n::Int,nr::Int)

transtion matrices TA and TI for GRS Master equation given matrix T
(introns are spliced)
returns T, TA, and TI
"""
function transition_rate_mat(T::Matrix,n,nr)

    nA = (n+1)*3^nr

    TA = zeros(nA,nA)
    TI = zeros(nA,nA)

    # R states are a base 3 number, e.g. 01200, which means pol2 no intron at step 2, pol2 and intron at step 3, and empty elsehwere
    for w=1:3^nr,ip=1:n+1,z=1:3^nr,i=1:n+1
        a = i + (n+1)*(z-1)
        b = ip + (n+1)*(w-1)
        zdigits = digits(z-1,base=3,pad=nr)
        wdigits = digits(w-1,base=3,pad=nr)
		TA[a,b] = T[a,b]*(any(wdigits.>1))
		TI[a,b] = T[a,b]*(~any(wdigits.>1))
    end
    return T,TA,TI
end
"""
transition_rate_mat_T(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)

transtion matrix T for GRS Master equation
(introns are spliced)
returns T
"""
function transition_rate_mat_T(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)

	T = transition_rate_mat_G(n,nr,gammap,gamman)
	for w=1:3^nr,z=1:3^nr
		zdigits = digits(z-1,base=3,pad=nr)
		wdigits = digits(w-1,base=3,pad=nr)
		z1 = zdigits[1]
		w1 = wdigits[1]
		zr = zdigits[nr]
		wr = wdigits[nr]
		zbar1 = zdigits[2:nr]
		wbar1 = wdigits[2:nr]
		zbarr = zdigits[1:nr-1]
		wbarr = wdigits[1:nr-1]
		B = nu[nr+1]*((zbarr==wbarr)*(((zr==0)-(zr==1))*(wr==1)+((zr==0)-(zr==2))*(wr==2)))
		C = eta[nr]*((zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2))
		T[(n+1)*z,(n+1)*w] += nu[1]*((zbar1==wbar1)*((z1==2)-(z1==0))*(w1==0))
		for i=1:n+1
			a = i + (n+1)*(z-1)
			b = i + (n+1)*(w-1)
			T[a,b] += B
			T[a,b] += C
			for j = 1:nr-1
				zbarj = zdigits[[1:j-1;j+2:nr]]
				wbarj = wdigits[[1:j-1;j+2:nr]]
				zbark = zdigits[[1:j-1;j+1:nr]]
				wbark = wdigits[[1:j-1;j+1:nr]]
				zj = zdigits[j]
				zj1 = zdigits[j+1]
				wj = wdigits[j]
				wj1 = wdigits[j+1]
				T[a,b]  +=  nu[j+1]*((zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0))) + eta[j]*((zbark==wbark)*((zj==1)-(zj==2))*(wj==2))
			end
		end
	end
	return T
end
