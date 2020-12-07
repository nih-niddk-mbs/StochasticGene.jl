"""
Transition rate matrices for G-R-S transcription model
n+1 G states, nr R steps, and nr S sites
"""

"""
transition_rate_mat(nu::Float64,nhist::Int)
Transition rate matrices for G=1 (n=0) transcription model
(i.e. Poisson model)
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
transition_rate_mat(n::Int,gammap::Array,gamman::Vector,nu::Float64,nhist::Int)
Transition rate matrices for G transcription model (i.e. classic telegraph model with n+1 states)
with mRNA decay time set to 1
"""
function transition_rate_mat(n::Int,gammap::Array,gamman::Vector,nu::Float64,nhist::Int)
    nT = n+1
	total = nhist + 2
	T = zeros(nT,nT)
	B = zeros(nT,nT)
	S = zeros(total,total)
	Sminus = similar(S)
	Splus = similar(S)
	# Generate T = A + B matrix and B matrix
	for ip=1:nT, i=1:nT
		T[i,ip] =  gammap[i]*(i-1==ip) - (gamman[i]+gammap[i+1])*(i==ip) + gamman[i+1]*(i+1==ip)
		B[i,ip] =  nu*(i==ip)*(i==n+1)
	end
	# Generate matrices for m transitions
	for m=1:total,mp=1:total
		S[m,mp] = -(mp-1)*(m==mp)
		Sminus[m,mp] = (m==mp+1)
		Splus[m,mp] = (mp-1)*(m==mp-1)
	end
	# Construct m truncated Master equation rate matrix
	M = kron(S,Matrix{Float64}(I,nT,nT)) + kron(Matrix{Float64}(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,Matrix{Float64}(I,nT,nT))
	M[end,end]+=nu  # boundary condition to ensure probability is conserved
	return M
end

"""
transition_rate_mat(n::Int,gammap::Array,gamman::Vector,nu::Float64,nhist::Int)
Transition rate matrices for G transcription model (i.e. classic telegraph model with n+1 states)
"""
function transition_rate_mat(n::Int,gammap::Array,gamman::Vector,nu::Float64,mu::Float64,nhist::Int)
    nT = n+1
	total = nhist + 2
	T = zeros(nT,nT)
	B = zeros(nT,nT)
	S = zeros(total,total)
	Sminus = similar(S)
	Splus = similar(S)
	# Generate T = A + B matrix and B matrix
	for ip=1:nT, i=1:nT
		T[i,ip] =  gammap[i]*(i-1==ip) - (gamman[i]+gammap[i+1])*(i==ip) + gamman[i+1]*(i+1==ip)
		B[i,ip] =  nu*(i==ip)*(i==n+1)
	end
	# Generate matrices for m transitions
	for m=1:total,mp=1:total
		S[m,mp] = -mu*(mp-1)*(m==mp)
		Sminus[m,mp] = (m==mp+1)
		Splus[m,mp] = mu*(mp-1)*(m==mp-1)
	end
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
    	T = zeros(nT,nT)
    	B = zeros(nT,nT)
    	# Generate T = A + B transition matrix and B transition matrix
        # R states are a Base 2 number, e.g. 0100, means presence of polymerase at step 2, other steps are empty
    	for w=1:2^nr, ip=1:n+1, z=1:2^nr, i=1:n+1
    		a = i + (n+1)*(z-1)
    		b = ip + (n+1)*(w-1)
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
    		T[a,b] =  (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)) + nu[nr+1]*float((i==ip)*(zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1))
    		B[a,b] =  nu[nr+1]*float((i==ip)*(zbarr==wbarr)*(zr==0)*(wr==1))
    		for j = 1:nr-1
    			zbarj = zdigits[[1:j-1;j+2:nr]]
    			wbarj = wdigits[[1:j-1;j+2:nr]]
    			zj = zdigits[j]
    			zj1 = zdigits[j+1]
    			wj = wdigits[j]
    			wj1 = wdigits[j+1]
    			T[a,b] += nu[j+1]*float((i==ip)*(zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
    		end
    	end
    	return T, B
    end

    """
    transition_rate_mat(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)
    transtion matrix for GRS Master equation
    where introns are spliced
    returns T, TA, and TI
    """
    function transition_rate_mat(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)
        # standard GR intron reporter model
        nA = (n+1)*3^nr
        nI = (n+1)*3^nr
        # Declare Transition matrices
        T = zeros(nA,nA)
        TA = zeros(nA,nA)
        TI = zeros(nI,nI)
        # Generate TA, TI, and T matrix
        # R states are a base 3 number, e.g. 01200, which means pol2 no intron at step 2, pol2 and intron at step 3, and empty elsehwere
        for w=1:3^nr,ip=1:n+1,z=1:3^nr,i=1:n+1
            a = i + (n+1)*(z-1)
            b = ip + (n+1)*(w-1)
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
            TA[a,b] = ((gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==2)-(z1==0))*(w1==0)) + nu[nr+1]*float((i==ip)*(zbarr==wbarr)*(((zr==0)-(zr==1))*(wr==1)+((zr==0)-(zr==2))*(wr==2))) + eta[nr]*float((i==ip)*(zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2)))*float(any(wdigits.>1))
            T[a,b]  =  (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==2)-(z1==0))*(w1==0)) + nu[nr+1]*float((i==ip)*(zbarr==wbarr)*(((zr==0)-(zr==1))*(wr==1)+((zr==0)-(zr==2))*(wr==2))) + eta[nr]*float((i==ip)*(zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2))
            TI[a,b] = ((gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==2)-(z1==0))*(w1==0)) + nu[nr+1]*float((i==ip)*(zbarr==wbarr)*(((zr==0)-(zr==1))*(wr==1)+((zr==0)-(zr==2))*(wr==2))) + eta[nr]*float((i==ip)*(zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2)))*float(~any(wdigits.>1))
            for j = 1:nr-1
                zbarj = zdigits[[1:j-1;j+2:nr]]
                wbarj = wdigits[[1:j-1;j+2:nr]]
                zbark = zdigits[[1:j-1;j+1:nr]]
                wbark = wdigits[[1:j-1;j+1:nr]]
                zj = zdigits[j]
                zj1 = zdigits[j+1]
                wj = wdigits[j]
                wj1 = wdigits[j+1]
                TA[a,b] += (nu[j+1]*float((i==ip)*(zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0))) + eta[j]*float((i==ip)*(zbark==wbark)*((zj==1)-(zj==2))*(wj==2)))*float(any(wdigits.>1))
                T[a,b]  +=  nu[j+1]*float((i==ip)*(zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0))) + eta[j]*float((i==ip)*(zbark==wbark)*((zj==1)-(zj==2))*(wj==2))
                TI[a,b] += (nu[j+1]*float((i==ip)*(zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+0*((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0))) + eta[j]*float((i==ip)*(zbark==wbark)*((zj==1)-(zj==2))*(wj==2)))*float(~any(wdigits.>1))
            end
        end
        return T,TA,TI
    end

	function transition_rate_mat_T(n::Int,nr::Int,gammap::Vector,gamman::Vector,nu::Vector,eta::Vector)
        # standard GR intron reporter model
        nA = (n+1)*3^nr
        # Declare Transition matrices
        T = zeros(nA,nA)
        # Generate TA, TI, and T matrix
        # R states are a base 3 number, e.g. 01200, which means pol2 no intron at step 2, pol2 and intron at step 3, and empty elsehwere
        for w=1:3^nr,ip=1:n+1,z=1:3^nr,i=1:n+1
            a = i + (n+1)*(z-1)
            b = ip + (n+1)*(w-1)
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
            T[a,b]  =  (gammap[i]*float(i-1==ip) - (gamman[i]+gammap[i+1])*float(i==ip) + gamman[i+1]*float(i+1==ip))*float(z==w) + nu[1]*float((i==ip)*(i==n+1)*(zbar1==wbar1)*((z1==2)-(z1==0))*(w1==0)) + nu[nr+1]*float((i==ip)*(zbarr==wbarr)*(((zr==0)-(zr==1))*(wr==1)+((zr==0)-(zr==2))*(wr==2))) + eta[nr]*float((i==ip)*(zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2))
            for j = 1:nr-1
                zbarj = zdigits[[1:j-1;j+2:nr]]
                wbarj = wdigits[[1:j-1;j+2:nr]]
                zbark = zdigits[[1:j-1;j+1:nr]]
                wbark = wdigits[[1:j-1;j+1:nr]]
                zj = zdigits[j]
                zj1 = zdigits[j+1]
                wj = wdigits[j]
                wj1 = wdigits[j+1]
                T[a,b]  +=  nu[j+1]*float((i==ip)*(zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0))) + eta[j]*float((i==ip)*(zbark==wbark)*((zj==1)-(zj==2))*(wj==2))
                end
        end
        return T
    end
