

"""
Element


"""
struct Element
	a::Int
	b::Int
	index::Int
	pm::Int
end


"""
MComponents

"""
struct MComponents
    elementsT::Vector{Element}
    elementsB::Vector{Element}
	nT::Int
	S::SparseMatrixCSC
	Sminus::SparseMatrixCSC
	Splus::SparseMatrixCSC
end

"""
TComponents

"""
struct TComponents
	nT::Int
    elementsT::Vector{Element}
    elementsTA::Vector{Element}
	elementsTI::Vector{Element}
end

"""
 MTComponents

"""
struct MTComponents
	mcomponents::MComponents
	tcomponents::TComponents
end
"""
Indices

"""
struct Indices
	gamma::Vector{Int}
	nu::Vector{Int}
	eta::Vector{Int}
	decay::Int
end

"""
set_indices(ntransitions,G,R,S)
set_indices(ntransitions,G,R)

"""
set_indices(ntransitions,G,R) = set_indices(ntransitions,G,R,R)
set_indices(ntransitions,G,R,S) = Indices(collect(1:ntransitions),collect(ntransitions+1:ntransitions + R + 1 ),collect(ntransitions + R + 2:ntransitions + R + S + 1),ntransitions + R + S + 2)


"""
set_elements_G!


"""
function set_elements_G!(elements,transitions,G,R,base,gamma)
	nT = G*base^R
	k = 1
	for j = 0 : G : nT-1
		set_elements_G!(elements,transitions,j,gamma)
	end
end
function set_elements_G!(elements,transitions,j=0,gamma=collect(1:length(transitions)))
	i = 1
	for t in transitions
		push!(elements,Element(t[1]+j, t[1]+j, gamma[i],-1))
		push!(elements,Element(t[2]+j, t[1]+j, gamma[i],1))
		i += 1
	end
end

"""
set_elements_RS!(elementsT,G,R,S,nu::Vector{Int},eta::Vector{Int},base::Int=3,type="")

"""
set_elements_R!(elementsT,G,R,indices::Indices) = set_elements_R!(elementsT,G,R,indices.nu)
set_elements_R!(elementsT,G,R,nu::Vector{Int}) = set_elements_RS!(elementsT,G,R,0,nu,Int[],2)
set_elements_RS!(elementT,G,R,indices::Indices) = set_elements_RS!(elementT,G,R,indices.nu,indices.eta)
set_elements_RS!(elementsT,G,R,nu::Vector{Int},eta::Vector{Int}) = set_elements_RS!(elementsT,G,R,R,nu,eta,3)


function set_elements_RS!(elementsT,G,R,S,nu::Vector{Int},eta::Vector{Int},base::Int=3,type="")
	if R > 0
		if type == "offeject"
			S = 0
			base = 2
		end
		for w=1:base^R,z=1:base^R
			zdigits = digits(z-1,base=base,pad=R)
			wdigits = digits(w-1,base=base,pad=R)
			z1 = zdigits[1]
			w1 = wdigits[1]
			zr = zdigits[R]
			wr = wdigits[R]
			zbar1 = zdigits[2:R]
			wbar1 = wdigits[2:R]
			zbarr = zdigits[1:R-1]
			wbarr = wdigits[1:R-1]
			sB = 0
			for l in 1:base-1
				sB += (zbarr==wbarr)*((zr==0)-(zr==l))*(wr==l)
			end
			if S > 0
				sC = (zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2)
			end
			for i=1:G
				a = i + G*(z-1)
				b = i + G*(w-1)
				if abs(sB) == 1
					push!(elementsT,Element(a,b,nu[R+1],sB))
				end
				if S > 0 && abs(sC) == 1
					push!(elementsT,Element(a,b,eta[R],sC))
				end
				if type == "offeject"
					s= (zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1)
					if  abs(s) == 1
						push!(elementsT,Element(a,b,eta[R],s))
					end
				end
				for j = 1:R-1
					zbarj = zdigits[[1:j-1;j+2:R]]
					wbarj = wdigits[[1:j-1;j+2:R]]
					zbark = zdigits[[1:j-1;j+1:R]]
					wbark = wdigits[[1:j-1;j+1:R]]
					zj = zdigits[j]
					zj1 = zdigits[j+1]
					wj = wdigits[j]
					wj1 = wdigits[j+1]
					s = 0
					for l in 1:base-1
						s += (zbarj==wbarj)*((zj==0)*(zj1==l)-(zj==l)*(zj1==0))*(wj==l)*(wj1==0)
					end
					if abs(s) == 1
						push!(elementsT,Element(a,b,nu[j+1],s))
					end
					if S > 0
						s = (zbark==wbark)*((zj==1)-(zj==2))*(wj==2)
						if abs(s) == 1
							push!(elementsT,Element(a,b,eta[j],s))
						end
					end
					if type == "offeject"
						s = (zbark==wbark)*((zj==0)-(zj==1))*(wj==1)
						if  abs(s) == 1
							push!(elementsT,Element(a,b,eta[j],s))
						end
					end
				end
			end
			s = (zbar1==wbar1)*((z1==base-1)-(z1==0))*(w1==0)
			if abs(s) == 1
				push!(elementsT,Element(G*z,G*w,nu[1],s))
			end
		end
	end
end

# """
# set_elements_R!
#
# returns T
# """
# set_elements_R!(elementsT,G,R,indices::Indices) = set_elements_R!(elementsT,G,R,indices.nu)
#
# function set_elements_R!(elementsT,G,R,nu::Vector)
# 	if R == 0
# 		set_elements_R!(elementsT,G,nu::Vector)
# 	else
# 		for w=1:2^R, z=1:2^R, i=1:G
# 			a = i + G*(z-1)
# 			b = i + G*(w-1)
# 			zdigits = digits(z-1,base=2,pad=R)
# 			wdigits = digits(w-1,base=2,pad=R)
# 			z1 = zdigits[1]
# 			w1 = wdigits[1]
# 			zr = zdigits[R]
# 			wr = wdigits[R]
# 			zbar1 = zdigits[2:R]
# 			wbar1 = wdigits[2:R]
# 			zbarr = zdigits[1:R-1]
# 			wbarr = wdigits[1:R-1]
# 			s = (i==G)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)
# 			if  abs(s) == 1
# 				push!(elementsT,Element(a,b,nu[1],s))
# 			end
# 			s = (zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1)
# 			if abs(s) == 1
# 				push!(elementsT,Element(a,b,nu[R+1],s))
# 			end
# 			for j = 1:R-1
# 				zbarj = zdigits[[1:j-1;j+2:R]]
# 				wbarj = wdigits[[1:j-1;j+2:R]]
# 				zj = zdigits[j]
# 				zj1 = zdigits[j+1]
# 				wj = wdigits[j]
# 				wj1 = wdigits[j+1]
# 				s = (zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)
# 				if abs(s) == 1
# 					push!(elementsT,Element(a,b,nu[j+1],s))
# 				end
# 			end
# 		end
# 	end
# end
#
# """
# set_elements_RS!(elementsT,G,R,ntransitions)
#
# return T
# """
#
# set_elements_RS!(elementT,G,R,indices::Indices) = set_elements_RS!(elementT,G,R,indices.nu,indices.eta)
#
# function set_elements_RS!(elementsT,G,R,nu::Vector{Int},eta::Vector{Int})
# 	for w=1:3^R,z=1:3^R
# 		zdigits = digits(z-1,base=3,pad=R)
# 		wdigits = digits(w-1,base=3,pad=R)
# 		z1 = zdigits[1]
# 		w1 = wdigits[1]
# 		zr = zdigits[R]
# 		wr = wdigits[R]
# 		zbar1 = zdigits[2:R]
# 		wbar1 = wdigits[2:R]
# 		zbarr = zdigits[1:R-1]
# 		wbarr = wdigits[1:R-1]
# 		sB = (zbarr==wbarr)*(((zr==0)-(zr==1))*(wr==1)+((zr==0)-(zr==2))*(wr==2))
# 		sC = (zbarr==wbarr)*((zr==1)-(zr==2))*(wr==2)
# 		s = (zbar1==wbar1)*((z1==2)-(z1==0))*(w1==0)
# 		if abs(s) == 1
# 			push!(elementsT,Element(G*z,G*w,nu[1],s))
# 		end
# 		for i=1:G
# 			a = i + G*(z-1)
# 			b = i + G*(w-1)
# 			if abs(sB) == 1
# 				push!(elementsT,Element(a,b,nu[R+1],sB))
# 			end
# 			if abs(sC) == 1
# 				push!(elementsT,Element(a,b,eta[R],sC))
# 			end
# 			# T[a,b] += B
# 			# T[a,b] += C
# 			for j = 1:R-1
# 				zbarj = zdigits[[1:j-1;j+2:R]]
# 				wbarj = wdigits[[1:j-1;j+2:R]]
# 				zbark = zdigits[[1:j-1;j+1:R]]
# 				wbark = wdigits[[1:j-1;j+1:R]]
# 				zj = zdigits[j]
# 				zj1 = zdigits[j+1]
# 				wj = wdigits[j]
# 				wj1 = wdigits[j+1]
# 				s = (zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0))
# 				if abs(s) == 1
# 					push!(elementsT,Element(a,b,nu[j+1],s))
# 				end
# 				s = (zbark==wbark)*((zj==1)-(zj==2))*(wj==2)
# 				if abs(s) == 1
# 					push!(elementsT,Element(a,b,eta[j],s))
# 				end
# 				# T[a,b] += nu[j+1]*((zbarj==wbarj)*(((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)+((zj==0)*(zj1==2)-(zj==2)*(zj1==0))*(wj==2)*(wj1==0)))
# 				# T[a,b] += eta[j]*((zbark==wbark)*((zj==1)-(zj==2))*(wj==2))
# 			end
# 		end
# 	end
# end

# """
# set_elements_RS_offeject!(G,R,ntransitions)
#
# """
# set_elements_R_offeject!(elementsT,G,R,indices::Indices) = set_elements_R_offeject!(elementsT,G,R,indices.nu,indices.eta)
#
#
# function set_elements_R_offeject!(elementsT,G,R,nu::Vector{Int},eta::Vector{Int})
# 	for w=1:2^R, z=1:2^R, i=1:G
# 		a = i + G*(z-1)
# 		b = i + G*(w-1)
# 		zdigits = digits(z-1,base=2,pad=R)
# 		wdigits = digits(w-1,base=2,pad=R)
# 		z1 = zdigits[1]
# 		w1 = wdigits[1]
# 		zr = zdigits[R]
# 		wr = wdigits[R]
# 		zbar1 = zdigits[2:R]
# 		wbar1 = wdigits[2:R]
# 		zbarr = zdigits[1:R-1]
# 		wbarr = wdigits[1:R-1]
# 		s = (i==G)*(zbar1==wbar1)*((z1==1)-(z1==0))*(w1==0)
# 		if  abs(s) == 1
# 			push!(elementsT,Element(a,b,nu[1],s))
# 		end
# 		# T[a,b] += nu[R+1]*(zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1)
# 		s = (zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1)
# 		if  abs(s) == 1
# 			push!(elementsT,Element(a,b,nu[R+1],s))
# 		end
# 		# T[a,b] += eta[nr]*(zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1)
# 		s= (zbarr==wbarr)*((zr==0)-(zr==1))*(wr==1)
# 		if  abs(s) == 1
# 			push!(elementsT,Element(a,b,eta[R],s))
# 		end
# 		for j = 1:R-1
# 			zbarj = zdigits[[1:j-1;j+2:R]]
# 			wbarj = wdigits[[1:j-1;j+2:R]]
# 			zbark = zdigits[[1:j-1;j+1:R]]
# 			wbark = wdigits[[1:j-1;j+1:R]]
# 			zj = zdigits[j]
# 			zj1 = zdigits[j+1]
# 			wj = wdigits[j]
# 			wj1 = wdigits[j+1]
# 			# T[a,b] += nu[j+1]*((zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0))
# 			s = (zbarj==wbarj)*((zj==0)*(zj1==1)-(zj==1)*(zj1==0))*(wj==1)*(wj1==0)
# 			if  abs(s) == 1
# 				push!(elementsT,Element(a,b,nu[j+1],s))
# 			end
# 			# T[a,b] += eta[j]*((zbark==wbark)*((zj==0)-(zj==1))*(wj==1))
# 			s = (zbark==wbark)*((zj==0)-(zj==1))*(wj==1)
# 			if  abs(s) == 1
# 				push!(elementsT,Element(a,b,eta[j],s))
# 			end
# 		end
# 	end
# end
"""
set_elements_TA!(elementsTA,elementsT)


"""
function set_elements_TA!(elementsTA,elementsT,G)
	for e in elementsT
		if e.b == G
			push!(elementsTA,e)
		end
	end
end
function set_elements_TA!(elementsTA,elementsT,G,R,base=3)
	for e in elementsT
		wdigits = digits(div(e.b-1,G),base=base,pad=R)
		if any(wdigits .> base-2)
			push!(elementsTA,e)
		end
	end
end
"""
set_elements_TI!(elementsTA,elementsT)


"""
function set_elements_TI!(elementsTI,elementsT,G)
	for e in elementsT
		if e.b < G
			push!(elementsTI,e)
		end
	end
end
function set_elements_TI!(elementsTI,elementsT,G,R,base=3)
	for e in elementsT
		wdigits = digits(div(e.b-1,G),base=base,pad=R)
		if ~any(wdigits .> base-2)
			push!(elementsTI,e)
		end
	end
end

"""
set_elements_T(transitions,gamma)
set_elements_T(transitions,G,R,base,f!,indices)

"""
function set_elements_T(transitions,gamma::Vector)
	elementsT = Vector{Element}(undef,0)
	set_elements_G!(elementsT,transitions,0,gamma)
	elementsT
end
function set_elements_T(transitions,G,R,base,f!,indices::Indices)
	elementsT = Vector{Element}(undef,0)
	set_elements_G!(elementsT,transitions,G,R,base,indices.gamma)
	f!(elementsT,G,R,indices)
	elementsT
end


"""
set_elements_B


"""
set_elements_B(G,ejectindex) = [Element(G,G,ejectindex,1)]

function set_elements_B(G,R,ejectindex,base=2)
	elementsB = Vector{Element}(undef,0)
	for w=1:base^R, z=1:base^R, i=1:G
		a = i + G*(z-1)
		b = i + G*(w-1)
		zdigits = digits(z-1,base=base,pad=R)
		wdigits = digits(w-1,base=base,pad=R)
		zr = zdigits[R]
		wr = wdigits[R]
		zbarr = zdigits[1:R-1]
		wbarr = wdigits[1:R-1]
		s = (zbarr==wbarr)*(zr==0)*(wr==1)
		if abs(s) == 1
			push!(elementsB,Element(a,b,ejectindex,s))
		end
	end
	return elementsB
end

"""
make_components_M

"""
function make_components_M(transitions,nT,total,decay)
	ntransitions = length(transitions)
	elementsT = set_elements_T(transitions,collect(1:ntransitions))
	elementsB = set_elements_B(nT,ntransitions+1)
	S,Sm,Sp = make_mat_S(total,decay)
	MComponents(elementsT,elementsB,nT,S,Sm,Sp)
end
function make_components_M(transitions,G,R,total,decay,base,indices)
	if R == 0
		return make_components_M(transitions,G,total,decay)
	else
		elementsT = set_elements_T(transitions,G,R,base,set_elements_R!,indices)
		elementsB = set_elements_B(G,R,indices.nu[R+1])
		S,Sm,Sp = make_mat_S(total,decay)
		return MComponents(elementsT,elementsB,G*base^R,S,Sm,Sp)
	end
end

"""
make_components_TAI

"""
function make_components_TAI(elementsT,G)
	elementsTA = Vector{Element}(undef,0)
	elementsTI = Vector{Element}(undef,0)
	set_elements_TA!(elementsTA,elementsT,G)
	set_elements_TI!(elementsTI,elementsT,G)
	TComponents(G,elementsT,elementsTA,elementsTI)
end
function make_components_TAI(elementsT,G,R,base)
	elementsTA = Vector{Element}(undef,0)
	elementsTI = Vector{Element}(undef,0)
	set_elements_TA!(elementsTA,elementsT,G,R,base)
	set_elements_TI!(elementsTI,elementsT,G,R,base)
	TComponents(G*base^R,elementsT,elementsTA,elementsTI)
end

"""
make_components(transitions,nT,r,data::RNAData)
make_components(transitions,G,R,r,data::RNALiveCellData,indices::Indices)
make_components(transitions,G,R,r,total,indices::Indices)
make_components(transitions,G,R,r,data::RNALiveCellData,type::String,indices::Indices)
make_components(transitions,G,R,r,total::Int,type::String,indices::Indices)

make the appropriate matrix components

"""
make_components(transitions,G,R,r,data::RNALiveCellData,indices::Indices) = make_components(transitions,G,R,r,data.nRNA,indices)
make_components(transitions,G,r,data::RNALiveCellData,indices::Indices) = make_components(transitions,G,r,data.nRNA+2,indices)
make_components(transitions,nT,r,data::RNAData) = make_components_M(transitions,nT,data.nRNA+2,r[2*nT])
make_components(transitions,G,R,r,data::RNALiveCellData,type::String,indices::Indices) = make_components(transitions,G,R,r,data.nRNA+2,type,indices)


function make_components(transitions,nT,r,data::RNAData{T1,T2}) where {T1 <: Array, T2 <: Array}
	c = Array{Mcomponents}(undef,0)
	for i in eachindex(data.nRNA)
		push!(c,make_components_M(transitions,nT,data.nRNA[i],r[2*nT*i]))
	end
	c
end

make_components(transitions,G,R,r,total,indices::Indices) = MTComponents(make_components_M(transitions,G,R,total,r[indices.decay],2,indices),make_components_TAI(set_elements_T(transitions,G,R,2,set_elements_R!,indices),G,R,2))


function make_components(transitions,G,r,total::Int,indices::Indices)
	elementsT = set_elements_T(transitions,indices.gamma)
	S,Sm,Sp = make_mat_S(total,r[indices.decay])
	MTComponents(MComponents(elementsT,set_elements_B(G,indices.nu[1]),G,S,Sm,Sp),make_components_TAI(elementsT,G))
end


function make_components(transitions,G,R,r,total::Int,type::String,indices::Indices)
	if type == "offeject"
		elementsT = set_elements_T(transitions,G,R,2,set_elements_R_offeject!,indices)
	else
		elementsT = set_elements_T(transitions,G,R,2,set_elements_R!,indices)
	end
	S,Sm,Sp = make_mat_S(total,r[indices.decay])
	MTComponents(MComponents(elementsT,set_elements_B(G,R,indices.nu[R+1]),G*2^R,S,Sm,Sp),make_components_TAI(set_elements_T(transitions,G,R,3,set_elements_RS!,indices),G,R,3))
end

"""
make_mat!(T::AbstractMatrix,elements::Vector,rates::Vector)

set matrix elements of T to those in vector argument elements
"""
function make_mat!(T::AbstractMatrix,elements::Vector,rates::Vector)
	for e in elements
		T[e.a,e.b] += e.pm * rates[e.index]
	end
end
"""
make_mat(elements,rates,nT)

create a matrix T of size nTxnT with matrix elements in vector elements

"""
function make_mat(elements::Vector,rates::Vector,nT::Int)
	T = spzeros(nT,nT)
	make_mat!(T,elements,rates)
	return T
end
"""
make_S_mat

"""
function make_mat_S(total::Int,mu::Float64)
	S = spzeros(total,total)
	Sminus = spzeros(total,total)
	Splus = spzeros(total,total)
	# Generate matrices for m transitions
	Splus[1,2] = mu
	for m=2:total-1
		S[m,m] = -mu*(m-1)
		Sminus[m,m-1] = 1
		Splus[m,m+1] = mu*m
	end
	S[total,total] = -mu*(total-1)
	Sminus[total,total-1] = 1
	return S,Sminus,Splus
end
"""
make_M_mat

return M matrix used to compute steady state RNA distribution
"""
function make_mat_M(T::SparseMatrixCSC,B::SparseMatrixCSC,mu::Float64,total::Int)
	S,Sminus,Splus = make_mat_S(total,mu)
	make_mat_M(T,B,S,Sminus,Splus)
end

function make_mat_M(T::SparseMatrixCSC,B::SparseMatrixCSC,S::SparseMatrixCSC,Sminus::SparseMatrixCSC,Splus::SparseMatrixCSC)
	nT = size(T,1)
	total = size(S,1)
	M = kron(S,sparse(I,nT,nT)) + kron(sparse(I,total,total),T-B) + kron(Sminus,B) + kron(Splus,sparse(I,nT,nT))
	M[end-size(B,1)+1:end,end-size(B,1)+1:end] .+= B  # boundary condition to ensure probability is conserved
	return M
end

function make_mat_M(components::MComponents,rates::Vector)
	T = make_mat(components.elementsT,rates,components.nT)
	B = make_mat(components.elementsB,rates,components.nT)
	make_mat_M(T,B,components.S,components.Sminus,components.Splus)
end

"""
make_T_mat


"""
make_mat_T(components::TComponents,rates) = make_mat(components.elementsT,rates,components.nT)

make_mat_TA(components::TComponents,rates) = make_mat(components.elementsTA,rates,components.nT)

make_mat_TI(components::TComponents,rates) = make_mat(components.elementsTI,rates,components.nT)
