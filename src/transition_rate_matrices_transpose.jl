# This file is part of StochasticGene.jl
#
# transition_rate_matrices_transpose.jl
#
# Functions to compute the transpose of the transition rate matrix for Stochastic transitions
# Matrices are the transpose to the convention for Kolmogorov forward equations
# i.e. reaction directions go from column to row,
# Notation:  M = transpose of transition rate matrix of full (countably infinite) stochastic process including mRNA kinetics
# Q = finite transition rate matrix for GRS continuous Markov process
# T = Q'
# transpose of Q is more convenient for solving chemical master equation
#

"""
	struct Element

structure for transition matrix elements
fields: 
`a`: row, 
`b`: column
`index`: rate vector index
`pm`: sign of elements

"""
struct Element
    a::Int
    b::Int
    index::Int
    pm::Int
end


"""
	struct MComponents

structure for matrix components for matrix M
"""
struct MComponents
    elementsT::Vector{Element}
    elementsB::Vector{Element}
    nT::Int
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end

"""
	AbstractTComponents

abstract type for components of transition matrices 
"""
abstract type AbstractTComponents end

"""
	struct TComponents

fields:
nT, elementsT

"""
struct TComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector
end

"""
	struct {elementType}

fields:
nT, elementsT, elementsTA, elementsTI

"""
struct TAIComponents{elementType} <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTA::Vector{elementType}
    elementsTI::Vector{elementType}
end
"""
	struct TDComponents

fields:
nT, elementsT, elementsTD:: Vector{Vector{Element}}

"""
struct TDComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTD::Vector{Vector{Element}}
end


"""
 	MTComponents

 structure for M and T components

fields:
mcomponents::MComponents
tcomponents::TComponents
"""
struct MTComponents
    mcomponents::MComponents
    tcomponents::TComponents
end
"""
 	MTComponents

 structure for M and TAI components

fields:
mcomponents::MComponents
tcomponents::TAIComponents
"""
struct MTAIComponents{elementType}
    mcomponents::MComponents
    tcomponents::TAIComponents{elementType}
end
"""
 	MTDComponents

 structure for M and TD components

fields:
mcomponents::MComponents
tcomponents::TDComponents
"""
struct MTDComponents
    mcomponents::MComponents
    tcomponents::TDComponents
end
"""
	struct Indices

index ranges for rates
"""
struct Indices
    gamma::Vector{Int}
    nu::Vector{Int}
    eta::Vector{Int}
    decay::Int
end


"""
    make_components_MTAI(transitions,G,R,S,insertstep,onstates,splicetype="")

make MTAI structure for GRS models (fitting mRNA, on time, and off time histograms)
"""
function make_components_MTAI(transitions, G, R, S, insertstep, onstates, nhist, decay, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    MTAIComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), make_components_TAI(elementsT, nT, onstates))
end


"""
    make_components_MT(transitions,G,R,S,insertstep,decay,splicetype="")

return MTcomponents structure for GRS models 

for fitting traces and mRNA histograms
"""
make_components_MT(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MTcomponents(make_components_M(transitions, G, R, nhist, decay, splicetype), make_components_T(transitions, G, R, S, insertstep, splicetype))

"""
    make_components_M(transitions, G, R, insertstep, decay, splicetype)

return Mcomponents structure  

matrix components for fitting mRNA histograms
"""
function make_components_M(transitions, G, R, nhist, decay, splicetype)
    indices = set_indices(length(transitions), R)
    elementsT, nT = set_elements_T(transitions, G, R, 0, 0, indices, splicetype)
    elementsB = set_elements_B(G, R, indices.nu[R+1])
    U, Um, Up = make_mat_U(nhist, decay)
    return MComponents(elementsT, elementsB, nT, U, Um, Up)
end

"""
    make_components_T(transitions, G, R, S, insertstep, splicetype)

return Tcomponent structure 

for fitting traces and also for creating TA and TI matrix components
"""
function make_components_T(transitions, G, R, S, insertstep, splicetype)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    TComponents(nT, elementsT)
end

"""
    make_components_TAI(elementsT, nT::Int, onstates::Vector)

return TAIComponent structure
"""
function make_components_TAI(elementsT, nT::Int, onstates::Vector)
    TAIComponents{Element}(nT, elementsT, set_elements_TA(elementsT, onstates), set_elements_TI(elementsT, onstates))
end

# function make_components_TAI(elementsT, G::Int, R::Int, base::Int)
#     TAIComponents{typeof(elementsT)}(nT, elementsT, set_elements_TA(elementsT, G, R, base), set_elements_TI(elementsT, G, R, base))
# end

function make_components_TAI(transitions, G, R, S, insertstep, onstates, splicetype::String="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    make_components_TAI(elementsT, nT, onstates)
end


"""
    state_index(G::Int,i,z)

return state index for state (i,z)
"""
state_index(G::Int, i, z) = i + G * (z - 1)



"""
    on_states(G, R, S, insertstep)

return vector of onstates for GRS model given reporter appears at insertstep
"""
function on_states(G, R, S, insertstep)
    base = S > 0 ? 3 : 2
    onstates = Int[]
    on_states!(onstates, G, R, insertstep, base)
    onstates
end

"""
    on_states!(onstates::Vector, G::Int, R::Int, insertstep, base)

return vector of on state indices for GR and GRS models
"""
function on_states!(onstates::Vector, G::Int, R::Int, insertstep, base)
    for i in 1:G, z in 1:base^R
        if any(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1)
            push!(onstates, state_index(G, i, z))
        end
    end
end


"""
    off_states(nT,onstates)

return barrier (off) states, complement of sojurn (on) states
"""
off_states(nT, onstates) = setdiff(collect(1:nT), onstates)

"""
    num_reporters_per_state(G::Int, R::Int, S::Int=0,insertstep=1,f=sum)

return number of a vector of the number reporters for each state index

if f = sum, returns total number of reporters
if f = any, returns 1 for presence of any reporter
"""
function num_reporters_per_state(G::Int, R::Int, S::Int=0, insertstep=1, f=sum)
    base = S > 0 ? 3 : 2
    reporters = Vector{Int}(undef, G * base^R)
    for i in 1:G, z in 1:base^R
        reporters[state_index(G, i, z)] = f(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1)
        # push!(reporters, f(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1))
    end
    reporters
end
function num_reporters_per_state(G::Int, onstates::Vector)
    reporters = Int[]
    for i in 1:G
        push!(reporters, i ∈ onstates ? 1 : 0)
    end
    reporters
end
"""
	set_indices(ntransitions, R, S, insertstep)
	set_indices(ntransitions,R)
	set_indices(ntransitions)

	return Indices structure
"""

function set_indices(ntransitions, R, S, insertstep)
    if insertstep > R > 0
        throw("insertstep>R")
    end
    if S > 0
        Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), collect(ntransitions+R+2:ntransitions+R+R-insertstep+2), ntransitions + R + R - insertstep + 3)
    elseif R > 0
        set_indices(ntransitions, R)
    else
        set_indices(ntransitions)
    end
end
set_indices(ntransitions, R) = Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), Int[], ntransitions + R + 2)
set_indices(ntransitions) = Indices(collect(1:ntransitions), [ntransitions + 1], Int[], ntransitions + 2)

"""
	set_elements_G!(elements,transitions,G,R,base,gamma)

	return G transition matrix elements

-`elements`: Vector of Element structures
-`transitions`: tupe of G state transitions
"""
function set_elements_G!(elements, transitions, G, gamma, nT)
    for j = 0:G:nT-1
        set_elements_G!(elements, transitions, gamma, j)
    end
end
function set_elements_G!(elements, transitions, gamma=collect(1:length(transitions)), j=0)
    i = 1
    for t in transitions
        push!(elements, Element(t[1] + j, t[1] + j, gamma[i], -1))
        push!(elements, Element(t[2] + j, t[1] + j, gamma[i], 1))
        i += 1
    end
end

"""
    set_elements_RS!(elementsT, G, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")

inplace update matrix elements in elementsT for GRS state transition matrix, do nothing if R == 0

"offeject" = pre-RNA is completely ejected when spliced
"""
function set_elements_RS!(elementsT, G, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")
    if R > 0
        if splicetype == "offeject"
            S = 0
            base = 2
        end
        if S == 0
            base = 2
        else
            base = 3
        end
        for w = 1:base^R, z = 1:base^R
            zdigits = digit_vector(z, base, R)
            wdigits = digit_vector(w, base, R)
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
                sB += (zbarr == wbarr) * ((zr == 0) - (zr == l)) * (wr == l)
            end
            if S > 0
                sC = (zbarr == wbarr) * ((zr == 1) - (zr == 2)) * (wr == 2)
            end
            for i = 1:G
                # a = i + G * (z - 1)
                # b = i + G * (w - 1)
                a = state_index(G, i, z)
                b = state_index(G, i, w)
                if abs(sB) == 1
                    push!(elementsT, Element(a, b, nu[R+1], sB))
                end
                if S > 0 && abs(sC) == 1
                    push!(elementsT, Element(a, b, eta[R-insertstep+1], sC))
                end
                if splicetype == "offeject"
                    s = (zbarr == wbarr) * ((zr == 0) - (zr == 1)) * (wr == 1)
                    if abs(s) == 1
                        push!(elementsT, Element(a, b, eta[R-insertstep+1], s))
                    end
                end
                for j = 1:R-1
                    zbarj = zdigits[[1:j-1; j+2:R]]
                    wbarj = wdigits[[1:j-1; j+2:R]]
                    zbark = zdigits[[1:j-1; j+1:R]]
                    wbark = wdigits[[1:j-1; j+1:R]]
                    zj = zdigits[j]
                    zj1 = zdigits[j+1]
                    wj = wdigits[j]
                    wj1 = wdigits[j+1]
                    s = 0
                    for l in 1:base-1
                        s += (zbarj == wbarj) * ((zj == 0) * (zj1 == l) - (zj == l) * (zj1 == 0)) * (wj == l) * (wj1 == 0)
                    end
                    if abs(s) == 1
                        push!(elementsT, Element(a, b, nu[j+1], s))
                    end
                    if S > 0 && j > insertstep - 1
                        s = (zbark == wbark) * ((zj == 1) - (zj == 2)) * (wj == 2)
                        if abs(s) == 1
                            push!(elementsT, Element(a, b, eta[j-insertstep+1], s))
                        end
                    end
                    if splicetype == "offeject" && j > insertstep - 1
                        s = (zbark == wbark) * ((zj == 0) - (zj == 1)) * (wj == 1)
                        if abs(s) == 1
                            push!(elementsT, Element(a, b, eta[j-insertstep+1], s))
                        end
                    end
                end
            end
            s = (zbar1 == wbar1) * ((z1 == base - 1) - (z1 == 0)) * (w1 == 0)
            if abs(s) == 1
                push!(elementsT, Element(G * z, G * w, nu[1], s))
            end
        end
    end
end

"""
    set_elements_TA(elementsT,onstates)

set onstate elements
"""
set_elements_TA(elementsT, onstates) = set_elements_TX(elementsT, onstates, set_elements_TA!)
"""
    set_elements_TA!(elementsTA, elementsT, onstates::Vector)

in place set onstate elements
"""
set_elements_TA!(elementsTA, elementsT, onstates) = set_elements_TX!(elementsTA, elementsT, onstates, ∈)

"""
    set_elements_TI!(elementsTI, elementsT, onstates::Vector)

set off state elements
"""
set_elements_TI(elementsT, onstates) = set_elements_TX(elementsT, onstates, set_elements_TI!)
"""
    set_elements_TI!(elementsTI, elementsT, onstates::Vector)

in place set off state elements
"""
set_elements_TI!(elementsTI, elementsT, onstates) = set_elements_TX!(elementsTI, elementsT, onstates, ∉)

"""
    set_elements_TX(elementsT, onstates::Vector{Int},f!)

set on or off state elements for onstates
"""
function set_elements_TX(elementsT, onstates::Vector{Int}, f!)
    elementsTX = Vector{Element}(undef, 0)
    f!(elementsTX, elementsT, onstates::Vector)
    elementsTX
end

"""
    set_elements_TX(elementsT, onstates::Vector{Vector{Int}},f!)

set on or off state elements for vector if onstates by calling f! (set_elements_TA! or set_elements_TI!)
"""
function set_elements_TX(elementsT, onstates::Vector{Vector{Int}}, f!)
    elementsTX = Vector{Vector{Element}}(undef, length(onstates))
    for i in eachindex(onstates)
        eTX = Vector{Element}(undef, 0)
        f!(eTX, elementsT, onstates::Vector)
        elementsTX[i] = eTX
    end
    elementsTX
end

"""
    set_elements_TX!(elementsTX, elementsT, onstates::Vector,f)

in place set on or off state elements by calling f (∈ or ∉)
"""
function set_elements_TX!(elementsTX, elementsT, onstates::Vector, f)
    for e in elementsT
        if f(e.b, onstates)
            push!(elementsTX, e)
        end
    end
end


# function set_elements_TA(elementsT, onstates::Vector{Int})
#     elementsTA = Vector{Element}(undef, 0)
#     set_elements_TA!(elementsTA, elementsT, onstates::Vector)
#     elementsTA
# end

# function set_elements_TA(elementsT, onstates::Vector{Vector{Int}})
#     elementsTA = Vector{Vector{Element}}(undef,length(onstates))
#     for i in eachindex(onstates)
#         eTA = Vector{Element}(undef, 0)
#         set_elements_TA!(eTA, elementsT, onstates::Vector)
#         elementsTA[i] = eTA
#     elementsTA
# end
"""
    set_elements_TA!(elementsTA, elementsT, onstates::Vector)


"""
# function set_elements_TA!(elementsTA, elementsT, onstates::Vector)
#     for e in elementsT
#         if e.b ∈ onstates
#             push!(elementsTA, e)
#         end
#     end
# end
# function set_elements_TA!(elementsTA, elementsT, G::Int, R::Int, insertstep=1, base::Int=3)
#     for e in elementsT
#         wdigits = digits(div(e.b - 1, G), base=base, pad=R)[insertstep:end]
#         if any(wdigits .> base - 2)
#             push!(elementsTA, e)
#         end
#     end
# end

"""
    set_elements_TI(elementsT,onstates)


"""
# function set_elements_TI(elementsT, onstates)
#     elementsTI = Vector{Element}(undef, 0)
#     set_elements_TI!(elementsTI, elementsT, onstates::Vector)
#     elementsTI
# end
# function set_elements_TI(elementsT, onstates::Vector{Vector{Int}})
#     elementsTI = Vector{Vector{Element}}(undef,length(onstates))
#     for i in eachindex(onstates)
#         eTi = Vector{Element}(undef, 0)
#         set_elements_TI!(eTI, elementsT, onstates::Vector)
#         elementsTI[i] = eTI
#     elementsTI
# end

# function set_elements_TI!(elementsTI, elementsT, onstates::Vector)
#     for e in elementsT
#         if e.b ∉ onstates
#             push!(elementsTI, e)
#         end
#     end
# end
# function set_elements_TI!(elementsTI, elementsT, G::Int, R::Int, base::Int=3)
#     for e in elementsT
#         wdigits = digits(div(e.b - 1, G), base=base, pad=R)
#         if ~any(wdigits .> base - 2)
#             push!(elementsTI, e)
#         end
#     end
# end


"""
    set_elements_T(transitions, G, R, S, indices::Indices, splicetype::String)

return T matrix elements (Element structure) and size of matrix
"""
function set_elements_T(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    if R > 0
        elementsT = Vector{Element}(undef, 0)
        base = S > 0 ? 3 : 2
        nT = G * base^R
        set_elements_G!(elementsT, transitions, G, indices.gamma, nT)
        set_elements_RS!(elementsT, G, R, S, insertstep, indices.nu, indices.eta, splicetype)
        return elementsT, nT
    else
        return set_elements_T(transitions, indices.gamma), G
    end
end
function set_elements_T(transitions, gamma::Vector)
    elementsT = Vector{Element}(undef, 0)
    set_elements_G!(elementsT, transitions, gamma)
    elementsT
end
"""
    set_elements_B(G, ejectindex)

return B matrix elements
"""
set_elements_B(G, ejectindex) = [Element(G, G, ejectindex, 1)]

function set_elements_B(G, R, ejectindex, base=2)
    if R > 0
        elementsB = Vector{Element}(undef, 0)
        for w = 1:base^R, z = 1:base^R, i = 1:G
            a = i + G * (z - 1)
            b = i + G * (w - 1)
            zdigits = digits(z - 1, base=base, pad=R)
            wdigits = digits(w - 1, base=base, pad=R)
            zr = zdigits[R]
            wr = wdigits[R]
            zbarr = zdigits[1:R-1]
            wbarr = wdigits[1:R-1]
            s = (zbarr == wbarr) * (zr == 0) * (wr == 1)
            if abs(s) == 1
                push!(elementsB, Element(a, b, ejectindex, s))
            end
        end
        return elementsB
    else
        set_elements_B(G, ejectindex)
    end
end
"""
make_mat!(T::AbstractMatrix,elements::Vector,rates::Vector)

set matrix elements of T to those in vector argument elements
"""
function make_mat!(T::AbstractMatrix, elements::Vector, rates::Vector)
    for e in elements
        T[e.a, e.b] += e.pm * rates[e.index]
    end
end
"""
make_mat(elements,rates,nT)

return an nTxnT sparse matrix T
"""
function make_mat(elements::Vector, rates::Vector, nT::Int)
    T = spzeros(nT, nT)
    make_mat!(T, elements, rates)
    return T
end
"""
make_S_mat

"""
function make_mat_U(total::Int, mu::Float64)
    U = spzeros(total, total)
    Uminus = spzeros(total, total)
    Uplus = spzeros(total, total)
    # Generate matrices for m transitions
    Uplus[1, 2] = mu
    for m = 2:total-1
        U[m, m] = -mu * (m - 1)
        Uminus[m, m-1] = 1
        Uplus[m, m+1] = mu * m
    end
    U[total, total] = -mu * (total - 1)
    Uminus[total, total-1] = 1
    return U, Uminus, Uplus
end
"""
make_M_mat

return M matrix used to compute steady state RNA distribution
"""
function make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, mu::Float64, total::Int)
    U, Uminus, Uplus = make_mat_U(total, mu)
    make_mat_M(T, B, U, Uminus, Uplus)
end
function make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, U::SparseMatrixCSC, Uminus::SparseMatrixCSC, Uplus::SparseMatrixCSC)
    nT = size(T, 1)
    total = size(U, 1)
    M = kron(U, sparse(I, nT, nT)) + kron(sparse(I, total, total), T - B) + kron(Uminus, B) + kron(Uplus, sparse(I, nT, nT))
    M[end-size(B, 1)+1:end, end-size(B, 1)+1:end] .+= B  # boundary condition to ensure probability is conserved
    return M
end

function make_mat_M(components::MComponents, rates::Vector)
    T = make_mat(components.elementsT, rates, components.nT)
    B = make_mat(components.elementsB, rates, components.nT)
    make_mat_M(T, B, components.U, components.Uminus, components.Uplus)
end

"""
make_T_mat


"""
make_mat_T(components::AbstractTComponents, rates) = make_mat(components.elementsT, rates, components.nT)

make_mat_TA(components::TAIComponents, rates) = make_mat(components.elementsTA, rates, components.nT)

make_mat_TI(components::TAIComponents, rates) = make_mat(components.elementsTI, rates, components.nT)

make_mat_TA(elementsTA, rates, nT) = make_mat(elementsTA, rates, nT)

make_mat_TI(elementsTI, rates, nT) = make_mat(elementsTI, rates, nT)