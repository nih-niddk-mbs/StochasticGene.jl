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

structure for T transition matrix elements
fields: 
- `a`: row, 
- `b`: column
- `index`: rate vector index
- `pm`: sign of elements

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
struct TAIComponents{elementType} <: AbstractTComponents

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
struct TDComponents <: AbstractTComponents

fields:
nT, elementsT, elementsTD:: Vector{Vector{Element}}

"""
struct TDComponents <: AbstractTComponents
    nT::Int
    elementsT::Vector{Element}
    elementsTG::Vector{Element}
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
gamma: G transitions
nu: R transitions
eta: splice transitions
decay: mRNA decay rate

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

function make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    elementsTG = set_elements_T(transitions, indices.gamma)
    c = Vector{Element}[]
    for i in eachindex(onstates)
        if dttype[i] == "ON"
            push!(c, set_elements_TA(elementsT, onstates[i]))
        elseif dttype[i] == "OFF"
            push!(c, set_elements_TI(elementsT, onstates[i]))
        elseif dttype[i] == "ONG"
            push!(c, set_elements_TA(elementsTG, onstates[i]))
        elseif dttype[i] == "OFFG"
            push!(c, set_elements_TI(elementsTG, onstates[i]))
        end
        # dttype[i] == "ON" ? push!(c, set_elements_TA(elementsT, onstates[i])) : push!(c, set_elements_TI(elementsT, onstates[i]))
    end
    MTDComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), TDComponents(nT, elementsT, elementsTG, c))
end

"""
    make_components_MT(transitions,G,R,S,insertstep,decay,splicetype="")

return MTcomponents structure for GRS models 

for fitting traces and mRNA histograms
"""
make_components_MT(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MTComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), make_components_T(transitions, G, R, S, insertstep, splicetype))

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
state_index(G::Int, g, z) = g + G * (z - 1)

function inverse_state(i::Int, G, R, S, insertstep::Int, f=sum)
    base = S > 0 ? 3 : 2
    g = mod(i - 1, G) + 1
    z = div(i - g, G) + 1
    zdigits = digit_vector(z, base, R)
    r = num_reporters_per_index(z, R, insertstep, base, f)
    return g, z, zdigits, r
end

function inverse_state(i::Vector{Int}, G, R, S, insertstep::Int)
    base = S > 0 ? 3 : 2
    z = Int[]
    g = Int[]
    zdigits = Vector[]
    r = Int[]
    for i in i
        gi, zi, zdi, ri = inverse_state(i, G, R, S, insertstep)
        push!(z, zi)
        push!(g, gi)
        push!(zdigits, zdi)
        push!(r, ri)
    end
    return g, z, zdigits, r
end

function inverse_state(i::Vector{Vector{Int}}, G, R, S, insertstep)
    z = Vector{Int}[]
    g = Vector{Int}[]
    zdigits = Vector{Vector}[]
    r = Vector{Int}[]
    for i in i
        gi, zi, zdi, ri = inverse_state(i, G, R, S, insertstep)
        push!(z, zi)
        push!(g, gi)
        push!(zdigits, zdi)
        push!(r, ri)
    end
    return g, z, zdigits, r
end

"""
    on_states(onstates, G, R, S, insertstep)

return vector of new onstates given onstates argument of fit function
"""
function on_states(onstates, G, R, S, insertstep)
    if isempty(onstates)
        return on_states(G, R, S, insertstep)
    else
        return on_states(onstates, G, R, S)
    end
end
"""
    on_states(G, R, S, insertstep)

return vector of onstates for GRS model given reporter appears at insertstep
"""
function on_states(G::Int, R, S, insertstep)
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
    (R == 0) && throw("Cannot use empty ON state [] for R = 0")
    for i in 1:G, z in 1:base^R
        if any(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1)
            push!(onstates, state_index(G, i, z))
        end
    end
end

function on_states(onstates::Vector, G, R, S)
    base = S > 0 ? 3 : 2
    o = Int[]
    for i in 1:G, z in 1:base^R
        if i ∈ onstates
            push!(o, state_index(G, i, z))
        end
    end
    o
end

function off_states(G::Int, R, S, insertstep)
    base = S > 0 ? 3 : 2
    nT = G * base^R
    off_states(nT, on_states(G, R, S, insertstep))
end
"""
    off_states(nT,onstates)

return barrier (off) states, complement of sojourn (on) states
"""
off_states(nT, onstates) = setdiff(collect(1:nT), onstates)


"""
    T_dimension(G, R, S)

Compute transition matrix dimension of GRS model
"""
function T_dimension(G, R, S)
    base = S > 0 ? 3 : 2
    G * base^R
end

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
        # reporters[state_index(G, i, z)] = f(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1)
        reporters[state_index(G, i, z)] = num_reporters_per_index(z, R, insertstep, base, f)
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

num_reporters_per_index(z, R, insertstep, base, f=sum) = f(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1)
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

function set_indices(ntransitions, R, S, insertstep, coupling)


end
"""
	set_elements_G!(elements,transitions,G,R,base,gamma)

	return G transition matrix elements

-`elements`: Vector of Element structures
-`transitions`: tupe of G state transitions
"""
function set_elements_G!(elements, transitions, G::Int, gamma, nT)
    for j = 0:G:nT-1
        set_elements_G!(elements, transitions, gamma, j)
    end
end
function set_elements_G!(elements, transitions, gamma::Vector=collect(1:length(transitions)), j=0)
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

add elements to elementsTX if column element satisfies f (∈ or ∉)
"""
function set_elements_TX!(elementsTX, elementsT, onstates::Vector, f)
    for e in elementsT
        if f(e.b, onstates)
            push!(elementsTX, e)
        end
    end
end

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
function make_mat_U(total::Int, decay::Float64)
    U = spzeros(total, total)
    Uminus = spzeros(total, total)
    Uplus = spzeros(total, total)
    # Generate matrices for m transitions
    Uplus[1, 2] = decay
    for m = 2:total-1
        U[m, m] = -decay * (m - 1)
        Uminus[m, m-1] = 1
        Uplus[m, m+1] = decay * m
    end
    U[total, total] = -decay * (total - 1)
    Uminus[total, total-1] = 1
    return U, Uminus, Uplus
end
"""
make_M_mat

return M matrix used to compute steady state RNA distribution
"""
function make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, decay::Float64, total::Int)
    U, Uminus, Uplus = make_mat_U(total, decay)
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

construct matrices from elements
"""
make_mat_T(components::AbstractTComponents, rates) = make_mat(components.elementsT, rates, components.nT)

make_mat_TA(components::TAIComponents, rates) = make_mat(components.elementsTA, rates, components.nT)

make_mat_TI(components::TAIComponents, rates) = make_mat(components.elementsTI, rates, components.nT)

make_mat_TA(elementsTA, rates, nT) = make_mat(elementsTA, rates, nT)

make_mat_TI(elementsTI, rates, nT) = make_mat(elementsTI, rates, nT)

function make_mat_GR(components, rates)
    nG = components.nG
    nR = components.nR
    G = make_mat(components.elementsG, rates, nG)
    R = make_mat(components.elementsR, rates, nR)
    RB = make_mat(components.elementsRB, rates, nR)
    G, R, RB
end

##### Coupling Under development

struct T2Components <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    elementsG::Vector
    elementsR::Vector
    elementsRG::Vector
end
struct TCSourceComponents <: AbstractTComponents
    nG::Int
    nR::Int
    elementsG::Vector
    elementsGC::Vector
    elementsRK::Vector
    elementsRKbar::Vector
end
struct TTargetComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    elementsGC::Vector
    elementsGCbar::Vector
    elementsRK::Vector
    elementsRKbar::Vector
end

struct TcoupledComponents <: AbstractTComponents
    nT::Int
    nG::Int
    nR::Int
    elementsG::Vector
    elementsGt::Vector
    elementsRG::Vector
    elementsRGbar::Vector
end

struct M2Components
    nG::Int
    nR::Int
    elementsG::Vector
    elementsRG::Vector
    elementsR::Vector
    elementsB::Vector{Element}
    U::SparseMatrixCSC
    Uminus::SparseMatrixCSC
    Uplus::SparseMatrixCSC
end
struct MT2Components
    mcomponents::M2Components
    tcomponents::T2Components
end

"""
    set_elements_RS2!(elementsRG, elementsR, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")

TBW
"""
function set_elements_RS2!(elementsRG, elementsR, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")
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
        nR = base^R
        for b = 1:nR, a = 1:nR
            zdigits = digit_vector(a, base, R)
            wdigits = digit_vector(b, base, R)
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
            if abs(sB) == 1
                push!(elementsR, Element(a, b, nu[R+1], sB))
            end
            if S > 0 && abs(sC) == 1
                push!(elementsR, Element(a, b, eta[R-insertstep+1], sC))
            end
            if splicetype == "offeject"
                s = (zbarr == wbarr) * ((zr == 0) - (zr == 1)) * (wr == 1)
                if abs(s) == 1
                    push!(elementsR, Element(a, b, eta[R-insertstep+1], s))
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
                    push!(elementsR, Element(a, b, nu[j+1], s))
                end
                if S > 0 && j > insertstep - 1
                    s = (zbark == wbark) * ((zj == 1) - (zj == 2)) * (wj == 2)
                    if abs(s) == 1
                        push!(elementsR, Element(a, b, eta[j-insertstep+1], s))
                    end
                end
                if splicetype == "offeject" && j > insertstep - 1
                    s = (zbark == wbark) * ((zj == 0) - (zj == 1)) * (wj == 1)
                    if abs(s) == 1
                        push!(elementsR, Element(a, b, eta[j-insertstep+1], s))
                    end
                end
            end
            s = (zbar1 == wbar1) * ((z1 == base - 1) - (z1 == 0)) * (w1 == 0)
            if abs(s) == 1
                push!(elementsRG, Element(a, b, nu[1], s))
            end
        end
    end
end


"""
    set_elements_G2!(elements, transitions, gamma::Vector=collect(1:length(transitions)), j=0)

TBW
"""
function set_elements_G2!(elements, transitions, gamma::Vector=collect(1:length(transitions)), j=0)
    i = 1
    for t in transitions
        push!(elements, Element(t[1] + j, t[1] + j, gamma[i], -1))
        push!(elements, Element(t[2] + j, t[1] + j, gamma[i], 1))
        i += 1
    end
end

function set_elements_Gt!(elements, transitions, target_transition=length(transitions), gamma::Vector=collect(1:length(transitions)), j=0)
    i = 1
    for t in transitions
        if i == target_transition
            push!(elements, Element(t[1] + j, t[1] + j, gamma[i], -1))
            push!(elements, Element(t[2] + j, t[1] + j, gamma[i], 1))
        end
        i += 1
    end
end

"""
    set_elements_B2(G, R, ejectindex, base=2)

TBW
"""
function set_elements_B2(G, R, ejectindex, base=2)
    if R > 0
        nR = base^R
        elementsB = Vector{Element}(undef, 0)
        for b = 1:nR, a = 1:nR
            zdigits = digits(a - 1, base=base, pad=R)
            wdigits = digits(b - 1, base=base, pad=R)
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
    set_elements_T2(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)

TBW
"""
function set_elements_T2(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    if R > 0
        elementsG = Vector{Element}(undef, 0)
        elementsR = Vector{Element}(undef, 0)
        elementsRG = Vector{Element}(undef, 0)
        base = S > 0 ? 3 : 2
        nR = base^R
        set_elements_G!(elementsG, transitions)
        set_elements_RS2!(elementsRG, elementsR, R, S, insertstep, indices.nu, indices.eta, splicetype)
        return elementsG, elementsRG, elementsR, nR, T_dimension(G, R, S)
    else
        return set_elements_G(transitions, indices.gamma), G
    end
end

function set_elements_Tcoupled(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    elementsG = Vector{Element}(undef, 0)
    elementsGt = Vector{Element}(undef, 0)
    set_elements_G2!(elementsG, transitions)
    set_elements_Gt!(elementsGt, transitions)
    if R > 0
        elementsRGbar = Vector{Element}(undef, 0)
        elementsRG = Vector{Element}(undef, 0)
        base = S > 0 ? 3 : 2
        nR = base^R
        set_elements_RS2!(elementsRG, elementsRGbar, R, S, insertstep, indices.nu, indices.eta, splicetype)
        return elementsG, elementsGt, elementsRG, elementsRGbar, nR, T_dimension(G, R, S)
    else
        return elementsG, elementsGt, G
    end
end

"""
    make_components_T2(transitions, G, R, S, insertstep, splicetype)

TBW
"""
function make_components_T2(transitions, G, R, S, insertstep, splicetype)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsRG, elementsR, nR, nT = set_elements_T2(transitions, G, R, S, insertstep, indices, splicetype)
    T2Components(nT, G, nR, elementsG, elementsR, elementsRG)
end

function make_components_Tcoupled(transitions, G, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsGt, elementsRG, elementsRGbar, nR, nT = set_elements_Tcoupled(transitions, G, R, S, insertstep, indices, splicetype)
    TcoupledComponents(nT, G, nR, elementsG, elementsGt, elementsRGbar, elementsRG)

end


function make_components_M2(transitions, G, R, nhist, decay)
    indices = set_indices(length(transitions), R)
    elementsG, elementsRG, elementsR, nR = set_elements_T2(transitions, G, R, 0, 1, indices, "")
    elementsB = set_elements_B2(G, R, indices.nu[R+1])
    U, Um, Up = make_mat_U(nhist, decay)
    M2Components(G, nR, elementsG, elementsRG, elementsR, elementsB, U, Um, Up)
end

function make_mat_GC(G, n, coupling=1.0)
    GR = spzeros(G, G)
    GR[n, n] = coupling
    return GR
end

make_components_MT2(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MT2Components(make_components_M2(transitions, G, R, nhist, decay), make_components_T2(transitions, G, R, S, insertstep, splicetype))

function make_mat_T2(components, rates)
    nG = components.nG
    nR = components.nR
    GR = make_mat_GC(nG, nG)
    G = make_mat(components.elementsG, rates, nG)
    R = make_mat(components.elementsR, rates, nR)
    RG = make_mat(components.elementsRG, rates, nR)
    make_mat_T2(G, GR, R, RG, nG, nR)
end

function make_mat_Tcoupled1(components, rates)
    nT = components.nT
    nG = components.nG
    nR = components.nR
    nS = nG
    GR = make_mat_GC(nG, nG)
    G = make_mat(components.elementsG, rates, nG)
    Gt = make_mat(components.elementsGt, rates, nG)
    Gs = make_mat_GC(nS, nS, rates[end])
    RGbar = make_mat(components.elementsRGbar, rates, nR)
    RG = make_mat(components.elementsRG, rates, nR)
    I2 = sparse(I, nT, nT)
    IR = sparse(I, nR, nR)
    T2 = make_mat_T2(G, GR, RGbar, RG, nG, nR)
    return T2, kron(I2, T2) + kron(T2, I2) + kron(kron(IR, Gs), kron(IR, Gt))
end

function make_mat_Tcoupled(components, rates)
    nT = components.nT
    nG = components.nG
    nR = components.nR
    nS = nG
    GR = make_mat_GC(nG, nG)
    G = make_mat(components.elementsG, rates, nG)
    RGbar = make_mat(components.elementsRGbar, rates, nR)
    RG = make_mat(components.elementsRG, rates, nR)
    T = make_mat_T2(G, GR, RGbar, RG, nG, nR)
    Gr = make_mat(components.elementsGt, rates, nG)
    Gs = make_mat_GC(nS, nS, rates[end])
    IR = sparse(I, nR, nR)
    return T, kron(IR,Gr), kron(IR,Gs), nT
end

function make_mat_Tc(ns,Tmats,gamma)
    n = prod(ns)
    Tc = spzeros(n,n)
    for i in eachindex(transcribers)
        n = prod(ns[1:i-1])
        m = prod(ns[i+1:end])
        Tc += kron(kron(sparse(I,n,n),[α[i]]),sparse(I,m,m))
    end
    for β in eachindex(transcribers)
        for α in 1:β-1
            Tc += gamma[α,β]*kron(kron(kron(kron(I1,S[α]),I2),R[β]),I3)
            Tc += gamma[β,α]*kron(kron(kron(kron(I1,R[α]),I2),S[β]),I3)
        end
    end
    return Tc
end

function make_mat_T2(G, GR, R, RG, nG, nR)
    kron(RG, GR) + kron(sparse(I, nR, nR), G) + kron(R, sparse(I, nG, nG))
end

function make_mat_B2(components, rates)
    RB = make_mat(components.elementsB, rates, components.nR)
    nG = components.nG
    kron(RB, sparse(I, nG, nG))
end

function make_mat_M2(components::M2Components, rates::Vector)
    T = make_mat_T2(components, rates)
    B = make_mat_B2(components, rates)
    make_mat_M(T, B, components.U, components.Uminus, components.Uplus)
end

function test_mat2(r, transitions, G, R, S, insertstep, nhist=20)
    components = make_components_MT2(transitions, G, R, S, insertstep, nhist, 1.01, "")
    T = make_mat_T2(components.tcomponents, r)
    M = make_mat_M2(components.mcomponents, r)
    return T, M, components
end

function test_mat(r, transitions, G, R, S, insertstep, nhist=20)
    components = make_components_MT(transitions, G, R, S, insertstep, nhist, 1.01, "")
    T = make_mat_T(components.tcomponents, r)
    M = make_mat_M(components.mcomponents, r)
    T, M, components
end

function test_mat_Tcoupled(r, transitions, G, R, S, insertstep)
    components = make_components_Tcoupled(transitions, G, R, S, insertstep)
    make_mat_Tcoupled(components, r)
end
