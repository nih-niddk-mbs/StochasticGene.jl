# This file is part of StochasticGene.jl
#
# transition_rate_set_elements.jl
#
# This file contains functions for creating and manipulating the individual elements
# that make up transition rate matrices in stochastic gene expression modeling.
# These functions are responsible for generating the specific matrix elements needed
# for different types of transition matrices used throughout the StochasticGene.jl package.
#
# Key functionality includes:
# - Gene state transition elements (G elements)
# - RNA/reporter state transition elements (R elements)
# - Splicing transition elements (S elements)
# - Combined GRS transition elements
# - On/off state filtering elements (TA/TI elements)
# - Dwell time analysis elements (TD elements)
# - Boundary condition elements (B elements)
# - Coupled system elements
# - Source and target state elements
#
# These functions work with the Element structure defined in transition_rate_structures.jl
# and are used by the matrix construction functions in transition_rate_make.jl to build
# the complete transition matrices needed for various types of stochastic gene expression analyses.

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
    set_elements_G(transitions, gamma::Vector)

Create gene state transition elements.

# Description
This function creates Element structures for gene state transitions based on the provided
transition tuples and gamma indices. Each transition tuple should contain [from_state, to_state].

# Arguments
- `transitions`: Tuple of gene state transitions, each containing [from_state, to_state]
- `gamma::Vector`: Vector of rate indices corresponding to each transition

# Returns
- `Vector{Element}`: Vector of Element structures representing gene state transitions
"""
function set_elements_G(transitions, gamma::Vector)
    elementsT = Vector{Element}(undef, 0)
    set_elements_G!(elementsT, transitions, gamma)
    elementsT
end

"""
    set_elements_RS!(elementsRGbar, elementsRG, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")

inplace update to RG and RGbar matrix elements
"""
function set_elements_RS!(elementsRGbar, elementsRG, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")
    if R > 0
        S, base = set_base(S, splicetype)
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
                push!(elementsRGbar, Element(a, b, nu[R+1], sB))
            end
            if S == R
                if S > insertstep - 1 && abs(sC) == 1
                    push!(elementsRGbar, Element(a, b, eta[S-insertstep+1], sC))
                end
            end
            if splicetype == "offeject"
                s = (zbarr == wbarr) * ((zr == 0) - (zr == 1)) * (wr == 1)
                if abs(s) == 1
                    push!(elementsRGbar, Element(a, b, eta[R-insertstep+1], s))
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
                    push!(elementsRGbar, Element(a, b, nu[j+1], s))
                end
                # if S > 0 && j > insertstep - 1
                if S >= j && j > insertstep - 1
                    s = (zbark == wbark) * ((zj == 1) - (zj == 2)) * (wj == 2)
                    if abs(s) == 1
                        push!(elementsRGbar, Element(a, b, eta[j-insertstep+1], s))
                    end
                end
                if splicetype == "offeject" && j > insertstep - 1
                    s = (zbark == wbark) * ((zj == 0) - (zj == 1)) * (wj == 1)
                    if abs(s) == 1
                        push!(elementsRGbar, Element(a, b, eta[j-insertstep+1], s))
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
    set_elements_GRS(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)

Create combined gene-RNA-splicing transition elements.

# Description
This function creates Element structures for the combined GRS (Gene-RNA-Splicing) model,
which includes gene state transitions, RNA/reporter state transitions, and splicing transitions.
The function handles both cases where R > 0 (with RNA states) and R == 0 (gene-only model).

# Arguments
- `transitions`: Tuple of gene state transitions
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `S`: Number of splicing sites
- `insertstep`: Insert step for RNA processing
- `indices::Indices`: Index structure containing rate indices
- `splicetype::String`: Type of splicing model ("", "offeject", etc.)

# Returns
- `Tuple`: (elementsG, elementsRGbar, elementsRG, nR, nT) for R > 0
- `Tuple`: (elementsG, G) for R == 0
"""
function set_elements_GRS(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    if R > 0
        elementsG = Vector{Element}(undef, 0)
        elementsRGbar = Vector{Element}(undef, 0)
        elementsRG = Vector{Element}(undef, 0)
        base = S > 0 ? 3 : 2
        nR = base^R
        set_elements_G!(elementsG, transitions)
        set_elements_RS!(elementsRGbar, elementsRG, R, S, insertstep, indices.nu, indices.eta, splicetype)
        return elementsG, elementsRGbar, elementsRG, nR, T_dimension(G, R, S)
    else
        return set_elements_G(transitions, indices.gamma), G
    end
end

function elements_TG!(elements, elementsG, G, nT)
    for e in elementsG
        for j in 0:G:nT-1
            push!(elements, Element(e.a + j, e.b + j, e.index, e.pm))
        end
    end

end

function elements_TR!(elements, elementsR, G::Int)
    for e in elementsR
        for j in 1:G
            push!(elements, Element(j + G * (e.a - 1), j + G * (e.b - 1), e.index, e.pm))
        end
    end
end

function elements_RG!(elements, elementsR, G::Int)
    for e in elementsR
        push!(elements, Element(G * e.a, G * e.b, e.index, e.pm))
    end
end

function set_elements_TGRS(elementsG::Vector, elementsRGbar, elementsRG, G, nT)
    elements = []
    elements_TG!(elements, elementsG, G, nT)
    elements_TR!(elements, elementsRGbar, G)
    elements_RG!(elements, elementsRG, G)
    elements, nT
end

function set_elements_TGRS(transitions::Tuple, G, R, S, insertstep, indices::Indices, splicetype::String)
    if R > 0
        elementsG, elementsRGbar, elementsRG, _, nT = set_elements_GRS(transitions, G, R, S, insertstep, indices, splicetype)
        return set_elements_TGRS(elementsG, elementsRGbar, elementsRG, G, nT)
    else
        return set_elements_G(transitions, indices.gamma), G
    end
end

"""
    set_elements_T!(elementsT, G, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")

inplace update matrix elements in elementsT for GRS state transition matrix, do nothing if R == 0

"offeject" = pre-RNA is completely ejected when spliced
"""
function set_elements_T!(elementsT, G, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")
    if R > 0
        S, base = set_base(S, splicetype)
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
                # if S > 0 && abs(sC) == 1
                #     push!(elementsT, Element(a, b, eta[R-insertstep+1], sC))
                # end
                if S == R
                    if S > insertstep - 1 && abs(sC) == 1
                        push!(elementsT, Element(a, b, eta[S-insertstep+1], sC))
                    end
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
                    # if S > 0 && j > insertstep - 1
                    if S >= j && j > insertstep - 1
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
    set_elements_T(transitions, G, R, S, indices::Indices, splicetype::String)

return T matrix elements (Element structure) and size of matrix
"""
function set_elements_T(transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    if R > 0
        elementsT = Vector{Element}(undef, 0)
        base = S > 0 ? 3 : 2
        nT = G * base^R
        set_elements_G!(elementsT, transitions, G, indices.gamma, nT)
        set_elements_T!(elementsT, G, R, S, insertstep, indices.nu, indices.eta, splicetype)
        return elementsT, nT
    else
        return set_elements_G(transitions, indices.gamma), G
    end
end

# function set_elements_T(transitions, gamma::Vector)
#     elementsT = Vector{Element}(undef, 0)
#     set_elements_G!(elementsT, transitions, gamma)
#     elementsT
# end

function set_elements_TD!(elements_TD, elementsT, sojourn::Vector)
    for e in elementsT
        if e.b ∈ sojourn
            push!(elements_TD, e)
        end
    end
end
function set_elements_TD(elementsT, sojourn::Vector)
    elementsTD = Vector{Element}(undef, 0)
    set_elements_TD!(elementsTD, elementsT, sojourn::Vector)
    elementsTD
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
    set_elements_TA!(elementsTA, elementsT, onstates::Vector)

in place set onstate elements
"""
set_elements_TA!(elementsTA, elementsT, onstates) = set_elements_TX!(elementsTA, elementsT, onstates, ∈)

"""
    set_elements_TA(elementsT,onstates)

set onstate elements
"""
set_elements_TA(elementsT, onstates) = set_elements_TX(elementsT, onstates, set_elements_TA!)
"""
    set_elements_TI!(elementsTI, elementsT, onstates::Vector)

in place set off state elements
"""
set_elements_TI!(elementsTI, elementsT, onstates) = set_elements_TX!(elementsTI, elementsT, onstates, ∉)

"""
    set_elements_TI!(elementsTI, elementsT, onstates::Vector)

set off state elements
"""
set_elements_TI(elementsT, onstates) = set_elements_TX(elementsT, onstates, set_elements_TI!)

"""
    set_elements_TDvec(elementsT, elementsG, onstates, dttype)

Create dwell time analysis elements for different data types.

# Description
This function creates Element structures for dwell time analysis based on the specified
data type. It filters transition elements according to whether they involve on states,
off states, or gene-specific states.

# Arguments
- `elementsT`: Vector of transition elements
- `elementsG`: Vector of gene state elements
- `onstates`: Vector of on state indices
- `dttype`: Data type string ("ON", "OFF", "ONG", "OFFG")

# Returns
- `Vector{Vector{Element}}`: Vector of filtered elements for each data type
"""
function set_elements_TDvec(elementsT, elementsG, onstates::Vector{Vector{Int}}, dttype::Vector)
    c = Vector{Element}[]
    for i in eachindex(onstates)
        if dttype[i] == "ON"
            push!(c, set_elements_TA(elementsT, onstates[i]))
        elseif dttype[i] == "OFF"
            push!(c, set_elements_TI(elementsT, onstates[i]))
        elseif dttype[i] == "ONG"
            push!(c, set_elements_TA(elementsG, onstates[i]))
        elseif dttype[i] == "OFFG"
            push!(c, set_elements_TI(elementsG, onstates[i]))
        end
    end
    c
end

function set_elements_TDvec(elementsT, elementsG, sojourn::Vector{Vector{Int}}, dttype, nT, nG)
    c = Vector{Element}[]
    d = Int[]
    for i in eachindex(sojourn)
        if occursin("G", dttype[i])
            push!(c, set_elements_TD(elementsG, sojourn[i]))
            push!(d, nG)
        else
            push!(c, set_elements_TD(elementsT, sojourn[i]))
            push!(d, nT)
        end
    end
    c, d
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
    set_elements_BRG(G, R, ejectindex, base=2)

Create boundary elements for gene-RNA systems.

# Description
This function creates Element structures for boundary conditions in gene-RNA systems.
These elements represent the ejection of RNA molecules from the system.

# Arguments
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `ejectindex`: Index for the ejection rate
- `base`: Base for RNA state representation (default: 2)

# Returns
- `Vector{Element}`: Vector of boundary elements for RNA ejection
"""
function set_elements_BRG(G, R, ejectindex, base=2)
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
    set_elements_Gt!(elements, transitions, target_transition, gamma, j)

Create target gene transition elements in-place.

# Description
This function creates Element structures for a specific target gene transition and adds
them to the provided elements vector in-place. Only the specified target transition
is processed.

# Arguments
- `elements`: Vector to append elements to (modified in-place)
- `transitions`: Tuple of gene state transitions
- `target_transition`: Index of the specific transition to process
- `gamma::Vector`: Vector of rate indices
- `j`: Offset for state indices (default: 0)

# Returns
- `Nothing`: Modifies the elements vector in-place
"""
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
    set_elements_Gt(transitions, target_transition, gamma)

Create target gene transition elements.

# Description
This function creates Element structures for a specific target gene transition.
It returns a new vector containing only the elements for the specified transition.

# Arguments
- `transitions`: Tuple of gene state transitions
- `target_transition`: Index of the specific transition to process
- `gamma`: Vector of rate indices

# Returns
- `Vector{Element}`: Vector of elements for the target transition
"""
function set_elements_Gt(transitions, target_transition, gamma)
    elementsGt = Vector{Element}(undef, 0)
    set_elements_Gt!(elementsGt, transitions, target_transition, gamma)
    return elementsGt
end

"""
    set_elements_Gs(nS)

Create source gene state elements.

# Description
This function creates Element structures for source gene states. These elements
represent the states from which coupling can occur in coupled gene systems.

# Arguments
- `nS`: Source state index or vector of source state indices

# Returns
- `Vector{Element}`: Vector of source state elements
"""
function set_elements_Gs(nS::Int)
    [Element(nS, nS, 0, 1)]
end

function set_elements_Gs(nS::Vector{Int})
    elementsGs = Vector{Element}(undef, 0)
    for i in nS
        push!(elementsGs, Element(i, i, 0, 1))
    end
    return elementsGs
end

function set_elements_Source(nS::Int)
    if nS > 0
        return [Element(nS, nS, 0, 1)]
    else
        return Element[]
    end
end

function set_elements_Source(nS::Vector{Int})
    elementsSource = Vector{Element}(undef, 0)
    for i in nS
        i > 0 && push!(elementsSource, Element(i, i, 0, 1))
    end
    return elementsSource
end

function set_elements_Rt!(elementsRGbar, target, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")
    if R > 0
        S, base = set_base(S, splicetype)
        nR = base^R
        for b = 1:nR, a = 1:nR
            zdigits = digit_vector(a, base, R)
            wdigits = digit_vector(b, base, R)
            zr = zdigits[R]
            wr = wdigits[R]
            zbarr = zdigits[1:R-1]
            wbarr = wdigits[1:R-1]
            sB = 0
            for l in 1:base-1
                sB += (zbarr == wbarr) * ((zr == 0) - (zr == l)) * (wr == l)
            end
            if S > 0
                sC = (zbarr == wbarr) * ((zr == 1) - (zr == 2)) * (wr == 2)
            end
            if abs(sB) == 1 && nu[R+1] ∈ target
                push!(elementsRGbar, Element(a, b, nu[R+1], sB))
            end
            if S == R && eta[S-insertstep+1] ∈ target
                if S > insertstep - 1 && abs(sC) == 1
                    push!(elementsRGbar, Element(a, b, eta[S-insertstep+1], sC))
                end
            end
            if splicetype == "offeject" && eta[S-insertstep+1] ∈ target
                s = (zbarr == wbarr) * ((zr == 0) - (zr == 1)) * (wr == 1)
                if abs(s) == 1
                    push!(elementsRGbar, Element(a, b, eta[R-insertstep+1], s))
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
                if abs(s) == 1 && nu[j+1] ∈ target
                    push!(elementsRGbar, Element(a, b, nu[j+1], s))
                end
                if S >= j && j > insertstep - 1 && eta[j-insertstep+1] ∈ target
                    s = (zbark == wbark) * ((zj == 1) - (zj == 2)) * (wj == 2)
                    if abs(s) == 1
                        push!(elementsRGbar, Element(a, b, eta[j-insertstep+1], s))
                    end
                end
                if splicetype == "offeject" && j > insertstep - 1 && eta[j-insertstep+1] ∈ target
                    s = (zbark == wbark) * ((zj == 0) - (zj == 1)) * (wj == 1)
                    if abs(s) == 1
                        push!(elementsRGbar, Element(a, b, eta[j-insertstep+1], s))
                    end
                end
            end
        end
    end
end

function set_elements_Rt(target, R, S, insertstep, nu::Vector{Int}, eta::Vector{Int}, splicetype="")
    elements = Element[]
    set_elements_Rt!(elements, target, R, S, insertstep, nu, eta, splicetype)
    elements
end

function set_elements_Init!(elementsRG, target, R, S, nu::Vector{Int}, splicetype="")
    if R > 0
        S, base = set_base(S, splicetype)
        nR = base^R
        for b = 1:nR, a = 1:nR
            zdigits = digit_vector(a, base, R)
            wdigits = digit_vector(b, base, R)
            z1 = zdigits[1]
            w1 = wdigits[1]
            zbar1 = zdigits[2:R]
            wbar1 = wdigits[2:R]
            s = (zbar1 == wbar1) * ((z1 == base - 1) - (z1 == 0)) * (w1 == 0)
            if abs(s) == 1 && nu[1] ∈ target
                push!(elementsRG, Element(a, b, nu[1], s))
            end
        end
    end
end

function set_elements_Init(target, R, S, nu::Vector{Int}, splicetype="")
    elements = Element[]
    set_elements_Init!(elements, target, R, S, nu, splicetype)
    elements
end

function set_elements_Source(sources, G::Int, R, S, splicetype="")
    elementsSource = Vector{Element}(undef, 0)
    _, base = set_base(S, splicetype)
    Gstates, Rstates = classify_states(sources, G, R, S, splicetype)
    for nS in Gstates
        elementsG = set_elements_Source(nS)
        elements_TG!(elementsSource, elementsG, G, G * base^R)
    end
    for nS in Rstates
        elementsR = set_elements_Source(nS)
        elements_TR!(elementsSource, elementsR, G)
    end
    elementsSource
end

"""
    set_elements_Target(target, Gtransitions, G, R, S, insertstep, indices, nT, splicetype)

Create target transition elements for coupled systems.

# Description
This function creates Element structures for target transitions in coupled gene systems.
It handles gene transitions, RNA initiation, and RNA processing transitions based on
the target specification.

# Arguments
- `target`: Target transition specification
- `Gtransitions`: Gene state transitions
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `S`: Number of splicing sites
- `insertstep`: Insert step for RNA processing
- `indices`: Index structure containing rate indices
- `nT`: Total number of states
- `splicetype`: Type of splicing model

# Returns
- `Vector{Element}`: Vector of target transition elements
"""
function set_elements_Target(target, Gtransitions, G, R, S, insertstep, indices, nT, splicetype)
    elementsTarget = Vector{Element}(undef, 0)
    Gts, Init, Rts = classify_transitions(target, indices)
    for t in Gts
        elements = set_elements_Gt(Gtransitions, t, indices.gamma)
        elements_TG!(elementsTarget, elements, G, nT)
    end
    for t in Init
        elements = set_elements_Init(t, R, S, indices.nu, splicetype)
        elements_RG!(elementsTarget, elements, G)
    end
    for t in Rts
        elements = set_elements_Rt(t, R, S, insertstep, indices.nu, indices.eta, splicetype)
        elements_TR!(elementsTarget, elements, G)
    end
    return elementsTarget
end

"""
    set_elements_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)

Create coupled gene-RNA transition elements.

# Description
This function creates Element structures for coupled gene-RNA systems, including
gene state transitions, target transitions, source states, and RNA processing transitions.

# Arguments
- `source_state`: Source state specification for coupling
- `target_transition`: Target transition specification
- `transitions`: Gene state transitions
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `S`: Number of splicing sites
- `insertstep`: Insert step for RNA processing
- `indices`: Index structure containing rate indices
- `splicetype`: Type of splicing model

# Returns
- `Tuple`: (elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT)
"""
function set_elements_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    elementsG, elementsRGbar, elementsRG, nR, nT = set_elements_GRS(transitions, G, R, S, insertstep, indices, splicetype)
    elementsGt = set_elements_Gt(transitions, target_transition, indices.gamma)
    elementsGs = set_elements_Gs(source_state)
    return elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT
end

function set_elements_TCoupledUnit(source_state, target_transition, transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
    elementsSource = set_elements_Source(source_state, G, R, S)
    elementsTarget = set_elements_Target(target_transition, transitions, G, R, S, insertstep, indices, nT, splicetype)
    return elementsT, elementsSource, elementsTarget, nT
end




#### Experimental

"""
    add_transition_element!(elements, state_a, state_b, rate, sign)

Helper function to add transition elements when condition is met.

# Description
This helper function adds a transition element to the elements vector if the sign
has magnitude 1. This is used to filter out zero-valued transitions.

# Arguments
- `elements`: Vector to append elements to (modified in-place)
- `state_a`: Source state index
- `state_b`: Target state index
- `rate`: Rate index
- `sign`: Sign of the transition (should be ±1 for valid transitions)

# Returns
- `Nothing`: Modifies the elements vector in-place
"""
function add_transition_element!(elements, state_a, state_b, rate, sign)
    if abs(sign) == 1
        push!(elements, Element(state_a, state_b, rate, sign))
    end
end

"""
    calculate_state_indices(zdigits, wdigits, R)

Calculate state indices for RNA transitions.

# Description
This helper function extracts various state indices from digit vectors for RNA
transitions. It provides convenient access to different parts of the state
representation needed for transition calculations.

# Arguments
- `zdigits`: Digit vector for source state
- `wdigits`: Digit vector for target state
- `R`: Number of RNA/reporter steps

# Returns
- `Tuple`: (zbar1, wbar1, zbarr, wbarr, z1, w1, zr, wr) - various state indices
"""
function calculate_state_indices(zdigits, wdigits, R)
    zbar1 = zdigits[2:R]
    wbar1 = wdigits[2:R]
    zbarr = zdigits[1:R-1]
    wbarr = wdigits[1:R-1]
    z1 = zdigits[1]
    w1 = wdigits[1]
    zr = zdigits[R]
    wr = wdigits[R]
    return zbar1, wbar1, zbarr, wbarr, z1, w1, zr, wr
end

"""
    process_rna_transitions!(elements, zdigits, wdigits, state_a, state_b, R, base, nu)

Process RNA transitions for both T and RS matrices.

# Description
This helper function processes RNA transitions by calculating transition elements
for both end transitions and internal transitions. It adds the appropriate
elements to the elements vector for RNA processing steps.

# Arguments
- `elements`: Vector to append elements to (modified in-place)
- `zdigits`: Digit vector for source state
- `wdigits`: Digit vector for target state
- `state_a`: Source state index
- `state_b`: Target state index
- `R`: Number of RNA/reporter steps
- `base`: Base for RNA state representation
- `nu`: Vector of RNA transition rate indices

# Returns
- `Nothing`: Modifies the elements vector in-place
"""
function process_rna_transitions!(elements, zdigits, wdigits, state_a, state_b, R, base, nu)
    # End transitions
    zbarr = zdigits[1:R-1]
    wbarr = wdigits[1:R-1]
    zr = zdigits[R]
    wr = wdigits[R]

    sB = 0
    for l in 1:base-1
        sB += (zbarr == wbarr) * ((zr == 0) - (zr == l)) * (wr == l)
    end
    add_transition_element!(elements, state_a, state_b, nu[R+1], sB)

    # Internal transitions
    for j = 1:R-1
        zbarj = zdigits[[1:j-1; j+2:R]]
        wbarj = wdigits[[1:j-1; j+2:R]]
        zj = zdigits[j]
        zj1 = zdigits[j+1]
        wj = wdigits[j]
        wj1 = wdigits[j+1]

        s = 0
        for l in 1:base-1
            s += (zbarj == wbarj) * ((zj == 0) * (zj1 == l) - (zj == l) * (zj1 == 0)) * (wj == l) * (wj1 == 0)
        end
        add_transition_element!(elements, state_a, state_b, nu[j+1], s)
    end
end

"""
    process_splicing!(elements, zdigits, wdigits, state_a, state_b, S, R, insertstep, eta, splicetype)

Process splicing transitions with fixed S==R bug.

# Description
This helper function processes splicing transitions by calculating transition elements
for both end splicing and internal splicing. It handles different splicing types
including "offeject" and standard splicing models.

# Arguments
- `elements`: Vector to append elements to (modified in-place)
- `zdigits`: Digit vector for source state
- `wdigits`: Digit vector for target state
- `state_a`: Source state index
- `state_b`: Target state index
- `S`: Number of splicing sites
- `R`: Number of RNA/reporter steps
- `insertstep`: Insert step for RNA processing
- `eta`: Vector of splicing transition rate indices
- `splicetype`: Type of splicing model

# Returns
- `Nothing`: Modifies the elements vector in-place
"""
function process_splicing!(elements, zdigits, wdigits, state_a, state_b, S, R, insertstep, eta, splicetype)
    zbarr = zdigits[1:R-1]
    wbarr = wdigits[1:R-1]
    zr = zdigits[R]
    wr = wdigits[R]

    if splicetype == "offeject"
        s = (zbarr == wbarr) * ((zr == 0) - (zr == 1)) * (wr == 1)
        add_transition_element!(elements, state_a, state_b, eta[R-insertstep+1], s)
    elseif S == R && S > insertstep - 1
        sC = (zbarr == wbarr) * ((zr == 1) - (zr == 2)) * (wr == 2)
        add_transition_element!(elements, state_a, state_b, eta[S-insertstep+1], sC)
    end

    process_internal_splicing!(elements, zdigits, wdigits, state_a, state_b, S, R, insertstep, eta, splicetype)
end

"""
    process_internal_splicing!(elements, zdigits, wdigits, state_a, state_b, S, R, insertstep, eta, splicetype)

Process internal splicing transitions.

# Description
This helper function processes internal splicing transitions by calculating
transition elements for splicing at internal positions in the RNA chain.
It handles both standard splicing and "offeject" splicing models.

# Arguments
- `elements`: Vector to append elements to (modified in-place)
- `zdigits`: Digit vector for source state
- `wdigits`: Digit vector for target state
- `state_a`: Source state index
- `state_b`: Target state index
- `S`: Number of splicing sites
- `R`: Number of RNA/reporter steps
- `insertstep`: Insert step for RNA processing
- `eta`: Vector of splicing transition rate indices
- `splicetype`: Type of splicing model

# Returns
- `Nothing`: Modifies the elements vector in-place
"""
function process_internal_splicing!(elements, zdigits, wdigits, state_a, state_b, S, R, insertstep, eta, splicetype)
    for j = 1:R-1
        zbark = zdigits[[1:j-1; j+1:R]]
        wbark = wdigits[[1:j-1; j+1:R]]
        zj = zdigits[j]
        wj = wdigits[j]

        if S >= j && j > insertstep - 1
            s = (zbark == wbark) * ((zj == 1) - (zj == 2)) * (wj == 2)
            add_transition_element!(elements, state_a, state_b, eta[j-insertstep+1], s)
        end

        if splicetype == "offeject" && j > insertstep - 1
            s = (zbark == wbark) * ((zj == 0) - (zj == 1)) * (wj == 1)
            add_transition_element!(elements, state_a, state_b, eta[j-insertstep+1], s)
        end
    end
end

"""
    process_states!(elements, G, R, S, insertstep, nu, eta, splicetype)

Shared state transition processing for T and RS matrices.

# Description
This helper function processes all state transitions by iterating through all
possible RNA states and gene states. It calls the RNA transition and splicing
processing functions to build the complete transition matrix elements.

# Arguments
- `elements`: Vector to append elements to (modified in-place)
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `S`: Number of splicing sites
- `insertstep`: Insert step for RNA processing
- `nu`: Vector of RNA transition rate indices
- `eta`: Vector of splicing transition rate indices
- `splicetype`: Type of splicing model

# Returns
- `Nothing`: Modifies the elements vector in-place
"""
function process_states!(elements, G, R, S, insertstep, nu, eta, splicetype)
    base = (splicetype == "offeject") ? 2 : 3
    for w = 1:base^R, z = 1:base^R
        zdigits = digit_vector(z, base, R)
        wdigits = digit_vector(w, base, R)

        for i = 1:G
            state_a = state_index(G, i, z)
            state_b = state_index(G, i, w)

            process_rna_transitions!(elements, zdigits, wdigits, state_a, state_b, R, base, nu)
            process_splicing!(elements, zdigits, wdigits, state_a, state_b, S, R, insertstep, eta, splicetype)
        end
    end
end

"""
    set_elements_T!(elementsT, G, R, S, insertstep, nu, eta, splicetype="")

Entry point for T matrix element creation.

# Description
This function serves as the entry point for creating T matrix elements. It calls
the shared state processing function to build all transition elements for the T matrix.

# Arguments
- `elementsT`: Vector to append elements to (modified in-place)
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `S`: Number of splicing sites
- `insertstep`: Insert step for RNA processing
- `nu`: Vector of RNA transition rate indices
- `eta`: Vector of splicing transition rate indices
- `splicetype`: Type of splicing model (default: "")

# Returns
- `Nothing`: Modifies the elementsT vector in-place
"""
function set_elements_T!(elementsT, G, R, S, insertstep, nu, eta, splicetype="")
    R > 0 && process_states!(elementsT, G, R, S, insertstep, nu, eta, splicetype)
end

"""
    set_elements_RS!(elementsRS, G, R, S, insertstep, nu, eta, splicetype="")

Entry point for RS matrix element creation.

# Description
This function serves as the entry point for creating RS matrix elements. It calls
the shared state processing function to build all transition elements for the RS matrix.

# Arguments
- `elementsRS`: Vector to append elements to (modified in-place)
- `G`: Number of gene states
- `R`: Number of RNA/reporter steps
- `S`: Number of splicing sites
- `insertstep`: Insert step for RNA processing
- `nu`: Vector of RNA transition rate indices
- `eta`: Vector of splicing transition rate indices
- `splicetype`: Type of splicing model (default: "")

# Returns
- `Nothing`: Modifies the elementsRS vector in-place
"""
function set_elements_RS!(elementsRS, G, R, S, insertstep, nu, eta, splicetype="")
    R > 0 && process_states!(elementsRS, G, R, S, insertstep, nu, eta, splicetype)
end

# function parse_sourcestring(s::String)
#     m = match(r"^([A-Za-z]+)(\d*)$", s)
#     if m === nothing
#         error("Input must consist of one or more letters optionally followed by digits: $s")
#     else
#         letters = m.captures[1]
#         digits = m.captures[2]
#         if isempty(digits)
#             return letters, nothing
#         else
#             return letters, parse(Int, digits)
#         end
#     end
# end

# function group_sources(strings::Vector{String}, G, R)
#     groups = Dict{String,Vector{Int}}()
#     for s in strings
#         letter, number = parse_sourcestring(s)
#         if !haskey(groups, letter)
#             groups[letter] = Vector{Int}()  # Initialize an empty vector for this letter.
#         end
#         if isnothing(number)
#             number = letter == "G" ? collect(1:G) : collect(1:R)
#             push!(groups[letter], number...)
#         else
#             push!(groups[letter], number)
#         end
#     end
#     return groups
# end