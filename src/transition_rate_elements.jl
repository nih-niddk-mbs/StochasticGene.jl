# This file is part of StochasticGene.jl
#
# transition_rate_set_elements.jl
#

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

TBW
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
                push!(elementsRGbar, Element(a, b, nu[R+1], sB))
            end
            if S > 0
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
                if S > 0 && S > j > insertstep - 1
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

TBW
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
                # if S > 0 && abs(sC) == 1
                #     push!(elementsT, Element(a, b, eta[R-insertstep+1], sC))
                # end
                if S > 0
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
                    if S > 0 && S > j > insertstep - 1
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
    set_elements_TDvec(elementsT, onstates, dttype)

TBW
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
    set_elements_Gt!(elements, transitions, target_transition=length(transitions), gamma::Vector=collect(1:length(transitions)), j=0)


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


"""
function set_elements_Gt(transitions, target_transition, gamma)
    elementsGt = Vector{Element}(undef, 0)
    set_elements_Gt!(elementsGt, transitions, target_transition, gamma)
    return elementsGt
end

"""
    set_elements_Gs(nS)


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

"""
    set_elements_Coupled(source_state, target_transition, transitions, G, R, S, insertstep, indices::Indices, splicetype::String)


"""
function set_elements_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    elementsG, elementsRGbar, elementsRG, nR, nT = set_elements_GRS(transitions, G, R, S, insertstep, indices, splicetype)
    elementsGt = set_elements_Gt(transitions, target_transition, indices.gamma)
    elementsGs = set_elements_Gs(source_state)
    return elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT
end

function set_elements_TCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    elementsT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
    elementsGt = set_elements_Gt(transitions, target_transition, indices.gamma)
    elementsGs = set_elements_Gs(source_state)
end

function set_elements_V!(elements, transitions, target_transition=length(transitions), gamma::Vector=collect(1:length(transitions)), j=0)
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
    set_elements_V(transitions, target_transition, gamma)


"""
function set_elements_V(transitions, target_transition, gamma)
    elementsV = Vector{Element}(undef, 0)
    set_elements_V!(elementsGt, transitions, target_transition, gamma)
    return elementsV
end

"""
    set_elements_U(nS::Int)


"""
function set_elements_U(nS::Int)
    [Element(nS, nS, 0, 1)]
end

function set_elements_U(nS::Vector{Int})
    elementsU = Vector{Element}(undef, 0)
    for i in nS
        push!(elementsU, Element(i, i, 0, 1))
    end
    return elementsU
end


function set_elements_TUVCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices::Indices, splicetype::String)
    elementsT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    elementsV = set_elements_V(transitions, target_transition, indices.gamma)
    elementsU = set_elements_U(source_state)
    return elementsT, elementsV, elementsU
end