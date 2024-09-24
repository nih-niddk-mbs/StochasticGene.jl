# This file is part of StochasticGene.jl
#
# transition_rate_functions.jl
#


"""
    state_index(G::Int,i,z)

return state index for state (i,z)
"""
state_index(G::Int, g, z) = g + G * (z - 1)


"""
    inverse_state(i::Int, G::Int, R, S, insertstep::Int, f=sum)


"""
function inverse_state(i::Int, G::Int, R, S, insertstep::Int, f=sum)
    base = S > 0 ? 3 : 2
    g = mod(i - 1, G) + 1
    z = div(i - g, G) + 1
    zdigits = digit_vector(z, base, R)
    r = num_reporters_per_index(z, R, insertstep, base, f)
    return g, z, zdigits, r
end

function unit_state(i::Int, G::Tuple, R, S, unit_model)
    nT = T_dimension(G, R, S, unit_model)
    rem = i
    unit = Int[]
    for j in unit_model
        push!(unit, mod(rem - 1, nT[j]) + 1)
        rem = div(rem - unit[j], nT[j]) + 1
    end
    unit
end

function inverse_state(i::Int, G::Tuple, R, S, insertstep, unit_model, f=sum)
    units = unit_state(i, G, R, S, unit_model)
    inverse_state(units, G, R, S, insertstep, f)
end

function inverse_state(units::Vector, G::Tuple, R, S, insertstep, f=sum)
    states = Tuple[]
    for i in eachindex(units)
        push!(states, inverse_state(units[i], G[i], R[i], S[i], insertstep[i], f))
    end
    states
end

function inverse_state(i::Vector{Int}, G::Int, R, S, insertstep::Int)
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

function inverse_state(i::Vector{Vector{Int}}, G::Int, R, S, insertstep)
    z = Vector{Int}[]
    g = Vector{Int}[]
    zdigits = Vector{Vector}[]
    r = Vector{Int}[]
    for i in i
        gi, zi, zdi, ri = inverse_state(i, G::Int, R, S, insertstep)
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

"""
    off_states(G::Int, R, S, insertstep)

return barrier (off) states, complement of sojourn (on) states
"""
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
    off_states(reporters)


"""
function off_states(reporters)
    offstates = Int[]
    for i in eachindex(reporters)
        (reporters[i]) < 1 && push!(offstates, i)
    end
    offstates
end

"""
    T_dimension(G, R, S)

Compute transition matrix dimension of GRS model
"""
function T_dimension(G, R, S)
    base = S > 0 ? 3 : 2
    G * base^R
end

function T_dimension(G::Tuple, R::Tuple, S::Tuple)
    nT = Int[]
    for m in eachindex(G)
        push!(nT, T_dimension(G[m], R[m], S[m]))
    end
    nT
end

function T_dimension(G::Tuple, R::Tuple, S::Tuple, unit_model)
    nT = Int[]
    for m in unit_model
        push!(nT, T_dimension(G[m], R[m], S[m]))
    end
    nT
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
function num_reporters_per_state_reduced(G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, sources, f=sum)
    reporters = Vector[]
    for i in eachindex(R)
        repi = num_reporters_per_state(G[i], R[i], S[i], insertstep[i], f)
        for j in 1:i-1
            (j ∈ sources[i]) && (repi = repeat(repi, outer=(G[j],)))
        end
        for j in i+1:length(R)
            (j ∈ sources[i]) && (repi = repeat(repi, inner=(G[j],)))
        end
        push!(reporters, repi)
    end
    reporters
end

function num_reporters_per_state(G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, unit_model, f=sum)
    nT = T_dimension.(G, R, S)
    reporters = Vector[]
    for m in unit_model
        rep = num_reporters_per_state(G[m], R[m], S[m], insertstep[m], f)
        for j in 1:m-1
            rep = repeat(rep, outer=(nT[unit_model[j]],))
        end
        for j in m+1:length(unit_model)
            rep = repeat(rep, inner=(nT[unit_model[j]],))
        end
        push!(reporters, rep)
    end
    reporters
end

"""
    num_reporters_per_index(z, R, insertstep, base, f=sum)

TBW
"""
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

"""
    set_indices(ntransitions, R, S, insertstep, offset)

"""
function set_indices(ntransitions, R, S, insertstep, offset)
    if insertstep > R > 0
        throw("insertstep>R")
    end
    if S > 0
        Indices(offset .+ collect(1:ntransitions), offset .+ collect(ntransitions+1:ntransitions+R+1), offset .+ collect(ntransitions+R+2:ntransitions+R+R-insertstep+2), offset + ntransitions + R + R - insertstep + 3)
    elseif R > 0
        Indices(offset .+ collect(1:ntransitions), offset .+ collect(ntransitions+1:ntransitions+R+1), Int[], offset .+ ntransitions + R + 2)
    else
        Indices(offset .+ collect(1:ntransitions), offset .+ [ntransitions + 1], Int[], offset + ntransitions + 2)
    end
end
