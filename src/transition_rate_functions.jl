# This file is part of StochasticGene.jl
#
# transition_rate_functions.jl
#


"""
    state_index(G::Int, g, z)

Return the state index for state `(g, z)`.

# Arguments
- `G::Int`: Total number of genes.
- `g`: Gene index.
- `z`: State index.

# Returns
- `Int`: The state index.
"""
state_index(G::Int, g, z) = g + G * (z - 1)

"""
    inverse_state(i::Int, G::Int, R, S, insertstep::Int, f=sum)
    inverse_state(i::Int, G::Tuple, R, S, insertstep, unit_model, f=sum)
    inverse_state(units::Vector, G::Tuple, R, S, insertstep, f=sum)
    inverse_state(i::Vector{Int}, G::Int, R, S, insertstep::Int)

Return the inverse state for a given state index or vector of units.

# Description
This function returns the inverse state (G and RS indices) for a given state index or vector of units. It can handle different types of inputs, including integers, tuples, and vectors.

# Arguments
- `i::Int`: State index.
- `G::Int`: Total number of genes.
- `G::Tuple`: Tuple of total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `unit_model`: Unit model.
- `units::Vector`: Vector of units.
- `f`: Function to apply (default is `sum`).

# Methods
- `inverse_state(i::Int, G::Int, R, S, insertstep::Int, f=sum)`: Handles `i` as an integer.
- `inverse_state(i::Int, G::Tuple, R, S, insertstep, unit_model, f=sum)`: Handles `i` as an integer with a tuple of genes.
- `inverse_state(units::Vector, G::Tuple, R, S, insertstep, f=sum)`: Handles `units` as a vector.
- `inverse_state(i::Vector{Int}, G::Int, R, S, insertstep::Int)`: Handles `i` as a vector of integers.

# Returns
- `Tuple`: Inverse state.
- `Vector{Tuple}`: Vector of inverse states.
"""
function inverse_state(i::Int, G::Int, R, S, insertstep::Int, f=sum)
    base = S > 0 ? 3 : 2
    g = mod(i - 1, G) + 1
    z = div(i - g, G) + 1
    zdigits = digit_vector(z, base, R)
    r = num_reporters_per_index(z, R, insertstep, base, f)
    return g, z, zdigits, r
end

function inverse_state(i::Int, G::Tuple, R, S, insertstep, unit_model, f=sum)
    units = unit_state(i, G, R, S, unit_model)
    inverse_state(units, G, R, S, insertstep, f)
end

function inverse_state(units, G::Tuple, R, S, insertstep, f=sum)
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
    unit_state(i::Int, G::Tuple, R, S, unit_model)

Return the unit state for a given state index.

# Arguments
- `i::Int`: State index.
- `G::Tuple`: Tuple of total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `unit_model`: Unit model.

# Returns
- `Tuple`: Unit state.
"""
function unit_state(i::Int, G::Tuple, R, S, unit_model)
    nT = T_dimension(G, R, S, unit_model)
    rem = i
    unit = Vector{Int}(undef, length(unit_model))
    for j in reverse(unit_model)
        unit[j] = mod(rem - 1, nT[j]) + 1
        rem = div(rem - unit[j], nT[j]) + 1
    end
    Tuple(unit)
end
"""
    on_states(onstates, G, R, S, insertstep)
    on_states(G::Int, R, S, insertstep)
    on_states(onstates::Vector, G, R, S)

Return vector of on state indices for GR and GRS models.

# Description
This function returns a vector of on state indices for GR and GRS models. It can handle both empty and non-empty `onstates`.

# Arguments
- `onstates`: Initial onstates. Can be a vector of integers or empty.
- `G::Int`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.

# Methods
- `on_states(onstates, G, R, S, insertstep)`: Handles both empty and non-empty `onstates`.
- `on_states(G::Int, R, S, insertstep)`: Generates `onstates` from scratch.
- `on_states(onstates::Vector, G, R, S)`: Generates `onstates` based on the provided vector.

# Returns
- `Vector{Int}`: Vector of new onstates.
"""
function on_states(onstates, G, R, S, insertstep)
    if isempty(onstates)
        return on_states(G, R, S, insertstep)
    else
        return on_states(onstates, G, R, S)
    end
end
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

function on_states(G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, unit_model)
    nT = T_dimension.(G, R, S)
    onstates = Vector[]
    for m in unit_model
        o = on_states(G[m], R[m], S[m], insertstep[m])
        println(o)
        nbefore = prod(nT[1:m-1])
        nafter = prod(nT[m+1:end])
        enlarged_indices = [(j - 1) * nbefore + 1:j * nbefore for j in o]
        enlarged_indices = vcat(enlarged_indices...)  # Flatten again into a single vector
        enlarged_indices = [i + (j-1)*nafter for i in enlarged_indices, j in 1:nafter]
        enlarged_indices = vcat(enlarged_indices...)  # Flatten into a single vector
        push!(onstates, enlarged_indices)
    end
    onstates
end

"""
    on_states!(onstates::Vector, G::Int, R::Int, insertstep, base)

In-place function to generate vector of on state indices for GR and GRS models.

# Description
This function modifies the `onstates` vector in place to generate a vector of on state indices for GR and GRS models.

# Arguments
- `onstates::Vector`: Initial onstates.
- `G::Int`: Total number of genes.
- `R::Int`: Number of reporters.
- `insertstep`: Insert step.
- `base`: Base value.

# Returns
- `Nothing`: This function modifies the `onstates` vector in place.
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
    off_states(nT, onstates)
    off_states(reporters)

Return barrier (off) states, complement of sojourn (on) states.

# Description
This function returns the barrier (off) states, which are the complement of the sojourn (on) states.

# Arguments
- `G::Int`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `nT`: Total number of states.
- `onstates`: Vector of on states.
- `reporters`: Vector of reporter states.

# Methods
- `off_states(G::Int, R, S, insertstep)`: Computes off states based on the given parameters.
- `off_states(nT, onstates)`: Computes off states as the complement of the given on states.
- `off_states(reporters)`: Computes off states based on the reporter states.

# Returns
- `Vector{Int}`: Vector of off states.
"""
function off_states(G::Int, R, S, insertstep)
    base = S > 0 ? 3 : 2
    nT = G * base^R
    off_states(nT, on_states(G, R, S, insertstep))
end

off_states(nT, onstates) = setdiff(collect(1:nT), onstates)

function off_states(reporters)
    offstates = Int[]
    for i in eachindex(reporters)
        (reporters[i]) < 1 && push!(offstates, i)
    end
    offstates
end

"""
    T_dimension(G, R, S)
    T_dimension(G::Tuple, R::Tuple, S::Tuple)
    T_dimension(G::Int, R::Int, S::Int)

Compute the transition matrix dimension of the GRS model.

# Description
This function computes the transition matrix dimension of the GRS model based on the provided parameters.

# Arguments
- `G`: Total number of genes or a tuple of total number of genes.
- `R`: Number of reporters or a tuple of number of reporters.
- `S`: Number of states or a tuple of number of states.

# Methods
- `T_dimension(G, R, S)`: Computes the dimension for given integers `G`, `R`, and `S`.
- `T_dimension(G::Tuple, R::Tuple, S::Tuple)`: Computes the dimension for given tuples `G`, `R`, and `S`.
- `T_dimension(G::Int, R::Int, S::Int)`: Computes the dimension for given integers `G`, `R`, and `S`.

# Returns
- `Int`: The dimension of the transition matrix.
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
    num_reporters_per_index(z, R, insertstep, base, f=sum)

Compute the number of reporters for a given state index.

# Description
This function computes the number of reporters for a given state index `z`. It uses the `digits` function to determine the number of reporters based on the specified `base` and `insertstep`.

# Arguments
- `z`: State index.
- `R`: Number of reporters.
- `insertstep`: Insert step.
- `base`: Base value.
- `f`: Function to apply (default is `sum`).

# Returns
- `Int`: The number of reporters for the given state index.
"""
num_reporters_per_index(z, R, insertstep, base, f=sum) = f(digits(z - 1, base=base, pad=R)[insertstep:end] .== base - 1)

"""
    num_reporters_per_state(G::Int, R::Int, S::Int=0, insertstep=1, f=sum)
    num_reporters_per_state(G::Int, onstates::Vector)
    num_reporters_per_state(G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, unit_model, f=sum)

Return the number of reporters for each state index.

# Description
This function returns the number of reporters for each state index. It can handle different types of inputs, including integers, tuples, and vectors.

# Arguments
- `G::Int`: Total number of genes.
- `R::Int`: Number of reporters.
- `S::Int`: Number of states (default is 0).
- `insertstep`: Insert step (default is 1).
- `f`: Function to apply (default is `sum`).
- `onstates::Vector`: Vector of on states.
- `G::Tuple`: Tuple of total number of genes.
- `R::Tuple`: Tuple of number of reporters.
- `S::Tuple`: Tuple of number of states.
- `insertstep::Tuple`: Tuple of insert steps.
- `unit_model`: Unit model.

# Methods
- `num_reporters_per_state(G::Int, R::Int, S::Int=0, insertstep=1, f=sum)`: Computes the number of reporters based on the given parameters.
- `num_reporters_per_state(G::Int, onstates::Vector)`: Computes the number of reporters based on the provided vector of on states.
- `num_reporters_per_state(G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, unit_model, f=sum)`: Computes the number of reporters based on the provided tuples and unit model.

# Returns
- `Vector{Int}`: Vector of the number of reporters for each state index.
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
# function num_reporters_per_state_reduced(G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, sources, f=sum)
#     reporters = Vector[]
#     for i in eachindex(R)
#         repi = num_reporters_per_state(G[i], R[i], S[i], insertstep[i], f)
#         for j in 1:i-1
#             (j ∈ sources[i]) && (repi = repeat(repi, outer=(G[j],)))
#         end
#         for j in i+1:length(R)
#             (j ∈ sources[i]) && (repi = repeat(repi, inner=(G[j],)))
#         end
#         push!(reporters, repi)
#     end
#     reporters
# end
"""
    set_indices(ntransitions, R, S, insertstep)
    set_indices(ntransitions, R)
    set_indices(ntransitions)
    set_indices(ntransitions, R, S, insertstep, offset)

Return an `Indices` structure based on the provided parameters.

# Description
This function returns an `Indices` structure based on the provided parameters. It can handle different combinations of `ntransitions`, `R`, `S`, `insertstep`, and `offset`.

# Arguments
- `ntransitions`: Number of transitions.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `offset`: Offset value (optional).

# Methods
- `set_indices(ntransitions, R, S, insertstep)`: Computes the `Indices` structure based on the given parameters.
- `set_indices(ntransitions, R)`: Computes the `Indices` structure based on `ntransitions` and `R`.
- `set_indices(ntransitions)`: Computes the `Indices` structure based on `ntransitions`.
- `set_indices(ntransitions, R, S, insertstep, offset)`: Computes the `Indices` structure based on the given parameters with an offset.

# Returns
- `Indices`: The computed `Indices` structure.
"""
function set_indices(ntransitions, R, S, insertstep)
    if insertstep > R > 0
        throw("insertstep>R")
    end
    if S > 0
        # Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), collect(ntransitions+R+2:ntransitions+R+R-insertstep+2), ntransitions + R + R - insertstep + 3)
        Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), collect(ntransitions+R+2:ntransitions+R+S-insertstep+2), ntransitions + R + S - insertstep + 3)
    elseif R > 0
        set_indices(ntransitions, R)
    else
        set_indices(ntransitions)
    end
end
set_indices(ntransitions, R) = Indices(collect(1:ntransitions), collect(ntransitions+1:ntransitions+R+1), Int[], ntransitions + R + 2)
set_indices(ntransitions) = Indices(collect(1:ntransitions), [ntransitions + 1], Int[], ntransitions + 2)

function set_indices(ntransitions, R, S, insertstep, offset)
    if insertstep > R > 0
        throw("insertstep>R")
    end
    if S > 0
        Indices(offset .+ collect(1:ntransitions), offset .+ collect(ntransitions+1:ntransitions+R+1), offset .+ collect(ntransitions+R+2:ntransitions+R+S-insertstep+2), offset + ntransitions + R + S - insertstep + 3)
    elseif R > 0
        Indices(offset .+ collect(1:ntransitions), offset .+ collect(ntransitions+1:ntransitions+R+1), Int[], offset + ntransitions + R + 2)
    else
        Indices(offset .+ collect(1:ntransitions), offset .+ [ntransitions + 1], Int[], offset + ntransitions + 2)
    end
end

"""
    kron_right(T, M::Vector, sources, unit_model, first, last)
    kron_right(T, M::Vector, unit_model, first, last)

Perform Kronecker right product operations.

# Description
This function performs a right Kronecker product operation on the matrix `T` using the identity matrices `I` and the `unit_model` indices. It can either include a `sources` vector to selectively apply the Kronecker product or apply it to all indices in the specified range.

# Arguments
- `T`: Initial matrix to be transformed.
- `M`: Vector of identity matrices.
- `sources`: Vector of source indices (optional).
- `unit_model`: Vector of unit model indices.
- `first`: Starting index for the operation.
- `last`: Ending index for the operation.

# Methods
- `kron_right(T, M, sources, unit_model, first, last)`: Applies the Kronecker product selectively based on the `sources` vector.
- `kron_right(T, M, unit_model, first, last)`: Applies the Kronecker product to all indices in the specified range.

# Returns
- `T`: The transformed matrix after applying the Kronecker product.
"""
function kron_right(T, M::Vector, sources, unit_model, first, last)
    for j in first:last
        if j ∈ sources
            T = kron(T, M[unit_model[j]])
        end
    end
    T
end

function kron_right(T, M::Vector, unit_model, first, last)
    for j in first:last
        T = kron(T, M[unit_model[j]])
    end
    T
end

"""
    kron_left(T, M::Vector, sources, unit_model, first, last)
    kron_left(T, M::Vectort unit_model, first, last)

Perform Kronecker left product operations.

# Description
This function performs a left Kronecker product operation on the matrix `T` using the identity matrices `M` and the `unit_model` indices. It can either include a `sources` vector to selectively apply the Kronecker product or apply it to all indices in the specified range.

# Arguments
- `T`: Initial matrix to be transformed.
- `M`: Vector of identity matrices.
- `sources`: Vector of source indices (optional).
- `unit_model`: Vector of unit model indices.
- `first`: Starting index for the operation.
- `last`: Ending index for the operation.

# Methods
- `kron_left(T, M, sources, unit_model, first, last)`: Applies the Kronecker product selectively based on the `sources` vector.
- `kron_left(T, M, unit_model, first, last)`: Applies the Kronecker product to all indices in the specified range.

# Returns
- `T`: The transformed matrix after applying the Kronecker product.
"""
function kron_left(T, M::Vector, sources, unit_model, first, last)
    for j in first:-1:last
        if j ∈ sources
            T = kron(M[unit_model[j]], T)
        end
    end
    T
end

function kron_left(T, M::Vector, unit_model, first, last)
    for j in first:-1:last
        T = kron(M[unit_model[j]], T)
    end
    T
end

function indices_kron_left(index::Vector, dimension::Int , right_dimension::Int)
    indices = Int[]
    total = dimension * right_dimension
    for i in 0:dimension:total-1
        append!(indices, index .+ i)
    end
    indices
end

function indices_kron_right(index::Vector, dimension::Int)
    indices = Int[]
    for i in 1:dimension
        append!(indices, i .+ dimension .*(index .- 1))
    end
    indices
end

