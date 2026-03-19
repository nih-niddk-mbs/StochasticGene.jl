# This file is part of StochasticGene.jl
#
# transition_rate_make.jl
#
# This file contains constructors and matrix building functions for creating the
# various types of transition rate matrices used in stochastic gene expression modeling.
# It serves as the main interface for constructing the different component structures
# and assembling them into complete transition matrices.
#
# Key functionality includes:
# - Component constructors for different types of transition matrices
# - Matrix construction functions that assemble elements into sparse matrices
# - Specialized matrix builders for coupled systems
# - mRNA distribution matrix constructors (M matrices)
# - Combined MT matrix constructors
# - Coupled system matrix constructors
# - Dwell time analysis matrix constructors
#
# This file ties together the structures from transition_rate_structures.jl,
# the element generators from transition_rate_elements.jl, and the utility functions
# from transition_rate_functions.jl to create the complete transition matrices
# needed for stochastic gene expression analysis.

### Component constructors

"""
    TComponents(transitions, G, R, S, insertstep, splicetype, f=set_elements_TGRS)

Constructor for TComponents structure

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of gene states.
- `R`: Number of reporter steps.
- `S`: Number of splice sites.
- `insertstep`: Insert step.
- `splicetype`: Splice type.

# Returns
- `TComponents`: The created TComponents structure.

"""
function TComponents(transitions::Tuple, G, R, S, insertstep, splicetype, f=set_elements_TGRS)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = f(transitions, G, R, S, insertstep, indices, splicetype)
    TComponents(nT, elementsT)
end


# """
#     TRGComponents(transitions, G, R, S, insertstep, splicetype)

# TBW
# """
# function TRGComponents(transitions::Tuple, G, R, S, insertstep, splicetype)
#     indices = set_indices(length(transitions), R, S, insertstep)
#     elementsG, elementsRGbar, elementsRG, nR, nT = set_elements_GRS(transitions, G, R, S, insertstep, indices, splicetype)
#     TRGComponents(nT, G, nR, elementsG, elementsRGbar, elementsRG)
# end

"""
    TAIComponents(elementsT::Vector, nT::Int, onstates::Vector)
    TAIComponents(transitions, G, R, S, insertstep, onstates, splicetype::String="")

Return TAIComponents structure.

# Description
This function returns a TAIComponents structure, which includes matrix components for fitting traces and creating TA and TI matrix components.

# Arguments
- `elementsT`: Transition elements.
- `nT::Int`: Number of transition elements.
- `onstates::Vector`: Vector of on states.
- `transitions`: Transition rates.
- `G`: Total number of gene states.
- `R`: Number of reporter steps.
- `S`: Number of splice sites.
- `insertstep`: Insert step.
- `splicetype`: Splice type (default is an empty string).

# Methods
- `TAIComponents(elementsT, nT::Int, onstates::Vector)`: Creates a TAIComponents structure from transition elements.
- `TAIComponents(transitions, G, R, S, insertstep, onstates, splicetype::String="")`: Creates a TAIComponents structure from transition rates and other parameters.

# Returns
- `TAIComponents`: The created TAIComponents structure.
"""
function TAIComponents(elementsT::Vector, nT::Int, onstates::Vector)
    TAIComponents{Element}(nT, elementsT, set_elements_TA(elementsT, onstates), set_elements_TI(elementsT, onstates))
end
function TAIComponents(transitions::Tuple, G, R, S, insertstep, onstates, splicetype::String="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
    TAIComponents(elementsT, nT, onstates)
end


# function make_components_TAI(transitions, G, R, S, insertstep, onstates, splicetype::String="")
#     indices = set_indices(length(transitions), R, S, insertstep)
#     elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
#     make_components_TAI(elementsT, nT, onstates)
# end
# function make_components_TAI(elementsT, nT::Int, onstates::Vector)
#     TAIComponents{Element}(nT, elementsT, set_elements_TA(elementsT, onstates), set_elements_TI(elementsT, onstates))
# end

# """
#     make_components_TD(transitions, G, R, S, insertstep, onstates, dttype, splicetype::String="")

# TBW
# """
# function make_components_TD(transitions, G, R, S, insertstep, onstates, dttype, splicetype::String="")
#     indices = set_indices(length(transitions), R, S, insertstep)
#     elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
#     elementsG = set_elements_G(transitions, indices.gamma)
#     c = set_elements_TDvec(elementsT, elementsG, onstates, dttype)
#     TDComponents(nT, G, elementsT, elementsG, c, [])
# end


function TDComponents(transitions::Tuple, G, R, S, insertstep, sojourn, dttype, splicetype::String="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
    elementsG = set_elements_G(transitions, indices.gamma)
    elementsTD, TDdims = set_elements_TDvec(elementsT, elementsG, sojourn, dttype, nT, G)
    TDComponents(nT, G, elementsT, elementsG, elementsTD, TDdims)
end

function TDCoupledUnitComponents(source_state, target_transition, transitions::Tuple, G, R, S, insertstep, sojourn, dttype, splicetype::String="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsSource = set_elements_Gs(source_state)
    elementsTarget = set_elements_Gt(transitions, target_transition, indices.gamma)
    elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
    elementsG = set_elements_G(transitions, indices.gamma)
    elementsTD, TDdims = set_elements_TDvec(elementsT, elementsG, sojourn, dttype, nT, G)
    TDCoupledUnitComponents(nT, G, elementsSource, elementsTarget, elementsT, elementsG, elementsTD, TDdims)
end

function TDCoupledComponents(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, sojourn, dttype, splicetype)
    comp = TDCoupledUnitComponents[]
    for i in eachindex(G)
        push!(comp, TDCoupledUnitComponents(source_state[i], target_transition[i], transitions[i], G[i], R[i], S[i], insertstep[i], sojourn[i], dttype[i], splicetype))
    end
    TCoupledComponents(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp, nothing)
end

# Build from (unit_model, connections). Source state and target transition are derived
# from connections via source_states_for_unit / target_transition_for_unit (single source of truth).
function TDCoupledComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, sojourn, dttype, splicetype="")
    unit_model, connections = coupling[1], coupling[2]
    n_units = length(unit_model)
    sources_vec = [Int[] for _ in 1:n_units]
    for (β, s, α, t) in connections
        push!(sources_vec[α], β)
    end
    sources = ntuple(i -> tuple(sources_vec[i]...), n_units)
    source_state = ntuple(i -> source_states_for_unit(connections, i), n_units)
    target_transition = ntuple(i -> target_transition_for_unit(connections, i), n_units)
    TDCoupledComponents(unit_model, sources, source_state, target_transition, transitions, G, R, S, insertstep, sojourn, dttype, splicetype)
end

"""
    make_components_TRGCoupledUnit(source_state, target_transition, transitions, G::Int, R, S, insertstep, splicetype="")

Create TRGCoupledUnitComponents structure for coupled models.

# Description
This function creates a TRGCoupledUnitComponents structure for coupled models, which includes
matrix components for fitting traces, mRNA histograms, and reporter gene data. It handles
the gene-RNA coupling in a single unit of a larger coupled system.

# Arguments
- `source_state`: Source state specification for coupling
- `target_transition`: Target transition specification
- `transitions`: Transition rates tuple
- `G::Int`: Total number of gene states
- `R`: Number of reporter steps
- `S`: Number of splice sites
- `insertstep`: Insert step for RNA processing
- `splicetype`: Splice type (default is an empty string)

# Returns
- `TRGCoupledUnitComponents`: The created TRGCoupledUnitComponents structure
"""
function make_components_TRGCoupledUnit(source_state, target_transition, transitions, G::Int, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT = set_elements_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)
    TRGCoupledUnitComponents(nT, G, nR, source_state, target_transition, elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG)
end

function TRGCoupledUnitComponents(source_state, target_transition, transitions, G::Int, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT = set_elements_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)
    TRGCoupledUnitComponents(nT, G, nR, source_state, target_transition, elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG)
end

"""
    make_components_TRGCoupled(coupling, transitions, G, R, S, insertstep, splicetype="")
    make_components_TRGCoupled(unit_model, sources, source_state, target_transition, transitions, G, R, S, insertstep, splicetype)

Build a `TCoupledComponents{Vector{TRGCoupledUnitComponents}}` for the RG-stack coupled path.

# Arguments
- `coupling::Tuple`: `(unit_model, connections)` where each connection is `(β, s, α, t)`.
- `transitions::Tuple`: per-model gene-state transition tuples.
- `G`, `R`, `S`, `insertstep`: per-model geometry (Tuples for the multi-unit form).
- `splicetype::String`: splicing model variant (default `""`).

# Returns
- `TCoupledComponents`: assembled coupled components.
"""

function make_components_TRGCoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)
    comp = TRGCoupledUnitComponents[]
    for i in eachindex(G)
        push!(comp, make_components_TRGCoupledUnit(source_state[i], target_transition[i], transitions[i], G[i], R[i], S[i], insertstep[i], splicetype))
    end
    TCoupledComponents{typeof(comp)}(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp, nothing)
end

# (unit_model, connections) format: derive per-unit data and dispatch to internal constructor.
function make_components_TRGCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")
    length(coupling) == 2 || return make_components_TRGCoupled(coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, splicetype)
    unit_model, connections = coupling[1], coupling[2]
    n_units = length(unit_model)
    sources_vec = [Int[] for _ in 1:n_units]
    for (β, s, α, t) in connections
        push!(sources_vec[α], β)
    end
    sources_t = ntuple(i -> tuple(sources_vec[i]...), n_units)
    source_state_t = ntuple(i -> source_states_for_unit(connections, i), n_units)
    target_transition_t = ntuple(i -> target_transition_for_unit(connections, i), n_units)
    make_components_TRGCoupled(unit_model, sources_t, source_state_t, target_transition_t, transitions, G, R, S, insertstep, splicetype)
end

TRGCoupledComponents(coupling, transitions, G, R, S, insertstep, splicetype="") = make_components_TRGCoupled(coupling, transitions, G, R, S, insertstep, splicetype)

function TCoupledUnitComponents(source_state, target_transition, transitions, G::Int, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, elementsSource, elementsTarget, nT = set_elements_TCoupledUnit(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)
    TCoupledUnitComponents{typeof(source_state),typeof(target_transition)}(nT, source_state, target_transition, elementsT, elementsSource, elementsTarget)
end

function TCoupledComponents(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)
    comp = TCoupledUnitComponents[]
    for i in eachindex(G)
        push!(comp, TCoupledUnitComponents(source_state[i], target_transition[i], transitions[i], G[i], R[i], S[i], insertstep[i], splicetype))
    end
    # Build per-connection records from per-unit components and model-level specification.
    connections = ConnectionRecord[]
    n_units = length(unit_model)
    for α in 1:n_units
        sα = sources[α]
        nconn = length(sα)
        nconn == 0 && continue
        comp_α = comp[α]
        ss = comp_α.sourceState
        tt = comp_α.targetTransition
        trans_α = transitions[α]
        G_α, R_α, S_α = G[α], R[α], S[α]
        nT_α = comp_α.nT
        indices_α = set_indices(length(trans_α), R_α, S_α, insertstep[α])
        for k in 1:nconn
            β = sα[k]
            s = (ss isa Tuple || ss isa AbstractVector) ? ss[k] : ss
            t = (tt isa Tuple || tt isa AbstractVector) ? tt[k] : tt
            (s isa Int && s == 0) && continue
            G_β, R_β, S_β = G[β], R[β], S[β]
            U_elements = set_elements_Source(s, G_β, R_β, S_β, splicetype)
            U = make_mat_S(U_elements, comp[β].nT)
            elementsTarget = set_elements_Target(t, trans_α, G_α, R_α, S_α, insertstep[α], indices_α, nT_α, splicetype)
            push!(connections, ConnectionRecord(β, α, U, elementsTarget, nT_α))
        end
    end
    # Always store connection list (can be empty). Tc = uncoupled T + sum over connections; empty loop does nothing.
    TCoupledComponents{typeof(comp)}(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp, connections)
end

"""
    TCoupledComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")

High-level constructor for `TCoupledComponents` using connection specifications.

Arguments:
- `coupling::Tuple`: `(unit_model, connections)` where
    - `unit_model`: tuple specifying the order of units
    - `connections`: vector of `ConnectionSpec` = `(β, s, α, t)` for each connection
- `transitions`, `G`, `R`, `S`, `insertstep`, `splicetype`: per-unit model specification.
"""
# Empty connections allowed: no interaction matrices are built; Tc = sum of uncoupled unit matrices.
# Build ConnectionRecords from the connections list here (not in the inner constructor), so we have
# correct (β, s, α, t) per connection; the inner constructor's loop used source_state[α] which is
# "states when α is a source" and would be empty for target-only units.
function TCoupledComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")
    unit_model, connections = coupling
    n_units = length(unit_model)
    sources_vec = [Int[] for _ in 1:n_units]
    for (β, s, α, t) in connections
        push!(sources_vec[α], β)
    end
    sources_t = ntuple(i -> tuple(sources_vec[i]...), n_units)
    source_state_t = ntuple(i -> source_states_for_unit(connections, i), n_units)
    target_transition_t = ntuple(i -> target_transition_for_unit(connections, i), n_units)
    comp = TCoupledUnitComponents[TCoupledUnitComponents(source_state_t[i], target_transition_t[i], transitions[i], G[i], R[i], S[i], insertstep[i], splicetype) for i in eachindex(G)]
    connection_data = ConnectionRecord[]
    for (β, s, α, t) in connections
        (s == 0) && continue
        comp_α = comp[α]
        trans_α = transitions[α]
        G_α, R_α, S_α = G[α], R[α], S[α]
        nT_α = comp_α.nT
        indices_α = set_indices(length(trans_α), R_α, S_α, insertstep[α])
        G_β, R_β, S_β = G[β], R[β], S[β]
        U_elements = set_elements_Source(s, G_β, R_β, S_β, splicetype)
        U = make_mat_S(U_elements, comp[β].nT)
        elementsTarget = set_elements_Target(t, trans_α, G_α, R_α, S_α, insertstep[α], indices_α, nT_α, splicetype)
        push!(connection_data, ConnectionRecord(β, α, U, elementsTarget, nT_α))
    end
    TCoupledComponents{typeof(comp)}(prod(T_dimension(G, R, S, unit_model)), unit_model, sources_t, comp, connection_data)
end

"""
    TCoupledFullComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")

Build full-matrix coupled components from the same arguments as TCoupledComponents.
Uncoupled part: expand each unit's elements to full (row,col) with rate indices in flat layout.
Coupling: expand each ConnectionRecord to full elements with coupling rate indices.
Legacy TCoupledComponents and make_mat_TC are unchanged.

# Building from hand-set full-space elements

You can also build TCoupledFullComponents from elements set directly in the full coupled space:

1. `nT = T_dimension(G, R, S, unit_model)` (or your slot dimensions); `N = prod(nT)`.
2. `rate_offset_per_unit = rate_offset_per_unit_from_lengths(n_rates_per_unit, n_coupling)`.
3. `elements = Element[]`; fill via [`set_elements_full_uncoupled!`](@ref) / [`set_elements_full_coupling!`](@ref) (from unit data) or
   `add_element_full!` for single entries (see [`full_state_index`](@ref) / [`full_state_from_index`](@ref)).
4. `TCoupledFullComponents(N, elements)`.
"""
function TCoupledFullComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")
    unit_model = coupling[1]
    nT = collect(T_dimension(G, R, S, unit_model))
    N = prod(nT)
    elements_base, elements_coupling = set_elements_TCoupledFull(coupling, transitions, G, R, S, insertstep, splicetype, unit_model)
    targets = [(unit_model[c[3]], c[4]) for c in coupling[2]]
    return TCoupledFullComponents(N, elements_base, elements_coupling, targets)
end

function TForcedComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype, f=set_elements_TGRS)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = f(transitions, G, R, S, insertstep, indices, splicetype)
    if length(coupling) == 2
        unit_model, connections = coupling
        n_units = length(unit_model)
        sources_vec = [Int[] for _ in 1:n_units]
        for (β, s, α, t) in connections
            push!(sources_vec[α], β)
        end
        sources_t = ntuple(i -> tuple(sources_vec[i]...), n_units)
        targets_out = targets((unit_model, sources_t))
        targets_t = ntuple(i -> tuple(targets_out[i]...), length(targets_out))
        TForcedComponents(nT, elementsT, targets_t)
    else
        TForcedComponents(nT, elementsT, coupling[4])
    end
end
"""
    MComponents(transitions::Tuple, G, R, nhist, decay, splicetype)

Return MComponents structure for fitting mRNA histograms.

# Description
This function returns an MComponents structure, which includes matrix components for fitting mRNA histograms.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of gene states.
- `R`: Number of reporter steps.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type.

# Returns
- `MComponents`: The created MComponents structure.
"""
function MComponents(transitions::Tuple, G, R, nhist, decay, splicetype, ejectnumber=1)
    indices = set_indices(length(transitions), R)
    elementsT, nT = set_elements_TGRS(transitions, G, R, 0, 1, indices, splicetype)
    elementsB = set_elements_B(G, R, indices.nu[R+1])
    U, Um, Up = make_mat_U(nhist, decay, ejectnumber)
    return MComponents(elementsT, elementsB, nT, U, Um, Up)
end

"""
    MTAIComponents(transitions::Tuple, G, R, S, insertstep, onstates, nhist, decay, splicetype="")

Make MTAI structure for GRS models (fitting mRNA, on time, and off time histograms).

# Description
This function creates an MTAI structure for GRS models, which is used for fitting mRNA, on time, and off time histograms.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of gene states.
- `R`: Number of reporter steps.
- `S`: Number of splice sites.
- `insertstep`: Insert step.
- `onstates`: Vector of on states.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTAIComponents`: The created MTAI structure.
"""
function MTAIComponents(transitions::Tuple, G, R, S, insertstep, onstates, nhist, decay, splicetype="", ejectnumber=1)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_TGRS(transitions, G, R, S, insertstep, indices, splicetype)
    MTAIComponents(MComponents(transitions, G, R, nhist, decay, splicetype, ejectnumber), TAIComponents(elementsT, nT, onstates))
end

"""
    MTDComponents(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype="", ejectnumber=1)

Construct `MTDComponents` for dwell-time GRS models: builds `MComponents` for the
mRNA matrix and `TDComponents` for the T and filtered TD matrices.

# Returns
- `MTDComponents`: combined M and TD components.
"""
function MTDComponents(transitions::Tuple, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="", ejectnumber=1)
    MTDComponents(MComponents(transitions, G, R, nhist, decay, splicetype, ejectnumber), TDComponents(transitions, G, R, S, insertstep, onstates, dttype, splicetype, ejectnumber))
end

"""
    MTComponents(transitions, G, R, S, insertstep, nhist, decay, splicetype="", ejectnumber=1)

Construct `MTComponents` for GRS models: combines `MComponents` (mRNA matrix)
with `TComponents` (transition matrix).

# Returns
- `MTComponents`: combined M and T components.
"""
MTComponents(transitions::Tuple, G, R, S, insertstep, nhist, decay, splicetype="", ejectnumber=1) = MTComponents(MComponents(transitions, G, R, nhist, decay, splicetype, ejectnumber), TComponents(transitions, G, R, S, insertstep, splicetype))

# """
#     make_components_MTRG(transitions, G, R, S, insertstep, nhist, decay, splicetype="")

# Return MTRGComponents structure for GRS models.

# # Description
# This function returns an MTRGComponents structure for GRS models, which is used for fitting traces, mRNA histograms, and reporter gene data.

# # Arguments
# - `transitions`: Transition rates.
# - `G`: Total number of gene states.
# - `R`: Number of reporter steps.
# - `S`: Number of splice sites.
# - `insertstep`: Insert step.
# - `nhist`: Number of histograms.
# - `decay`: Decay rates.
# - `splicetype`: Splice type (default is an empty string).

# # Returns
# - `MTRGComponents`: The created MTRGComponents structure.
# """
# make_components_MTRG(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MTRGComponents(make_components_MRG(transitions, G, R, nhist, decay), make_components_TRG(transitions, G, R, S, insertstep, splicetype))

# MTRGComponents(transitions::Tuple, G, R, S, insertstep, nhist, decay, splicetype="") = MTRGComponents(MComponents(transitions, G, R, nhist, decay, splicetype), TRGComponents(transitions, G, R, S, insertstep, splicetype))

### matrix constructors

"""
    make_mat!(T::AbstractMatrix, elements::Vector, rates::Vector)

Set matrix elements of T to those in vector argument elements.

# Description
This function sets the matrix elements of T to those in the vector argument elements, using the provided rates.

# Arguments
- `T::AbstractMatrix`: The matrix to be modified.
- `elements::Vector`: Vector of elements to set in the matrix.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Returns
- `Nothing`: This function modifies the matrix T in place.
"""
function make_mat!(T::AbstractMatrix, elements::Vector, rates::Vector)
    for e in elements
        T[e.a, e.b] += e.pm * rates[e.index]
    end
end

"""
    make_mat(elements::Vector, rates::Vector, nT::Int)

Return an nT x nT sparse matrix T.

# Description
This function returns an nT x nT sparse matrix T, with elements set according to the provided elements and rates.

# Arguments
- `elements::Vector`: Vector of elements to set in the matrix.
- `rates::Vector`: Vector of rates corresponding to the elements.
- `nT::Int`: Size of the matrix (nT x nT).

# Returns
- `SparseMatrixCSC`: The created sparse matrix T.
"""
function make_mat(elements::Vector, rates::Vector, nT::Int)
    T = spzeros(nT, nT)
    make_mat!(T, elements, rates)
    return T
end

"""
    make_mat_U(total::Int, decay::Float64)

Return matrices U, Uminus, and Uplus for m transitions.

# Description
This function returns matrices U, Uminus, and Uplus for m transitions, using the provided total number of transitions and decay rate.

# Arguments
- `total::Int`: Total number of transitions.
- `decay::Float64`: Decay rate.

# Returns
- `Tuple{SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC}`: The created matrices U, Uminus, and Uplus.
"""
function make_mat_U(total::Int, decay::Float64, ejectnumber::Int)
    U = spzeros(total, total)
    Uminus = spzeros(total, total)
    Uplus = spzeros(total, total)
    # Generate matrices for m transitions
    Uplus[1, 2] = decay
    for m = 2:total-1
        U[m, m] = -decay * (m - 1)
        m - ejectnumber > 0 && (Uminus[m, m-ejectnumber] = 1)
        Uplus[m, m+1] = decay * m
    end
    U[total, total] = -decay * (total - 1)
    total - ejectnumber > 0 && (Uminus[total, total-ejectnumber] = 1)
    return U, Uminus, Uplus
end

"""
    make_mat_U(total::Int, decay::Float64, mean_eject::Float64)

Create U matrix with Poisson-distributed ejection

# Arguments
- `total::Int`: Total number of mRNA states
- `decay::Float64`: mRNA decay rate
- `mean_eject::Float64`: Mean number of mRNAs ejected

# Returns
- `U::SparseMatrixCSC`: Transition matrix
- `Uminus::SparseMatrixCSC`: Ejection matrix
- `Uplus::SparseMatrixCSC`: Birth matrix
"""

function make_mat_U(total::Int, decay::Float64, ejectnumber::Tuple)
    U = spzeros(total, total)
    Uminus = spzeros(total, total)
    Uplus = spzeros(total, total)

    # Generate matrices for m transitions
    Uplus[1, 2] = decay
    for m = 2:total-1
        U[m, m] = -decay * (m - 1)
        Uplus[m, m+1] = decay * m
    end
    U[total, total] = -decay * (total - 1)

    d = NegativeBinomial(ejectnumber[1], ejectnumber[2])
    max_eject = ceil(Int, mean(d) + 3 * std(d))  # 3 standard deviations
    # Set ejection terms using Poisson distribution
    for m = 1:total
        # Calculate Poisson probabilities
        for k = 0:min(m, max_eject)
            if m - k >= 1
                Uminus[m, m-k] = pdf(d, k)
            end
        end
    end

    return U, Uminus, Uplus
end


"""
    make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, decay::Float64, total::Int)
    make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, U::SparseMatrixCSC, Uminus::SparseMatrixCSC, Uplus::SparseMatrixCSC)
    make_mat_M(components::MComponents, rates::Vector)

Return M matrix used to compute steady state RNA distribution.

# Description
This function returns the M matrix used to compute the steady state RNA distribution. It can either generate the U, Uminus, and Uplus matrices internally, accept them as arguments, or use components to create the matrix.

# Arguments
- `T::SparseMatrixCSC`: Transition matrix.
- `B::SparseMatrixCSC`: Boundary matrix.
- `decay::Float64`: Decay rate.
- `total::Int`: Total number of transitions.
- `U::SparseMatrixCSC`: Precomputed U matrix.
- `Uminus::SparseMatrixCSC`: Precomputed Uminus matrix.
- `Uplus::SparseMatrixCSC`: Precomputed Uplus matrix.
- `components::MComponents`: Components structure containing elements and matrices.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Methods
- `make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, decay::Float64, total::Int)`: Generates the U, Uminus, and Uplus matrices internally.
- `make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, U::SparseMatrixCSC, Uminus::SparseMatrixCSC, Uplus::SparseMatrixCSC)`: Accepts precomputed U, Uminus, and Uplus matrices.
- `make_mat_M(components::MComponents, rates::Vector)`: Uses components to create the matrix.

# Returns
- `SparseMatrixCSC`: The created M matrix.
"""
function make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, decay::Float64, total::Int, ejectnumber=1)
    U, Uminus, Uplus = make_mat_U(total, decay, ejectnumber)
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

# """
#     make_mat_MRG(components::MRGComponents, rates::Vector)

# Return MRG matrix used to compute steady state RNA distribution for GRS models.

# # Description
# This function returns the MRG matrix used to compute the steady state RNA distribution for GRS models. It uses the provided components and rates to create the matrix.

# # Arguments
# - `components::MRGComponents`: Components structure containing elements and matrices for GRS models.
# - `rates::Vector`: Vector of rates corresponding to the elements.

# # Returns
# - `SparseMatrixCSC`: The created MRG matrix.
# """
# function make_mat_MRG(components::MRGComponents, rates::Vector)
#     T = make_mat_TRG(components, rates)
#     B = make_mat_B2(components, rates)
#     make_mat_M(T, B, components.U, components.Uminus, components.Uplus)
# end

"""
    make_mat_B2(components, rates)

Create boundary matrix B2 for GRS models.

# Description
This function creates the boundary matrix B2 for GRS models, which represents
the boundary conditions for RNA processing. It constructs the matrix using
the provided components and rates.

# Arguments
- `components`: Components structure containing elements and matrices for GRS models
- `rates::Vector`: Vector of rates corresponding to the elements

# Returns
- `SparseMatrixCSC`: The created boundary matrix B2
"""
function make_mat_B2(components, rates)
    RB = make_mat(components.elementsB, rates, components.nR)
    nG = components.nG
    kron(RB, sparse(I, nG, nG))
end

make_mat_G(components, rates) = make_mat(components.elementsG, rates, components.nG)


"""
    make_mat_T(components::AbstractTComponents, rates)

Construct matrices from elements.

# Description
This function constructs matrices from elements using the provided components and rates.

# Arguments
- `components::AbstractTComponents`: Components structure containing elements and matrices.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Returns
- `SparseMatrixCSC`: The created matrix T.
"""
make_mat_T(components::AbstractTComponents, rates) = make_mat(components.elementsT, rates, components.nT)

"""
    make_mat_TRG(G, GR, RGbar, RG, nG, nR)
    make_mat_TRG(components, rates)

Return TRG matrix used to compute steady state RNA distribution for GRS models.

# Description
This function returns the TRG matrix used to compute the steady state RNA distribution for GRS models. It can either accept precomputed matrices or use components and rates to create the matrix.

# Arguments
- `G`: Gene matrix.
- `GR`: Gene reporter matrix.
- `RGbar`: Reporter gene boundary matrix.
- `RG`: Reporter gene matrix.
- `nG`: Number of genes.
- `nR`: Number of reporters.
- `components`: Components structure containing elements and matrices for GRS models.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Methods
- `make_mat_TRG(G, GR, RGbar, RG, nG, nR)`: Accepts precomputed matrices.
- `make_mat_TRG(components, rates)`: Uses components and rates to create the matrix.

# Returns
- `SparseMatrixCSC`: The created TRG matrix.
"""
function make_mat_TRG(G, GR, RGbar, RG, nG, nR)
    kron(RG, GR) + kron(sparse(I, nR, nR), G) + kron(RGbar, sparse(I, nG, nG))
end

function make_mat_TRG(components, rates)
    nG = components.nG
    nR = components.nR
    GR = make_mat_GR(nG)
    G = make_mat(components.elementsG, rates, nG)
    RGbar = make_mat(components.elementsRGbar, rates, nR)
    RG = make_mat(components.elementsRG, rates, nR)
    make_mat_TRG(G, GR, RGbar, RG, nG, nR)
end

"""
    make_mat_TA(elementsTA, rates, nT)
    make_mat_TA(components::TAIComponents, rates)

Return TA matrix used to compute steady state RNA distribution.

# Description
This function returns the TA matrix used to compute the steady state RNA distribution. It can either accept elements directly or use components and rates to create the matrix.

# Arguments
- `elementsTA`: Transition elements for TA.
- `rates::Vector`: Vector of rates corresponding to the elements.
- `nT::Int`: Number of transition elements.
- `components::TAIComponents`: Components structure containing elements and matrices for TA.

# Methods
- `make_mat_TA(elementsTA, rates, nT)`: Accepts elements directly.
- `make_mat_TA(components::TAIComponents, rates)`: Uses components and rates to create the matrix.

# Returns
- `SparseMatrixCSC`: The created TA matrix.
"""
make_mat_TA(elementsTA, rates, nT) = make_mat(elementsTA, rates, nT)
make_mat_TA(components::TAIComponents, rates) = make_mat(components.elementsTA, rates, components.nT)

"""
    make_mat_TI(components::TAIComponents, rates)
    make_mat_TI(elementsTI, rates, nT)

Return TI matrix used to compute steady state RNA distribution.

# Description
This function returns the TI matrix used to compute the steady state RNA distribution. It can either accept elements directly or use components and rates to create the matrix.

# Arguments
- `components::TAIComponents`: Components structure containing elements and matrices for TI.
- `rates::Vector`: Vector of rates corresponding to the elements.
- `elementsTI`: Transition elements for TI.
- `nT::Int`: Number of transition elements.

# Methods
- `make_mat_TI(components::TAIComponents, rates)`: Uses components and rates to create the matrix.
- `make_mat_TI(elementsTI, rates, nT)`: Accepts elements directly.

# Returns
- `SparseMatrixCSC`: The created TI matrix.
"""
make_mat_TI(components::TAIComponents, rates) = make_mat(components.elementsTI, rates, components.nT)

make_mat_TI(elementsTI, rates, nT) = make_mat(elementsTI, rates, nT)

"""
    make_mat_GR(components, rates)
    make_mat_GR(G)

Return GR matrix used to compute steady state RNA distribution.

# Description
This function returns the GR matrix used to compute the steady state RNA distribution. It can either accept components and rates to create the matrix or use a given size to create an identity matrix.

# Arguments
- `components`: Components structure containing elements and matrices for GR.
- `rates::Vector`: Vector of rates corresponding to the elements.
- `G::Int`: Size of the GR matrix.

# Methods
- `make_mat_GR(components, rates)`: Uses components and rates to create the matrix.
- `make_mat_GR(G)`: Creates an identity matrix of size G.

# Returns
- `Tuple{SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC}`: The created GR matrices (G, R, RB) for the first method.
- `SparseMatrixCSC`: The created GR matrix for the second method.
"""
function make_mat_GR(components, rates)
    nG = components.nG
    nR = components.nR
    G = make_mat(components.elementsG, rates, nG)
    R = make_mat(components.elementsRGbar, rates, nR)
    RB = make_mat(components.elementsRB, rates, nR)
    G, R, RB
end
function make_mat_GR(G)
    GR = spzeros(G, G)
    GR[G, G] = 1
    return GR
end

"""
    make_mat_S(elements, nG)

Create source matrix Gs for coupled systems.

# Description
This function creates the Gs matrix used in coupled systems to represent
source states. It constructs a sparse matrix from the provided elements.

# Arguments
- `elements`: Transition elements for Gs
- `nG::Int`: Number of gene states

# Returns
- `SparseMatrixCSC`: The created Gs matrix
"""
function make_mat_S(elements, nG)
    G = spzeros(nG, nG)
    for e in elements
        G[e.a, e.b] += 1.0
    end
    return G
end

# Only build when there is something to build. Empty elements => no source coupling for this unit.
function _make_mat_source(elements, nG)
    isempty(elements) && return spzeros(0, 0)
    return make_mat_S(elements, nG)
end

"""
    make_mat_C(components, rates)

Create matrices for coupled system analysis.

# Description
This function creates the matrices needed for coupled system analysis, including
transition matrices, gene matrices, target matrices, source matrices, and identity
matrices. The specific matrices returned depend on the component type.

# Arguments
- `components`: Components structure containing elements and matrices for coupled models
- `rates::Vector`: Vector of rates corresponding to the elements

# Returns
- `Tuple`: Various matrices depending on component type (T, G, Gt, Gs, I_G, I_R, I_T for TRGCoupledUnitComponents)
- `Tuple`: (T, Source, Target, IT) for TCoupledUnitComponents
- `Tuple`: (T, TD, G, Gt, Gs, IG, IR, IT) for TDCoupledUnitComponents
"""
function make_mat_C(components::TRGCoupledUnitComponents, rates)
    nT = components.nT
    nG = components.nG
    nR = components.nR
    GR = make_mat_GR(nG)
    G = make_mat(components.elementsG, rates, nG)
    RGbar = make_mat(components.elementsRGbar, rates, nR)
    RG = make_mat(components.elementsRG, rates, nR)
    T = make_mat_TRG(G, GR, RGbar, RG, nG, nR)
    Gt = make_mat(components.elementsGt, rates, nG)
    Gs = _make_mat_source(components.elementsGs, nG)
    return T, G, Gt, Gs, sparse(I, nG, nG), sparse(I, nR, nR), sparse(I, nT, nT)
end

function make_mat_C(components::TCoupledUnitComponents, rates)
    nT = components.nT
    T = make_mat_T(components, rates)
    Target = make_mat(components.elementsTarget, rates, nT)
    Source = make_mat_S(components.elementsSource, nT)
    return T, Source, Target, sparse(I, nT, nT)
end

function make_mat_CD(components::TDCoupledUnitComponents, rates)
    nT = components.nT
    T = make_mat_T(components, rates)
    Target = make_mat(components.elementsTarget, rates, nT)
    Source = make_mat_S(components.elementsSource, nT)
    TD = SparseMatrixCSC[]
    for i in eachindex(components.elementsTD)
        push!(TD, make_mat(components.elementsTD[i], rates, components.TDdims[i]))
    end
    return T, TD, Source, Target, sparse(I, nT, nT)
end

function make_mat_C(components::TDCoupledUnitComponents, rates)
    nT = components.nT
    nG = components.nG
    nR = div(nT, nG)
    T = make_mat(components.elementsT, rates, nT)
    TD = SparseMatrixCSC[]
    for i in eachindex(components.elementsTD)
        push!(TD, make_mat(components.elementsTD[i], rates, components.TDdims[i]))
    end
    G = make_mat(components.elementsG, rates, nG)
    Gt = make_mat(components.elementsTarget, rates, nG)
    Gs = _make_mat_source(components.elementsSource, nG)
    return T, TD, G, Gt, Gs, sparse(I, nG, nG), sparse(I, nR, nR), sparse(I, nT, nT)
end


"""
    make_matvec_C(components, rates)

Create matrix vectors for coupled system analysis.

# Description
This function creates vectors of matrices for coupled system analysis. It processes
each unit in the coupled system and creates the appropriate matrices for that unit.
The specific matrices returned depend on the component type.

# Arguments
- `components`: Components structure containing elements and matrices for coupled models
- `rates::Vector`: Vector of rates corresponding to the elements

# Returns
- `Tuple`: Vectors of matrices depending on component type
"""
function make_matvec_C(components::TCoupledComponents{Vector{TRGCoupledUnitComponents}}, rates)
    n = length(components.model)
    T = Vector{SparseMatrixCSC}(undef, n)
    G = Vector{SparseMatrixCSC}(undef, n)
    Gt = Vector{SparseMatrixCSC}(undef, n)
    Gs = Vector{SparseMatrixCSC}(undef, n)
    IG = Vector{SparseMatrixCSC}(undef, n)
    IR = Vector{SparseMatrixCSC}(undef, n)
    IT = Vector{SparseMatrixCSC}(undef, n)
    for i in eachindex(components.model)
        T[i], G[i], Gt[i], Gs[i], IG[i], IR[i], IT[i] = make_mat_C(components.modelcomponents[i], rates[i])
    end
    return T, G, Gt, Gs, IG, IR, IT
end

function make_matvec_C(components::TCoupledComponents{Vector{TCoupledUnitComponents}}, rates)
    n = length(components.model)
    T = Vector{SparseMatrixCSC}(undef, n)
    Source = Vector{SparseMatrixCSC}(undef, n)
    Target = Vector{SparseMatrixCSC}(undef, n)
    IT = Vector{SparseMatrixCSC}(undef, n)
    for i in eachindex(components.model)
        T[i], Source[i], Target[i], IT[i] = make_mat_C(components.modelcomponents[i], rates[i])
    end
    return T, Source, Target, IT
end

function make_matvec_C(components::TCoupledComponents{Vector{TDCoupledUnitComponents}}, rates)
    n = length(components.model)
    T = Vector{SparseMatrixCSC}(undef, n)
    TD = Vector{Vector{SparseMatrixCSC}}(undef, n)
    G = Vector{SparseMatrixCSC}(undef, n)
    Gt = Vector{SparseMatrixCSC}(undef, n)
    Gs = Vector{SparseMatrixCSC}(undef, n)
    IG = Vector{SparseMatrixCSC}(undef, n)
    IR = Vector{SparseMatrixCSC}(undef, n)
    IT = Vector{SparseMatrixCSC}(undef, n)
    for i in eachindex(rates)
        T[i], TD[i], G[i], Gt[i], Gs[i], IG[i], IR[i], IT[i] = make_mat_C(components.modelcomponents[i], rates[i])
    end
    return T, TD, G, Gt, Gs, IG, IR, IT
end

"""
    make_mat_TC(coupling_strength, T, U, V, IT, sources, unit_model)
    make_mat_TC(components, rates, coupling_strength)

Create coupled transition matrix TC.

# Description
This function creates the coupled transition matrix TC used for coupled system analysis.
It combines individual unit transition matrices with coupling terms to represent
the full coupled system dynamics.

# Arguments
- `coupling_strength`: Strength of the coupling between units
- `T`: Vector of individual unit transition matrices
- `U`: Vector of source matrices
- `V`: Vector of target matrices
- `IT`: Vector of identity matrices
- `sources`: Vector of source indices for each unit
- `unit_model`: Vector specifying the order of units
- `components`: Components structure containing elements and matrices for coupled models
- `rates::Vector`: Vector of rates corresponding to the elements

# Methods
- `make_mat_TC(coupling_strength, T, U, V, IT, sources, unit_model)`: Uses the provided matrices and indices directly
- `make_mat_TC(components, rates, coupling_strength)`: Uses components and rates to create the matrix

# Returns
- `SparseMatrixCSC`: The created coupled transition matrix TC (uncoupled T + sum over connections; empty connection list adds nothing).
"""
function make_mat_TC(coupling_strength, T, U, V, IT, sources, unit_model)
    n = length(unit_model)
    N = prod(size.(IT, 2))
    # Canonical order of (source, target) pairs: same as coupling_strength vector order
    coupling_pairs = Tuple{Int,Int}[(β, α) for α in 1:n for β in sources[α]]
    Tc = spzeros(N, N)
    for α in 1:n
        Tα = T[unit_model[α]]
        Tα = kron_left(Tα, IT, unit_model, α - 1, 1)
        Tα = kron_right(Tα, IT, unit_model, α + 1, n)
        for β in 1:α-1
            if β ∈ sources[α]
                Vβ = V[unit_model[α]]
                Vβ = kron_right(Vβ, IT, unit_model, α + 1, n)
                Vβ = kron_left(Vβ, IT, unit_model, α - 1, β + 1)
                Vβ = kron(U[unit_model[β]], Vβ)
                Vβ = kron_left(Vβ, IT, unit_model, β - 1, 1)
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Tα += coupling_strength[idx] * Vβ
            end
        end
        for β in α+1:n
            if β ∈ sources[α]
                # Coupling term must have same state order (1..n) as Tα. Build V[α]⊗I⊗..⊗U[β]⊗I..
                # (Original order kron_left(V[α],..,β-1,1); kron(U[β],Vβ) gave (β,α,α) → size n_β·n_α² ≠ N when n_α≠1.)
                Vβ = V[unit_model[α]]
                Vβ = kron_right(Vβ, IT, unit_model, α + 1, β - 1)
                Vβ = kron(Vβ, U[unit_model[β]])
                Vβ = kron_right(Vβ, IT, unit_model, β + 1, n)
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Tα += coupling_strength[idx] * Vβ
            end
        end
        Tc += Tα
    end
    return Tc
end

function make_mat_TC_with_self_coupling(coupling_strength, T, U, V, IT, sources, unit_model)
    n = length(unit_model)
    N = prod(size.(IT, 2))
    coupling_pairs = Tuple{Int,Int}[(β, α) for α in 1:n for β in sources[α]]
    Tc = spzeros(N, N)
    for α in 1:n
        Tα = T[unit_model[α]]
        Tα = kron_left(Tα, IT, unit_model, α - 1, 1)
        Tα = kron_right(Tα, IT, unit_model, α + 1, n)
        
        # Add self-coupling case: (U[α]*V[α]) in α slot, same state order (1..n) as Tα
        if α ∈ sources[α]  # Check if unit α can couple to itself
            Vα = U[unit_model[α]] * V[unit_model[α]]  # matrix product, n_α×n_α
            Vα = kron_left(Vα, IT, unit_model, α - 1, 1)
            Vα = kron_right(Vα, IT, unit_model, α + 1, n)
            idx = findfirst(isequal((α, α)), coupling_pairs)
            Tα += coupling_strength[idx] * Vα
        end
        
        # Coupling from previous units
        for β in 1:α-1
            if β ∈ sources[α]
                Vβ = V[unit_model[α]]
                Vβ = kron_right(Vβ, IT, unit_model, α + 1, n)
                Vβ = kron_left(Vβ, IT, unit_model, α - 1, β + 1)
                Vβ = kron(U[unit_model[β]], Vβ)
                Vβ = kron_left(Vβ, IT, unit_model, β - 1, 1)
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Tα += coupling_strength[idx] * Vβ
            end
        end
        
        # Coupling from later units
        for β in α+1:n
            if β ∈ sources[α]
                Vβ = V[unit_model[α]]
                Vβ = kron_right(Vβ, IT, unit_model, α + 1, β - 1)
                Vβ = kron(Vβ, U[unit_model[β]])
                Vβ = kron_right(Vβ, IT, unit_model, β + 1, n)
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Tα += coupling_strength[idx] * Vβ
            end
        end

        Tc += Tα
    end
    return Tc
end

function make_mat_TC(components::TCoupledComponents{Vector{TRGCoupledUnitComponents}}, rates, coupling_strength)
    T, _, Gt, Gs, _, IR, IT = make_matvec_C(components, rates)
    make_mat_TC(coupling_strength, T, kron.(IR, Gs), kron.(IR, Gt), IT, components.sources, components.model)
end

function make_mat_TC(components::TCoupledComponents{Vector{TCoupledUnitComponents}}, rates, coupling_strength)
    T, Source, Target, IT = make_matvec_C(components, rates)
    # Build from connection list: Tc = uncoupled T + sum over connections (loop over empty does nothing).
    connections = components.connections
    if connections !== nothing && length(connections) == length(coupling_strength)
        return _make_mat_TC_connection_data(T, IT, components.model, connections, coupling_strength, rates)
    end
    make_mat_TC(coupling_strength, T, Source, Target, IT, components.sources, components.model)
end

# Build coupling term matrix for one connection: U at position β, V at position α (same state order 1..n as Tα).
function _coupling_term(U, V, β::Int, α::Int, IT, unit_model, n::Int)
    if β < α
        Vβ = kron_right(V, IT, unit_model, α + 1, n)
        Vβ = kron_left(Vβ, IT, unit_model, α - 1, β + 1)
        Vβ = kron(U, Vβ)
        Vβ = kron_left(Vβ, IT, unit_model, β - 1, 1)
    else
        Vβ = kron_right(V, IT, unit_model, α + 1, β - 1)
        Vβ = kron(Vβ, U)
        Vβ = kron_right(Vβ, IT, unit_model, β + 1, n)
        if α > 1
            Vβ = kron_left(Vβ, IT, unit_model, α - 1, 1)
        end
    end
    Vβ
end

# connection_data: list of ConnectionRecord (may be empty; loop adds nothing and Tc = uncoupled T). Build V_k from rates in loop.
function _make_mat_TC_connection_data(T, IT, unit_model, connection_data, coupling_strength, rates)
    n = length(unit_model)
    N = prod(size.(IT, 2))
    Tc = spzeros(eltype(T[1]), N, N)
    for α in 1:n
        Tα = T[unit_model[α]]
        Tα = kron_left(Tα, IT, unit_model, α - 1, 1)
        Tα = kron_right(Tα, IT, unit_model, α + 1, n)
        Tc += Tα
    end
    for k in eachindex(connection_data)
        rec = connection_data[k]
        U_k = rec.U
        elementsTarget_k = rec.elementsTarget
        nT_α_k = rec.nTα
        β_k = rec.β
        α_k = rec.α
        V_k = make_mat(elementsTarget_k, rates[α_k], nT_α_k)
        Tc += coupling_strength[k] * _coupling_term(U_k, V_k, β_k, α_k, IT, unit_model, n)
    end
    Tc
end

function make_mat_T(components::TCoupledComponents{Vector{TCoupledUnitComponents}}, r::Tuple)
    rates, coupling_strength = r
    make_mat_TC(components, rates, coupling_strength)
end

"""
    make_mat_TC(components::TCoupledFullComponents, rates, coupling_rates)

Build the full N×N coupled transition matrix from `TCoupledFullComponents` elements.
Iterates `elements_base` (indexed by `rates[model][localindex]`) and
`elements_coupling` (indexed by `coupling_rates[localindex]`).

# Arguments
- `rates::Vector{Vector{Float64}}`: per-model rate vectors.
- `coupling_rates::Vector{Float64}`: precomputed `γ_k * base_rate_k` for each coupling.

# Returns
- `SparseMatrixCSC{Float64}`: the assembled T matrix.
"""
function make_mat_TC(components::TCoupledFullComponents,
                     rates::Vector{Vector{Float64}},
                     coupling_rates::Vector{Float64})
    T = spzeros(Float64, components.N, components.N)
    for e in components.elements_base
        T[e.a, e.b] += e.pm * rates[e.idx.model][e.idx.localindex]
    end
    for e in components.elements_coupling
        T[e.a, e.b] += e.pm * coupling_rates[e.idx.localindex]
    end
    return T
end

"""
    make_mat_TC(components::TDCoupledFullComponents, rates, coupling_rates)

Build the full T matrix from a `TDCoupledFullComponents` (same element loop as `TCoupledFullComponents`).
"""
function make_mat_TC(components::TDCoupledFullComponents,
                     rates::Vector{Vector{Float64}},
                     coupling_rates::Vector{Float64})
    T = spzeros(Float64, components.N, components.N)
    for e in components.elements_base
        T[e.a, e.b] += e.pm * rates[e.idx.model][e.idx.localindex]
    end
    for e in components.elements_coupling
        T[e.a, e.b] += e.pm * coupling_rates[e.idx.localindex]
    end
    return T
end

"""
    make_mat_TCD(components::TDCoupledFullComponents, unit, dtype, rates, coupling_rates)

Build the TD matrix for `unit` and `dtype` index from pre-filtered elements. The TD matrix
is T restricted to sojourn-state columns (transitions to non-sojourn states are absent).
"""
function make_mat_TCD(components::TDCoupledFullComponents, unit::Int, dtype::Int,
                      rates::Vector{Vector{Float64}}, coupling_rates::Vector{Float64})
    TD = spzeros(Float64, components.N, components.N)
    for e in components.elementsTD_base[unit][dtype]
        TD[e.a, e.b] += e.pm * rates[e.idx.model][e.idx.localindex]
    end
    for e in components.elementsTD_coupling[unit][dtype]
        TD[e.a, e.b] += e.pm * coupling_rates[e.idx.localindex]
    end
    return TD
end

"""
    make_mat_TC_gene_full(coupling_strength, Gm, Gs, Gt, IG, sources, model)

Gene-level coupled transition matrix (G-only, no R): builds N_gene×N_gene with correct
Kronecker ordering for mixed nG per unit. Use for dt=true (e.g. ONG, OFFG) when units have
different G sizes or hidden units. Does not call legacy make_mat_TC.
"""
function make_mat_TC_gene_full(coupling_strength, Gm, Gs, Gt, IG, sources, model)
    unit_model = model
    n = length(unit_model)
    N = prod(size.(IG, 2))
    Tc = spzeros(N, N)
    coupling_pairs = Tuple{Int,Int}[(β, α) for α in 1:n for β in sources[α]]
    function kron_slots(Ms::Vector)
        out = Ms[1]
        for j in 2:n
            out = kron(out, Ms[j])
        end
        out
    end
    for α in 1:n
        Ms = [j == α ? Gm[unit_model[α]] : IG[unit_model[j]] for j in 1:n]
        Tα = kron_slots(Ms)
        for β in 1:α-1
            if β ∈ sources[α]
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Ms_c = [j == β ? Gs[unit_model[β]] : (j == α ? Gt[unit_model[α]] : IG[unit_model[j]]) for j in 1:n]
                Tα += coupling_strength[idx] * kron_slots(Ms_c)
            end
        end
        for β in α+1:n
            if β ∈ sources[α]
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Ms_c = [j == β ? Gs[unit_model[β]] : (j == α ? Gt[unit_model[α]] : IG[unit_model[j]]) for j in 1:n]
                Tα += coupling_strength[idx] * kron_slots(Ms_c)
            end
        end
        Tc += Tα
    end
    return Tc
end

"""
    make_mat_TC_full(unit, Tin, T_full, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)

General coupled transition matrix for dwell times: always uses full per-unit state dimensions
(N = ∏ n_α) and full U/V/IT—no marginalization shortcuts. Use this path for hidden units,
mixed R per unit, or when legacy make_mat_TC(unit, ...) would give dimension mismatches.
Implements its own N×N build so Kronecker ordering is correct for mixed dimensions; does not
call the shared make_mat_TC (legacy path unchanged).
"""
function make_mat_TC_full(unit::Int, Tin::SparseMatrixCSC, T_full, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
    unit_model = model
    n = length(unit_model)
    # Full-size U, V for every unit
    U = Vector{SparseMatrixCSC}(undef, n)
    V = Vector{SparseMatrixCSC}(undef, n)
    for i in 1:n
        if size(IR[i], 1) > 0
            U[i] = kron(IR[i], Gs[i])
            V[i] = kron(IR[i], Gt[i])
        else
            U[i] = Gs[i]
            V[i] = Gt[i]
        end
    end
    T = copy(T_full)
    T[unit] = Tin
    N = prod(size.(IT, 2))
    Tc = spzeros(N, N)
    coupling_pairs = Tuple{Int,Int}[(β, α) for α in 1:n for β in sources[α]]
    # Helper: Kronecker product in state order (1..n) with given per-slot matrices
    function kron_slots(Ms::Vector)
        out = Ms[1]
        for j in 2:n
            out = kron(out, Ms[j])
        end
        out
    end
    for α in 1:n
        # Uncoupled block: slot α = T[α], others = IT
        Ms = [j == α ? T[unit_model[α]] : IT[unit_model[j]] for j in 1:n]
        Tα = kron_slots(Ms)
        for β in 1:α-1
            if β ∈ sources[α]
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Ms_c = [j == β ? U[unit_model[β]] : (j == α ? V[unit_model[α]] : IT[unit_model[j]]) for j in 1:n]
                Tα += coupling_strength[idx] * kron_slots(Ms_c)
            end
        end
        for β in α+1:n
            if β ∈ sources[α]
                idx = findfirst(isequal((β, α)), coupling_pairs)
                Ms_c = [j == β ? U[unit_model[β]] : (j == α ? V[unit_model[α]] : IT[unit_model[j]]) for j in 1:n]
                Tα += coupling_strength[idx] * kron_slots(Ms_c)
            end
        end
        Tc += Tα
    end
    return Tc
end

# Legacy: uses G-only matrices for non-current units (marginalization shortcut). Do not use for hidden units or mixed R.
function make_mat_TC(unit::Int, Tin::SparseMatrixCSC, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
    U = copy(Gs)
    V = copy(Gt)
    ITC = copy(IG)
    T = copy(Gm)
    T[unit] = Tin
    if size(T[unit]) == size(IT[unit])
        ITC[unit] = IT[unit]
        U[unit] = kron(IR[unit], Gs[unit])
        V[unit] = kron(IR[unit], Gt[unit])
    end
    make_mat_TC(coupling_strength, T, U, V, ITC, sources, model)
end

"""
    make_mat_TCD!(TCD::SparseMatrixCSC, sojourn::Vector{Int})

Remove transitions to non-sojourn states from TCD matrix in-place.

# Description
This function modifies the TCD matrix in-place by setting to zero all transitions
that lead to states not in the sojourn set. This is used for dwell time analysis
to focus only on transitions within the sojourn states.

# Arguments
- `TCD::SparseMatrixCSC`: Transition matrix to modify (modified in-place)
- `sojourn::Vector{Int}`: Vector of sojourn state indices

# Returns
- `Nothing`: Modifies the TCD matrix in-place
"""
function make_mat_TCD!(TCD::SparseMatrixCSC, sojourn::Vector{Int})
    rows, cols, vals = findnz(TCD)
    for i in eachindex(cols)
        a = rows[i]
        b = cols[i]
        if b ∉ sojourn
            TCD[a, b] = 0
        end
    end
end

"""
    make_mat_TCD(TC::SparseMatrixCSC, sojourn::Vector{Int})

Create dwell time transition matrix TCD.

# Description
This function creates a dwell time transition matrix TCD by copying the input
transition matrix and removing transitions to non-sojourn states. The resulting
matrix is used for dwell time analysis.

# Arguments
- `TC::SparseMatrixCSC`: Input transition matrix
- `sojourn::Vector{Int}`: Vector of sojourn state indices

# Returns
- `SparseMatrixCSC`: Dwell time transition matrix with zeros removed
"""
function make_mat_TCD(TC::SparseMatrixCSC, sojourn::Vector{Int})
    TCD = copy(TC)
    make_mat_TCD!(TCD, sojourn)
    return dropzeros(TCD)
end

"""
    make_mat_TCD(unit::Int, TD::Vector, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn::Vector{Vector{Int}})

Legacy: create dwell time transition matrices using G-marginalization path (make_mat_TC(unit, ...)).
"""
function make_mat_TCD(unit::Int, TD::Vector, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn::Vector{Vector{Int}})
    TCD = Vector{SparseMatrixCSC}(undef, length(TD))
    for i in eachindex(TD)
        TCD[i] = make_mat_TC(unit, TD[i], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
        make_mat_TCD!(TCD[i], sojourn[i])
    end
    return dropzeros.(TCD)
end

"""
    make_mat_TCD_full(unit::Int, TD::Vector, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn::Vector{Vector{Int}})

General path: dwell time TCD matrices using full dimensions (make_mat_TC_full). Use for hidden units / mixed R.
"""
function make_mat_TCD_full(unit::Int, TD::Vector, T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn::Vector{Vector{Int}})
    TCD = Vector{SparseMatrixCSC}(undef, length(TD))
    for i in eachindex(TD)
        TCD[i] = make_mat_TC_full(unit, TD[i], T, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
        make_mat_TCD!(TCD[i], sojourn[i])
    end
    return dropzeros.(TCD)
end

# ---- TCoupledFullComponents: full-matrix element expansion (parallel to legacy TCoupledComponents) ----

"""
    expand_unit_elements_to_full(α::Int, elements::Vector{Element}, nT::Vector{Int}, rate_base::Int)

Compute full N×N transition elements for one uncoupled unit from its local elements.

For each local element (a,b) in the unit and each state of the other units, produces
one full-space element with row/col in 1:N and rate index offset by rate_base.

# Arguments
- `α`: Slot index (1-based).
- `elements`: Unit's local `Vector{Element}` (a,b in 1:nT[α]).
- `nT`: Slot dimensions (length n).
- `rate_base`: First rate index for this unit in the flat rate vector.

# Returns
- `Vector{Element}` with (a,b) in full space 1:N.
"""
function expand_unit_elements_to_full(α::Int, elements::Vector{Element}, nT::Vector{Int}, rate_base::Int)
    n = length(nT)
    N = prod(nT)
    block = N ÷ nT[α]
    out = Element[]
    positions = [1:(α-1); (α+1):n]
    for e in elements
        rate_idx = rate_base + e.index - 1
        for other_k in 1:block
            state = zeros(Int, n)
            r = other_k - 1
            # Decode other_k-1 so slot n varies fastest (match full_state_index convention)
            for pos in reverse(positions)
                d = nT[pos]
                state[pos] = mod(r, d) + 1
                r = div(r, d)
            end
            state[α] = e.a
            row = full_state_index(state, nT)
            state[α] = e.b
            col = full_state_index(state, nT)
            push!(out, Element(row, col, rate_idx, e.pm))
        end
    end
    return out
end

"""
    expand_coupling_to_full(U, elementsTarget, β, α, nT, gamma_index, uncoupled_elements, uncoupled_pos_index)

Compute full N×N transition elements for one coupling (source operator U at slot β,
target elements at slot α). Each element stores `index = gamma_index`, pointing
into the coupling part of the flat rate vector; the associated target rate index
is tracked separately (per-connection) in `components.target_rates`.

# Arguments
- `U`: Source unit's operator matrix (nT[β]×nT[β]).
- `elementsTarget`: Target unit's transition elements (local indices).
- `β`, `α`: Source and target slot indices (1-based).
- `nT`: Slot dimensions.
- `gamma_index`: Index of the coupling parameter γ_k in the flat rate vector.
- `uncoupled_elements`: Full uncoupled elements (used to import `pm` by position).
- `uncoupled_pos_index`: Fast lookup map `(row,col) -> index` into `uncoupled_elements`.

# Returns
- `Vector{Element}` with (a,b) in full space 1:N; index = gamma_index for all elements.
"""
function expand_coupling_to_full(U::SparseMatrixCSC,
                                 elementsTarget::Vector{Element},
                                 β::Int, α::Int,
                                 nT::Vector{Int},
                                 gamma_index::Int,
                                 uncoupled_elements::Vector{ElementCoupledFull},
                                 uncoupled_pos_index::Dict{Tuple{Int,Int},Int})
    n = length(nT)
    N = prod(nT)
    block_other = N ÷ (nT[β] * nT[α])
    positions = [1:(β-1); (β+1):(α-1); (α+1):n]
    other_dims = nT[positions]
    out = ElementCoupledFull[]
    Urows, Ucols, Uvals = findnz(U)
    for o in 1:block_other
        state = zeros(Int, n)
        r = o - 1
        for (ii, pos) in enumerate(reverse(positions))
            state[pos] = mod(r, other_dims[ii]) + 1
            r = div(r, other_dims[ii])
        end
        for k in eachindex(Urows)
            iβ, jβ = Urows[k], Ucols[k]
            u_sign = sign(Uvals[k])
            u_sign == 0 && continue
            for e in elementsTarget
                row_state = copy(state)
                row_state[β] = jβ
                row_state[α] = e.a
                col_state = copy(state)
                col_state[β] = iβ
                col_state[α] = e.b
                row = full_state_index(row_state, nT)
                col = full_state_index(col_state, nT)
                # Import sign from uncoupled element at same matrix position.
                # This enforces "coupling modifies existing transition rate with same sign".
                pos = (row, col)
                haskey(uncoupled_pos_index, pos) || continue
                base_idx = uncoupled_pos_index[pos]
                pm = uncoupled_elements[base_idx].pm
                idx = IndexCoupledFull(0, gamma_index)
                push!(out, ElementCoupledFull(row, col, idx, Int8(pm)))
            end
        end
    end
    return out
end

"""
    rate_offset_per_unit_from_lengths(n_rates_per_unit::Vector{Int}, n_coupling::Int)

Build `rate_offset_per_unit` for TCoupledFullComponents from per-unit rate counts.

# Arguments
- `n_rates_per_unit`: Length n_units; n_rates_per_unit[i] = number of rate indices for slot i.
- `n_coupling`: Number of coupling strength parameters.

# Returns
- `Vector{Int}` of length n_units+1: rate_offset_per_unit[i] = last index of slot i-1
  (with rate_offset_per_unit[1] = 0). Base-rate indices for slot i are then
  (rate_offset_per_unit[i] + 1) : rate_offset_per_unit[i] + n_rates_per_unit[i],
  and the first coupling index is rate_offset_per_unit[end] + 1.
"""
function rate_offset_per_unit_from_lengths(n_rates_per_unit::Vector{Int}, n_coupling::Int)
    offsets = Int[0]
    for r in n_rates_per_unit
        push!(offsets, offsets[end] + r)
    end
    return offsets
end

"""
    TDCoupledFullComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, sojourn, dttype)

Build dwell-time coupled components in full state space. Shares element construction
with `TCoupledFullComponents`; sojourn sets are expanded to full-state indices.
"""
function TDCoupledFullComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, sojourn, dttype)
    unit_model = coupling[1]
    nT_vec = collect(T_dimension(G, R, S, unit_model))
    n_units = length(unit_model)
    elements_base, elements_coupling = set_elements_TCoupledFull(coupling, transitions, G, R, S, insertstep, "", unit_model)
    N = prod(nT_vec)
    elementsTD_base = Vector{Vector{Vector{ElementCoupledFull}}}(undef, n_units)
    elementsTD_coupling = Vector{Vector{Vector{ElementCoupledFull}}}(undef, n_units)
    for k in 1:n_units
        α = unit_model[k]
        n_dtype = length(sojourn[α])
        elementsTD_base[k] = Vector{Vector{ElementCoupledFull}}(undef, n_dtype)
        elementsTD_coupling[k] = Vector{Vector{ElementCoupledFull}}(undef, n_dtype)
        for i in 1:n_dtype
            soj_unit = sojourn[α][i]
            if occursin("G", dttype[α][i])
                soj_unit = g_sojourn_to_T_sojourn(soj_unit, G[α], R[α], S[α])
            end
            soj = full_state_indices_for_unit_sojourn(k, soj_unit, nT_vec)
            elementsTD_base[k][i] = filter(e -> e.b ∈ soj, elements_base)
            elementsTD_coupling[k][i] = filter(e -> e.b ∈ soj, elements_coupling)
        end
    end
    targets = [(unit_model[c[3]], c[4]) for c in coupling[2]]
    return TDCoupledFullComponents(N, elements_base, elements_coupling, elementsTD_base, elementsTD_coupling, targets)
end
