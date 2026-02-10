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
    build_spec = (transitions, G, R, S, insertstep, splicetype)
    TCoupledComponents(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp, build_spec)
end

TDCoupledComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, sojourn, dttype, splicetype="") = TDCoupledComponents(coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, sojourn, dttype, splicetype)

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
    make_components_TCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")
    make_components_TCoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)

Create TCoupledComponents structure for coupled models.

# Description
This function creates a TCoupledComponents structure for coupled models, which includes
matrix components for fitting traces, mRNA histograms, and reporter gene data. It handles
the construction of coupled systems with multiple interacting units.

# Arguments
- `coupling::Tuple`: Tuple of coupling parameters
- `transitions::Tuple`: Tuple of transition rates
- `G`: Total number of gene states
- `R`: Number of reporter steps
- `S`: Number of splice sites
- `insertstep`: Insert step for RNA processing
- `splicetype`: Splice type (default is an empty string)
- `unit_model`: Unit model specification
- `sources`: Source units for each unit
- `source_state`: Source state specification
- `target_transition`: Target transition specification

# Methods
- `make_components_TCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")`: Creates a TCoupledComponents structure from coupling parameters and transition rates.
- `make_components_TCoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)`: Creates a TCoupledComponents structure from unit model and other parameters.

# Returns
- `TCoupledComponents`: The created TCoupledComponents structure
"""

function make_components_TRGCoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)
    comp = TRGCoupledUnitComponents[]
    for i in eachindex(G)
        push!(comp, make_components_TRGCoupledUnit(source_state[i], target_transition[i], transitions[i], G[i], R[i], S[i], insertstep[i], splicetype))
    end
    build_spec = (transitions, G, R, S, insertstep, splicetype)
    TCoupledComponents{typeof(comp)}(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp, build_spec)
end

make_components_TRGCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="") = make_components_TRGCoupled(coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, splicetype)

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
    build_spec = (transitions, G, R, S, insertstep, splicetype)
    TCoupledComponents{typeof(comp)}(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp, build_spec)
end

TCoupledComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="") = TCoupledComponents(coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, splicetype)


function TForcedComponents(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype, f=set_elements_TGRS)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = f(transitions, G, R, S, insertstep, indices, splicetype)
    TForcedComponents(nT, elementsT, coupling[4])
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
    make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="")

Make MTD structure for GRS models.

# Description
This function creates an MTD structure for GRS models, which is used for various types of data fitting.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of gene states.
- `R`: Number of reporter steps.
- `S`: Number of splice sites.
- `insertstep`: Insert step.
- `onstates`: Vector of on states.
- `dttype`: Data type.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTDComponents`: The created MTD structure.
"""

function MTDComponents(transitions::Tuple, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="", ejectnumber=1)
    MTDComponents(MComponents(transitions, G, R, nhist, decay, splicetype, ejectnumber), TDComponents(transitions, G, R, S, insertstep, onstates, dttype, splicetype, ejectnumber))
end

"""
    make_components_MT(transitions, G, R, S, insertstep, nhist, decay, splicetype="")

Return MTComponents structure for GRS models.

# Description
This function returns an MTComponents structure for GRS models, which is used for fitting traces and mRNA histograms.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of gene states.
- `R`: Number of reporter steps.
- `S`: Number of splice sites.
- `insertstep`: Insert step.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTComponents`: The created MTComponents structure.
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
    if components.elementsGs[1].a < 1 || components.elementsGs[1].b < 1
        Gs = spzeros(0)
    else
        # Gs = make_mat_GC(nS, nS, couplingStrength)
        Gs = make_mat_S(components.elementsGs, nG)
    end
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
    if components.elementsSource[1].a < 1 || components.elementsSource[1].b < 1
        Gs = spzeros(0)
    else
        Gs = make_mat_S(components.elementsSource, nG)
    end
    Gt = make_mat(components.elementsTarget, rates, nG)
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

# Build per-connection (U, elementsTarget, nT_α, β, α) from (s,t) lists in components using build_spec; same order as to_connections.
function _connection_data_from_components(components::TCoupledComponents{Vector{TCoupledUnitComponents}}, build_spec)
    transitions, G, R, S, insertstep, splicetype = build_spec
    comp = components.modelcomponents
    ncoupling = sum(length.(components.sources))
    ncoupling == 0 && return nothing
    source_state = tuple([c.sourceState for c in comp]...)
    target_transition = tuple([c.targetTransition for c in comp]...)
    conns = to_connections((components.model, components.sources, source_state, target_transition, ncoupling))
    data = Any[]
    for (β, α, s, t) in conns
        (s isa Int && s == 0) && continue
        trans_α = transitions[α]
        G_α, R_α, S_α = G[α], R[α], S[α]
        nT_α = comp[α].nT
        indices_α = set_indices(length(trans_α), R_α, S_α, insertstep[α])
        U = make_mat_S(set_elements_Source(s, G[β], R[β], S[β], splicetype), comp[β].nT)
        elementsTarget = set_elements_Target(t, trans_α, G_α, R_α, S_α, insertstep[α], indices_α, nT_α, splicetype)
        push!(data, (U, elementsTarget, nT_α, β, α))
    end
    length(data) == ncoupling ? data : nothing
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
    connection_UV = (components.build_spec !== nothing ? _connection_data_from_components(components, components.build_spec) : nothing)
    return T, Source, Target, IT, connection_UV
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
- `SparseMatrixCSC`: The created coupled transition matrix TC
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
    T, Source, Target, IT, connection_UV = make_matvec_C(components, rates)
    # Generalized coupling: per-connection (U,V) built from (s,t) lists when build_spec set
    if connection_UV !== nothing && length(connection_UV) == length(coupling_strength)
        return _make_mat_TC_connection_data(T, IT, components.model, connection_UV, coupling_strength, rates)
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

# connection_data: list of (U_k, elementsTarget_k, nT_α_k, β_k, α_k); build V_k from rates in loop.
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
        U_k, elementsTarget_k, nT_α_k, β_k, α_k = connection_data[k]
        V_k = make_mat(elementsTarget_k, rates[α_k], nT_α_k)
        Tc += coupling_strength[k] * _coupling_term(U_k, V_k, β_k, α_k, IT, unit_model, n)
    end
    Tc
end

function make_mat_T(components::TCoupledComponents{Vector{TCoupledUnitComponents}}, r::Tuple)
    rates, coupling_strength = r
    make_mat_TC(components, rates, coupling_strength)
end

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

Create dwell time transition matrices for coupled system.

# Description
This function creates dwell time transition matrices for a coupled system by
applying the coupled transition matrix construction and then filtering for
sojourn states. It handles multiple dwell time matrices for different units.

# Arguments
- `unit::Int`: Unit index in the coupled system
- `TD::Vector`: Vector of dwell time matrices
- `Gm`: Vector of gene matrices
- `Gs`: Vector of source matrices
- `Gt`: Vector of target matrices
- `IT`: Vector of identity matrices
- `IG`: Vector of gene identity matrices
- `IR`: Vector of RNA identity matrices
- `coupling_strength`: Coupling strength parameters
- `sources`: Source specifications
- `model`: Model specification
- `sojourn::Vector{Vector{Int}}`: Vector of sojourn state vectors for each unit

# Returns
- `Vector{SparseMatrixCSC}`: Vector of dwell time transition matrices with zeros removed
"""
function make_mat_TCD(unit::Int, TD::Vector, Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model, sojourn::Vector{Vector{Int}})
    TCD = Vector{SparseMatrixCSC}(undef, length(TD))
    for i in eachindex(TD)
        TCD[i] = make_mat_TC(unit, TD[i], Gm, Gs, Gt, IT, IG, IR, coupling_strength, sources, model)
        make_mat_TCD!(TCD[i], sojourn[i])
    end
    return dropzeros.(TCD)
end
