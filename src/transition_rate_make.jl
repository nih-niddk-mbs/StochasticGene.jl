# This file is part of StochasticGene.jl
#
# transition_rate_make.jl
#

"""
    make_components_T(transitions, G, R, S, insertstep, splicetype)

Return TComponents structure for fitting traces and creating TA and TI matrix components.

# Description
This function returns a TComponents structure, which is used for fitting traces and creating TA and TI matrix components.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `splicetype`: Splice type.

# Returns
- `TComponents`: The created TComponents structure.
"""
function make_components_T(transitions, G, R, S, insertstep, splicetype)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    TComponents(nT, elementsT)
end
"""
    make_components_TRG(transitions, G, R, S, insertstep, splicetype)

Return TRGComponents structure for GRS models.

# Description
This function returns a TRGComponents structure for GRS models, which includes matrix components for fitting traces, mRNA histograms, and reporter gene data.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `splicetype`: Splice type.

# Returns
- `TRGComponents`: The created TRGComponents structure.
"""
function make_components_TRG(transitions, G, R, S, insertstep, splicetype)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsRGbar, elementsRG, nR, nT = set_elements_GRS(transitions, G, R, S, insertstep, indices, splicetype)
    TRGComponents(nT, G, nR, elementsG, elementsRGbar, elementsRG)
end

function make_components_TD(transitions, G, R, S, insertstep, onstates, dttype, splicetype::String="")
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
    end
    TDComponents(nT, elementsT, elementsTG, c)
end
"""
    make_components_M(transitions, G, R, nhist, decay, splicetype)

Return MComponents structure for fitting mRNA histograms.

# Description
This function returns an MComponents structure, which includes matrix components for fitting mRNA histograms.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type.

# Returns
- `MComponents`: The created MComponents structure.
"""
function make_components_M(transitions, G, R, nhist, decay, splicetype)
    indices = set_indices(length(transitions), R)
    elementsT, nT = set_elements_T(transitions, G, R, 0, 1, indices, splicetype)
    elementsB = set_elements_B(G, R, indices.nu[R+1])
    U, Um, Up = make_mat_U(nhist, decay)
    return MComponents(elementsT, elementsB, nT, U, Um, Up)
end

"""
    make_components_MRG(transitions, G, R, nhist, decay)

Return MRGComponents structure for GRS models.

# Description
This function returns an MRGComponents structure for GRS models, which includes matrix components for fitting mRNA histograms and reporter gene data.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `nhist`: Number of histograms.
- `decay`: Decay rates.

# Returns
- `MRGComponents`: The created MRGComponents structure.
"""
function make_components_MRG(transitions, G, R, nhist, decay)
    indices = set_indices(length(transitions), R)
    elementsG, elementsRGbar, elementsRG, nR = set_elements_GRS(transitions, G, R, 0, 1, indices, "")
    elementsB = set_elements_BRG(G, R, indices.nu[R+1])
    U, Um, Up = make_mat_U(nhist, decay)
    MRGComponents(G, nR, elementsG, elementsRGbar, elementsRG, elementsB, U, Um, Up)
end



"""
    make_components_MTAI(transitions, G, R, S, insertstep, onstates, nhist, decay, splicetype="")

Make MTAI structure for GRS models (fitting mRNA, on time, and off time histograms).

# Description
This function creates an MTAI structure for GRS models, which is used for fitting mRNA, on time, and off time histograms.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `onstates`: Vector of on states.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTAIComponents`: The created MTAI structure.
"""
function make_components_MTAI(transitions, G, R, S, insertstep, onstates, nhist, decay, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
    MTAIComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), make_components_TAI(elementsT, nT, onstates))
end


"""
    make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="")

Make MTD structure for GRS models.

# Description
This function creates an MTD structure for GRS models, which is used for various types of data fitting.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `onstates`: Vector of on states.
- `dttype`: Data type.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTDComponents`: The created MTD structure.
"""
# function make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="")
#     indices = set_indices(length(transitions), R, S, insertstep)
#     elementsT, nT = set_elements_T(transitions, G, R, S, insertstep, indices, splicetype)
#     elementsTG = set_elements_T(transitions, indices.gamma)
#     c = Vector{Element}[]
#     for i in eachindex(onstates)
#         if dttype[i] == "ON"
#             push!(c, set_elements_TA(elementsT, onstates[i]))
#         elseif dttype[i] == "OFF"
#             push!(c, set_elements_TI(elementsT, onstates[i]))
#         elseif dttype[i] == "ONG"
#             push!(c, set_elements_TA(elementsTG, onstates[i]))
#         elseif dttype[i] == "OFFG"
#             push!(c, set_elements_TI(elementsTG, onstates[i]))
#         end
#         # dttype[i] == "ON" ? push!(c, set_elements_TA(elementsT, onstates[i])) : push!(c, set_elements_TI(elementsT, onstates[i]))
#     end
#     MTDComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), TDComponents(nT, elementsT, elementsTG, c))
# end

function make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="")
    MTDComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), make_components_TD(transitions, G, R, S, insertstep, onstates, dttype, splicetype))
end



"""
    make_components_MT(transitions, G, R, S, insertstep, nhist, decay, splicetype="")

Return MTComponents structure for GRS models.

# Description
This function returns an MTComponents structure for GRS models, which is used for fitting traces and mRNA histograms.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTComponents`: The created MTComponents structure.
"""
make_components_MT(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MTComponents(make_components_M(transitions, G, R, nhist, decay, splicetype), make_components_T(transitions, G, R, S, insertstep, splicetype))

"""
    make_components_MTRG(transitions, G, R, S, insertstep, nhist, decay, splicetype="")

Return MTRGComponents structure for GRS models.

# Description
This function returns an MTRGComponents structure for GRS models, which is used for fitting traces, mRNA histograms, and reporter gene data.

# Arguments
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `nhist`: Number of histograms.
- `decay`: Decay rates.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `MTRGComponents`: The created MTRGComponents structure.
"""
make_components_MTRG(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MTRGComponents(make_components_MRG(transitions, G, R, nhist, decay), make_components_TRG(transitions, G, R, S, insertstep, splicetype))



"""
    make_components_TAI(elementsT, nT::Int, onstates::Vector)
    make_components_TAI(transitions, G, R, S, insertstep, onstates, splicetype::String="")

Return TAIComponents structure.

# Description
This function returns a TAIComponents structure, which includes matrix components for fitting traces and creating TA and TI matrix components.

# Arguments
- `elementsT`: Transition elements.
- `nT::Int`: Number of transition elements.
- `onstates::Vector`: Vector of on states.
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `splicetype`: Splice type (default is an empty string).

# Methods
- `make_components_TAI(elementsT, nT::Int, onstates::Vector)`: Creates a TAIComponents structure from transition elements.
- `make_components_TAI(transitions, G, R, S, insertstep, onstates, splicetype::String="")`: Creates a TAIComponents structure from transition rates and other parameters.

# Returns
- `TAIComponents`: The created TAIComponents structure.
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
    make_components_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, splicetype="")

Return TRGCoupledComponents structure for coupled models.

# Description
This function returns a TRGCoupledComponents structure for coupled models, which includes matrix components for fitting traces, mRNA histograms, and reporter gene data.

# Arguments
- `source_state`: Source state.
- `target_transition`: Target transition.
- `transitions`: Transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `splicetype`: Splice type (default is an empty string).

# Returns
- `TRGCoupledComponents`: The created TRGCoupledComponents structure.
"""
function make_components_TRGCoupledUnit(source_state, target_transition, transitions, G::Int, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT = set_elements_TRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)
    TRGCoupledUnitComponents(nT, G, nR, source_state, target_transition, elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG)
end

function make_components_TDRG(source_state, target_transition, transitions, G::Int, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, elementsGD, elementsRD, nR, nT = set_elements_TDRGCoupled(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)
    TDRGCoupledUnitComponents(nT, G, nR, source_state, target_transition, elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, elementsGD, elementsRD)
end

"""
    make_components_TCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")
    make_components_TCoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)

Return TCoupledComponents structure for coupled models.

# Description
This function returns a TCoupledComponents structure for coupled models, which includes matrix components for fitting traces, mRNA histograms, and reporter gene data.

# Arguments
- `coupling::Tuple`: Tuple of coupling parameters.
- `transitions::Tuple`: Tuple of transition rates.
- `G`: Total number of genes.
- `R`: Number of reporters.
- `S`: Number of states.
- `insertstep`: Insert step.
- `splicetype`: Splice type (default is an empty string).
- `unit_model`: Unit model.
- `sources`: Source units for each unit.
- `source_state`: Source state.
- `target_transition`: Target transition.

# Methods
- `make_components_TCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")`: Creates a TCoupledComponents structure from coupling parameters and transition rates.
- `make_components_TCoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)`: Creates a TCoupledComponents structure from unit model and other parameters.

# Returns
- `TCoupledComponents`: The created TCoupledComponents structure.
"""
m
function make_components_TCoupled(f, comp, unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)
    for i in eachindex(G)
        push!(comp, f(source_state[i], target_transition[i], transitions[i], G[i], R[i], S[i], insertstep[i], splicetype))
    end
    TCoupledComponents{typeof(comp)}(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp)
end

make_components_TRGCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="") = make_components_TCoupled(make_components_TRGCoupledUnit, TRGCoupledUnitComponents[], coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, splicetype)

make_components_TDRGCoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="") = make_components_TCoupled(make_components_TDRGCoupledUnit, TDRGCoupledUnitComponents[], coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, splicetype)

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
    make_mat_MRG(components::MRGComponents, rates::Vector)

Return MRG matrix used to compute steady state RNA distribution for GRS models.

# Description
This function returns the MRG matrix used to compute the steady state RNA distribution for GRS models. It uses the provided components and rates to create the matrix.

# Arguments
- `components::MRGComponents`: Components structure containing elements and matrices for GRS models.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Returns
- `SparseMatrixCSC`: The created MRG matrix.
"""
function make_mat_MRG(components::MRGComponents, rates::Vector)
    T = make_mat_TRG(components, rates)
    B = make_mat_B2(components, rates)
    make_mat_M(T, B, components.U, components.Uminus, components.Uplus)
end

"""
    make_mat_B2(components, rates)

Return boundary matrix B2 for GRS models.

# Description
This function returns the boundary matrix B2 for GRS models, using the provided components and rates.

# Arguments
- `components`: Components structure containing elements and matrices for GRS models.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Returns
- `SparseMatrixCSC`: The created boundary matrix B2.
"""
function make_mat_B2(components, rates)
    RB = make_mat(components.elementsB, rates, components.nR)
    nG = components.nG
    kron(RB, sparse(I, nG, nG))
end

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
    make_mat_Gs(elements, nG)

Return Gs matrix used to compute steady state RNA distribution.

# Description
This function returns the Gs matrix used to compute the steady state RNA distribution. It accepts transition elements directly and constructs the matrix.

# Arguments
- `elements`: Transition elements for Gs.
- `nG::Int`: Number of gene elements.

# Returns
- `SparseMatrixCSC`: The created Gs matrix.
"""
function make_mat_Gs(elements, nG)
    G = spzeros(nG, nG)
    for e in elements
        G[e.a, e.b] += 1.0
    end
    return G
end

"""
    make_mat_C(components, rates)

Return matrices used to compute steady state RNA distribution for coupled models.

# Description
This function returns the matrices used to compute the steady state RNA distribution for coupled models. It uses the provided components and rates to create the matrices.

# Arguments
- `components`: Components structure containing elements and matrices for coupled models.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Returns
- `Tuple{SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC}`: The created matrices (T, G, Gt, Gs, I_G, I_R, I_T).
"""
function make_mat_C(components, rates)
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
        Gs = make_mat_Gs(components.elementsGs, nG)
    end
    return T, G, Gt, Gs, sparse(I, nG, nG), sparse(I, nR, nR), sparse(I, nT, nT)
end

function make_mat_TDC(components, rates, onstates)
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
        Gs = make_mat_Gs(components.elementsGs, nG)
    end
    return TD
end

"""
    make_matvec_C(components, rates)

Return matrix-vector product used to compute steady state RNA distribution for coupled models.

# Description
This function returns the matrix-vector product used to compute the steady state RNA distribution for coupled models. It uses the provided components and rates to create the product.

# Arguments
- `components`: Components structure containing elements and matrices for coupled models.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Returns
- `SparseVector`: The created matrix-vector product.
"""
function make_matvec_C(components, rates)
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

function make_matvec_DC(components, rates)
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


"""
    make_mat_TC(coupling_strength, T, U, V, IT, sources, unit_model)
    make_mat_TC(components, rates, coupling_strength)

Return TC matrix used to compute steady state RNA distribution for coupled models.

# Description
This function returns the TC matrix used to compute the steady state RNA distribution for coupled models. It can either use the provided coupling strength, matrices, identity matrices, sources, and unit model indices directly, or use components and rates to create the matrix.

# Arguments
- `coupling_strength`: Strength of the coupling.
- `T`: Initial matrix.
- `U`: Matrix U.
- `V`: Matrix V.
- `IT`: Identity matrix for T.
- `sources`: Vector of source indices.
- `unit_model`: Vector of unit model indices.
- `components`: Components structure containing elements and matrices for coupled models.
- `rates::Vector`: Vector of rates corresponding to the elements.

# Methods
- `make_mat_TC(coupling_strength, T, U, V, IT, sources, unit_model)`: Uses the provided matrices and indices directly.
- `make_mat_TC(components, rates, coupling_strength)`: Uses components and rates to create the matrix.

# Returns
- `SparseMatrixCSC`: The created TC matrix.
"""
function make_mat_TC(coupling_strength, T, U, V, IT, sources, unit_model)
    n = length(unit_model)
    N = prod(size.(IT, 2))
    Tc = spzeros(N, N)
    for α in 1:n
        Tα = T[unit_model[α]]
        Tα = kron_backward(Tα, IT, unit_model, α - 1, 1)
        Tα = kron_forward(Tα, IT, unit_model, α + 1, n)
        for β in 1:α-1
            if β ∈ sources[α]
                Vβ = V[unit_model[α]]
                Vβ = kron_forward(Vβ, IT, unit_model, α + 1, n)
                Vβ = kron_backward(Vβ, IT, unit_model, α - 1, β + 1)
                Vβ = kron(U[unit_model[β]], Vβ)
                Vβ = kron_backward(Vβ, IT, unit_model, β - 1, 1)
                Tα += coupling_strength[sources[α][β]] * Vβ
            end
        end
        for β in α+1:n
            if β ∈ sources[α]
                Vβ = V[unit_model[α]]
                Vβ = kron_backward(Vβ, IT, unit_model, β - 1, 1)
                Vβ = kron_forward(Vβ, IT, unit_model, α + 1, β - 1)
                Vβ = kron(U[unit_model[β]], Vβ)
                Vβ = kron_forward(Vβ, IT, unit_model, β + 1, n)
                Tα += coupling_strength[sources[α][β]] * Vβ
            end
        end
        Tc += Tα
    end
    return Tc
end

function make_mat_TCD(coupling_strength, T, G, Gs, V, IG, IT, sources, model)
    n = length(model)
    Tc = SparseMatrixCSC[]
    for α in 1:n
        Tα = T[model[α]]
        Tα = kron_backward(Tα, IG, sources[α], model, α - 1, 1)
        Tα = kron_forward(Tα, IG, sources[α], model, α + 1, n)
        for β in 1:α-1
            if β ∈ sources[α]
                Gβ = IT[model[α]]
                Gβ = kron_forward(Gβ, IG, sources[α], model, α + 1, n)
                Gβ = kron_backward(Gβ, IG, sources[α], model, α - 1, β + 1)
                Gβ = kron(G[model[β]], Gβ)
                Gβ = kron_backward(Gβ, IG, sources[α], model, β - 1, 1)
                Vβ = V[model[α]]
                Vβ = kron_forward(Vβ, IG, sources[α], model, α + 1, n)
                Vβ = kron_backward(Vβ, IG, sources[α], model, α - 1, β + 1)
                Vβ = kron(Gs[model[β]], Vβ)
                Vβ = kron_backward(Vβ, IG, sources[α], model, β - 1, 1)
                Tα += Gβ + coupling_strength[sources[α][β]] * Vβ
            end
        end
        for β in α+1:n
            if β ∈ sources[α]
                Gβ = IT[model[α]]
                Gβ = kron_backward(Gβ, IG, sources[α], model, α - 1, 1)
                Gβ = kron_forward(Gβ, IG, sources[α], model, α + 1, β - 1)
                Gβ = kron(Gβ, G[model[β]])
                Gβ = kron_forward(Gβ, IG, sources[α], model, β + 1, n)
                Vβ = V[model[α]]
                Vβ = kron_backward(Vβ, IG, sources[α], model, β - 1, 1)
                Vβ = kron_forward(Vβ, IG, sources[α], model, α + 1, β - 1)
                Vβ = kron(Vβ, Gs[model[β]])
                Vβ = kron_forward(Vβ, IG, sources[α], model, β + 1, n)
                Tα += Gβ + coupling_strength[sources[α][β]] * Vβ
            end
        end
        push!(Tc, Tα)
    end
    return Tc
end

function make_mat_TC(components, rates, coupling_strength)
    T, _, Gt, Gs, _, IR, IT = make_matvec_C(components, rates)
    make_mat_TC(coupling_strength, T, kron.(IR, Gs), kron.(IR, Gt), IT, components.sources, components.model)
end


function make_mat_TCD(components::TCoupledComponents{Vector{TDRGCoupledUnitComponents}}, rates, coupling_strength)
    T, G, Gt, Gs, IG, IR, IT = make_mat_C(components, rates)
    make_mat_TA
    make_mat_TCD(coupling_strength, T, G, Gs, V, IG, IT, components.sources, components.model)

end

