# This file is part of StochasticGene.jl
#
# transition_rate_make.jl
#

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
    make_components_MTD(transitions, G, R, S, insertstep, onstates, dttype, nhist, decay, splicetype::String="")

TBW
"""
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
    make_components_MTRG(transitions, G, R, S, insertstep, nhist, decay, splicetype="")

TBW
"""
make_components_MTRG(transitions, G, R, S, insertstep, nhist, decay, splicetype="") = MTRGComponents(make_components_MRG(transitions, G, R, nhist, decay), make_components_TRG(transitions, G, R, S, insertstep, splicetype))



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

function make_components_MRG(transitions, G, R, nhist, decay)
    indices = set_indices(length(transitions), R)
    elementsG, elementsRGbar, elementsRG, nR = set_elements_GRS(transitions, G, R, 0, 1, indices, "")
    elementsB = set_elements_BRG(G, R, indices.nu[R+1])
    U, Um, Up = make_mat_U(nhist, decay)
    MRGComponents(G, nR, elementsG, elementsRGbar, elementsRG, elementsB, U, Um, Up)
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
    make_components_TRG(transitions, G, R, S, insertstep, splicetype)

TBW
"""
function make_components_TRG(transitions, G, R, S, insertstep, splicetype)
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsRGbar, elementsRG, nR, nT = set_elements_GRS(transitions, G, R, S, insertstep, indices, splicetype)
    TRGComponents(nT, G, nR, elementsG, elementsRGbar, elementsRG)
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
    make_components_ModelCoupled(source_state, target_transition, transitions, G, R, S, insertstep, splicetype="")


"""
function make_components_ModelCoupled(source_state, target_transition, transitions, G, R, S, insertstep, splicetype="")
    indices = set_indices(length(transitions), R, S, insertstep)
    elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG, nR, nT = set_elements_Coupled(source_state, target_transition, transitions, G, R, S, insertstep, indices, splicetype)
    ModelCoupledComponents(nT, G, nR, source_state, target_transition, elementsG, elementsGt, elementsGs, elementsRGbar, elementsRG)
end

"""
    make_components_Tcoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="")


"""
make_components_Tcoupled(coupling::Tuple, transitions::Tuple, G, R, S, insertstep, splicetype="") = make_components_Tcoupled(coupling[1], coupling[2], coupling[3], coupling[4], transitions, G, R, S, insertstep, splicetype)

"""
    make_components_Tcoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)


"""
function make_components_Tcoupled(unit_model, sources, source_state, target_transition, transitions, G::Tuple, R::Tuple, S::Tuple, insertstep::Tuple, splicetype)
    comp = ModelCoupledComponents[]
    for i in eachindex(G)
        push!(comp, make_components_ModelCoupled(source_state[i], target_transition[i], transitions[i], G[i], R[i], S[i], insertstep[i], splicetype))
    end
    TCoupledComponents(prod(T_dimension(G, R, S, unit_model)), unit_model, sources, comp)
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
    make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, decay::Float64, total::Int)

return M matrix used to compute steady state RNA distribution
"""
function make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, decay::Float64, total::Int)
    U, Uminus, Uplus = make_mat_U(total, decay)
    make_mat_M(T, B, U, Uminus, Uplus)
end
"""
    make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, U::SparseMatrixCSC, Uminus::SparseMatrixCSC, Uplus::SparseMatrixCSC)

TBW
"""
function make_mat_M(T::SparseMatrixCSC, B::SparseMatrixCSC, U::SparseMatrixCSC, Uminus::SparseMatrixCSC, Uplus::SparseMatrixCSC)
    nT = size(T, 1)
    total = size(U, 1)
    M = kron(U, sparse(I, nT, nT)) + kron(sparse(I, total, total), T - B) + kron(Uminus, B) + kron(Uplus, sparse(I, nT, nT))
    M[end-size(B, 1)+1:end, end-size(B, 1)+1:end] .+= B  # boundary condition to ensure probability is conserved
    return M
end

"""
    make_mat_M(components::MComponents, rates::Vector)

TBW
"""
function make_mat_M(components::MComponents, rates::Vector)
    T = make_mat(components.elementsT, rates, components.nT)
    B = make_mat(components.elementsB, rates, components.nT)
    make_mat_M(T, B, components.U, components.Uminus, components.Uplus)
end

function make_mat_B2(components, rates)
    RB = make_mat(components.elementsB, rates, components.nR)
    nG = components.nG
    kron(RB, sparse(I, nG, nG))
end

function make_mat_MRG(components::MRGComponents, rates::Vector)
    T = make_mat_TRG(components, rates)
    B = make_mat_B2(components, rates)
    make_mat_M(T, B, components.U, components.Uminus, components.Uplus)
end


"""
    make_mat_T(components::AbstractTComponents, rates)

construct matrices from elements
"""
make_mat_T(components::AbstractTComponents, rates) = make_mat(components.elementsT, rates, components.nT)

"""
    make_mat_TRG(G, GR, R, RG, nG, nR)

TBW
"""
function make_mat_TRG(G, GR, RGbar, RG, nG, nR)
    kron(RG, GR) + kron(sparse(I, nR, nR), G) + kron(RGbar, sparse(I, nG, nG))
end

"""
    make_mat_TRG(components, rates)

TBW
"""
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
    make_mat_TA(components::TAIComponents, rates)

TBW
"""
make_mat_TA(components::TAIComponents, rates) = make_mat(components.elementsTA, rates, components.nT)

"""
    make_mat_TI(components::TAIComponents, rates)

TBW
"""
make_mat_TI(components::TAIComponents, rates) = make_mat(components.elementsTI, rates, components.nT)

"""
    make_mat_TA(elementsTA, rates, nT)

TBW
"""
make_mat_TA(elementsTA, rates, nT) = make_mat(elementsTA, rates, nT)

"""
    make_mat_TI(elementsTI, rates, nT)

TBW
"""
make_mat_TI(elementsTI, rates, nT) = make_mat(elementsTI, rates, nT)

"""
    make_mat_GR(components, rates)

TBW
"""
function make_mat_GR(components, rates)
    nG = components.nG
    nR = components.nR
    G = make_mat(components.elementsG, rates, nG)
    R = make_mat(components.elementsRGbar, rates, nR)
    RB = make_mat(components.elementsRB, rates, nR)
    G, R, RB
end

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

function make_mat_GR(G)
    GR = spzeros(G, G)
    GR[G, G] = 1
    return GR
end

function make_mat_Gs(elements, nG)
    G = spzeros(nG, nG)
    for e in elements
        G[e.a, e.b] += 1.0
    end
    return G
end

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

function make_mat_TCr(components, rates, coupling_strength)
    T, G, Gt, Gs, IG, IR, IT = make_matvec_C(components, rates)
    make_mat_TCr2(coupling_strength, T, G, Gs, kron.(IR, Gt), IG, IT, components.sources, components.model), make_mat_TC(coupling_strength, G, Gs, Gt, IG, components.sources, components.model)
end

function make_mat_TC(components, rates, coupling_strength)
    T, _, Gt, Gs, _, IR, IT = make_matvec_C(components, rates)
    make_mat_TC(coupling_strength, T, kron.(IR, Gs), kron.(IR, Gt), IT, components.sources, components.model)
end

function make_mat_Tkron(T::Matrix, IT, sources, unit_model)
    Tc = spzeros(size(unit_model))
    for α in eachindex(unit_model)
        Tα = T[unit_model[α]]
        for j in α-1:-1:1
            if j ∈ sources[α]
                Tα = kron(IT[unit_model[j]], Tα)
            end
        end
        for j in α+1:n
            if j ∈ sources[α]
                Tα = kron(Tα, IT[unit_model[j]])
            end
        end
        Tc += Tα
    end
    Tc
end

function kron_forward(T, I, sources, unit_model, first, last)
    for j in first:last
        if j ∈ sources
            T = kron(T, I[unit_model[j]])
        end
    end
    T
end

function kron_forward(T, I, unit_model, first, last)
    for j in first:last
        T = kron(T, I[unit_model[j]])
    end
    T
end

function kron_backward(T, I, sources, unit_model, first, last)
    for j in first:-1:last
        if j ∈ sources
            T = kron(I[unit_model[j]], T)
        end
    end
    T
end

function kron_backward(T, I, unit_model, first, last)
    for j in first:-1:last
        T = kron(I[unit_model[j]], T)
    end
    T
end

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



######



# function make_mat_TCr(coupling_strength, T, G, Gs, V, IG, IT, sources, model)
#     n = length(model)
#     Tc = SparseMatrixCSC[]
#     for α in 1:n
#         Tα = T[model[α]]
#         Tα = kron_backward(Tα, IG, model, α - 1, 1)
#         Tα = kron_forward(Tα, IG, model, α + 1, n)
#         for β in 1:α-1
#             Gβ = IT[model[α]]
#             Gβ = kron_forward(Gβ, IG, model, α + 1, n)
#             Gβ = kron_backward(Gβ, IG, model, α - 1, β + 1)
#             Gβ = kron(G[model[β]], Gβ)
#             Gβ = kron_backward(Gβ, IG, model, β - 1, 1)
#             Tα += Gβ
#         end
#         for β in α+1:n
#             Gβ = IT[model[α]]
#             Gβ = kron_backward(Gβ, IG, model, α - 1, 1)
#             Gβ = kron_forward(Gβ, IG, model, α + 1, β - 1)
#             Gβ = kron(Gβ, G[model[β]])
#             Gβ = kron_forward(Gβ, IG, model, β + 1, n)
#             Tα += Gβ
#         end
#         for β in 1:α-1
#             Gβ = IT[model[α]]
#             Gβ = kron_forward(Gβ, IG, model, α + 1, n)
#             Gβ = kron_backward(Gβ, IG, model, α - 1, β + 1)
#             Gβ = kron(G[model[β]], Gβ)
#             Gβ = kron_backward(Gβ, IG, model, β - 1, 1)
#             Tα += Gβ
#             if β ∈ sources[α]
#                 Vβ = V[model[α]]
#                 Vβ = kron_forward(Vβ, IG, model, α + 1, n)
#                 Vβ = kron_backward(Vβ, IG, model, α - 1, β + 1)
#                 Vβ = kron(Gs[model[β]], Vβ)
#                 Vβ = kron_backward(Vβ, IG, model, β - 1, 1)
#                 Tα += coupling_strength[sources[α][β]] * Vβ
#             end
#         end
#         for β in α+1:n
#             Gβ = IT[model[α]]
#             Gβ = kron_backward(Gβ, IG, model, α - 1, 1)
#             Gβ = kron_forward(Gβ, IG, model, α + 1, β - 1)
#             Gβ = kron(Gβ, G[model[β]])
#             Gβ = kron_forward(Gβ, IG, model, β + 1, n)
#             Tα += Gβ
#             if β ∈ sources[α]
#                 Vβ = V[model[α]]
#                 Vβ = kron_backward(Vβ, IG, model, α - 1, 1)
#                 Vβ = kron_forward(Vβ, IG, model, α + 1, β - 1)
#                 Vβ = kron(Vβ, Gs[model[β]])
#                 Vβ = kron_forward(Vβ, IG, model, β + 1, n)
#                 Tα += coupling_strength[sources[α][β]] * Vβ
#             end
#         end
#         push!(Tc, Tα)
#     end
#     return Tc
# end


# function make_mat_TCr2(coupling_strength, T, G, Gs, V, IG, IT, sources, model)
#     n = length(model)
#     Tc = SparseMatrixCSC[]
#     for α in 1:n
#         Tα = T[model[α]]
#         Tα = kron_backward(Tα, IG, sources[α], model, α - 1, 1)
#         Tα = kron_forward(Tα, IG, sources[α], model, α + 1, n)
#         for β in 1:α-1
#             if β ∈ sources[α]
#                 Gβ = IT[model[α]]
#                 Gβ = kron_forward(Gβ, IG, sources[α], model, α + 1, n)
#                 Gβ = kron_backward(Gβ, IG, sources[α], model, α - 1, β + 1)
#                 Gβ = kron(G[model[β]], Gβ)
#                 Gβ = kron_backward(Gβ, IG, sources[α], model, β - 1, 1)
#                 Vβ = V[model[α]]
#                 Vβ = kron_forward(Vβ, IG, sources[α], model, α + 1, n)
#                 Vβ = kron_backward(Vβ, IG, sources[α], model, α - 1, β + 1)
#                 Vβ = kron(Gs[model[β]], Vβ)
#                 Vβ = kron_backward(Vβ, IG, sources[α], model, β - 1, 1)
#                 Tα += Gβ + coupling_strength[sources[α][β]] * Vβ
#             end
#         end
#         for β in α+1:n
#             if β ∈ sources[α]
#                 Gβ = IT[model[α]]
#                 Gβ = kron_backward(Gβ, IG, sources[α], model, α - 1, 1)
#                 Gβ = kron_forward(Gβ, IG, sources[α], model, α + 1, β - 1)
#                 Gβ = kron(Gβ, G[model[β]])
#                 Gβ = kron_forward(Gβ, IG, sources[α], model, β + 1, n)
#                 Vβ = V[model[α]]
#                 Vβ = kron_backward(Vβ, IG, sources[α], model, β - 1, 1)
#                 Vβ = kron_forward(Vβ, IG, sources[α], model, α + 1, β - 1)
#                 Vβ = kron(Vβ, Gs[model[β]])
#                 Vβ = kron_forward(Vβ, IG, sources[α], model, β + 1, n)
#                 Tα += Gβ + coupling_strength[sources[α][β]] * Vβ
#             end
#         end
#         push!(Tc, Tα)
#     end
#     return Tc
# end

# function make_mat_TCreduced(coupling_strength, T, G, Gs, V, IG, IT, sources, model)
#     n = length(model)
#     Tc = SparseMatrixCSC[]
#     for α in 1:n
#         Tα = T[model[α]]
#         for j in α-1:-1:1
#             if j ∈ sources[α]
#                 Tα = kron(IG[model[j]], Tα)
#             end
#         end
#         for j in α+1:n
#             if j ∈ sources[α]
#                 Tα = kron(Tα, IG[model[j]])
#             end
#         end
#         for β in 1:α-1
#             if β ∈ sources[α]
#                 Gβ = IT[model[α]]
#                 Vβ = V[model[α]]
#                 for j in α+1:n
#                     if j ∈ sources[α]
#                         Gβ = kron(Gβ, IG[model[j]])
#                         Vβ = kron(Vβ, IG[model[j]])
#                     end
#                 end
#                 for j in α-1:-1:β+1
#                     if j ∈ sources[α]
#                         Gβ = kron(IG[model[j]], Gβ)
#                         Vβ = kron(IG[model[j]], Vβ)
#                     end
#                 end
#                 Gβ = kron(G[model[β]], Gβ)
#                 Vβ = kron(Gs[model[β]], Vβ)
#                 for j in β-1:-1:1
#                     if j ∈ sources[α]
#                         Gβ = kron(G[model[j]], Gβ)
#                         Vβ = kron(G[model[j]], Vβ)
#                     end
#                 end
#                 Tα += Gβ
#                 println(Vβ)
#                 Tα += coupling_strength[sources[α][β]] * Vβ
#             end
#         end
#         for β in α+1:n
#             if β ∈ sources[α]
#                 Gβ = IT[model[α]]
#                 Vβ = V[model[α]]
#                 for j in α-1:-1:1
#                     if j ∈ sources[α]
#                         Gβ = kron(G[model[β]], Gβ)
#                         Vβ = kron(G[model[β]], Vβ)
#                     end
#                 end
#                 for j in α+1:β-1
#                     if j ∈ sources[α]
#                         Gβ = kron(Gβ, IG[model[j]])
#                         Vβ = kron(Vβ, IG[model[j]])
#                     end
#                 end
#                 Gβ = kron(Gβ, G[model[β]])
#                 Vβ = kron(Vβ, Gs[model[β]])
#                 for j in β+1:n
#                     if j ∈ sources[α]
#                         Gβ = kron(Gβ, IG[model[j]])
#                         Vβ = kron(Vβ, IG[model[j]])
#                     end
#                 end
#                 Tα += Gβ
#                 println(coupling_strength[sources[α][β]])
#                 Tα += coupling_strength[sources[α][β]] * Vβ
#             end
#         end
#         push!(Tc, Tα)
#     end
#     return Tc
# end