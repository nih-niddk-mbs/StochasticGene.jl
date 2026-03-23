# Optional include for makescriptcoupled.jl workflows (copy or merge into your script).
# Usage: include("contrib/makeswarm_latent_helpers.jl")  # from StochasticGene repo root
#
# Requires: using StochasticGene

"""Number of units with observed traces (hidden unit 3 with R[end]==0 → 2)."""
function n_observed_trace_units(coupling, R::Tuple)
    length(coupling[1]) < 3 && return length(coupling[1])
    length(R) >= 3 && R[end] == 0 && return 2
    return length(coupling[1])
end

"""Tie all ncoupling strengths to one parameter (group shares index totalrates+1)."""
function fixedeffects_rsum_coupling_identical(ncoupling::Int, totalrates::Int)
    ncoupling < 2 && return tuple()
    return (collect(totrates .+ (1:ncoupling)),)
end

"""Fitted params: targets + only first coupling index when all couplings are tied."""
function fittedparam_rsum_tied(target_inds::Vector{Int}, totalrates::Int, ncoupling::Int)
    ncoupling < 2 && return vcat(target_inds, collect(totrates .+ (1:ncoupling)))
    return vcat(target_inds, [totalrates + 1])
end

"""Hidden G=3: r12=r32, r21=r23 (matches `_transitions_hidden_g3_all_pairs` order)."""
function fixedeffects_hidden_g3_symmetry(unit3_rate_start::Int)
    b = unit3_rate_start - 1
    return ([b + 1, b + 4], [b + 2, b + 3])
end

"""Fitted indices for hidden block under that symmetry."""
function fittedparam_hidden_g3_symmetry(unit3_rate_start::Int)
    return [
        unit3_rate_start, unit3_rate_start + 1, unit3_rate_start + 4,
        unit3_rate_start + 5, unit3_rate_start + 6, unit3_rate_start + 7,
    ]
end
