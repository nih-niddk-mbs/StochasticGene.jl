# This file is part of StochasticGene.jl

"""
    Likelihood functions for composable trait models
"""

"""
    likelihoodfn(model::ComposableModel, data::ComposableData)

Compute likelihood function based on model and data traits.
"""
function likelihoodfn(model::ComposableModel, data::ComposableData)
    # Base case - RNA data
    if has_data_trait(data, RNATrait)
        rna_trait = get_data_trait(data, RNATrait)
        return likelihood_rna(model, rna_trait.nRNA, rna_trait.histRNA)
    end
    
    # ON/OFF data
    if has_data_trait(data, OnOffTrait)
        onoff_trait = get_data_trait(data, OnOffTrait)
        return likelihood_onoff(model, onoff_trait.bins, onoff_trait.ON, onoff_trait.OFF)
    end
    
    # Dwell time data
    if has_data_trait(data, DwellTimeTrait)
        dwell_trait = get_data_trait(data, DwellTimeTrait)
        return likelihood_dwelltime(model, dwell_trait.bins, dwell_trait.dwell_times, dwell_trait.dt_types)
    end
    
    # Trace data
    if has_data_trait(data, TraceTrait)
        trace_trait = get_data_trait(data, TraceTrait)
        return likelihood_trace(model, trace_trait.interval, trace_trait.trace)
    end
    
    error("No matching likelihood function found for given model and data traits")
end

"""
    likelihood_rna(model::ComposableModel, nRNA::Int, histRNA::Vector)

RNA histogram likelihood for composable models.
"""
function likelihood_rna(model::ComposableModel, nRNA::Int, histRNA::Vector)
    # Handle different model traits
    if has_trait(model, CouplingTrait)
        return likelihood_rna_coupled(model, nRNA, histRNA)
    elseif has_trait(model, GridTrait)
        return likelihood_rna_grid(model, nRNA, histRNA)
    elseif has_trait(model, HierarchicalTrait)
        return likelihood_rna_hierarchical(model, nRNA, histRNA)
    else
        # Base GRSM case
        return likelihood_rna_base(model, nRNA, histRNA)
    end
end

"""
    likelihood_onoff(model::ComposableModel, bins::Vector, ON::Vector, OFF::Vector)

ON/OFF time likelihood for composable models.
"""
function likelihood_onoff(model::ComposableModel, bins::Vector, ON::Vector, OFF::Vector)
    if has_trait(model, CouplingTrait)
        return likelihood_onoff_coupled(model, bins, ON, OFF)
    elseif has_trait(model, GridTrait)
        return likelihood_onoff_grid(model, bins, ON, OFF)
    elseif has_trait(model, HierarchicalTrait)
        return likelihood_onoff_hierarchical(model, bins, ON, OFF)
    else
        return likelihood_onoff_base(model, bins, ON, OFF)
    end
end

"""
    likelihood_dwelltime(model::ComposableModel, bins::Vector{Vector}, dwell_times::Vector{Vector}, dt_types::Vector)

Dwell time likelihood for composable models.
"""
function likelihood_dwelltime(model::ComposableModel, bins::Vector{Vector}, dwell_times::Vector{Vector}, dt_types::Vector)
    if has_trait(model, CouplingTrait)
        return likelihood_dwelltime_coupled(model, bins, dwell_times, dt_types)
    elseif has_trait(model, GridTrait)
        return likelihood_dwelltime_grid(model, bins, dwell_times, dt_types)
    elseif has_trait(model, HierarchicalTrait)
        return likelihood_dwelltime_hierarchical(model, bins, dwell_times, dt_types)
    else
        return likelihood_dwelltime_base(model, bins, dwell_times, dt_types)
    end
end

"""
    likelihood_trace(model::ComposableModel, interval::Float64, trace::Any)

Trace likelihood for composable models.
"""
function likelihood_trace(model::ComposableModel, interval::Float64, trace::Any)
    if has_trait(model, CouplingTrait)
        return likelihood_trace_coupled(model, interval, trace)
    elseif has_trait(model, GridTrait)
        return likelihood_trace_grid(model, interval, trace)
    elseif has_trait(model, HierarchicalTrait)
        return likelihood_trace_hierarchical(model, interval, trace)
    else
        return likelihood_trace_base(model, interval, trace)
    end
end

# Base implementations
function likelihood_rna_base(model::ComposableModel, nRNA::Int, histRNA::Vector)
    # Implementation for base GRSM model
    # This would be similar to the original implementation but using model.traits
end

function likelihood_onoff_base(model::ComposableModel, bins::Vector, ON::Vector, OFF::Vector)
    # Base implementation for ON/OFF likelihood
end

function likelihood_dwelltime_base(model::ComposableModel, bins::Vector{Vector}, dwell_times::Vector{Vector}, dt_types::Vector)
    # Base implementation for dwell time likelihood
end

function likelihood_trace_base(model::ComposableModel, interval::Float64, trace::Any)
    # Base implementation for trace likelihood
end

# Coupled implementations
function likelihood_rna_coupled(model::ComposableModel, nRNA::Int, histRNA::Vector)
    coupling_trait = get_trait(model, CouplingTrait)
    # Implementation for coupled models using coupling_trait information
end

# Grid implementations
function likelihood_rna_grid(model::ComposableModel, nRNA::Int, histRNA::Vector)
    grid_trait = get_trait(model, GridTrait)
    # Implementation for grid models using grid_trait information
end

# Hierarchical implementations
function likelihood_rna_hierarchical(model::ComposableModel, nRNA::Int, histRNA::Vector)
    hier_trait = get_trait(model, HierarchicalTrait)
    # Implementation for hierarchical models using hier_trait information
end

# Helper functions
"""
    prepare_likelihood_components(model::ComposableModel)

Prepare necessary components for likelihood calculation based on model traits.
"""
function prepare_likelihood_components(model::ComposableModel)
    components = Dict()
    
    if has_trait(model, CouplingTrait)
        coupling_trait = get_trait(model, CouplingTrait)
        components[:coupling] = prepare_coupling_components(coupling_trait)
    end
    
    if has_trait(model, GridTrait)
        grid_trait = get_trait(model, GridTrait)
        components[:grid] = prepare_grid_components(grid_trait)
    end
    
    if has_trait(model, HierarchicalTrait)
        hier_trait = get_trait(model, HierarchicalTrait)
        components[:hierarchical] = prepare_hierarchical_components(hier_trait)
    end
    
    components
end

# Add other helper functions as needed 