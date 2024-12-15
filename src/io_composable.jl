# This file is part of StochasticGene.jl

"""
    IO functions for composable trait models and data
"""

using CSV
using DataFrames

"""
    load_data(datapath::String, label::String, gene::String; kwargs...)

Load experimental data and create a ComposableData object with appropriate traits.
Automatically detects data type from file contents and structure.

# Arguments
- `datapath::String`: Path to data file or directory
- `label::String`: Label for the dataset
- `gene::String`: Gene name
- `kwargs...`: Additional keyword arguments for specific data types
"""
function load_data(datapath::String, label::String, gene::String; kwargs...)
    traits = DataTrait[]
    
    # Detect and load RNA data
    if has_rna_data(datapath)
        nRNA, histRNA = load_rna_data(datapath)
        push!(traits, RNATrait(nRNA, histRNA))
    end
    
    # Detect and load ON/OFF data
    if has_onoff_data(datapath)
        bins, ON, OFF = load_onoff_data(datapath)
        push!(traits, OnOffTrait(bins, ON, OFF))
    end
    
    # Detect and load dwell time data
    if has_dwelltime_data(datapath)
        bins, dwell_times, dt_types = load_dwelltime_data(datapath)
        push!(traits, DwellTimeTrait(bins, dwell_times, dt_types))
    end
    
    # Detect and load trace data
    if has_trace_data(datapath)
        interval, trace = load_trace_data(datapath)
        push!(traits, TraceTrait(interval, trace))
    end
    
    isempty(traits) && error("No valid data found in $datapath")
    
    ComposableData(Tuple(traits), label, gene)
end

"""
    save_data(data::ComposableData, outpath::String)

Save ComposableData to files based on its traits.
"""
function save_data(data::ComposableData, outpath::String)
    # Save RNA data
    if has_data_trait(data, RNATrait)
        rna_trait = get_data_trait(data, RNATrait)
        save_rna_data(rna_trait, outpath, data.label, data.gene)
    end
    
    # Save ON/OFF data
    if has_data_trait(data, OnOffTrait)
        onoff_trait = get_data_trait(data, OnOffTrait)
        save_onoff_data(onoff_trait, outpath, data.label, data.gene)
    end
    
    # Save dwell time data
    if has_data_trait(data, DwellTimeTrait)
        dwell_trait = get_data_trait(data, DwellTimeTrait)
        save_dwelltime_data(dwell_trait, outpath, data.label, data.gene)
    end
    
    # Save trace data
    if has_data_trait(data, TraceTrait)
        trace_trait = get_data_trait(data, TraceTrait)
        save_trace_data(trace_trait, outpath, data.label, data.gene)
    end
end

"""
    save_model(model::ComposableModel, outpath::String)

Save ComposableModel to file, including all trait information.
"""
function save_model(model::ComposableModel, outpath::String)
    # Create model info dictionary
    model_info = Dict{String,Any}(
        "type" => "ComposableModel",
        "traits" => String[],
        "parameters" => Dict{String,Any}()
    )
    
    # Add trait information
    if has_trait(model, CouplingTrait)
        coupling = get_trait(model, CouplingTrait)
        push!(model_info["traits"], "coupling")
        model_info["parameters"]["coupling"] = Dict(
            "source_models" => coupling.source_models,
            "source_units" => coupling.source_units,
            "source_states" => coupling.source_states,
            "target_transitions" => coupling.target_transitions,
            "n_coupling_params" => coupling.n_coupling_params
        )
    end
    
    if has_trait(model, GridTrait)
        grid = get_trait(model, GridTrait)
        push!(model_info["traits"], "grid")
        model_info["parameters"]["grid"] = Dict(
            "rate_range" => [grid.rate_range.start, grid.rate_range.stop],
            "noise_range" => [grid.noise_range.start, grid.noise_range.stop],
            "grid_range" => [grid.grid_range.start, grid.grid_range.stop],
            "n_grid" => grid.n_grid
        )
    end
    
    if has_trait(model, HierarchicalTrait)
        hier = get_trait(model, HierarchicalTrait)
        push!(model_info["traits"], "hierarchical")
        model_info["parameters"]["hierarchical"] = Dict(
            "pool" => Dict(
                "nhyper" => hier.pool.nhyper,
                "nrates" => hier.pool.nrates,
                "nparams" => hier.pool.nparams,
                "nindividuals" => hier.pool.nindividuals,
                "ratestart" => hier.pool.ratestart,
                "paramstart" => hier.pool.paramstart,
                "hyperindices" => hier.pool.hyperindices
            )
        )
    end
    
    # Add base model parameters
    model_info["parameters"]["base"] = Dict(
        "rates" => model.rates,
        "Gtransitions" => model.Gtransitions,
        "G" => model.G,
        "R" => model.R,
        "S" => model.S,
        "insertstep" => model.insertstep,
        "nalleles" => model.nalleles,
        "splicetype" => model.splicetype,
        "rateprior" => model.rateprior,
        "proposal" => model.proposal,
        "fittedparam" => model.fittedparam,
        "fixedeffects" => model.fixedeffects,
        "method" => model.method
    )
    
    # Save to JSON file
    open(outpath, "w") do io
        JSON.print(io, model_info, 4)
    end
end

"""
    load_model(filepath::String)

Load ComposableModel from file.
"""
function load_model(filepath::String)
    model_info = JSON.parsefile(filepath)
    
    # Verify it's a composable model
    model_info["type"] == "ComposableModel" || error("Not a ComposableModel file")
    
    # Create traits
    traits = [GRSMTrait()]  # Base trait
    
    # Add additional traits based on saved information
    if "coupling" in model_info["traits"]
        cp = model_info["parameters"]["coupling"]
        push!(traits, CouplingTrait(
            cp["source_models"],
            cp["source_units"],
            cp["source_states"],
            cp["target_transitions"],
            cp["n_coupling_params"]
        ))
    end
    
    if "grid" in model_info["traits"]
        gp = model_info["parameters"]["grid"]
        push!(traits, GridTrait(
            gp["rate_range"][1]:gp["rate_range"][2],
            gp["noise_range"][1]:gp["noise_range"][2],
            gp["grid_range"][1]:gp["grid_range"][2],
            gp["n_grid"]
        ))
    end
    
    if "hierarchical" in model_info["traits"]
        hp = model_info["parameters"]["hierarchical"]["pool"]
        pool = Pool(
            hp["nhyper"],
            hp["nrates"],
            hp["nparams"],
            hp["nindividuals"],
            hp["ratestart"],
            hp["paramstart"],
            hp["hyperindices"]
        )
        push!(traits, HierarchicalTrait(pool))
    end
    
    # Get base parameters
    bp = model_info["parameters"]["base"]
    
    # Create model
    ComposableModel(
        Tuple(traits),
        bp["rates"],
        bp["Gtransitions"],
        bp["G"],
        bp["R"],
        bp["S"],
        bp["insertstep"],
        bp["nalleles"],
        bp["splicetype"],
        bp["rateprior"],
        bp["proposal"],
        bp["fittedparam"],
        bp["fixedeffects"],
        bp["method"],
        make_components(traits),  # This would need to be implemented
        make_reporter(traits)     # This would need to be implemented
    )
end

# Helper functions for data detection and loading
function has_rna_data(path::String)
    # Implementation to detect RNA data files
end

function has_onoff_data(path::String)
    # Implementation to detect ON/OFF data files
end

function has_dwelltime_data(path::String)
    # Implementation to detect dwell time data files
end

function has_trace_data(path::String)
    # Implementation to detect trace data files
end

function load_rna_data(path::String)
    # Implementation to load RNA data
end

function load_onoff_data(path::String)
    # Implementation to load ON/OFF data
end

function load_dwelltime_data(path::String)
    # Implementation to load dwell time data
end

function load_trace_data(path::String)
    # Implementation to load trace data
end

# Helper functions for data saving
function save_rna_data(trait::RNATrait, outpath::String, label::String, gene::String)
    # Implementation to save RNA data
end

function save_onoff_data(trait::OnOffTrait, outpath::String, label::String, gene::String)
    # Implementation to save ON/OFF data
end

function save_dwelltime_data(trait::DwellTimeTrait, outpath::String, label::String, gene::String)
    # Implementation to save dwell time data
end

function save_trace_data(trait::TraceTrait, outpath::String, label::String, gene::String)
    # Implementation to save trace data
end

# Results handling functions
"""
    save_results(model::ComposableModel, data::ComposableData, results::Dict, outpath::String)

Save fitting or simulation results for composable models.
"""
function save_results(model::ComposableModel, data::ComposableData, results::Dict, outpath::String)
    # Create results directory if it doesn't exist
    mkpath(outpath)
    
    # Save model
    save_model(joinpath(outpath, "model.json"), model)
    
    # Save data
    save_data(data, outpath)
    
    # Save results based on model and data traits
    if has_trait(model, HierarchicalTrait)
        save_hierarchical_results(model, results, outpath)
    end
    
    if has_trait(model, GridTrait)
        save_grid_results(model, results, outpath)
    end
    
    if has_trait(model, CouplingTrait)
        save_coupling_results(model, results, outpath)
    end
    
    # Save basic results
    save_basic_results(model, results, outpath)
end

function save_hierarchical_results(model::ComposableModel, results::Dict, outpath::String)
    # Implementation for hierarchical results
end

function save_grid_results(model::ComposableModel, results::Dict, outpath::String)
    # Implementation for grid results
end

function save_coupling_results(model::ComposableModel, results::Dict, outpath::String)
    # Implementation for coupling results
end

function save_basic_results(model::ComposableModel, results::Dict, outpath::String)
    # Implementation for basic results
end

# Add other IO-related functions as needed 